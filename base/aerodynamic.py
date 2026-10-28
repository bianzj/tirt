import math

from .control import KAPPA
from .eb import CP, air_density

BARE_SOIL_ROUGHNESS = 0.01
# Below this LAI the canopy exchange network blends smoothly into bare soil.
SPARSE_CANOPY_LAI = 0.5


def _set_stability_reference(state, length, reference_height):
    state['stability_reference_height'] = reference_height
    # _stability_functions caps the stable reference-height zeta at 1.
    state['stability_limited'] = bool(length > 0.0 and reference_height / length > 1.0)


def _stability_functions(length, reference_height=None):
    """Businger--Dyer unstable and Webb stable integrated similarities.

    The reference-height stability is limited once, for the whole profile.
    Clipping each height independently would cancel both roughness corrections
    in strong stability and incorrectly restore neutral heat transfer.
    """
    stable_inverse = 1.0 / length if length > 0.0 else 0.0
    if reference_height is not None:
        stable_inverse = min(stable_inverse, 1.0 / max(float(reference_height),1e-6))

    def psim(value):
        if length < 0:
            x = max(1.0 - 16.0 * value / length, 1e-8) ** 0.25
            return 2.0 * math.log((1.0 + x) / 2.0) + math.log((1.0 + x * x) / 2.0) - 2.0 * math.atan(x) + math.pi / 2.0
        if length > 0:
            return -5.0 * max(value, 0.0) * stable_inverse
        return 0.0

    def psih(value):
        if length < 0:
            x = max(1.0 - 16.0 * value / length, 1e-8) ** 0.25
            return 2.0 * math.log((1.0 + x * x) / 2.0)
        if length > 0:
            return -5.0 * max(value, 0.0) * stable_inverse
        return 0.0

    def phstar(value, reference, displacement):
        if length < 0:
            x = max(1.0 - 16.0 * value / length, 1e-8) ** 0.25
            return ((value - displacement) / (reference - displacement)) * ((x * x - 1.0) / (x * x + 1.0))
        if length > 0:
            return -5.0 * max(value, 0.0) * stable_inverse
        return 0.0

    return psim, psih, phstar


def _bulk_resistance(reference_height, roughness_length, wind_speed, length, psim, psih):
    """Bulk heat-transfer resistance from a local wind reference to a surface."""
    reference_height = max(float(reference_height), 1.01)
    roughness_length = min(max(float(roughness_length), 1.0e-4), reference_height / 2.0)
    wind_speed = max(float(wind_speed), 0.05)
    heat_roughness = max(roughness_length * 0.1, 1.0e-5)
    momentum_log = (
        math.log(reference_height / roughness_length)
        - psim(reference_height)
        + psim(roughness_length)
    )
    heat_log = (
        math.log(reference_height / heat_roughness)
        - psih(reference_height)
        + psih(heat_roughness)
    )
    return max(momentum_log * heat_log / (KAPPA * KAPPA * wind_speed), 1.0e-3)


def aerodynamic_resistance(settings, meteo, state, urban=False):
    """Complete equivalent of old/aerodynamic.cpp::aeresist and aeresist_urban."""
    lai = max(float(settings["lai"]), 0.0)
    hc = max(float(settings["height"]), 0.1)
    wind = max(float(meteo["wind_speed"]), float(settings.get("min_wind_speed", 0.05)), 0.01)
    # Convective gustiness keeps free-convective daytime exchange finite.
    buoyancy = max(float(state.get("scene_sensible_heat_W_m2", 0.0)), 0.0)
    ta = float(meteo["air_temperature_c"]) + 273.15
    wstar = (9.81 * buoyancy * 1000.0 / (air_density(meteo) * CP * ta)) ** (1.0 / 3.0) if buoyancy > 0.0 else 0.0
    wind = math.hypot(wind, wstar)
    z = max(float(settings["z"]), hc + 0.1)
    length = float(state.get("L", 1e9))
    profile_reference = z + _canopy_weight(lai) * (max(z, 2.5 * hc) - z)
    _set_stability_reference(state, length, profile_reference)
    psim, psih, phstar = _stability_functions(length, profile_reference)

    if urban:
        hc = max(float(settings.get("urban_building_height", 20.0)), 0.1)
        width = max(float(settings.get("urban_building_width", 10.0)), 0.1)
        building_length = max(float(settings.get("urban_building_length", 10.0)), 0.1)
        density = max(float(settings.get("urban_density", 0.0015)), 0.0)
        scene_height = max(hc, float(settings.get("height", 0.0)))
        z = max(float(settings.get("z", 15.0)), scene_height + 0.1)
        profile_reference = z
        _set_stability_reference(state, length, profile_reference)
        psim, psih, phstar = _stability_functions(length, profile_reference)
        state["z"] = z
        bai = 1.0
        pixel_area = 1e6
        building_area = building_length * width
        count = max(density * pixel_area, 1e-6)
        lambda_p = min(max(count * building_area / pixel_area, 0.0), 0.99)
        displacement = (1.0 + 4.0 ** (-lambda_p) * (lambda_p - 1.0)) * hc
        roughness = max(float(settings.get("urban_roughness_length", 0.15)), 0.01)
        numerator = math.log(max((z - displacement) / roughness, 1.01)) - psim(z - displacement) + psim(roughness)
        ustar = max(0.001, KAPPA * wind / max(numerator, 1e-6))
        zr = 2.5 * hc
        pm_h = psim(hc - displacement)
        uh = max(ustar / KAPPA * (math.log(max((hc - displacement) / roughness, 1.01)) - pm_h + psim(roughness)), 0.01)
        town_z0 = min(5.0, hc * 0.1)
        k = 0.35
        r_ratio = 0.74
        if length > 0:
            psi_m1 = -4.7 * z / length
            psi_m2 = -4.7 * town_z0 / length / r_ratio
        elif length < 0:
            x = max(1.0 - 15.0 * z / length, 1e-8) ** 0.25
            psi_m1 = math.log(((1.0 + x) / 2.0) ** 2 * (1.0 + x * x) / 2.0) - 2.0 * math.atan(x) + math.pi / 2.0
            y = max(1.0 - 9.0 * town_z0 / length, 1e-8) ** 0.5
            psi_m2 = math.log(((1.0 + y) / 2.0) ** 2)
        else:
            psi_m1 = psi_m2 = 0.0
        denominator = max(math.log(max(z / roughness, 1.01)) - (psi_m1 - psi_m2), 1e-6)
        cd = ((k * k) / r_ratio) / denominator**4
        w_can = math.sqrt(max(cd, 0.0)) * abs(uh)
        exp_term = math.exp(-0.25 * hc / max(width, 0.1))
        log1 = math.log(max(hc / 3.0 / town_z0, 1.01))
        log2 = math.log(max(hc / 3.0 / town_z0, 1.01))
        u_can = (2.0 / math.pi) * exp_term * log1 / max(log2, 1e-6) * abs(uh)
        ustar = max(0.001, KAPPA * wind / max(math.log(max((z - displacement) / roughness, 1.01)) - psim(z - displacement) + psim(roughness), 1e-6))

        # Use the same momentum/heat bulk-transfer formulation for every
        # urban component. Roof exchange is referenced to the above-canopy
        # wind; street and wall exchange use the slower canyon wind and a
        # geometry-derived canyon roughness length.
        urban_rbs = max(float(settings.get("urban_rbs", 10.0)), 0.0)
        cap = _resistance_cap(settings)
        roof_reference = max(z - displacement, 1.01)
        roof_bulk = min(_bulk_resistance(
            roof_reference, roughness, wind, length, psim, psih
        ), cap)
        roof = roof_bulk + urban_rbs

        canyon_scale = max(
            float(settings.get("urban_canyon_roughness_ratio", 0.05))
            * math.sqrt(width * building_length),
            0.05,
        )
        canyon_reference = max(2.0 * canyon_scale, 1.01)  # in-canyon scale, matches u_can
        canyon_bulk = min(_bulk_resistance(
            canyon_reference, canyon_scale, u_can, length, psim, psih
        ), cap)
        canyon_ratio = max(float(settings.get("urban_canyon_min_ratio", 1.05)), 1.0)
        canyon = max(canyon_bulk + urban_rbs, roof * canyon_ratio)

        lai_value = max(float(settings.get("lai", 0.0)), 0.0)
        if lai_value > 0.0:
            leaf_decay = max(float(settings.get("urban_leaf_wind_decay", 0.5)), 0.0)
            leaf_wind = max(u_can * math.exp(-leaf_decay * lai_value), 0.05)
            leaf_boundary_factor = max(
                float(settings.get("urban_leaf_boundary_factor", 70.0)), 0.0
            )
            leaf_boundary = min(leaf_boundary_factor / lai_value * math.sqrt(
                max(float(settings.get("leaf_width", 0.2)) / leaf_wind, 1.0e-8)
            ), cap)
            leaf = canyon + leaf_boundary
        else:
            leaf_wind = 0.0
            leaf_boundary = 0.0
            leaf = 0.0

        state.update(
            ustar=ustar,
            u_canopy=max(u_can, 0.0),
            u_leaf=max(leaf_wind, 0.0),
            d_urban=displacement,
            z0_urban=roughness,
            z0_canyon=canyon_scale,
            r_bulk_roof=roof_bulk,
            r_bulk_canyon=canyon_bulk,
            r_boundary_urban=urban_rbs,
            r_boundary_leaf=leaf_boundary,
            raa_roof=roof,
            raa_street=canyon,
            raa_wall=canyon,
            raa_leaf=leaf,
            raa_soil=canyon,
        )
        return state

    if lai <= 0.0:
        ustar, raa_soil = _bare_soil(z, wind, length, psim, psih, _resistance_cap(settings))
        state.update(ustar=ustar, raa_leaf=0.0, raa_soil=raa_soil,
                     aerodynamic_reference_height=z, aerodynamic_displacement=0.0,
                     aerodynamic_roughness_length=BARE_SOIL_ROUGHNESS, aerodynamic_canopy_weight=0.0)
        return state

    ustar, raa_leaf, raa_soil, displacement, roughness = _canopy(settings, lai, hc, z, wind, length,
                                                                 profile_reference, psim, psih, phstar)
    weight = _canopy_weight(lai)
    if weight < 1.0:
        bare_ustar, bare_soil = _bare_soil(z, wind, length, psim, psih, _resistance_cap(settings))
        ustar = weight * ustar + (1.0 - weight) * bare_ustar
        raa_soil = weight * raa_soil + (1.0 - weight) * bare_soil
        displacement *= weight
        roughness = math.exp(weight * math.log(roughness) + (1.0 - weight) * math.log(BARE_SOIL_ROUGHNESS))
    state.update(ustar=ustar, raa_leaf=raa_leaf, raa_soil=raa_soil,
                 aerodynamic_reference_height=z, aerodynamic_displacement=displacement,
                 aerodynamic_roughness_length=roughness, aerodynamic_canopy_weight=weight)
    return state


def _canopy_weight(lai):
    """Smoothstep weight of the canopy exchange network; bare-soil exchange as LAI -> 0."""
    x = min(max(float(lai) / SPARSE_CANOPY_LAI, 0.0), 1.0)
    return x * x * (3.0 - 2.0 * x)


def _bare_soil(z, wind, length, psim, psih, cap):
    roughness = BARE_SOIL_ROUGHNESS
    ustar = max(0.001, KAPPA * wind / max(math.log(max(z / roughness, 1.01)) - psim(z) + psim(roughness), 1e-6))
    return ustar, min(_bulk_resistance(z, roughness, wind, length, psim, psih), cap)


def _resistance_cap(settings):
    # SCOPE (resistances.m) limits each aerodynamic resistance segment to 400 s m-1.
    cap = float(settings.get("max_aerodynamic_resistance", 400.0))
    return cap if cap > 0.0 else math.inf


def _canopy(settings, lai, hc, z, wind, length, profile_reference, psim, psih, phstar):
    # Raupach (1994) d and z0 as in SCOPE: U/u*(h) = max(3.3, (Cs + Cr*LAI/2)^-1/2),
    # so z0 decreases toward the substrate-drag limit as LAI -> 0.
    sq = math.sqrt(max(float(settings["aero_cd1"]) * lai / 2.0, 1e-8))
    g1 = max(3.3, max(float(settings["aero_cssoil"]) + float(settings["aero_cr"]) * lai / 2.0, 1e-8) ** -0.5)
    n = float(settings["aero_cd"]) * lai / (2.0 * KAPPA * KAPPA)
    zr = 2.5 * hc
    displacement = hc * (1.0 - (1.0 - math.exp(-sq)) / sq)
    roughness = max((hc - displacement) * math.exp(-KAPPA * g1 + float(settings["aero_psicor"])), 1e-4)
    if hc > z:
        ustar_temp = max(0.001, KAPPA * wind / max(math.log(max((z - displacement) / roughness, 1.01)) - psim(z - displacement), 1e-6))
        wind = max(0.01, ustar_temp / KAPPA * (math.log(max((hc - displacement) / roughness, 1.01)) - psim(hc - displacement)))
    ustar = max(0.001, KAPPA * wind / max(math.log(max((z - displacement) / roughness, 1.01)) - psim(z - displacement) + psim(roughness), 1e-6))
    kh = KAPPA * ustar * max(zr - displacement, 1e-6)
    if length < 0:
        kh = KAPPA * ustar * max(zr - displacement, 1e-6) * math.sqrt(max(1.0 - 16.0 * (hc - displacement) / length, 1e-6))
    elif length > 0:
        kh = KAPPA * ustar * max(zr - displacement, 1e-6) / (1.0 + 5.0 * (hc-displacement)*min(1.0/length,1.0/profile_reference))
    canopy_wind = max(ustar / KAPPA * (math.log(max((hc - displacement) / roughness, 1.01)) - psim(hc - displacement)), 0.01)
    rai = 1.0 / (KAPPA * ustar) * (math.log(max((z - displacement) / max(zr - displacement, 1e-6), 1.01)) - psih(z - displacement) + psih(zr - displacement)) if z > zr else 0.0
    rar = 1.0 / (KAPPA * ustar) * ((zr - hc) / max(zr - displacement, 1e-6)) - phstar(zr, zr, displacement) + phstar(hc, zr, displacement)

    def log_ratio(argument):
        argument = max(argument, 1e-8)
        if argument > 50.0:
            return math.log1p(-2.0 * math.exp(-argument))
        e = math.exp(argument)
        return math.log((e - 1.0) / (e + 1.0))

    kh = max(kh, 1e-6)
    rac = hc * math.sinh(n) / max(n * kh, 1e-6) * (log_ratio(n) - log_ratio(n * (roughness + displacement) / hc)) if n else 0.0
    wind_leaf = max(canopy_wind * math.exp(n * ((roughness + displacement) / hc - 1.0)), 0.01)
    boundary = 70.0 / lai * math.sqrt(max(float(settings["leaf_width"]) / wind_leaf, 1e-8))
    soil = hc * math.sinh(n) / max(n * kh, 1e-6) * (log_ratio(n * (roughness + displacement) / hc) - log_ratio(n * BARE_SOIL_ROUGHNESS / hc)) if n else 0.0
    cap = _resistance_cap(settings)
    above = min(rai + rar + rac, cap)
    return (ustar, max(above + min(boundary, cap), 1e-3),
            max(above + min(soil + float(settings["aero_rbs"]), cap), 1e-3), displacement, roughness)


__all__ = ["aerodynamic_resistance"]
