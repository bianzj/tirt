import math
from datetime import date, timedelta

from .aerodynamic import GRAVITY, aerodynamic_resistance
from .biochemical import suresist, suresist_urban
from .control import KAPPA
from .eb import CP, RHOA, rebalance, rebalance_urban
from .evapo import evapotranspiration, urban_evapotranspiration, update_thermal_history
from .radiation import (
    flux_weights,
    merge_radiation,
    prepare_longwave,
    prepare_shortwave,
    radiative_transfer_diagnostics,
    urban_canopy_fraction,
)
from .solar import slope_direction, solar_angles
from .structs import FluxStructure

# Convergence tolerance of the Monin-Obukhov iteration in z/L.
OBUKHOV_ZETA_TOLERANCE = 0.02


def run_node(settings, meteo, node, year, doy, latitude, longitude, dem=0.0, vza=0.0, vaa=0.0, initial_state=None, max_iter=100, utc_offset_hours=0.0):
    # Meteorological time is local time; solar geometry is evaluated in UTC.
    local_time = float(meteo["time"])
    utc_time = local_time - float(utc_offset_hours) / 24.0
    local_date = date(int(year), 1, 1) + timedelta(days=int(doy) - 1 + math.floor(local_time) - int(doy))
    utc_date = local_date + timedelta(days=math.floor(utc_time) - math.floor(local_time))
    utc_doy = (utc_date - date(utc_date.year, 1, 1)).days + 1
    global_sza, global_saa = solar_angles(
        utc_date.year, utc_doy, utc_time, latitude, longitude, timezone=0.0
    )
    global_vza, global_vaa = float(vza), float(vaa)
    slope, aspect = float(settings.get("slope_angle", 0.0)), float(settings.get("slope_aspect", 0.0))
    if settings.get("surface_model") == "slope":
        sza, saa, sunlit_surface = slope_direction(global_sza, global_saa, slope, aspect)
        vza, vaa, _ = slope_direction(global_vza, global_vaa, slope, aspect)
    else:
        sza, saa = global_sza, global_saa
        sunlit_surface = global_sza < 90.0 if settings.get("surface_model") == "terrain" else True
    is_urban = settings.get("surface_model") == "urban"
    soil_initial = float(settings.get("soil_temperature_c", 0.02)) + 273.15
    urban_initial = float(settings.get("urban_temperature_c", settings.get("soil_temperature_c", 0.02))) + 273.15
    wall_layers = max(int(settings.get("urban_wall_layers", 3)), 1)
    wall_boundary = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
    wall_profile = [wall_boundary] * (wall_layers + 1)
    surface_initial = urban_initial if is_urban else soil_initial
    state = {
        "Tss": soil_initial, "Tsh": soil_initial, "Tls": surface_initial, "Tlh": surface_initial,
        "Troof_sunlit": surface_initial, "Troof_shaded": surface_initial,
        "Twall_sunlit": surface_initial, "Twall_shaded": surface_initial,
        "L": 1e9, "H_scene": 0.0,
        "dt_seconds": float(settings.get("dtime", 1800.0)),
        "surface_history_sunlit": [soil_initial], "surface_history_shaded": [soil_initial],
        "surface_history_roof_sunlit": [urban_initial], "surface_history_roof_shaded": [urban_initial],
        "surface_history_wall_sunlit": [urban_initial], "surface_history_wall_shaded": [urban_initial],
        "surface_history_street_sunlit": [soil_initial], "surface_history_street_shaded": [soil_initial],
        "depth_profile_sunlit": [soil_initial] * 4, "depth_profile_shaded": [soil_initial] * 4,
        "depth_profile_roof_sunlit": [urban_initial] * 4, "depth_profile_roof_shaded": [urban_initial] * 4,
        "depth_profile_wall_sunlit": wall_profile[:], "depth_profile_wall_shaded": wall_profile[:],
        "depth_profile_street_sunlit": [soil_initial] * 4, "depth_profile_street_shaded": [soil_initial] * 4,
        "cs": settings.get("ca", 380.0), "ci": settings.get("ca", 380.0),
    }
    if initial_state:
        for key in ("Tss", "Tsh", "Tls", "Tlh", "Troof_sunlit", "Troof_shaded", "Twall_sunlit", "Twall_shaded", "L", "H_scene", "dt_seconds",
                    "surface_history_sunlit", "surface_history_shaded",
                    "surface_history_roof_sunlit", "surface_history_roof_shaded",
                    "surface_history_wall_sunlit", "surface_history_wall_shaded",
                    "surface_history_street_sunlit", "surface_history_street_shaded",
                    "depth_profile_sunlit", "depth_profile_shaded",
                    "depth_profile_roof_sunlit", "depth_profile_roof_shaded",
                    "depth_profile_wall_sunlit", "depth_profile_wall_shaded",
                    "depth_profile_street_sunlit", "depth_profile_street_shaded"):
            if key in initial_state:
                state[key] = list(initial_state[key]) if isinstance(initial_state[key], (list, tuple)) else initial_state[key]
        if "time" in initial_state:
            state["dt_seconds"] = max(abs(float(meteo["time"]) - float(initial_state["time"])) * 86400.0, 1.0)
    # Keep the wall profile shape and its indoor boundary fixed when carrying state forward.
    for profile_key in ("depth_profile_wall_sunlit", "depth_profile_wall_shaded"):
        profile = list(state.get(profile_key, ()))
        state[profile_key] = (profile[:wall_layers + 1] + wall_profile)[:wall_layers + 1]
        state[profile_key][-1] = wall_boundary
    # Wall illumination must be available during the first radiation iteration.
    state.update(sza=sza, saa=saa, vza=vza, vaa=vaa)
    shortwave = prepare_shortwave(
        settings, meteo, sza, saa, vza, vaa, sunlit=sunlit_surface
    )
    iterations = 0
    closed = False
    netrad = {}
    # Per-ground weights of the soil (sunlit, shaded) and leaf (sunlit, shaded) component fluxes.
    weights = flux_weights(settings, sza, saa, sunlit_surface)
    scene_weights = _urban_plan_weights(settings) if is_urban else weights

    def evaluate_flux(rebalance_state=True):
        nonlocal netrad
        longwave = prepare_longwave(settings, meteo, state, sza, saa, sunlit_surface)
        netrad = merge_radiation(shortwave, longwave)
        aerodynamic_resistance(settings, meteo, state, is_urban)
        if is_urban:
            # Urban aerodynamic_resistance computes separate canyon and leaf
            # transfer paths; only the asphalt street shares its resistance
            # with the legacy soil state used by the radiation/flux helpers.
            state["raa_soil"] = state["raa_street"]
            suresist_urban(settings, meteo, state, netrad)
            _ensure_urban_state(state)
            urban_evapotranspiration(settings, meteo, state, netrad)
            balanced = rebalance_urban(settings, meteo, state) if rebalance_state else False
        else:
            suresist(settings, meteo, state, netrad)
            evapotranspiration(settings, meteo, state, netrad)
            balanced = rebalance(settings, meteo, state) if rebalance_state else False
        return balanced

    obukhov = _ObukhovIteration(state, settings.get("z", 15.0))
    # The input surface option is authoritative; surface_type=13 is kept for
    # Keep the input surface option authoritative for all scenarios.
    for kiter in range(max_iter):
        balanced = evaluate_flux()
        if is_urban:
            _ensure_urban_state(state)
        iterations = kiter
        # Outer Monin-Obukhov iteration: move L only once the energy balance has closed for the current L.
        if balanced and obukhov.update(_scene_sensible_heat(state, scene_weights, is_urban), meteo):
            closed = True
            break
    # Report fluxes evaluated at the solved temperatures rather than at the last Newton iterate.
    evaluate_flux(rebalance_state=False)
    state["H_scene"] = _scene_sensible_heat(state, scene_weights, is_urban)
    for key in [key for key in state if key.startswith("_eb_")]:
        del state[key]

    if is_urban:
        _ensure_urban_state(state)
    update_thermal_history(settings, state)
    state.update(
        time=meteo["time"], sza=sza, saa=saa, vza=vza, vaa=vaa,
        sza_global=global_sza, saa_global=global_saa,
        vza_global=global_vza, vaa_global=global_vaa,
        slope_angle=slope, slope_aspect=aspect,
        sunlit_surface=sunlit_surface,
        iterations=iterations, closed=closed,
        soil_sunlit_fraction=weights[0], soil_shaded_fraction=weights[1],
        lai_sunlit=weights[2], lai_shaded=weights[3],
    )
    state.update(radiative_transfer_diagnostics(settings, state, vza, vaa))
    state.update(
        shortwave_direct=netrad.get("shortwave_direct", 0.0),
        shortwave_diffuse=netrad.get("shortwave_diffuse", 0.0),
        shortwave_global=netrad.get("shortwave_global", 0.0),
        shortwave_direct_incident=netrad.get("shortwave_direct_incident", 0.0),
        shortwave_diffuse_incident=netrad.get("shortwave_diffuse_incident", 0.0),
        shortwave_diffuse_fraction=netrad.get("shortwave_diffuse_fraction", 0.0),
        longwave_down=netrad.get("longwave_down", 0.0),
        tir_scatter_soil=netrad.get("tir_scatter_soil", 0.0),
        tir_scatter_leaf=netrad.get("tir_scatter_leaf", 0.0),
    )
    if not is_urban:
        for component in ("soil", "leaf"):
            for suffix, sunlit in (("sunlit", True), ("shaded", False)):
                direct = float(netrad.get(f"directVrad_{component}", 0.0)) if sunlit else 0.0
                diffuse = float(netrad.get(f"diffuseVrad_{component}", 0.0))
                longwave = float(netrad.get(f"diffuseTrad_{component}", 0.0))
                if sunlit:
                    longwave += float(netrad.get(f"directTrad_{component}", 0.0))
                state[f"SW_direct_{component}_{suffix}"] = direct
                state[f"SW_diffuse_{component}_{suffix}"] = diffuse
                state[f"LW_input_{component}_{suffix}"] = longwave
                state[f"radiation_input_{component}_{suffix}"] = direct + diffuse + longwave
    return FluxStructure.from_mapping(state)


def _urban_plan_weights(settings):
    """Surface area per unit plan area for roof, wall, open street and street canopy."""
    density = max(float(settings.get("urban_density", 0.0015)), 0.0)
    length = max(float(settings.get("urban_building_length", 10.0)), 0.0)
    width = max(float(settings.get("urban_building_width", 10.0)), 0.0)
    height = max(float(settings.get("urban_building_height", 20.0)), 0.0)
    roof = min(density * length * width, 0.99)
    wall = density * 2.0 * height * (length + width)
    canopy = urban_canopy_fraction(settings, 0.0)
    street = 1.0 - roof
    return {"roof": roof, "wall": wall, "street": street * (1.0 - canopy), "leaf": street * canopy}


def _scene_sensible_heat(state, weights, is_urban):
    if is_urban:
        return sum(
            weight * 0.5 * (float(state.get(f"H_{component}sunlit", 0.0)) + float(state.get(f"H_{component}shaded", 0.0)))
            for component, weight in weights.items()
        )
    keys = ("H_soilsunlit", "H_soilshaded", "H_leafsunlit", "H_leafshaded")
    return sum(float(weight) * float(state.get(key, 0.0)) for weight, key in zip(weights, keys))


class _ObukhovIteration:
    """Solve x = 1/L = F(x) by Illinois regula falsi; |x| <= 1 keeps the root bracketed in [-1, 1]."""

    def __init__(self, state, reference):
        self.state = state
        self.reference = max(float(reference), 1.0)
        self.low = (-1.0, None)
        self.high = (1.0, None)
        self.side = 0

    def update(self, sensible, meteo):
        state = self.state
        state["H_scene"] = sensible
        ta = float(meteo["air_temperature_c"]) + 273.15
        ustar = max(float(state.get("ustar", 0.0)), 1.0e-3)
        # Obukhov (1946): L = -rho*cp*u*^3*Ta/(kappa*g*H), limited to |L| >= 1 m.
        target = min(max(-KAPPA * GRAVITY * sensible / (RHOA * CP * ustar ** 3 * ta), -1.0), 1.0)
        x = 1.0 / float(state.get("L", 1e9))
        g = target - x
        side = -1 if g >= 0.0 else 1
        if side < 0:
            self.low = (x, g)
        else:
            self.high = (x, g)
        tolerance = OBUKHOV_ZETA_TOLERANCE / self.reference
        if abs(g) <= tolerance or self.high[0] - self.low[0] <= tolerance:
            return True
        (x_low, g_low), (x_high, g_high) = self.low, self.high
        if g_low is not None and g_high is not None:
            if side == self.side:
                if side < 0:
                    g_high *= 0.5
                    self.high = (x_high, g_high)
                else:
                    g_low *= 0.5
                    self.low = (x_low, g_low)
            new = (x_low * g_high - x_high * g_low) / (g_high - g_low)
        else:
            # Fixed-point step; it always lands inside (low, high] and so builds the bracket.
            new = target
        self.side = side
        state["L"] = 1.0 / new if abs(new) > 1.0e-9 else 1e9
        # Residual signs recorded for the temperature damping refer to the previous L.
        for key in [key for key in state if key.startswith("_eb_res_")]:
            del state[key]
        return False


def _ensure_urban_state(state):
    state.setdefault("Troof_sunlit", state["Tls"])
    state.setdefault("Troof_shaded", state["Tlh"])
    state.setdefault("Twall_sunlit", state["Tls"])
    state.setdefault("Twall_shaded", state["Tlh"])
    state["Trs"] = state["Troof_sunlit"]
    state["Trh"] = state["Troof_shaded"]
    state["Tws"] = state["Twall_sunlit"]
    state["Twh"] = state["Twall_shaded"]
    state["Tstreet_sunlit"] = state["Tss"]
    state["Tstreet_shaded"] = state["Tsh"]
