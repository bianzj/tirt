import math
import numpy as np

from .proportion import proportion_bidirectional_crown_one, proportion_bidirectional_hom_one, proportion_bidirectional_row_one
from .terrain_rt import compound_terrain_proportions
from .soilopt import hapke_reflectance

SOLAR_CONSTANT = 1361.0


def relative_azimuth(saa, vaa):
    """Return the circular view-sun azimuth separation in [0, 180] degrees."""
    difference = np.mod(np.asarray(vaa, dtype=float) - np.asarray(saa, dtype=float), 360.0)
    return np.minimum(difference, 360.0 - difference)


def _vegetation_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True, row_azimuth=None):
    lai = max(float(settings["lai"]), 0.0)
    if lai <= 0.0:
        # Bare ground still has to honour the illumination state: a shaded facet is not sunlit soil.
        return np.asarray([0.0, 1.0, 0.0, 0.0]) if not sunlit else np.asarray([1.0, 0.0, 0.0, 0.0])
    sza_ = np.asarray([min(max(float(sza), 0.0), 89.9)])
    vza_ = np.asarray([min(max(float(vza), 0.0), 89.9)])
    saa_ = np.asarray([float(saa)])
    vaa_ = np.asarray([float(vaa)])
    vsa_ = relative_azimuth(saa_, vaa_)
    vegetation_model = str(settings.get("vegetation_model", "")).strip().lower()
    if vegetation_model in {"crown", "forest"} or settings["surface_type"] in (1, 2, 3, 4, 5):
        values = proportion_bidirectional_crown_one(lai, max(settings["stand_density"], 1e-6), settings["hspot"], max(settings["height"], 1e-6), max(settings["width"], 1e-6), vza_, sza_, vsa_)
    elif vegetation_model in {"row", "raw"} or settings["distribution"] == 2 or settings["surface_type"] in (12, 14):
        values = proportion_bidirectional_row_one(
            lai,
            settings["hspot"],
            settings["row_width"],
            settings["row_blank"],
            settings["row_height"],
            vza_,
            vaa_,
            sza_,
            saa_,
            np.asarray([settings["row_azimuth"] if row_azimuth is None else row_azimuth]),
        )
    else:
        values = proportion_bidirectional_hom_one(lai, settings["hspot"], vza_, sza_, vsa_)
    values = np.maximum(np.asarray(values, dtype=float).reshape(4), 0.0)
    if not sunlit:
        values = np.asarray([0.0, values[0] + values[1], 0.0, values[2] + values[3]])
    return values / values.sum() if values.sum() else np.asarray([1.0, 0.0, 0.0, 0.0])


def _component_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    if (
        settings.get("surface_model") == "terrain"
        and int(settings.get("terrain_type", 2)) == 2
    ):
        return compound_terrain_proportions(
            settings,
            sza,
            saa,
            vza,
            vaa,
            lambda local_sza, local_saa, local_vza, local_vaa, local_sunlit, local_row_azimuth=None: _vegetation_proportions(
                settings,
                local_sza,
                local_saa,
                local_vza,
                local_vaa,
                local_sunlit,
                local_row_azimuth,
            ),
            sunlit,
        )
    return _vegetation_proportions(settings, sza, saa, vza, vaa, sunlit)


def optical_properties(settings):
    rho = min(max(float(settings["leaf_reflectance"]), 0.0), 0.99)
    tau = min(max(float(settings["leaf_transmittance"]), 0.0), 0.99 - rho)
    soil = min(max(float(settings["soil_reflectance"]), 0.0), 0.99)
    return rho, tau, soil, 1.0 - rho - tau, 1.0 - soil


def directional_reflectance(
    settings, sza, saa, vza, vaa, leaf_rho, leaf_tau, soil_rho,
    sunlit=True, use_hapke=None, diffuse_fraction=None,
):
    """Optical BRF/HDRF with direct shadows and a diffuse transport field.

    Mixed illumination defaults to 80% direct + 20% diffuse, matching the shared
    shortwave convention; pass ``diffuse_fraction`` explicitly to override it.
    """
    from .vnir import canopy_reflectance
    model = str(settings.get('vegetation_model','hom')).lower()
    if diffuse_fraction is None:
        configured = float(settings.get('shortwave_diffuse_fraction', 0.2))
        diffuse_fraction = 0.2 if configured < 0.0 else configured
    if use_hapke is None:
        use_hapke = int(settings.get('hapke_enable',settings.get('hapke_optical_enable',0))) == 1
    is_crown = (model in ('crown','forest') or settings.get('surface_type',0) in (1,2,3,4,5))
    from . import canopy_closure
    bare_terrain = (canopy_closure.closure(settings) != 'legacy' and not use_hapke
                    and settings.get('surface_model','plane') == 'terrain'
                    and int(settings.get('terrain_type', 2)) == 2
                    and not is_crown and float(settings.get('lai',0)) <= 0.0)
    if (canopy_closure.closure(settings) != 'legacy' and settings.get('surface_model','plane') == 'terrain'
            and int(settings.get('terrain_type', 2)) == 2 and not use_hapke and model in ('hom','row','raw')
            and float(settings.get('lai',0)) > 0):
        from .optical_terrain_standard import reflectance
        return reflectance(dict(settings, vegetation_model='row' if model in ('row','raw') else 'hom'),
                           float(sza),float(saa),float(vza),float(vaa),
                           float(leaf_rho),float(leaf_tau),float(soil_rho),diffuse_fraction,sunlit)
    if (is_crown and settings.get('surface_model','plane') == 'terrain'
            and int(settings.get('terrain_type', 2)) == 2 and not use_hapke):
        method = settings.get('optical_terrain_method', 'facet_exchange')
        common=settings.get('forest_transport')
        if common not in (None,'legacy','reciprocal_exchange'):
            raise ValueError('Unknown common forest transport: '+str(common))
        if common=='reciprocal_exchange':
            method='reciprocal_exchange'
        if method not in ('facet_exchange', 'effective_slab', 'finite_crown_ports', 'reciprocal_exchange'):
            raise ValueError('Unknown optical terrain method: '+str(method))
        from . import canopy_closure
        if canopy_closure.closure(settings) != 'legacy' and method != 'effective_slab':
            from .optical_terrain_standard import reflectance
            return reflectance(settings,float(sza),float(saa),float(vza),float(vaa),
                               float(leaf_rho),float(leaf_tau),float(soil_rho),diffuse_fraction,sunlit)
        if method == 'reciprocal_exchange':
            from .optical_terrain_reciprocal import reflectance
            return reflectance(settings,float(sza),float(saa),float(vza),float(vaa),
                               float(leaf_rho),float(leaf_tau),float(soil_rho),diffuse_fraction,sunlit)
        if method == 'finite_crown_ports':
            from .optical_terrain_ports import reflectance
            return reflectance(settings,float(sza),float(saa),float(vza),float(vaa),
                               float(leaf_rho),float(leaf_tau),float(soil_rho),diffuse_fraction,sunlit)
        if method == 'facet_exchange':
            from .optical_terrain import terrain_crown_reflectance
            return terrain_crown_reflectance(settings, float(sza), float(saa), float(vza), float(vaa),
                                              float(leaf_rho), float(leaf_tau), float(soil_rho),
                                              diffuse_fraction, sunlit)
    # The generic path below computes illumination from the joint integral;
    # the thermal two-layer hotspot approximation is not used for VNIR.
    homogeneous = (model not in ('row','raw','crown','forest')
                   and settings.get('surface_type',0) not in (1,2,3,4,5,12,14)
                   and settings.get('distribution',0) != 2
                   and not (settings.get('surface_model') == 'terrain'
                            and int(settings.get('terrain_type', 2)) == 2))
    if bare_terrain:
        # Bare compound terrain: the terrain already supplies the sunlit/shaded soil fractions,
        # so the Lambertian soil reflects the beam only where the surface is both seen and lit.
        parts = np.asarray(_component_proportions(settings, sza, saa, vza, vaa, sunlit), float).reshape(4)
        if not sunlit:
            parts = np.asarray([0.0, parts[0] + parts[1], 0.0, parts[2] + parts[3]])
        visible_soil = float(parts[0] + parts[1])
        sunlit_soil = float(parts[0]) / visible_soil if visible_soil > 1.0e-12 else 0.0
        direct = (0.0 if diffuse_fraction >= 1.0 else (1.0 - diffuse_fraction) * soil_rho * sunlit_soil)
        return float(direct + diffuse_fraction * soil_rho)
    # Optical-only dispatch. Thermal observations call the separate legacy
    # function above; shared component and emissivity functions are unchanged.
    is_crown = (model in ('crown','forest') or settings.get('surface_type',0) in (1,2,3,4,5))
    is_row = (model in ('row','raw') or settings.get('distribution',0) == 2
              or settings.get('surface_type',0) in (12,14)) and not is_crown
    closure = str(settings.get('row_closure') or '').strip().lower()
    if is_crown and settings.get('surface_model','plane') == 'plane' and not use_hapke:
        method = settings.get('optical_crown_method','component_exchange')
        if method not in ('component_exchange','effective_slab'):
            raise ValueError('Unknown optical crown method: '+str(method))
        if method == 'component_exchange':
            from .optical_crown import crown_reflectance
            return crown_reflectance(settings,float(sza),float(vza),float(relative_azimuth(saa,vaa)),
                                     float(leaf_rho),float(leaf_tau),float(soil_rho),
                                     diffuse_fraction,sunlit)
    soil_direct = (hapke_reflectance(vza,settings.get('hapke_K',1),soil_rho)
                   if use_hapke else soil_rho)
    if is_row and closure and settings.get('surface_model','plane') == 'plane':
        from .row_closure import optical_reflectance
        return optical_reflectance(closure, closure, settings, sza, saa, vza, vaa, float(leaf_rho),
                                   float(leaf_tau), float(soil_rho), diffuse_fraction, sunlit, soil_direct)
    from . import canopy_closure
    if (canopy_closure.closure(settings) != 'legacy' and settings.get('surface_model','plane') == 'plane'
            and not is_crown and (homogeneous or is_row)):
        kind = 'row' if is_row else 'hom'
        shape = tuple(float(settings[key]) for key in ('row_width','row_blank','row_height')) if is_row else None
        return canopy_closure.reflectance(
            kind, float(settings.get('lai',0)), float(settings.get('hspot',0.1)), shape, float(leaf_rho),
            float(leaf_tau), float(soil_rho), float(vza), float(vaa), float(sza), float(saa),
            float(settings.get('row_azimuth',0.0)) if is_row else 0.0, None, diffuse_fraction, sunlit, float(soil_direct),
            canopy_closure.closure(settings))
    if homogeneous:
        view_gap = math.exp(-0.5*float(settings.get('lai',0))/max(math.cos(math.radians(vza)),1e-6))
        fractions = np.asarray([view_gap,0,1-view_gap,0])
        solar = fractions
    else:
        fractions = _component_proportions(settings,sza,saa,vza,vaa,sunlit)
        solar = _component_proportions(settings,sza,saa,sza,saa,True)
    value = canopy_reflectance(
        float(settings.get('lai',0)),float(settings.get('hspot',0.1)),
        float(leaf_rho),float(leaf_tau),float(soil_rho),float(sza),float(vza),
        float(relative_azimuth(saa,vaa)),fractions,float(solar[0]+solar[1]),
        homogeneous,diffuse_fraction,sunlit,soil_direct,
    )
    if is_row and not use_hapke and sunlit and diffuse_fraction < 1:
        method = settings.get('optical_row_method', 'conditioned_corridor')
        if method not in ('conditioned_corridor', 'effective_slab'):
            raise ValueError('Unknown optical row method: '+str(method))
        if method == 'conditioned_corridor':
            from .optical_row import primary_correction
            value += (1-diffuse_fraction)*primary_correction(
                settings, sza, saa, vza, vaa, leaf_rho, leaf_tau, soil_rho,
                float(solar[0]+solar[1]), fractions)
    return float(value)


def _requested_shortwave_diffuse_fraction(settings, global_shortwave, sza):
    """Partition horizontal global shortwave with Erbs when no fraction is supplied."""
    configured = float(settings.get("shortwave_diffuse_fraction", -1.0))
    if 0.0 <= configured <= 1.0:
        return configured
    cos_sza = max(math.cos(math.radians(float(sza))), 1.0e-6)
    # Same top-of-atmosphere reference (with the orbital factor) as clearness_index.
    extraterrestrial_horizontal = SOLAR_CONSTANT * _orbital_factor(settings) * cos_sza
    clearness = min(max(float(global_shortwave) / extraterrestrial_horizontal, 0.0), 1.2)
    if clearness <= 0.22:
        diffuse_fraction = 1.0 - 0.09 * clearness
    elif clearness <= 0.80:
        diffuse_fraction = (
            0.9511 - 0.1604 * clearness + 4.388 * clearness ** 2
            - 16.638 * clearness ** 3 + 12.336 * clearness ** 4
        )
    else:
        diffuse_fraction = 0.165
    return min(max(diffuse_fraction, 0.0), 1.0)


def _orbital_factor(settings):
    day = settings.get('solar_day_of_year')
    # Without a date use the perihelion value, the largest admissible factor.
    return 1.0+0.033*math.cos(2.0*math.pi*float(day)/365.0) if day is not None else 1.033


def _solar_direct_normal_limit(settings, sza):
    """Clear-sky direct normal irradiance [W m-2] used as the beam upper bound."""
    sza = float(sza)
    if sza >= 90.0:
        return 0.0
    # Meinel & Meinel (1976) DNI = 1367*E0*0.7**(AM**0.678); air mass from Kasten & Young (1989).
    air_mass = 1.0/(math.cos(math.radians(sza))+0.50572*(96.07995-sza)**-1.6364)
    return 1367.0*_orbital_factor(settings)*0.7**(air_mass**0.678)


def clearness_index(settings, global_shortwave, sza):
    """Global horizontal over top-of-atmosphere horizontal irradiance; nan with the sun down."""
    cos_sza = math.cos(math.radians(float(sza)))
    if cos_sza <= 0.0:
        return float('nan')
    return float(global_shortwave)/(SOLAR_CONSTANT*_orbital_factor(settings)*cos_sza)


def clearness_index_warning(states):
    """Summary of steps whose measured shortwave exceeds the top-of-atmosphere value."""
    exceeded = [(float(state['time']), float(state['clearness_index'])) for state in states
                if state.get('clearness_index_exceeds_one')]
    if not exceeded:
        return None
    listed = '; '.join(f'{time:.4f} (kt={value:.2f})' for time, value in exceeded)
    return (f'Clearness index kt>1 at {len(exceeded)}/{len(states)} time steps: {listed}. '
            'Meteorological values are used as instantaneous values at their timestamps, including '
            'period means; near sunrise/sunset a period mean can exceed the instantaneous top-of-atmosphere '
            'irradiance. Total shortwave is kept and direct power above the clear-sky limit is treated as diffuse.')


def _shortwave_diffuse_fraction(settings, global_shortwave, sza):
    """Admissible horizontal partition; excess direct power stays as diffuse.

    This preserves measured total shortwave. The limit can be disabled for
    deliberately artificial forcing with shortwave_enforce_solar_limit=0.
    """
    requested = _requested_shortwave_diffuse_fraction(settings, global_shortwave, sza)
    if not int(settings.get('shortwave_enforce_solar_limit',1)):
        return requested
    total = max(float(global_shortwave),0.0)
    if total <= 0.0:
        return requested
    horizontal_limit = _solar_direct_normal_limit(settings, sza)*max(math.cos(math.radians(sza)),0.0)
    return max(requested,1.0-min(horizontal_limit/total,1.0))


__all__ = ["relative_azimuth", "optical_properties", "directional_reflectance",
           "clearness_index", "clearness_index_warning"]
