import math
from functools import lru_cache
import numpy as np

from .gap import gap_probability_hom_analytical
from .hotspot import hotspot_analytical, hotspot_layer, hotspot_vegetation_volume
from .proportion import proportion_bidirectional_crown_one, proportion_bidirectional_hom_one, proportion_bidirectional_row_one
from .rt_transfer import (
    canopy_longwave,
    canopy_shortwave,
    directional_single_scatter,
    tirt_direct,
    tirt_scatter,
)
from .terrain_rt import compound_terrain_proportions
from .soilopt import hapke_reflectance


# Leaf projection of the spherical leaf angle distribution assumed by the TiRT canopy models.
LEAF_PROJECTION = 0.5


def _shade_components(values):
    components = np.asarray(values, dtype=float).reshape(4, -1)
    return np.vstack((
        np.zeros_like(components[0]),
        components[0] + components[1],
        np.zeros_like(components[0]),
        components[2] + components[3],
    ))


def relative_azimuth(saa, vaa):
    """Return the circular view-sun azimuth separation in [0, 180] degrees."""
    difference = np.mod(np.asarray(vaa, dtype=float) - np.asarray(saa, dtype=float), 360.0)
    return np.minimum(difference, 360.0 - difference)


_PROPORTION_KEYS = (
    "lai", "vegetation_model", "hspot",
    "crown_rad_a", "crown_rad_b", "crown_height", "crown_density",
    "row_width", "row_blank", "row_height", "row_azimuth",
)


def _proportion_key(settings):
    # The key must list every settings entry read by _vegetation_proportions_uncached.
    return tuple((name, settings[name]) for name in _PROPORTION_KEYS if name in settings)


def _vegetation_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True, key=None):
    key = _proportion_key(settings) if key is None else key
    return np.array(_vegetation_proportions_cached(key, float(sza), float(saa), float(vza), float(vaa), bool(sunlit)))


@lru_cache(maxsize=65536)
def _vegetation_proportions_cached(key, sza, saa, vza, vaa, sunlit):
    return tuple(_vegetation_proportions_uncached(dict(key), sza, saa, vza, vaa, sunlit))


def _vegetation_proportions_uncached(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    lai = max(float(settings["lai"]), 0.0)
    if lai <= 0.0:
        return np.asarray([1.0, 0.0, 0.0, 0.0])
    sza_ = np.asarray([min(max(float(sza), 0.0), 89.9)])
    vza_ = np.asarray([min(max(float(vza), 0.0), 89.9)])
    saa_ = np.asarray([float(saa)])
    vaa_ = np.asarray([float(vaa)])
    vsa_ = relative_azimuth(saa_, vaa_)
    vegetation_model = str(settings.get("vegetation_model", "")).strip().lower()
    if vegetation_model == "crown":
        values = proportion_bidirectional_crown_one(
            lai, settings["crown_density"], settings["hspot"], settings["crown_rad_b"], settings["crown_rad_a"],
            vza_, sza_, vsa_, hcb=settings["crown_height"],
        )
    elif vegetation_model == "row":
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
            np.asarray([settings["row_azimuth"]]),
        )
    else:
        values = proportion_bidirectional_hom_one(lai, settings["hspot"], vza_, sza_, vsa_)
    values = np.maximum(np.asarray(values, dtype=float).reshape(4), 0.0)
    if not sunlit:
        values = np.asarray([0.0, values[0] + values[1], 0.0, values[2] + values[3]])
    return values / values.sum() if values.sum() else np.asarray([1.0, 0.0, 0.0, 0.0])


def _component_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    if settings.get("surface_model") == "terrain":
        key = _proportion_key(settings)
        return compound_terrain_proportions(
            settings,
            sza,
            saa,
            vza,
            vaa,
            lambda local_sza, local_saa, local_vza, local_vaa, local_sunlit: _vegetation_proportions(
                settings,
                local_sza,
                local_saa,
                local_vza,
                local_vaa,
                local_sunlit,
                key,
            ),
            sunlit,
        )
    return _vegetation_proportions(settings, sza, saa, vza, vaa, sunlit)


def component_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    return _component_proportions(settings, sza, saa, vza, vaa, sunlit)


def urban_canopy_fraction(settings, vza=0.0):
    """Street fraction hidden by the street crowns along the view direction (TiRT crown gap probability)."""
    if float(settings.get("lai", 0.0)) <= 0.0:
        return 0.0
    vza = min(max(float(vza), 0.0), 89.0)
    fractions = _vegetation_proportions(settings, vza, 0.0, vza, 0.0, True)
    return float(fractions[2] + fractions[3])


def sun_extinction(sza):
    """Direct-beam extinction per unit LAI, G / cos(sza); the zenith is capped at 85 degrees as in SCOPE."""
    return LEAF_PROJECTION / math.cos(math.radians(min(float(sza), 85.0)))


def flux_weights(settings, sza, saa, sunlit=True):
    """Per unit ground: sunlit and shaded soil fractions and sunlit and shaded leaf area indices.

    Component fluxes are per unit soil or leaf area; these weights give the scene fluxes per unit
    ground, as in SCOPE ebal. The sunlit leaf area follows from the beam the canopy intercepts:
    (1 - gap seen from the sun) = LAI_sunlit * G / cos(sza).
    """
    lai = max(float(settings.get("lai", 0.0)), 0.0)
    if not sunlit or float(sza) >= 90.0:
        return 0.0, 1.0, 0.0, lai
    if lai <= 0.0:
        return 1.0, 0.0, 0.0, 0.0
    fractions = _vegetation_proportions(settings, min(float(sza), 89.0), saa, min(float(sza), 89.0), saa, True)
    gap = float(fractions[0] + fractions[1])
    lai_sunlit = min((1.0 - gap) / sun_extinction(sza), lai)
    return gap, 1.0 - gap, lai_sunlit, lai - lai_sunlit


def optical_properties(settings):
    rho = min(max(float(settings["leaf_reflectance"]), 0.0), 0.99)
    tau = min(max(float(settings["leaf_transmittance"]), 0.0), 0.99 - rho)
    soil = min(max(float(settings["soil_reflectance"]), 0.0), 0.99)
    return rho, tau, soil, 1.0 - rho - tau, 1.0 - soil


def directional_reflectance(
    settings, sza, saa, vza, vaa, leaf_rho, leaf_tau, soil_rho,
    sunlit=True, use_hapke=None,
):
    """Calculate a view-direction reflectance with the same RT angular terms."""
    proportions = _component_proportions(settings, sza, saa, vza, vaa, sunlit)
    soil_fraction = float(proportions[0] + proportions[1])
    leaf_fraction = float(proportions[2] + proportions[3])
    if use_hapke is None:
        use_hapke = int(settings.get("hapke_enable", settings.get("hapke_optical_enable", 0))) == 1
    if use_hapke:
        soil_direct = hapke_reflectance(
            vza,
            settings.get("hapke_K", 1.0),
            soil_rho,
        )
    else:
        soil_direct = float(soil_rho)
    base = soil_fraction * soil_direct + leaf_fraction * float(leaf_rho)
    if leaf_fraction <= 1.0e-12 or float(settings.get("lai", 0.0)) <= 0.0:
        return float(np.clip(base, 0.0, 1.0))
    raa = float(relative_azimuth(saa, vaa))
    directional = directional_single_scatter(
        settings["lai"], leaf_rho, leaf_tau, sza, vza, raa,
    )
    reference = directional_single_scatter(
        settings["lai"], leaf_rho, leaf_tau, sza, vza, 90.0,
    )
    return float(np.clip(base + directional - reference, 0.0, 1.0))


def _empty_shortwave():
    return {
        "directVrad_incident": 0.0, "directVrad_incident_spectrum": [],
        "diffuseVrad_incident": 0.0, "diffuseVrad_incident_spectrum": [],
        "directVrad_leaf": 0.0, "diffuseVrad_leaf": 0.0,
        "directVrad_soil": 0.0, "diffuseVrad_soil": 0.0,
        "diffuseVrad_soil_spectrum": [], "diffuseVrad_soil_incident_spectrum": [],
        "directPrad_leaf": 0.0, "diffusePrad_leaf": 0.0,
        "shortwave_global": 0.0, "shortwave_direct_incident": 0.0,
        "shortwave_diffuse_incident": 0.0, "shortwave_diffuse_fraction": 0.0,
        "fractions": np.zeros(4),
    }


def _shortwave_diffuse_fraction(settings, global_shortwave, sza):
    """Partition horizontal global shortwave with Erbs when no fraction is supplied."""
    configured = float(settings.get("shortwave_diffuse_fraction", -1.0))
    if 0.0 <= configured <= 1.0:
        return configured
    cos_sza = max(math.cos(math.radians(float(sza))), 1.0e-6)
    extraterrestrial_horizontal = 1361.0 * cos_sza
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


def prepare_shortwave(settings, meteo, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    result = _empty_shortwave()
    rin = float(meteo["shortwave"])
    if rin < 20.0 or sza >= 89.9:
        return result
    diffuse_fraction = _shortwave_diffuse_fraction(settings, rin, sza)
    direct_horizontal = rin * (1.0 - diffuse_fraction)
    diffuse_horizontal = rin * diffuse_fraction
    result["shortwave_global"] = rin
    result["shortwave_direct_incident"] = direct_horizontal
    result["shortwave_diffuse_incident"] = diffuse_horizontal
    result["shortwave_diffuse_fraction"] = diffuse_fraction
    rho, tau, soil, leaf_abs, soil_abs = optical_properties(settings)
    leaf_rho_spectrum = settings.get("leaf_reflectance_spectrum", [])
    leaf_tau_spectrum = settings.get("leaf_transmittance_spectrum", [])
    soil_spectrum = settings.get("soil_reflectance_spectrum", [])
    street_spectrum = settings.get("urban_reflectance_street_spectrum", [])
    urban_street = str(settings.get("surface_model", "")).strip().lower() == "urban"
    wavelengths = settings.get("wavelengths", [])[:2001]
    fesun = settings.get("fesun", [])[:2001]
    fesky = settings.get("fesky", [])[:2001]
    if len(fesun) != 2001 or len(fesky) != 2001:
        fesun = np.full(2001, 1.0 / 2001.0)
        fesky = np.full(2001, 1.0 / 2001.0)
        wavelengths = np.arange(400.0, 2401.0)
    # A sunlit leaf intercepts the beam with its projection G per unit leaf area: horizontal irradiance * G / cos(sza).
    leaf_projection = sun_extinction(sza)
    for index in range(2001):
        lrho = leaf_rho_spectrum[index] if len(leaf_rho_spectrum) == 2001 else rho
        ltau = leaf_tau_spectrum[index] if len(leaf_tau_spectrum) == 2001 else tau
        rs = (street_spectrum[index] if urban_street and len(street_spectrum) == 2001 else soil_spectrum[index] if len(soil_spectrum) == 2001 else soil)
        labs = 1.0 - lrho - ltau
        sabs = 1.0 - rs
        esun = direct_horizontal * float(fesun[index]) * 0.001
        esky = diffuse_horizontal * float(fesky[index]) * 0.001
        try:
            if settings["lai"] > 0.0:
                diffuse_leaf, diffuse_soil, diffuse_soil_incident = canopy_shortwave(
                    settings["lai"], lrho, ltau, rs, sza, esun if sunlit else 0.0, esky,
                )
            else:
                diffuse_leaf, diffuse_soil, diffuse_soil_incident = 0.0, esky * sabs, esky
        except (FloatingPointError, IndexError, ValueError, ZeroDivisionError, OverflowError):
            diffuse_leaf, diffuse_soil, diffuse_soil_incident = 0.0, 0.0, 0.0
        # The meteorological shortwave driver is the incoming horizontal
        # irradiance.  Converting it back to normal irradiance with 1/cos(SZA)
        # creates an artificial sunrise/sunset singularity near 90 degrees.
        result["directVrad_incident"] += esun if sunlit else 0.0
        result["directVrad_incident_spectrum"].append(esun if sunlit else 0.0)
        result["diffuseVrad_incident"] += esky
        result["diffuseVrad_incident_spectrum"].append(esky)
        direct_leaf = esun * leaf_projection * labs if sunlit else 0.0
        direct_soil = esun * sabs if sunlit else 0.0
        result["directVrad_leaf"] += direct_leaf
        result["directVrad_soil"] += direct_soil
        result["diffuseVrad_leaf"] += diffuse_leaf
        result["diffuseVrad_soil"] += diffuse_soil
        result["diffuseVrad_soil_spectrum"].append(diffuse_soil)
        result["diffuseVrad_soil_incident_spectrum"].append(diffuse_soil_incident)
        wl = wavelengths[index]
        if 400.0 <= wl <= 700.0:
            result["directPrad_leaf"] += direct_leaf * wl * 1e-3 / (6.02214e23 * 6.6262e-34 * 299792458.0)
            result["diffusePrad_leaf"] += diffuse_leaf * wl * 1e-3 / (6.02214e23 * 6.6262e-34 * 299792458.0)
    result["fractions"] = _component_proportions(settings, sza, saa, vza, vaa, sunlit)
    result.update({"rho_leaf": rho, "tau_leaf": tau, "rho_soil": soil})
    return result


def prepare_longwave(settings, meteo, state, sza=0.0, saa=0.0, sunlit=True):
    """Absorbed longwave per unit area: leaves (both faces; sunlit = diffuse + direct entry) and soil."""
    sky = float(meteo["longwave"])
    lai = float(settings["lai"])
    soil_emissivity = max(float(settings["soil_emissivity"]), 0.0)
    if lai > 0.0:
        soil_sunlit_fraction = flux_weights(settings, sza, saa, sunlit)[0]
        extinction = sun_extinction(sza) if sunlit and float(sza) < 90.0 else 0.0
        leaf_sunlit, leaf_shaded, soil = canopy_longwave(
            lai, max(float(settings["leaf_emissivity"]), 0.0), soil_emissivity, sky,
            state["Tls"], state["Tlh"], state["Tss"], state["Tsh"], extinction, soil_sunlit_fraction,
        )
    else:
        leaf_sunlit = leaf_shaded = 0.0
        soil = soil_emissivity * sky
    result = {
        "directTrad_leaf": leaf_sunlit - leaf_shaded, "diffuseTrad_leaf": leaf_shaded,
        "directTrad_soil": 0.0, "diffuseTrad_soil": soil,
    }
    result["longwave_down"] = sky
    result["tir_scatter_soil"] = soil
    result["tir_scatter_leaf"] = leaf_shaded
    return result


def merge_radiation(shortwave, longwave):
    result = dict(shortwave)
    result.update(longwave)
    result["shortwave_direct"] = shortwave.get("directVrad_leaf", 0.0) + shortwave.get("directVrad_soil", 0.0)
    result["shortwave_diffuse"] = shortwave.get("diffuseVrad_leaf", 0.0) + shortwave.get("diffuseVrad_soil", 0.0)
    for component in ("soil", "leaf"):
        diffuse = shortwave.get(f"diffuseVrad_{component}", 0.0) + longwave.get(f"diffuseTrad_{component}", 0.0)
        direct = shortwave.get(f"directVrad_{component}", 0.0) + longwave.get(f"directTrad_{component}", 0.0)
        result[f"{component}_sunlit"] = direct + diffuse
        result[f"{component}_shaded"] = diffuse
    return result


def radiative_transfer_diagnostics(settings, state, vza=0.0, vaa=0.0):
    lai = max(float(settings["lai"]), 0.0)
    if lai <= 0.0:
        return {"gap_view": 1.0, "gap_sun": 1.0, "hotspot_leaf": 0.0, "hotspot_soil": 1.0, "tirt_direct_soil_sun": 1.0, "tirt_direct_soil_shade": 0.0, "tirt_direct_leaf_sun": 0.0, "tirt_direct_leaf_shade": 0.0, "tirt_scatter_soil_sun": 0.0, "tirt_scatter_soil_shade": 0.0, "tirt_scatter_leaf_sun": 0.0, "tirt_scatter_leaf_shade": 0.0}
    sza, saa = state["sza"], state["saa"]
    if sza >= 89.0:
        return {"gap_view": float(gap_probability_hom_analytical(lai, vza)), "gap_sun": 0.0, "hotspot_leaf": 0.0, "hotspot_soil": 0.0, "sunlit_volume": 0.0, "tirt_direct_soil_sun": 0.0, "tirt_direct_soil_shade": 0.0, "tirt_direct_leaf_sun": 0.0, "tirt_direct_leaf_shade": 0.0, "tirt_scatter_soil_sun": 0.0, "tirt_scatter_soil_shade": 0.0, "tirt_scatter_leaf_sun": 0.0, "tirt_scatter_leaf_shade": 0.0}
    direct = tirt_direct(lai, settings["hspot"], sza, saa, vza, vaa)
    scatter = tirt_scatter(lai, settings["soil_reflectance"], settings["leaf_reflectance"], sza)
    if not state.get("sunlit_surface", True):
        direct = _shade_components(direct)
        scatter = _shade_components(scatter)
    return {
        "gap_view": float(gap_probability_hom_analytical(lai, vza)),
        "gap_sun": float(gap_probability_hom_analytical(lai, sza)),
        "hotspot_leaf": float(np.asarray(hotspot_analytical(lai, settings["hspot"], np.asarray([vza]), np.asarray([sza]), relative_azimuth(saa, vaa)).reshape(-1)[0])),
        "hotspot_soil": float(np.asarray(hotspot_layer(lai, settings["hspot"], np.asarray([vza]), np.asarray([sza]), relative_azimuth(saa, vaa)).reshape(-1)[0])),
        "sunlit_volume": float(np.asarray(hotspot_vegetation_volume(lai, np.asarray([sza]))).reshape(-1)[0]),
        "tirt_direct_soil_sun": direct[0], "tirt_direct_soil_shade": direct[1],
        "tirt_direct_leaf_sun": direct[2], "tirt_direct_leaf_shade": direct[3],
        "tirt_scatter_soil_sun": scatter[0], "tirt_scatter_soil_shade": scatter[1],
        "tirt_scatter_leaf_sun": scatter[2], "tirt_scatter_leaf_shade": scatter[3],
    }


__all__ = ["relative_azimuth", "component_proportions", "urban_canopy_fraction", "sun_extinction", "flux_weights", "optical_properties", "directional_reflectance", "prepare_shortwave", "prepare_longwave", "merge_radiation", "radiative_transfer_diagnostics"]
