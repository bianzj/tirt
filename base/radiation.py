import math
import numpy as np

from .gap import gap_probability_hom_analytical
from .hotspot import hotspot_analytical, hotspot_layer, hotspot_vegetation_volume
from .proportion import proportion_bidirectional_crown_one, proportion_bidirectional_hom_one, proportion_bidirectional_row_one
from .rt_transfer import (
    diffuse_scatter_tir,
    diffuse_scatter_vnir,
    directional_single_scatter,
    tirt_direct,
    tirt_scatter,
)
from .terrain_rt import compound_terrain_proportions
from .soilopt import hapke_reflectance


SIGMA = 5.67e-8


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


def _vegetation_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    lai = max(float(settings["lai"]), 0.0)
    if lai <= 0.0:
        return np.asarray([1.0, 0.0, 0.0, 0.0])
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
            np.asarray([settings["row_azimuth"]]),
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
        and int(settings.get("terrain_type", 1)) == 2
    ):
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
            ),
            sunlit,
        )
    return _vegetation_proportions(settings, sza, saa, vza, vaa, sunlit)


def component_proportions(settings, sza, saa, vza=0.0, vaa=0.0, sunlit=True):
    return _component_proportions(settings, sza, saa, vza, vaa, sunlit)


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
                diffuse_leaf, diffuse_soil, diffuse_soil_incident = diffuse_scatter_vnir(
                    settings["lai"], lrho, ltau, rs, sza, esun, esky, vza,
                    float(relative_azimuth(saa, vaa)), return_incident=True,
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
        direct_leaf = esun * labs if sunlit else 0.0
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


def prepare_longwave(settings, meteo, state, sza=0.0, saa=0.0, vza=0.0, vaa=0.0):
    result = {"directTrad_leaf": 0.0, "diffuseTrad_leaf": 0.0, "directTrad_soil": 0.0, "diffuseTrad_soil": 0.0}
    sky = float(meteo["longwave"])
    lai = float(settings["lai"])
    shortwave = max(float(meteo.get("shortwave", 0.0)), 0.0)
    leaf_emissivity = max(1.0 - float(settings["leaf_reflectance"]), 0.0)
    # The upper leaf side sees the sky, while the lower side sees the
    # vegetation/soil boundary.  Treating both sides as sky was the main
    # reason for unrealistically cold leaves at night.
    lower_surface = 0.5 * SIGMA * (state["Tss"] ** 4 + state["Tsh"] ** 4)
    sky_leaf = leaf_emissivity * (sky + lower_surface)
    sky_soil = sky * (1.0 - settings["soil_reflectance"])
    # The TIR transfer approximation uses the solar direction as its source
    # direction and becomes ill-conditioned at the horizon.  Fade it out
    # before that limit and use the sky boundary when the shortwave driver
    # says that the sun is effectively absent.
    tir_weight = max(0.0, min(1.0, (80.0 - float(sza)) / 10.0))
    if lai > 0.0 and shortwave >= 20.0 and sza < 89.0 and tir_weight > 0.0:
        lai_rt = lai
        if settings["surface_type"] in (1, 2, 3, 4, 5):
            lai_rt = lai
        # diffuse_scatter_tir returns (soil, leaf); keep the component order
        # explicit here so the hotter leaf emission cannot heat shaded soil.
        tir_soil, tir_leaf = diffuse_scatter_tir(lai_rt, settings["leaf_reflectance"], 0.0, settings["soil_reflectance"], sza, sky, state["Tss"], state["Tsh"], state["Tls"], state["Tlh"], vza, float(relative_azimuth(saa, vaa)))
        result["diffuseTrad_leaf"] = tir_weight * max(float(tir_leaf), 0.0) + (1.0 - tir_weight) * sky_leaf
        result["diffuseTrad_soil"] = tir_weight * max(float(tir_soil), 0.0) + (1.0 - tir_weight) * sky_soil
    elif lai > 0.0:
        # At night (or with the sun below the local horizon), the canopy still
        # exchanges longwave radiation with the sky.  Returning zero here made
        # every leaf lose two-sided emission without an atmospheric input.
        result["diffuseTrad_leaf"] = sky_leaf
        result["diffuseTrad_soil"] = sky_soil
    else:
        result["diffuseTrad_soil"] = sky_soil
    result["longwave_down"] = sky
    result["tir_scatter_soil"] = result["diffuseTrad_soil"]
    result["tir_scatter_leaf"] = result["diffuseTrad_leaf"]
    return result


def merge_radiation(shortwave, longwave):
    result = dict(shortwave)
    result.update(longwave)
    result["shortwave_direct"] = shortwave.get("directVrad_leaf", 0.0) + shortwave.get("directVrad_soil", 0.0)
    result["shortwave_diffuse"] = shortwave.get("diffuseVrad_leaf", 0.0) + shortwave.get("diffuseVrad_soil", 0.0)
    result["soil_sunlit"] = shortwave.get("directVrad_soil", 0.0) + shortwave.get("diffuseVrad_soil", 0.0) + longwave.get("diffuseTrad_soil", 0.0)
    result["soil_shaded"] = shortwave.get("diffuseVrad_soil", 0.0) + longwave.get("diffuseTrad_soil", 0.0)
    result["leaf_sunlit"] = shortwave.get("directVrad_leaf", 0.0) + shortwave.get("diffuseVrad_leaf", 0.0) + longwave.get("diffuseTrad_leaf", 0.0)
    result["leaf_shaded"] = shortwave.get("diffuseVrad_leaf", 0.0) + longwave.get("diffuseTrad_leaf", 0.0)
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


__all__ = ["relative_azimuth", "component_proportions", "optical_properties", "directional_reflectance", "prepare_shortwave", "prepare_longwave", "merge_radiation", "radiative_transfer_diagnostics"]
