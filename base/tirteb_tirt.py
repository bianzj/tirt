"""TIRTEB thermal observations computed with the TiRT directional model (run._surface_result).

Maps the energy-balance settings and component temperatures onto the TiRT input contract, so that
TIRTEB and TiRT give the same directional radiance for the same scene and temperatures.
"""
import numpy as np

from .brightness import _urban_spectral_value

_VEGETATION = {"bare": "bare", "bare_soil": "bare", "soil": "bare", "barren": "bare",
               "hom": "hom", "homogeneous": "hom", "row": "row", "raw": "row",
               "crown": "crown", "forest": "crown"}


def _surface(settings):
    surface = str(settings.get("surface_model", "plane")).strip().lower()
    if surface == "urban" or int(settings.get("surface_type", -1)) == 13:
        return "urban"
    if surface in ("terrain", "slope"):
        return "slope" if int(settings.get("terrain_type", 1)) == 1 or surface == "slope" else "terrain"
    return "plane"


def tirt_config(settings, state):
    """TiRT surface, vegetation and config for the current TIRTEB state."""
    surface = _surface(settings)
    vegetation = _VEGETATION.get(str(settings.get("vegetation_model", "hom")).strip().lower(), "hom")
    if float(settings.get("lai", 0.0)) <= 0.0:
        vegetation = "bare"
    height, width = float(settings.get("height", 3.0)), float(settings.get("width", 3.0))
    config = {
        "lai": settings.get("lai", 0.0), "hspot": settings.get("hspot", 0.2),
        "row_width": settings.get("row_width", 0.5), "row_blank": settings.get("row_blank", 0.5),
        "row_height": settings.get("row_height", 1.0), "row_azimuth": settings.get("row_azimuth", 0.0),
        # TIRTEB crowns: height = vertical semiaxis, width = horizontal radius, stand_density = crowns per m2
        "crown_rad_a": width, "crown_rad_b": height,
        "crown_density": settings.get("stand_density", 0.02),
        "crown_height": settings.get("crown_base_height", height),
        "slope_angle": settings.get("terrain_slope", 0.0), "slope_aspect": settings.get("terrain_aspect", 0.0),
        "terrain_height": settings.get("terrain_height", 20.0), "terrain_radius": settings.get("terrain_radius", 10.0),
        "terrain_density": settings.get("terrain_density", 0.001),
        "building_density": settings.get("urban_density", 0.0015),
        "building_length": settings.get("urban_building_length", 10.0),
        "building_width": settings.get("urban_building_width", 10.0),
        "building_height": settings.get("urban_building_height", 20.0),
        "building_azimuth": settings.get("urban_building_azimuth", 0.0),
        "hapke_enable": settings.get("hapke_enable", 0), "hapke_k": settings.get("hapke_K", 1.0),
        "soil_reflectance": settings.get("soil_reflectance"),
    }
    for key in ("thermal_crown_geometry", "thermal_row_geometry", "forest_transport",
                "thermal_hotspot_inner", "thermal_hotspot_outer", "sky_temperature_k"):
        if settings.get(key) not in (None, ""):
            config[key] = settings[key]

    sunlit = bool(state.get("sunlit_surface", True)) and float(state.get("sza_global", state["sza"])) < 90.0

    def pair(sunlit_key, shaded_key, fallback_sunlit, fallback_shaded):
        shaded = float(state.get(shaded_key, state[fallback_shaded]))
        return (float(state.get(sunlit_key, state[fallback_sunlit])) if sunlit else shaded), shaded

    temperatures = {"soil": pair("Tss", "Tsh", "Tss", "Tsh"), "leaf": pair("Tls", "Tlh", "Tls", "Tlh")}
    temperatures["terrain"] = temperatures["soil"]
    temperatures["roof"] = pair("Troof_sunlit", "Troof_shaded", "Tls", "Tlh")
    temperatures["wall"] = pair("Twall_sunlit", "Twall_shaded", "Tls", "Tlh")
    temperatures["street"] = pair("Tstreet_sunlit", "Tstreet_shaded", "Tss", "Tsh")
    for name, (sun_k, shade_k) in temperatures.items():
        config[f"{name}_sunlit_k"] = str(sun_k)
        config[f"{name}_shaded_k"] = str(shade_k)
    return surface, vegetation, {key: item for key, item in config.items() if item is not None}


def tir_spectral_value(settings, key, scalar_key, wavelength):
    """Thermal-band value of a TIRTEB spectrum, falling back to the scalar setting."""
    spectrum = settings.get(key, [])
    wavelengths = settings.get("tir_wavelengths", [])
    if len(spectrum) == len(wavelengths) and wavelengths:
        index = int(np.argmin(np.abs(np.asarray(wavelengths, dtype=float) - float(wavelength))))
        return float(spectrum[index])
    return float(settings.get(scalar_key, 0.0))


def tir_band_settings(settings, wavelength):
    """Settings with leaf and soil reflectance/emissivity of one thermal band (emissivity = 1 - rho - tau)."""
    band = dict(settings)
    leaf_rho = tir_spectral_value(settings, "leaf_reflectance_tir_spectrum", "leaf_reflectance", wavelength)
    leaf_tau = tir_spectral_value(settings, "leaf_transmittance_tir_spectrum", "leaf_transmittance", wavelength)
    soil_rho = tir_spectral_value(settings, "soil_reflectance_tir_spectrum", "soil_reflectance", wavelength)
    band.update(leaf_reflectance=leaf_rho, leaf_transmittance=leaf_tau,
                leaf_emissivity=max(1.0 - leaf_rho - leaf_tau, 1.0e-6),
                soil_reflectance=soil_rho, soil_emissivity=max(1.0 - soil_rho, 1.0e-6))
    return band


def band_emissivities(settings, wavelength):
    """TiRT spectra dict for one thermal band from the TIRTEB band settings."""
    soil = float(settings["soil_emissivity"])
    return {
        "leaf_e": float(settings["leaf_emissivity"]), "soil_e": soil, "terrain_e": soil,
        "roof_e": _urban_spectral_value(settings, "urban_emissivity_roof_spectrum", wavelength, "urban_emissivity_roof"),
        "wall_e": _urban_spectral_value(settings, "urban_emissivity_wall_spectrum", wavelength, "urban_emissivity_wall"),
        "street_e": _urban_spectral_value(settings, "urban_emissivity_street_spectrum", wavelength, "urban_emissivity_street"),
    }


def tirt_brightness(settings, state, wavelength, vza, vaa):
    """Directional radiance and brightness temperature from TiRT; NaN for slope back faces.

    Urban vegetation covers the fraction `urban_vegetation_fraction` of the street (as in the TIRTEB
    energy balance and optical urban path); the rest is bare street, mixed linearly in radiance.
    """
    import run
    from .physicsF import inv_planck
    surface, vegetation, config = tirt_config(settings, state)
    # Below the horizon every component is shaded, so the sun only has to stay at a finite zenith.
    sun = (min(float(state.get("sza_global", state["sza"])), 89.0), float(state.get("saa_global", state["saa"])))
    geometry = (float(vza), float(vaa), *sun)
    if surface == "slope" and run._slope_backface(config, geometry):
        return float("nan"), float("nan")
    spectra = band_emissivities(settings, wavelength)
    sky = config.get("sky_temperature_k")

    def patch(kind):
        radiance, temperature = run._surface_result(config, surface, kind, float(wavelength), geometry, spectra)
        if sky not in (None, "") and np.isfinite(radiance):
            radiance, temperature = run._with_sky_reflection(config, surface, kind, float(wavelength), geometry,
                                                             spectra, radiance, float(sky))
        return float(radiance), float(temperature)

    if surface == "urban" and vegetation != "bare":
        cover = min(max(float(settings.get("urban_vegetation_fraction", 1.0)), 0.0), 1.0)
        if cover < 1.0:
            radiance = cover * patch(vegetation)[0] + (1.0 - cover) * patch("bare")[0]
            return radiance, float(np.asarray(inv_planck(float(wavelength), radiance)).reshape(-1)[0])
    return patch(vegetation)
