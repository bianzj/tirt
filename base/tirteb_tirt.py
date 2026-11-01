"""TIRTEB thermal observations computed with the TiRT directional model (run._surface_result).

Maps the energy-balance settings and component temperatures onto the TiRT input contract, so that
TIRTEB and TiRT give the same directional radiance for the same scene and temperatures.
"""
import numpy as np

from .brightness import _urban_spectral_value, _vegetation_spectral_value

_VEGETATION = {"bare": "bare", "bare_soil": "bare", "soil": "bare", "barren": "bare",
               "hom": "hom", "homogeneous": "hom", "row": "row", "raw": "row",
               "crown": "crown", "forest": "crown"}


def _surface(settings):
    surface = str(settings.get("surface_model", "plane")).strip().lower()
    if surface == "urban" or int(settings.get("surface_type", -1)) == 13:
        return "urban"
    if surface in ("terrain", "slope"):
        return "slope" if int(settings.get("terrain_type", 2)) == 1 or surface == "slope" else "terrain"
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
        # Input crown_rad_b -> height (vertical semiaxis), crown_rad_a -> width (horizontal radius),
        # crown_density -> stand_density (crowns per m2), crown_height -> crown_base_height (crown base above ground)
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
        "canopy_closure": "standard",
    }
    sky = settings.get("sky_temperature_k")
    if str(sky).strip().lower() == "auto":
        # Equivalent sky temperature from the measured broadband downwelling longwave
        # of this time step. Exact for an overcast (grey) sky; for a clear sky it
        # overestimates the 8-14 um window sky radiance.
        down = float(state.get("longwave_down", float("nan")))
        sky = (down / 5.67e-8) ** 0.25 if np.isfinite(down) and down > 0.0 else None
    if sky not in (None, ""):
        config["sky_temperature_k"] = float(sky)

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


def shortwave_optics(settings, wavelength):
    """Shortwave properties from the PROSPECT/BSM (or measured) spectra and the urban spectra."""
    w = float(wavelength)
    return {
        "leaf_reflectance": _vegetation_spectral_value(settings, "leaf_reflectance_spectrum", "leaf_reflectance", w),
        "leaf_transmittance": _vegetation_spectral_value(settings, "leaf_transmittance_spectrum", "leaf_transmittance", w),
        "soil_reflectance": _vegetation_spectral_value(settings, "soil_reflectance_spectrum", "soil_reflectance", w),
        "street_reflectance": _urban_spectral_value(settings, "urban_reflectance_street_spectrum", w, "urban_reflectance_street"),
        "wall_reflectance": _urban_spectral_value(settings, "urban_reflectance_wall_spectrum", w, "urban_reflectance_wall"),
        "roof_reflectance": _urban_spectral_value(settings, "urban_reflectance_roof_spectrum", w, "urban_reflectance_roof"),
    }


def band_optics(settings, wavelength):
    """Optical properties of one band: the resolved band table, else the model spectra or energy-balance means."""
    w = float(wavelength)
    for band_wavelength, band in zip(settings.get("observation_wavelengths", ()), settings.get("observation_bands", ())):
        if abs(float(band_wavelength) - w) < 1.0e-9:
            return band
    if w < 2.5:
        return shortwave_optics(settings, w)
    return {
        "leaf_reflectance": float(settings["leaf_reflectance"]),
        "leaf_transmittance": float(settings["leaf_transmittance"]),
        "soil_reflectance": float(settings["soil_reflectance"]),
        "street_reflectance": 1.0 - float(settings["urban_emissivity_street"]),
        "wall_reflectance": 1.0 - float(settings["urban_emissivity_wall"]),
        "roof_reflectance": 1.0 - float(settings["urban_emissivity_roof"]),
    }


def tir_band_settings(settings, wavelength):
    """Settings carrying one thermal band's Kirchhoff emissivities (leaf 1 - rho - tau, others 1 - rho)."""
    optics = band_optics(settings, wavelength)
    band = dict(settings)
    leaf_rho, leaf_tau, soil_rho = optics["leaf_reflectance"], optics["leaf_transmittance"], optics["soil_reflectance"]
    band.update(leaf_reflectance=leaf_rho, leaf_transmittance=leaf_tau,
                leaf_emissivity=max(1.0 - leaf_rho - leaf_tau, 1.0e-6),
                soil_reflectance=soil_rho, soil_emissivity=max(1.0 - soil_rho, 1.0e-6),
                urban_emissivity_roof=1.0 - optics["roof_reflectance"],
                urban_emissivity_wall=1.0 - optics["wall_reflectance"],
                urban_emissivity_street=1.0 - optics["street_reflectance"])
    return band


def band_emissivities(settings):
    """TiRT spectra dict for one thermal band from tir_band_settings()."""
    soil = float(settings["soil_emissivity"])
    return {
        "leaf_e": float(settings["leaf_emissivity"]), "soil_e": soil, "terrain_e": soil,
        "roof_e": float(settings["urban_emissivity_roof"]),
        "wall_e": float(settings["urban_emissivity_wall"]),
        "street_e": float(settings["urban_emissivity_street"]),
    }


def tirt_brightness(settings, state, wavelength, vza, vaa):
    """Directional radiance and brightness temperature from TiRT; NaN for slope back faces.

    Urban vegetation is part of the TiRT urban model itself (homogeneous, row or
    crown); no extra street coverage fraction is applied.
    """
    import run
    surface, vegetation, config = tirt_config(settings, state)
    # Below the horizon every component is shaded, so the sun only has to stay at a finite zenith.
    sun = (min(float(state.get("sza_global", state["sza"])), 89.0), float(state.get("saa_global", state["saa"])))
    geometry = (float(vza), float(vaa), *sun)
    if surface == "slope" and run._slope_backface(config, geometry):
        return float("nan"), float("nan")
    spectra = band_emissivities(settings)
    sky = config.get("sky_temperature_k")

    def patch(kind):
        radiance, temperature = run._surface_result(config, surface, kind, float(wavelength), geometry, spectra)
        if sky not in (None, "") and np.isfinite(radiance):
            radiance, temperature = run._with_sky_reflection(config, surface, kind, float(wavelength), geometry,
                                                             spectra, radiance, float(sky))
        return float(radiance), float(temperature)

    return patch(vegetation)


_TEMPERATURE_KEYS = ("Tss", "Tsh", "Tls", "Tlh", "Troof_sunlit", "Troof_shaded", "Twall_sunlit", "Twall_shaded",
                     "Tstreet_sunlit", "Tstreet_shaded")


def tirt_emissivity(settings, state, wavelength, vza, vaa):
    """Directional emissivity (black sky, all scattering orders) of the TiRT model; NaN for slope back faces."""
    import run
    from .physicsF import planck
    temperature = 300.0
    surface, vegetation, config = tirt_config(settings, dict(state, **{key: temperature for key in _TEMPERATURE_KEYS}))
    config.pop("sky_temperature_k", None)
    sun = (min(float(state.get("sza_global", state["sza"])), 89.0), float(state.get("saa_global", state["saa"])))
    geometry = (float(vza), float(vaa), *sun)
    if surface == "slope" and run._slope_backface(config, geometry):
        return float("nan")
    radiance, _ = run._surface_result(config, surface, vegetation, float(wavelength), geometry, band_emissivities(settings))
    return float(radiance) / float(planck(float(wavelength), temperature))
