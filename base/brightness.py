import math

import numpy as np

from .emissivity import emissivity_direct, emissivity_scatter_analytical
from .physicsF import inv_planck, planck
from .radiance import radiance_direct, radiance_scatter
from .radiation import component_proportions, directional_reflectance, legacy_directional_reflectance, relative_azimuth
from .soilopt import hapke_thermal_emissivity
from .urban import Urban
from .thermal_crown import geometry_mode, is_forest, thermal_emission


def _urban_spectral_value(settings, key, wavelength, default):
    spectrum = settings.get(key, [])
    wavelengths = settings.get("wavelengths", [])
    # The settings grid is in nm; the public observation wavelength is in um.
    index = int(np.argmin(np.abs(np.asarray(wavelengths, dtype=float) - float(wavelength) * 1000.0))) if wavelengths else -1
    if len(spectrum) == len(wavelengths) and 0 <= index < len(spectrum):
        return float(spectrum[index])
    return float(settings.get(default, default if isinstance(default, (int, float)) else 0.9))


def _vegetation_spectral_value(settings, key, default, wavelength):
    spectrum = settings.get(key, [])
    wavelengths = settings.get("wavelengths", [])
    # PROSPECT/BSM optical spectra contain only the first 2001 grid samples.
    if len(spectrum) > 0 and (len(spectrum) == len(wavelengths) or (len(spectrum) == 2001 and len(wavelengths) >= 2001)):
        grid_nm = np.asarray(wavelengths[:len(spectrum)], dtype=float)
        index = int(np.argmin(np.abs(grid_nm - float(wavelength) * 1000.0)))
        return float(spectrum[index])
    return float(settings.get(default, 0.0))


def _urban_shapes(settings):
    """Build the urban geometry expected by the TIRT urban projection model."""
    length = max(float(settings.get("urban_building_length", 10.0)), 0.0)
    width = max(float(settings.get("urban_building_width", 10.0)), 0.0)
    height = max(float(settings.get("urban_building_height", 20.0)), 0.0)
    density = max(float(settings.get("urban_density", 0.0015)), 0.0)
    return np.asarray([[length, width, height, density, 0.0, 90.0]], dtype=float)


def urban_component_fractions(settings, vza, vaa, sza, saa):
    """Return urban wall/street/roof sunlit-shaded view fractions."""
    urban = Urban()
    urban.set_strcutural_input(_urban_shapes(settings))
    urban.set_angular_input(
        np.asarray([vza], dtype=float),
        np.asarray([vaa], dtype=float),
        float(sza),
        float(saa),
    )
    values = np.asarray(urban.calculate_direct_emissivity(0, 3), dtype=float).reshape(-1)
    values = np.maximum(values, 0.0)
    total = float(values.sum())
    return values / total if total > 0.0 else np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


def urban_canopy_fraction(settings, vza=0.0):
    """Fraction of the visible street covered by the selected canopy."""
    lai = max(float(settings.get("lai", 0.0)), 0.0)
    footprint = min(max(float(settings.get("urban_vegetation_fraction", 0.30)), 0.0), 1.0)
    if lai <= 0.0 or footprint <= 0.0:
        return 0.0
    mu = max(math.cos(math.radians(min(max(float(vza), 0.0), 89.0))), 0.1)
    return float(np.clip(footprint * (1.0 - math.exp(-0.5 * lai / mu)), 0.0, 0.95))


def urban_directional_reflectance(settings, wavelength, vza, vaa, sza, saa):
    """Optical urban shadows, facet incidence and wall/street/crown exchange."""
    rho_roof = _urban_spectral_value(settings, "urban_reflectance_roof_spectrum", wavelength, "urban_reflectance_roof")
    rho_wall = _urban_spectral_value(settings, "urban_reflectance_wall_spectrum", wavelength, "urban_reflectance_wall")
    rho_street = _urban_spectral_value(settings, "urban_reflectance_street_spectrum", wavelength, "urban_reflectance_street")
    leaf_rho = _vegetation_spectral_value(settings, "leaf_reflectance_spectrum", "leaf_reflectance", wavelength)
    leaf_tau = _vegetation_spectral_value(settings, "leaf_transmittance_spectrum", "leaf_transmittance", wavelength)
    from .optical_urban import urban_reflectance
    return urban_reflectance(settings, float(vza), float(vaa), float(sza), float(saa),
                             rho_roof, rho_wall, rho_street, leaf_rho, leaf_tau)


def vegetation_brightness(settings, state, wavelength, vza=0.0, vaa=0.0):
    tss, tsh, tls, tlh = state["Tss"], state["Tsh"], state["Tls"], state["Tlh"]
    vza_ = np.asarray([vza], dtype=float)
    sza_ = np.asarray([min(max(state["sza"], 0.0), 89.0)], dtype=float)
    vaa_ = np.asarray([vaa], dtype=float)
    saa_ = np.asarray([state["saa"]], dtype=float)
    vsa_ = relative_azimuth(saa_, vaa_)
    lai = settings["lai"]
    soil_emissivity = float(settings["soil_emissivity"])
    if int(settings.get("hapke_enable", settings.get("hapke_thermal_enable", 0))) == 1:
        soil_emissivity = float(hapke_thermal_emissivity(
            vza,
            settings.get("hapke_K", 1.0),
            settings["soil_reflectance"],
            settings.get("soil_moisture", 0.0),
            settings.get("soil_tir_water_effect", 0),
            settings.get("soil_tir_moisture_dry", 0.05),
            settings.get("soil_tir_moisture_wet", settings.get("satwater", 0.45)),
            settings.get("soil_tir_emissivity_dry", settings["soil_emissivity"]),
            settings.get("soil_tir_emissivity_wet", 0.99),
        ))
    if geometry_mode(settings) != "legacy" and is_forest(settings):
        if settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 1)) == 2:
            raise ValueError("Dual thermal crown geometry currently requires a plane or one local slope.")
        radiance, temperature, _ = thermal_emission(
            dict(settings, soil_emissivity=soil_emissivity), [tss, tsh, tls, tlh],
            wavelength, state["sza"], vza, relative_azimuth(state["saa"], vaa),
            sunlit=state.get("sunlit_surface", True))
        return float(radiance[0]), float(temperature[0])
    pcom = component_proportions(
        settings, state["sza"], state["saa"], vza, vaa,
        sunlit=state.get("sunlit_surface", True))
    ecom = np.asarray([[soil_emissivity], [soil_emissivity], [settings["leaf_emissivity"]], [settings["leaf_emissivity"]]])
    rcom = planck(wavelength, np.asarray([tss, tsh, tls, tlh], dtype=float))
    direct, _ = emissivity_direct(pcom, ecom)
    if lai <= 0.0:
        radiance = float(np.asarray(radiance_direct(direct, rcom)).reshape(-1)[0])
        return radiance, float(np.asarray(inv_planck(wavelength, radiance)).reshape(-1)[0])
    scatter, _ = emissivity_scatter_analytical(lai, np.asarray([settings["leaf_emissivity"]]), np.asarray([settings["soil_emissivity"]]), vza_, sza_)
    radiance = float(np.asarray(radiance_direct(direct, rcom) + radiance_scatter(scatter, rcom)).reshape(-1)[0])
    if wavelength >= 8.0:
        directional_rho = legacy_directional_reflectance(
            settings, state["sza"], state["saa"], vza, vaa,
            settings["leaf_reflectance"], settings["leaf_transmittance"],
            settings["soil_reflectance"],
            sunlit=state.get("sunlit_surface", True),
            use_hapke=False,
        )
        reference_rho = legacy_directional_reflectance(
            settings, state["sza"], state["saa"], vza, state["saa"] + 90.0,
            settings["leaf_reflectance"], settings["leaf_transmittance"],
            settings["soil_reflectance"],
            sunlit=state.get("sunlit_surface", True),
            use_hapke=False,
        )
        effective_temperature = float(np.mean([tss, tsh, tls, tlh]))
        radiance -= (directional_rho - reference_rho) * float(
            np.asarray(planck(wavelength, effective_temperature)).reshape(-1)[0]
        )
        radiance = max(radiance, 1e-12)
    return radiance, float(np.asarray(inv_planck(wavelength, radiance)).reshape(-1)[0])


def urban_brightness(settings, state, wavelength, vza=0.0, vaa=0.0):
    urban = Urban()
    urban.set_strcutural_input(_urban_shapes(settings))
    roof_sun = state.get("Troof_sunlit", state["Tls"])
    roof_shade = state.get("Troof_shaded", state["Tlh"])
    wall_sun = state.get("Twall_sunlit", state["Tls"])
    wall_shade = state.get("Twall_shaded", state["Tlh"])
    urban.set_thermal_input(roof_sun, roof_shade, wall_sun, wall_shade, state["Tss"], state["Tsh"])
    urban.set_angular_input(np.asarray([vza]), np.asarray([vaa]), state["sza"], state["saa"])
    urban.set_spectral_input(
        _urban_spectral_value(settings, "urban_emissivity_street_spectrum", wavelength, "urban_emissivity_street"),
        _urban_spectral_value(settings, "urban_emissivity_wall_spectrum", wavelength, "urban_emissivity_wall"),
        _urban_spectral_value(settings, "urban_emissivity_roof_spectrum", wavelength, "urban_emissivity_roof"),
    )
    e = np.asarray(urban.calculate_effective_component_emissivity(2), dtype=float).reshape(-1)
    temperatures = np.asarray([wall_sun, wall_shade, state["Tss"], state["Tsh"], roof_sun, roof_shade])
    radiance = float(np.sum(e * planck(wavelength, temperatures)))
    if float(settings.get("lai", 0.0)) > 0.0:
        canopy_radiance, _ = vegetation_brightness(settings, state, wavelength, vza, vaa)
        fractions = urban_component_fractions(settings, vza, vaa, state["sza"], state["saa"])
        street = float(fractions[2] + fractions[3])
        street_radiance = float(np.sum(e[2:4] * planck(wavelength, np.asarray([state["Tss"], state["Tsh"]]))))
        radiance += urban_canopy_fraction(settings, vza) * (street * canopy_radiance - street_radiance)
    radiance = max(radiance, 1.0e-12)
    return radiance, float(np.asarray(inv_planck(wavelength, radiance)).reshape(-1)[0])


__all__ = ["urban_canopy_fraction", "urban_component_fractions", "urban_directional_reflectance", "vegetation_brightness", "urban_brightness"]
