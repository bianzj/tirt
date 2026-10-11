from functools import lru_cache

import numpy as np

from .radiation import directional_reflectance, urban_canopy_fraction
from .urban import Urban


_WAVELENGTH_GRIDS = {}


def nearest_wavelength_index(wavelengths, wavelength):
    """Return int(np.argmin(np.abs(np.asarray(wavelengths, dtype=float) - float(wavelength))))."""
    if not isinstance(wavelengths, list):
        return int(np.argmin(np.abs(np.asarray(wavelengths, dtype=float) - float(wavelength))))
    entry = _WAVELENGTH_GRIDS.get(id(wavelengths))
    if entry is None or entry[0] is not wavelengths or entry[1] != wavelengths:
        if len(_WAVELENGTH_GRIDS) >= 8:
            _WAVELENGTH_GRIDS.clear()
        entry = (wavelengths, list(wavelengths), np.asarray(wavelengths, dtype=float), {})
        _WAVELENGTH_GRIDS[id(wavelengths)] = entry
    target = float(wavelength)
    index = entry[3].get(target)
    if index is None:
        index = int(np.argmin(np.abs(entry[2] - target)))
        entry[3][target] = index
    return index


def _urban_spectral_value(settings, key, wavelength, default):
    spectrum = settings.get(key, [])
    wavelengths = settings.get("wavelengths", [])
    index = nearest_wavelength_index(wavelengths, wavelength) if wavelengths else -1
    if len(spectrum) == len(wavelengths) and 0 <= index < len(spectrum):
        return float(spectrum[index])
    return float(settings.get(default, default if isinstance(default, (int, float)) else 0.9))


def _vegetation_spectral_value(settings, key, default, wavelength):
    spectrum = settings.get(key, [])
    wavelengths = settings.get("wavelengths", [])
    if len(spectrum) == len(wavelengths) and wavelengths:
        index = nearest_wavelength_index(wavelengths, wavelength)
        return float(spectrum[index])
    return float(settings.get(default, 0.0))


def _urban_shapes(settings):
    """Build the urban geometry expected by the TIRT urban projection model."""
    length = max(float(settings.get("urban_building_length", 10.0)), 0.0)
    width = max(float(settings.get("urban_building_width", 10.0)), 0.0)
    height = max(float(settings.get("urban_building_height", 20.0)), 0.0)
    density = max(float(settings.get("urban_density", 0.0015)), 0.0)
    azimuth = float(settings.get("urban_building_azimuth", 0.0))
    return np.asarray([[length, width, height, density, azimuth, azimuth + 90.0]], dtype=float)


def urban_component_fractions(settings, vza, vaa, sza, saa):
    """Return urban wall/street/roof sunlit-shaded view fractions."""
    shape = tuple(_urban_shapes(settings)[0])
    return np.array(_urban_component_fractions(shape, float(vza), float(vaa), float(sza), float(saa)))


@lru_cache(maxsize=65536)
def _urban_component_fractions(shape, vza, vaa, sza, saa):
    urban = Urban()
    urban.set_strcutural_input(np.asarray([shape], dtype=float))
    urban.set_angular_input(
        np.asarray([vza], dtype=float),
        np.asarray([vaa], dtype=float),
        float(sza),
        float(saa),
    )
    values = np.asarray(urban.calculate_direct_emissivity(0, 3), dtype=float).reshape(-1)
    values = np.maximum(values, 0.0)
    total = float(values.sum())
    return tuple(values / total if total > 0.0 else np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))


def urban_directional_reflectance(settings, wavelength, vza, vaa, sza, saa):
    """Mix urban surfaces with a canopy layer over the visible street."""
    fractions = urban_component_fractions(settings, vza, vaa, sza, saa)
    wall_fraction = float(fractions[0] + fractions[1])
    street_fraction = float(fractions[2] + fractions[3])
    roof_fraction = float(fractions[4] + fractions[5])
    rho_roof = _urban_spectral_value(settings, "urban_reflectance_roof_spectrum", wavelength, "urban_reflectance_roof")
    rho_wall = _urban_spectral_value(settings, "urban_reflectance_wall_spectrum", wavelength, "urban_reflectance_wall")
    rho_street = _urban_spectral_value(settings, "urban_reflectance_street_spectrum", wavelength, "urban_reflectance_street")
    city = wall_fraction * rho_wall + street_fraction * rho_street + roof_fraction * rho_roof
    if float(settings.get("lai", 0.0)) <= 0.0 or street_fraction <= 0.0:
        return float(np.clip(city, 0.0, 1.0))
    leaf_rho = _vegetation_spectral_value(settings, "leaf_reflectance_spectrum", "leaf_reflectance", wavelength)
    leaf_tau = _vegetation_spectral_value(settings, "leaf_transmittance_spectrum", "leaf_transmittance", wavelength)
    soil_rho = _vegetation_spectral_value(settings, "soil_reflectance_spectrum", "soil_reflectance", wavelength)
    canopy = directional_reflectance(
        settings, sza, saa, vza, vaa, leaf_rho, leaf_tau, soil_rho,
        sunlit=True, use_hapke=False,
    )
    mixed = city + street_fraction * urban_canopy_fraction(settings, vza) * (canopy - rho_street)
    return float(np.clip(mixed, 0.0, 1.0))


__all__ = ["nearest_wavelength_index", "urban_canopy_fraction", "urban_component_fractions", "urban_directional_reflectance"]
