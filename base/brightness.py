import numpy as np

from .urban import Urban


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
    bearing = float(settings.get("urban_building_azimuth", settings.get("building_azimuth", 0.0)))
    return np.asarray([[length, width, height, density, bearing, bearing+90.0]], dtype=float)


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


def urban_canopy_fraction(settings, vza=0.0, vaa=0.0):
    """Visible leaf fraction of the street canopy along a view: 1 - view gap of the standard closure."""
    from .tirteb_radiation import view_fractions
    parts = view_fractions(dict(settings, surface_model="plane"), vza, vaa, vza, vaa)
    return float(np.clip(parts[2] + parts[3], 0.0, 1.0))


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


__all__ = ["urban_canopy_fraction", "urban_component_fractions", "urban_directional_reflectance"]
