"""Explicit VNIR interface with separated direct/diffuse canopy transport.

Forests use crown-volume first scattering and conservative component exchange;
rows retain effective-slab diffuse transport and urban surfaces their mixture.
These approximations require scene-specific validation.
"""
import math

from .brightness import urban_directional_reflectance
from .radiation import directional_reflectance


def reflectance(settings, wavelength_nm, vza, vaa, sza, saa):
    wave = float(wavelength_nm)
    if not math.isfinite(wave) or not 400 <= wave <= 2500:
        raise ValueError("Optical wavelength_nm must be in [400, 2500].")
    angles = (vza, vaa, sza, saa)
    if not all(math.isfinite(float(a)) for a in angles):
        raise ValueError("Observation angles must be finite.")
    if not 0 <= vza < 90 or not 0 <= sza < 90:
        raise ValueError("VNIR requires above-horizon sun and view directions.")
    surface = settings.get("surface_model", "plane")
    if surface not in ("plane", "slope", "terrain", "urban"):
        raise ValueError("Optical surface must be plane, slope, terrain or urban.")
    if surface in ('slope', 'terrain') and settings.get('vegetation_model') not in ('crown', 'forest'):
        raise ValueError('The slope/terrain optical adapter requires a forest.')
    if surface == 'terrain' and int(settings.get('terrain_type', 1)) != 2:
        raise ValueError('The optical terrain adapter requires compound terrain_type=2.')
    if surface == "slope":
        from .solar import slope_direction
        slope, aspect = float(settings.get('slope_angle', 0.)), float(settings.get('slope_aspect', 0.))
        local_view, local_vaa, visible = slope_direction(vza, vaa, slope, aspect)
        local_sun, local_saa, lit = slope_direction(sza, saa, slope, aspect)
        if not visible:
            raise ValueError('Optical view is behind the slope surface.')
        if not lit:
            local_sun = min(local_sun, 89.9)
        return reflectance(dict(settings, surface_model='plane',
                                shortwave_sunlit_surface=lit and settings.get('shortwave_sunlit_surface', True)),
                           wave, local_view, local_vaa, local_sun, local_saa)
    grid = list(settings["wavelengths"])
    if not grid or wave < min(grid) or wave > max(grid):
        raise ValueError("Requested wavelength is outside the supplied optical grid.")
    # Interpolate explicitly in nm, never reinterpret thermal emissivity as VNIR.
    import numpy as np
    order = np.argsort(grid)
    if len(set(grid)) != len(grid) or not np.isfinite(grid).all():
        raise ValueError("Optical grid must contain unique finite wavelengths.")
    local = dict(settings)
    for key in ("leaf_reflectance", "leaf_transmittance", "soil_reflectance",
                "urban_reflectance_roof", "urban_reflectance_wall", "urban_reflectance_street"):
        values = settings.get(key + "_spectrum")
        if values is not None:
            if len(values) != len(grid) or not np.isfinite(values).all():
                raise ValueError(f"Invalid {key}_spectrum.")
            if np.any(np.asarray(values) < 0) or np.any(np.asarray(values) > 1):
                raise ValueError(f"{key} must be in [0, 1].")
            local[key] = float(np.interp(wave, np.asarray(grid)[order], np.asarray(values)[order]))
            local.pop(key + "_spectrum", None)
        if key in local and (not math.isfinite(float(local[key])) or not 0 <= float(local[key]) <= 1):
            raise ValueError(f"{key} must be finite and in [0, 1].")
    if local["leaf_reflectance"] + local["leaf_transmittance"] > 1 + 1e-12:
        raise ValueError("Leaf reflectance + transmittance must not exceed one.")
    if not math.isfinite(float(local["lai"])) or local["lai"] < 0:
        raise ValueError("LAI must be finite and nonnegative.")
    if surface == "urban":
        return urban_directional_reflectance(local, wave / 1000, vza, vaa, sza, saa)
    return directional_reflectance(local, sza, saa, vza, vaa,
                                  local["leaf_reflectance"], local["leaf_transmittance"],
                                  local["soil_reflectance"], use_hapke=False,
                                  sunlit=bool(local.get('shortwave_sunlit_surface', True)))
