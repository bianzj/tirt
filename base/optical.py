"""Explicit VNIR interface with separated direct/diffuse canopy transport.

Forests use crown-volume first scattering and conservative component exchange;
rows retain effective-slab diffuse transport and urban surfaces their mixture.
These approximations require scene-specific validation.
"""
import math

from .brightness import urban_directional_reflectance
from .radiation import directional_reflectance


class InvisibleDirection(ValueError):
    """The sun or view direction is below the horizon or behind the slope."""


def reflectance(settings, wavelength_nm, vza, vaa, sza, saa):
    wave = float(wavelength_nm)
    if not math.isfinite(wave) or not 400 <= wave <= 2500:
        raise ValueError("Optical wavelength_nm must be in [400, 2500].")
    angles = (vza, vaa, sza, saa)
    if not all(math.isfinite(float(a)) for a in angles):
        raise ValueError("Observation angles must be finite.")
    if not 0 <= vza < 90 or not 0 <= sza < 90:
        raise InvisibleDirection("VNIR requires above-horizon sun and view directions.")
    surface = settings.get("surface_model", "plane")
    if surface not in ("plane", "slope", "terrain", "urban"):
        raise ValueError("Optical surface must be plane, slope, terrain or urban.")
    if (surface == 'terrain' and settings.get('vegetation_model') in ('crown', 'forest')
            and int(settings.get('terrain_type', 2)) != 2):
        raise ValueError('The optical terrain forest adapter requires compound terrain_type=2.')
    if surface == "slope":
        # Single-slope BRF is referred to the horizontal global irradiance, as in the slope
        # energy balance: local direct = max(cos i, 0)/cos(sza), local sky = (1+cos(slope))/2
        # (isotropic sky, no reflection from surrounding terrain).
        from .solar import slope_direction
        slope, aspect = float(settings.get('slope_angle', 0.)), float(settings.get('slope_aspect', 0.))
        local_view, local_vaa, visible = slope_direction(vza, vaa, slope, aspect)
        local_sun, local_saa, _ = slope_direction(sza, saa, slope, aspect)
        if not visible:
            raise InvisibleDirection('Optical view is behind the slope surface.')
        configured = float(settings.get('shortwave_diffuse_fraction', 0.2))
        fraction = 0.2 if configured < 0.0 else configured
        b, z = math.radians(slope), math.radians(float(sza))
        incidence = (math.cos(b) * math.cos(z)
                     + math.sin(b) * math.sin(z) * math.cos(math.radians(float(saa) - aspect)))
        sunlit = bool(settings.get('shortwave_sunlit_surface', True))
        direct = (1.0 - fraction) * max(incidence, 0.0) / math.cos(z) if sunlit else 0.0
        sky = fraction * (1.0 + math.cos(b)) / 2.0
        local = dict(settings, surface_model='plane')
        from .canopy_closure import closure
        if closure(settings) != 'legacy' and str(settings.get('vegetation_model', '')).lower() in ('row', 'raw'):
            local['row_azimuth'] = slope_direction(90.0, float(settings.get('row_azimuth', 0.)), slope, aspect)[1]
        sun = (min(local_sun, 89.9), local_saa)
        value = 0.0
        if direct > 0.0:
            value += direct * reflectance(dict(local, shortwave_diffuse_fraction=0.0, shortwave_sunlit_surface=True),
                                          wave, local_view, local_vaa, *sun)
        if sky > 0.0:
            value += sky * reflectance(dict(local, shortwave_diffuse_fraction=1.0),
                                       wave, local_view, local_vaa, *sun)
        return value
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
    if local.get("lai") in (None, ""):
        raise ValueError("Optical settings require 'lai' (leaf area index).")
    if not math.isfinite(float(local["lai"])) or local["lai"] < 0:
        raise ValueError("LAI must be finite and nonnegative.")
    if surface == "urban":
        return urban_directional_reflectance(local, wave / 1000, vza, vaa, sza, saa)
    if surface == "terrain":
        local.setdefault("terrain_type", 2)
        local.setdefault("surface_type", 0)
        local.setdefault("distribution", 2 if str(local.get("vegetation_model", "")).lower() in ("row", "raw") else 0)
    return directional_reflectance(local, sza, saa, vza, vaa,
                                  local["leaf_reflectance"], local["leaf_transmittance"],
                                  local["soil_reflectance"], use_hapke=False,
                                  sunlit=bool(local.get('shortwave_sunlit_surface', True)))
