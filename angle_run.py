import argparse
import csv
import itertools
from pathlib import Path

import numpy as np

from base.crown import Crown
from base.hom import Hom
from base.inputio import (
    interpolate,
    numbers,
    read_geometry,
    read_input,
    read_spectrum,
    read_wavelengths,
    repeat_to,
    resolve_path,
    value,
    write_input_snapshot,
)
from base.physicsF import inv_planck, planck, slope1
from base.row import Row
from base.soilopt import hapke_thermal_emissivity


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input_tirteb.csv"
OUTPUT_FIELDS = (
    "index", "scenario_id", "scenario_label", "surface_model", "vegetation_model", "model", "wavelength_um",
    "vza", "vaa", "sza", "saa", "raa", "emissivity_leaf", "emissivity_soil",
    "emissivity_roof", "emissivity_wall", "emissivity_street",
    "radiance", "brightness_temperature_K", "brightness_temperature_C",
)

SCENARIO_PARAMETERS = (
    ("lai", ("lai", "vegetation_lai"), "LAI"),
    ("hspot", ("hspot", "vegetation_hspot"), "Hotspot"),
    ("sza", ("sza", "solar_zenith"), "SZA"),
    ("saa", ("saa", "solar_azimuth"), "SAA"),
    ("slope_angle", ("slope_angle", "terrain_slope"), "Slope angle"),
    ("slope_aspect", ("slope_aspect", "terrain_aspect"), "Slope aspect"),
    ("row_width", ("row_width",), "Row width"),
    ("row_blank", ("row_blank",), "Row gap"),
    ("row_height", ("row_height",), "Row height"),
    ("row_azimuth", ("row_azimuth",), "Row azimuth"),
    ("crown_rad_a", ("crown_rad_a", "crown_radius"), "Crown rad-a"),
    ("crown_rad_b", ("crown_rad_b", "crown_width"), "Crown rad-b"),
    ("crown_height", ("crown_height",), "Crown height"),
    ("crown_density", ("crown_density",), "Crown density"),
)


def _float(config, *keys, default=0.0):
    return float(value(config, *keys, default=default))


def _temperature(config, name, default_c):
    kelvin = value(config, f"{name}_k")
    if kelvin not in (None, ""):
        return float(kelvin)
    return _float(config, f"{name}_c", default=default_c) + 273.15


def _stepped_angles(maximum, step):
    angles = list(np.arange(0.0, maximum + 1.0e-12, step))
    if not angles or not np.isclose(angles[-1], maximum):
        angles.append(maximum)
    return tuple(angles)


def _canonical(config):
    surface = str(value(config, "surface_model", "surface", default="plane")).strip().lower()
    surface = {"flat": "plane", "bare": "plane", "single_slope": "slope", "sloped": "slope", "building": "urban", "city": "urban"}.get(surface, surface)
    vegetation = str(value(config, "vegetation_model", "vegetation", default="hom")).strip().lower()
    vegetation = {"homogeneous": "hom", "raw": "row", "row_crop": "row", "forest": "crown", "bare_soil": "bare"}.get(vegetation, vegetation)
    if surface not in {"plane", "slope", "terrain", "urban"}:
        raise ValueError("surface_model 必须是 plane、slope、terrain 或 urban")
    if vegetation not in {"bare", "hom", "row", "crown"}:
        raise ValueError("vegetation_model 必须是 bare、hom、row 或 crown")
    return surface, vegetation


def _geometry(config, input_path):
    from run import _geometry as shared_geometry
    return shared_geometry(config, input_path)


def _wavelengths(config, input_path):
    source = str(value(config, "spectral_source", "wavelength_source", default="direct")).strip().lower()
    if source in {"file", "txt", "read"}:
        path = resolve_path(input_path, value(config, "spectrals_file", "wavelengths_file", default="wavelengths.txt"))
        return np.asarray(read_wavelengths(path), dtype=float)
    raw = value(config, "wavelengths", "observation_wavelengths", "wavelength", default=10.5)
    result = np.asarray(numbers(raw), dtype=float)
    if result.size == 0 or np.any(result <= 0.0):
        raise ValueError("wavelengths 必须包含正数波长")
    return result


def _typical_band_values(config, keys, wavelengths, defaults, label):
    raw = value(config, *keys)
    if raw in (None, ""):
        return np.interp(wavelengths, np.asarray([8.0, 10.5, 12.0]), np.asarray(defaults, dtype=float))
    result = np.asarray(numbers(raw), dtype=float)
    if result.size == 1:
        return np.repeat(result.item(), len(wavelengths))
    if result.size != len(wavelengths):
        raise ValueError(f"{label} 的数目必须与 wavelengths 的数目一致: {result.size} != {len(wavelengths)}")
    return result


def _spectra(config, input_path, wavelengths):
    size = len(wavelengths)
    file_mode = {"file", "txt", "read"}
    leaf_from_file = str(value(config, "vegetation_spectrum_source", default="direct")).strip().lower() in file_mode
    soil_from_file = str(value(config, "soil_spectrum_source", default="direct")).strip().lower() in file_mode
    terrain_from_file = str(value(config, "terrain_spectrum_source", default="direct")).strip().lower() in file_mode
    urban_from_file = str(value(config, "urban_spectrum_source", default="direct")).strip().lower() in file_mode
    roof_from_file = str(value(config, "roof_spectrum_source", default="direct")).strip().lower() in file_mode
    wall_from_file = str(value(config, "wall_spectrum_source", default="direct")).strip().lower() in file_mode
    street_from_file = str(value(config, "street_spectrum_source", default="direct")).strip().lower() in file_mode
    leaf_default = (0.975, 0.985, 0.980)
    soil_default = (0.920, 0.955, 0.965)
    terrain_default = (0.920, 0.950, 0.960)
    roof_default = (0.920, 0.950, 0.940)
    wall_default = (0.900, 0.930, 0.920)
    street_default = (0.940, 0.955, 0.950)
    leaf_e = np.interp(wavelengths, [8.0, 10.5, 12.0], leaf_default) if leaf_from_file else _typical_band_values(config, ("emissivity_leaf", "leaf_emissivity"), wavelengths, leaf_default, "emissivity_leaf")
    soil_e = np.interp(wavelengths, [8.0, 10.5, 12.0], soil_default) if soil_from_file else _typical_band_values(config, ("emissivity_soil", "soil_emissivity"), wavelengths, soil_default, "emissivity_soil")
    terrain_e = np.interp(wavelengths, [8.0, 10.5, 12.0], terrain_default) if terrain_from_file else _typical_band_values(config, ("emissivity_terrain", "terrain_emissivity"), wavelengths, terrain_default, "emissivity_terrain")
    roof_e = np.interp(wavelengths, [8.0, 10.5, 12.0], roof_default) if urban_from_file or roof_from_file else _typical_band_values(config, ("emissivity_roof", "urban_emissivity_roof"), wavelengths, roof_default, "emissivity_roof")
    wall_e = np.interp(wavelengths, [8.0, 10.5, 12.0], wall_default) if urban_from_file or wall_from_file else _typical_band_values(config, ("emissivity_wall", "urban_emissivity_wall"), wavelengths, wall_default, "emissivity_wall")
    street_e = np.interp(wavelengths, [8.0, 10.5, 12.0], street_default) if urban_from_file or street_from_file else _typical_band_values(config, ("emissivity_street", "urban_emissivity_street"), wavelengths, street_default, "emissivity_street")

    if leaf_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "vegetation_spectrum_file")))
        leaf_e = interpolate(table, wavelengths, 1)

    if soil_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "soil_spectrum_file")))
        soil_e = interpolate(table, wavelengths, 1)

    if terrain_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "terrain_spectrum_file")))
        terrain_e = interpolate(table, wavelengths, 1)

    if urban_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "urban_spectrum_file")))
        if table.shape[1] >= 7:
            roof_e, wall_e, street_e = (interpolate(table, wavelengths, col) for col in (4, 5, 6))
        elif table.shape[1] >= 4:
            roof_e, wall_e, street_e = (interpolate(table, wavelengths, col) for col in (1, 2, 3))
        else:
            raise ValueError("urban_spectrum_file 至少需要 wavelength roof wall street 四列")

    for component, enabled, key in (
        ("roof", roof_from_file, "roof_spectrum_file"),
        ("wall", wall_from_file, "wall_spectrum_file"),
        ("street", street_from_file, "street_spectrum_file"),
    ):
        if enabled:
            table = read_spectrum(resolve_path(input_path, value(config, key)))
            values = interpolate(table, wavelengths, 1)
            if component == "roof": roof_e = values
            elif component == "wall": wall_e = values
            else: street_e = values

    result = {
        "leaf_e": np.asarray(leaf_e), "soil_e": np.asarray(soil_e), "terrain_e": np.asarray(terrain_e),
        "roof_e": np.asarray(roof_e), "wall_e": np.asarray(wall_e), "street_e": np.asarray(street_e),
    }
    for component, emissivities in result.items():
        if not np.all(np.isfinite(emissivities)) or np.any(emissivities < 0.0) or np.any(emissivities > 1.0):
            raise ValueError(f"{component} emissivity must be finite and between 0 and 1")
    return result


def _vegetation_temperatures(config):
    return (
        _temperature(config, "soil_sunlit", 45.0), _temperature(config, "soil_shaded", 30.0),
        _temperature(config, "leaf_sunlit", 33.0), _temperature(config, "leaf_shaded", 30.0),
    )


def _soil_emissivity(config, view_zenith, base_emissivity):
    enabled = str(value(config, "hapke_enable", "soil_hapke_enable", default="0")).strip().lower()
    if enabled in {"0", "false", "no", "off"}:
        return float(base_emissivity)
    albedo = _float(config, "soil_reflectance", "soil_albedo", default=1.0 - float(base_emissivity))
    K = _float(config, "hapke_K", "soil_hapke_K", default=1.0)
    result = hapke_thermal_emissivity(view_zenith, K, albedo)
    return float(np.asarray(result).reshape(-1)[0])


def _shape_values(config, key, default):
    return numbers(value(config, key), (default,))


def _terrain_shapes(config):
    heights = _shape_values(config, "terrain_height", 20.0)
    radii = _shape_values(config, "terrain_radius", 10.0)
    densities = _shape_values(config, "terrain_density", 0.001)
    size = max(len(heights), len(radii), len(densities))
    shapes = np.column_stack([repeat_to(heights, size, "terrain_height"), repeat_to(radii, size, "terrain_radius"), repeat_to(densities, size, "terrain_density")])
    occupied = float(np.sum(shapes[:, 2] * np.pi * shapes[:, 1] ** 2))
    if occupied >= 1.0:
        # TIRTEB's compound-terrain density is a relative scene density,
        # while the legacy TIRT terrain kernel expects an areal number
        # density. Normalize only the angular compatibility path so the
        # temporal energy-balance input keeps its original meaning.
        shapes[:, 2] *= 0.8 / max(occupied, 1.0e-12)
    return shapes


def _urban_shapes(config):
    lengths = _shape_values(config, "building_length", 10.0)
    widths = _shape_values(config, "building_width", 10.0)
    heights = _shape_values(config, "building_height", 20.0)
    densities = _shape_values(config, "building_density", 0.0015)
    azimuths = _shape_values(config, "building_azimuth", 0.0)
    size = max(map(len, (lengths, widths, heights, densities, azimuths)))
    azimuths = repeat_to(azimuths, size, "building_azimuth")
    return np.column_stack([
        repeat_to(lengths, size, "building_length"), repeat_to(widths, size, "building_width"),
        repeat_to(heights, size, "building_height"), repeat_to(densities, size, "building_density"),
        azimuths, azimuths + 90.0,
    ])


def _urban_model_geometry(geometry):
    """Use a defined internal azimuth for nadir urban observations.

    At VZA=0 the view azimuth is physically undefined. The urban scattering
    approximation uses relative azimuth to split wall emission into sunlit and
    shaded terms, so retaining an arbitrary input VAA creates a false 0/180
    difference. Keep the reported geometry unchanged and use the solar azimuth
    only inside the urban model at nadir.
    """
    vza, vaa, sza, saa = geometry
    if abs(float(vza)) <= 1.0e-8:
        return float(vza), float(saa), float(sza), float(saa)
    return geometry


def _checked_radiance(radiance, wavelength, model, geometry):
    radiance = float(np.asarray(radiance).reshape(-1)[0])
    if not np.isfinite(radiance) or radiance < 0.0:
        raise FloatingPointError(
            f"Invalid radiance for {model}, wavelength={wavelength}, geometry={geometry}: {radiance}")
    brightness = 0.0 if radiance == 0.0 else float(inv_planck(wavelength, radiance))
    if not np.isfinite(brightness):
        raise FloatingPointError(f"Invalid brightness temperature for {model}, geometry={geometry}")
    return radiance, brightness


def _checked_components(emissivity, model, geometry, wavelength):
    emissivity = np.asarray(emissivity, float).reshape(-1)
    if (not np.all(np.isfinite(emissivity)) or np.any(emissivity < -1.0e-10) or
            float(np.sum(emissivity)) > 1.0 + 1.0e-7):
        raise FloatingPointError(
            f"Invalid component emissivities for {model}, wavelength={wavelength}, geometry={geometry}: {emissivity.tolist()}")
    return np.maximum(emissivity, 0.0)


def _set_canopy_geometry(model, config, vegetation):
    model.set_vegetation_input(
        vegetation, _float(config, "row_width", default=0.5),
        _float(config, "row_blank", default=0.5),
        _float(config, "row_height", default=1.0),
        _float(config, "row_azimuth", default=0.0))


def _vegetation_result(config, vegetation, wavelength, geometry, spectra):
    vza, vaa, sza, saa = geometry
    # The analytical crown and hotspot equations are singular exactly at the
    # horizon. Keep the reported geometry unchanged, but evaluate the model at
    # the nearest finite angle.
    model_vza = min(abs(float(vza)), 89.999)
    model_sza = min(abs(float(sza)), 89.999)
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    if vegetation == "bare":
        soil_emissivity = _soil_emissivity(config, model_vza, spectra["soil_e"])
        radiance = soil_emissivity * planck(wavelength, soil_sun)
        return _checked_radiance(radiance, wavelength, "plane+bare", geometry)
    lai = _float(config, "lai", "vegetation_lai", default=1.5)
    hspot = _float(config, "hspot", "vegetation_hspot", default=0.15)
    if lai == 0.0 or (vegetation == "crown" and _float(config, "crown_density", default=0.02) == 0.0):
        return _vegetation_result(config, "bare", wavelength, geometry, spectra)
    if vegetation == "hom":
        model = Hom()
        model.set_structure(lai, hspot)
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    elif vegetation == "row":
        model = Row()
        model.set_structure(lai, hspot, _float(config, "row_width", default=0.5), _float(config, "row_blank", default=0.5), _float(config, "row_height", default=1.0))
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([vaa]), np.asarray([saa]), np.asarray([_float(config, "row_azimuth", default=0.0)]))
    else:
        model = Crown()
        model.set_hotspot_geometry(config.get("thermal_crown_geometry", "legacy"),
                                   config.get("thermal_hotspot_inner"), config.get("thermal_hotspot_outer"))
        rad_a = max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6)
        rad_b = max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6)
        model.set_structure(lai, hspot, _float(config, "crown_density", default=0.02), rad_b, rad_a)
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance = float(np.asarray(model.run(ifradiance=1)).reshape(-1)[0])
    return _checked_radiance(radiance, wavelength, f"plane+{vegetation}", geometry)


def _slope_local_geometry(config, geometry):
    """Transform global view and solar directions into slope-local angles."""
    vza, vaa, sza, saa = geometry
    slope_angle = _float(config, "slope_angle", "terrain_slope", default=0.0)
    slope_aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)

    def local_direction(zenith, azimuth):
        with np.errstate(divide="ignore", invalid="ignore"):
            local_zenith, local_azimuth = slope1(
                np.asarray([abs(float(zenith))], dtype=float),
                np.asarray([float(azimuth)], dtype=float),
                slope_angle,
                slope_aspect,
            )
        return float(local_zenith[0]), float(local_azimuth[0])

    local_vza, local_vaa = local_direction(vza, vaa)
    local_sza, local_saa = local_direction(sza, saa)
    return local_vza, local_vaa, local_sza, local_saa


def _slope_result(config, vegetation, wavelength, geometry, spectra):
    if vegetation == "bare":
        slope = _float(config, "slope_angle", "terrain_slope", default=0.0)
        aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)
        cosine_sun = (np.cos(np.deg2rad(slope)) * np.cos(np.deg2rad(geometry[2])) +
                      np.sin(np.deg2rad(slope)) * np.sin(np.deg2rad(geometry[2])) *
                      np.cos(np.deg2rad(geometry[3] - aspect)))
        local_config = config.copy()
        if cosine_sun <= 0.0:
            local_config["soil_sunlit_k"] = str(_temperature(config, "soil_shaded", 30.0))
        return _vegetation_result(local_config, vegetation, wavelength, _slope_local_geometry(config, geometry), spectra)

    from base.slope_veg import Slope_Veg

    vza, vaa, sza, saa = geometry
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    lai = _float(config, "lai", "vegetation_lai", default=1.5)
    hspot = _float(config, "hspot", "vegetation_hspot", default=0.15)
    stand_density = _float(config, "slope_density", "crown_density", default=0.02)
    crown_rad_a = max(_float(config, "slope_crown_rad_a", "crown_rad_a", "slope_crown_radius", "crown_radius", default=3.0), 1.0e-6)
    crown_rad_b = max(_float(config, "slope_crown_rad_b", "crown_rad_b", "slope_crown_width", "crown_width", default=3.0), 1.0e-6)
    slope_angle = _float(config, "slope_angle", "terrain_slope", default=0.0)
    slope_aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)

    model = Slope_Veg()
    model.set_angle(np.asarray([abs(vza)]), np.asarray([sza]), np.asarray([vaa]), np.asarray([saa]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_structure(lai, hspot, stand_density, crown_rad_b, crown_rad_a)
    _set_canopy_geometry(model, config, vegetation)
    model.set_slope(slope_angle, slope_aspect)
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    return _checked_radiance(model.run(ifradiance=1), wavelength, f"slope+{vegetation}", geometry)


def _slope_backface(config, geometry):
    """Return whether the slope-corrected view zenith exceeds 89 degrees."""
    vza, vaa = geometry[0], geometry[1]
    slope_angle = _float(config, "slope_angle", "terrain_slope", default=0.0)
    slope_aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)
    # Dot product with the slope normal retains the backside information that
    # slope1 folds into the upper hemisphere for vegetation calculations.
    cosine = (
        np.cos(np.deg2rad(slope_angle)) * np.cos(np.deg2rad(float(vza)))
        + np.sin(np.deg2rad(slope_angle)) * np.sin(np.deg2rad(float(vza)))
        * np.cos(np.deg2rad(float(vaa) - slope_aspect))
    )
    corrected_vza = np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return corrected_vza > 89.0 + 1.0e-9


def _surface_result(config, surface, vegetation, wavelength, geometry, spectra):
    if vegetation == "crown" and str(config.get("thermal_crown_geometry", "legacy")).strip().lower() != "legacy" and surface != "plane":
        raise ValueError("New thermal crown geometry in angle_run.py requires surface=plane.")
    vza, vaa, sza, saa = geometry
    if surface == "plane":
        return _vegetation_result(config, vegetation, wavelength, geometry, spectra)
    if surface == "slope":
        if _slope_backface(config, geometry):
            return float("nan"), float("nan")
        return _slope_result(config, vegetation, wavelength, geometry, spectra)
    if surface == "terrain":
        from base.terrain import Terrain
        from base.terrain_veg import Terrain_Veg

        model = Terrain() if vegetation == "bare" else Terrain_Veg()
        model.set_strcutural_input(_terrain_shapes(config))
        terrain_vza = min(max(abs(vza), 1.0e-3), 89.999)
        terrain_vaa = saa if abs(float(vza)) <= 1.0e-8 else vaa
        model.set_angular_input(np.asarray([terrain_vza]), np.asarray([terrain_vaa]), sza, saa)
        soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
        terrain_sun = _temperature(config, "terrain_sunlit", leaf_sun - 273.15)
        terrain_shade = _temperature(config, "terrain_shaded", leaf_shade - 273.15)
        if vegetation == "bare":
            soil_emissivity = _soil_emissivity(config, terrain_vza, spectra["soil_e"])
            model.set_spectral_input(soil_emissivity, spectra["terrain_e"])
            model.set_thermal_input(soil_sun, soil_shade, terrain_sun, terrain_shade)
            component_temperatures = (soil_sun, soil_shade, terrain_sun, terrain_shade)
        else:
            model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
            _set_canopy_geometry(model, config, vegetation)
            model.set_spectral_input(spectra["soil_e"], spectra["leaf_e"])
            model.set_thermal_input(soil_sun, soil_shade, leaf_sun, leaf_shade)
            component_temperatures = (soil_sun, soil_shade, leaf_sun, leaf_shade)
        emissivity = _checked_components(
            model.calculate_effective_component_emissivity(2), f"terrain+{vegetation}", geometry, wavelength)
        radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures))))
        return _checked_radiance(radiance, wavelength, f"terrain+{vegetation}", geometry)

    from base.urban import Urban
    from base.urban_veg import Urban_Veg

    model = Urban() if vegetation == "bare" else Urban_Veg()
    model.set_strcutural_input(_urban_shapes(config))
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    roof_sun = _temperature(config, "roof_sunlit", 45.0)
    roof_shade = _temperature(config, "roof_shaded", 30.0)
    wall_sun = _temperature(config, "wall_sunlit", 45.0)
    wall_shade = _temperature(config, "wall_shaded", 30.0)
    street_sun = _temperature(config, "street_sunlit", 45.0)
    street_shade = _temperature(config, "street_shaded", 30.0)
    model_geometry = _urban_model_geometry(geometry)
    model_vza, model_vaa, model_sza, model_saa = model_geometry
    # Evaluate the analytical urban kernel at the finite one-sided horizon.
    model_vza = min(abs(model_vza), 89.999)
    model.set_angular_input(np.asarray([model_vza]), np.asarray([model_vaa]), model_sza, model_saa)
    model.set_thermal_input(roof_sun, roof_shade, wall_sun, wall_shade, street_sun, street_shade)
    if vegetation == "bare":
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade)
    else:
        model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
        _set_canopy_geometry(model, config, vegetation)
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"], spectra["leaf_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade, leaf_sun, leaf_shade)
    emissivity = _checked_components(
        model.calculate_effective_component_emissivity(2), f"urban+{vegetation}", geometry, wavelength)
    radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures[:len(emissivity)]))))
    return _checked_radiance(radiance, wavelength, f"urban+{vegetation}", geometry)


def _scenario_configs(config):
    """Expand semicolon-separated scenario parameters into a Cartesian batch."""
    entries = []
    for key, aliases, label in SCENARIO_PARAMETERS:
        raw = value(config, *aliases)
        if raw in (None, ""):
            continue
        values = tuple(float(item) for item in numbers(raw))
        if not values:
            continue
        entries.append((key, aliases, label, values))

    if not entries:
        return [(1, "default", config.copy())]

    varying = [entry for entry in entries if len(entry[3]) > 1]
    scenario_configs = []
    products = itertools.product(*(entry[3] for entry in entries))
    for scenario_id, combination in enumerate(products, 1):
        scenario_config = config.copy()
        label_parts = []
        for (key, aliases, label, values), selected in zip(entries, combination):
            scenario_config[key] = str(selected)
            for alias in aliases:
                scenario_config[alias] = str(selected)
        for entry, selected in zip(entries, combination):
            if len(entry[3]) > 1:
                label_parts.append(f"{entry[2]}={selected:g}")
        scenario_label = "; ".join(label_parts) if varying else "default"
        scenario_configs.append((scenario_id, scenario_label, scenario_config))
    return scenario_configs


def _simulate_rows(config, input_path, surface, vegetation, wavelengths, spectra, scenario_id, scenario_label):
    geometries = _geometry(config, input_path)
    rows = []
    index = 0
    for geometry in geometries:
        vza, vaa, sza, saa = geometry
        raa = abs(vaa - saa) % 360.0
        raa = min(raa, 360.0 - raa)
        for band_index, wavelength in enumerate(wavelengths):
            band_spectra = {key: values[band_index] for key, values in spectra.items()}
            try:
                radiance, brightness_temperature = _surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)
            except FloatingPointError as exc:
                raise FloatingPointError(
                    f"{surface}+{vegetation}, wavelength={wavelength}, geometry={geometry}: {exc}") from exc
            rows.append({
                "index": index, "scenario_id": scenario_id, "scenario_label": scenario_label,
                "surface_model": surface, "vegetation_model": vegetation,
                "model": f"{surface}+{vegetation}", "wavelength_um": float(wavelength),
                "vza": float(vza), "vaa": float(vaa), "sza": float(sza), "saa": float(saa), "raa": float(raa),
                "emissivity_leaf": float(spectra["leaf_e"][band_index]), "emissivity_soil": float(spectra["soil_e"][band_index]),
                "emissivity_roof": float(spectra["roof_e"][band_index]), "emissivity_wall": float(spectra["wall_e"][band_index]), "emissivity_street": float(spectra["street_e"][band_index]),
                "radiance": radiance, "brightness_temperature_K": brightness_temperature, "brightness_temperature_C": brightness_temperature - 273.15,
            })
            index += 1
    return rows


def run_angle(input_path=DEFAULT_INPUT, output_path=None, model=None, **overrides):
    input_path = Path(input_path).resolve()
    config = read_input(input_path)
    config.update({key: str(item) for key, item in overrides.items() if item is not None})
    surface, vegetation = _canonical(config)
    if model is not None:
        vegetation = {"homogeneous": "hom", "raw": "row", "forest": "crown"}.get(str(model).lower(), str(model).lower())
        if vegetation not in {"bare", "hom", "row", "crown"}:
            raise ValueError("model must be bare, hom, row, or crown")
    config["surface_model"], config["vegetation_model"] = surface, vegetation
    wavelengths = _wavelengths(config, input_path)
    spectra = _spectra(config, input_path, wavelengths)
    if output_path is None:
        case_name = str(value(config, "case_name", default=f"{surface}_{vegetation}"))
        result_dir = value(config, "result_dir")
        if result_dir in (None, ""):
            result_dir = Path(value(config, "cases_dir", default="cases")) / case_name
        output_path = Path(result_dir) / "output.csv"
    else:
        output_path = Path(output_path)
    if not output_path.is_absolute():
        output_path = input_path.parent / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    snapshot = output_path.parent / "input.csv"
    if snapshot.resolve() == input_path:
        snapshot = output_path.parent / "input_used.csv"

    rows = []
    for scenario_id, scenario_label, scenario_config in _scenario_configs(config):
        rows.extend(_simulate_rows(
            scenario_config, input_path, surface, vegetation, wavelengths, spectra,
            scenario_id, scenario_label,
        ))
    # Do not expose invalid radiometric records to the GUI or plot readers.
    rows = [
        row for row in rows
        if np.isfinite(float(row["radiance"])) and np.isfinite(float(row["brightness_temperature_K"]))
    ]
    if not rows:
        raise ValueError("未找到可见的观测方向；请检查坡面朝向和观测角度")
    for output_index, row in enumerate(rows):
        row["index"] = output_index
    write_input_snapshot(config, input_path, snapshot, output_dir=output_path.parent)
    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--model", choices=("bare", "hom", "row", "crown"))
    args = parser.parse_args()
    print(run_angle(args.input, args.output, args.model))


if __name__ == "__main__":
    main()
