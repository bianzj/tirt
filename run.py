import argparse
import csv
import itertools
import shutil
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
)
from base.physicsF import inv_planck, planck, slope1
from base.row import Row
from base.soilopt import hapke_thermal_emissivity


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input.csv"
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
    sza_values = numbers(value(config, "sza", "solar_zenith"), (30.0,))
    saa_values = numbers(value(config, "saa", "solar_azimuth"), (0.0,))
    sza_default, saa_default = sza_values[0], saa_values[0]
    raw_mode = value(config, "geometry_mode", "observation_mode", "direction_mode")
    if raw_mode not in (None, ""):
        try:
            mode = int(float(raw_mode))
        except ValueError as exc:
            raise ValueError("geometry_mode 必须是 0、1、2 或 3") from exc
        if mode not in {0, 1, 2, 3}:
            raise ValueError("geometry_mode 必须是 0、1、2 或 3")
    else:
        source = str(value(config, "geometry_source", "direction_source", default="direct")).strip().lower()
        mode = 2 if source in {"file", "txt", "read"} else 0

    if mode == 2:
        path = resolve_path(input_path, value(config, "geometry_file", "directions_file", default="geometry.txt"))
        vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
        return list(zip(vza, vaa, sza, saa))

    if mode == 1:
        plane = str(value(config, "geometry_view", "principal_plane", default="both")).strip().lower()
        relative_azimuths = {
            "parallel": (0.0, 180.0),
            "solar_principal": (0.0, 180.0),
            "perpendicular": (90.0, 270.0),
            "vertical": (90.0, 270.0),
            "both": (0.0, 90.0, 180.0, 270.0),
        }.get(plane)
        if relative_azimuths is None:
            raise ValueError("geometry_view 必须是 parallel、perpendicular 或 both")
        raw_vza = value(config, "principal_vza", "principal_plane_vza")
        if raw_vza not in (None, ""):
            principal_vza = numbers(raw_vza)
        else:
            step = _float(config, "principal_vza_step", default=10.0)
            maximum = _float(config, "principal_vza_max", default=75.0)
            if step <= 0.0 or maximum < 0.0:
                raise ValueError("principal_vza_step 必须大于 0，principal_vza_max 不能小于 0")
            principal_vza = _stepped_angles(maximum, step)
        if not principal_vza or any(angle < 0.0 or angle > 90.0 for angle in principal_vza):
            raise ValueError("principal_vza 必须位于 0 到 90 度之间")
        return [
            (view_zenith, (solar_azimuth + relative_azimuth) % 360.0, solar_zenith, solar_azimuth)
            for solar_zenith in sza_values
            for solar_azimuth in saa_values
            for relative_azimuth in relative_azimuths
            for view_zenith in principal_vza
        ]

    if mode == 3:
        raw_step = _float(config, "principal_vza_step", default=10.0)
        maximum = _float(config, "principal_vza_max", default=75.0)
        if raw_step <= 0.0 or maximum < 0.0:
            raise ValueError("principal_vza_step 必须大于 0，principal_vza_max 不能小于 0")
        principal_vza = _stepped_angles(maximum, raw_step)
        include_principal = str(value(config, "include_principal", default="1")).strip().lower() not in {"0", "false", "no", "off"}
        include_hemisphere = str(value(config, "include_hemisphere", default="0")).strip().lower() not in {"0", "false", "no", "off"}
        principal = []
        if include_principal:
            principal = [
                (view_zenith, (solar_azimuth + relative_azimuth) % 360.0, solar_zenith, solar_azimuth)
                for solar_zenith in sza_values
                for solar_azimuth in saa_values
                for relative_azimuth in (0.0, 180.0, 90.0, 270.0)
                for view_zenith in principal_vza
            ]
        source = str(value(config, "additional_geometry_source", "observation_source", default="none")).strip().lower()
        additional = []
        if include_hemisphere:
            path = resolve_path(input_path, value(config, "additional_geometry_file", "observation_geometry_file", default="data/directions/hemisphere_directions.txt"))
            vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
            additional = [item for item in zip(vza, vaa, sza, saa) if item[0] <= maximum + 1.0e-9]
        if not include_hemisphere and source in {"file", "txt", "read"}:
            raw_file = value(config, "additional_geometry_file", "observation_geometry_file")
            if raw_file not in (None, ""):
                path = resolve_path(input_path, raw_file)
                vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
                additional.extend(zip(vza, vaa, sza, saa))
        if source in {"direct", "manual", "angles"}:
            vza = numbers(value(config, "observation_vza", "additional_vza"))
            vaa = numbers(value(config, "observation_vaa", "additional_vaa"))
            if len(vza) != len(vaa):
                raise ValueError("手动观测角度必须一一对应：observation_vza 和 observation_vaa 数量不同")
            additional.extend((view_zenith, view_azimuth, sza_default, saa_default)
                              for view_zenith, view_azimuth in zip(vza, vaa))
        combined = []
        seen = set()
        for item in principal + additional:
            key = tuple(round(float(value), 8) for value in item)
            if key not in seen:
                seen.add(key)
                combined.append(item)
        return combined

    vza = numbers(value(config, "vza", "observation_vza"), (0.0, 30.0, 60.0))
    vaa = numbers(value(config, "vaa", "observation_vaa"), (0.0, 90.0, 180.0))
    pair_mode = str(value(config, "geometry_pair_mode", default="cross")).strip().lower()
    if pair_mode in {"zip", "pair", "paired"}:
        size = max(len(vza), len(vaa), len(sza_values), len(saa_values))
        return list(zip(
            repeat_to(vza, size, "vza"), repeat_to(vaa, size, "vaa"),
            repeat_to(sza_values, size, "sza"), repeat_to(saa_values, size, "saa"),
        ))
    return [(view_zenith, view_azimuth, solar_zenith, solar_azimuth)
            for view_zenith in vza for view_azimuth in vaa
            for solar_zenith in sza_values for solar_azimuth in saa_values]


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
    if result.size != len(wavelengths):
        raise ValueError(f"{label} 的数目必须与 wavelengths 的数目一致: {result.size} != {len(wavelengths)}")
    return result


def _spectra(config, input_path, wavelengths):
    size = len(wavelengths)
    file_mode = {"file", "txt", "read"}
    leaf_from_file = str(value(config, "vegetation_spectrum_source", default="direct")).lower() in file_mode
    soil_from_file = str(value(config, "soil_spectrum_source", default="direct")).lower() in file_mode
    terrain_from_file = str(value(config, "terrain_spectrum_source", default="direct")).lower() in file_mode
    urban_from_file = str(value(config, "urban_spectrum_source", default="direct")).lower() in file_mode
    roof_from_file = str(value(config, "roof_spectrum_source", default="direct")).lower() in file_mode
    wall_from_file = str(value(config, "wall_spectrum_source", default="direct")).lower() in file_mode
    street_from_file = str(value(config, "street_spectrum_source", default="direct")).lower() in file_mode
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

    return {
        "leaf_e": np.asarray(leaf_e), "soil_e": np.asarray(soil_e), "terrain_e": np.asarray(terrain_e),
        "roof_e": np.asarray(roof_e), "wall_e": np.asarray(wall_e), "street_e": np.asarray(street_e),
    }


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
        raise ValueError("terrain_density * pi * terrain_radius^2 的总占据率必须小于 1")
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


def _vegetation_result(config, vegetation, wavelength, geometry, spectra):
    vza, vaa, sza, saa = geometry
    # The analytical crown and hotspot equations are singular exactly at the
    # horizon. Keep the reported geometry unchanged, but evaluate the model at
    # the nearest finite angle.
    model_vza = min(max(abs(float(vza)), 1.0e-3), 89.999)
    model_sza = min(max(abs(float(sza)), 1.0e-3), 89.999)
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    if vegetation == "bare":
        soil_emissivity = _soil_emissivity(config, model_vza, spectra["soil_e"])
        radiance = 0.5 * soil_emissivity * (planck(wavelength, soil_sun) + planck(wavelength, soil_shade))
        return float(radiance), float(inv_planck(wavelength, radiance))
    lai = _float(config, "lai", "vegetation_lai", default=1.5)
    hspot = _float(config, "hspot", "vegetation_hspot", default=0.15)
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
        rad_a = max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6)
        rad_b = max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6)
        model.set_structure(lai, hspot, _float(config, "crown_density", default=0.02), rad_b, rad_a)
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance = float(np.asarray(model.run(ifradiance=1)).reshape(-1)[0])
    if not np.isfinite(radiance):
        radiance = float(0.5 * spectra["soil_e"] * (planck(wavelength, soil_sun) + planck(wavelength, soil_shade)))
    return radiance, float(inv_planck(wavelength, radiance))


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
        return _vegetation_result(config, vegetation, wavelength, _slope_local_geometry(config, geometry), spectra)

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
    model.set_slope(slope_angle, slope_aspect)
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance_values = np.asarray(model.run(ifradiance=1)).reshape(-1)
    if not np.all(np.isfinite(radiance_values)):
        return _vegetation_result(config, vegetation, wavelength, _slope_local_geometry(config, geometry), spectra)
    radiance = float(radiance_values[0])
    return radiance, float(inv_planck(wavelength, radiance))


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
        model.set_structural_input(_terrain_shapes(config))
        terrain_vza = min(max(abs(vza), 1.0e-3), 89.999)
        model.set_angular_input(np.asarray([terrain_vza]), np.asarray([vaa]), sza, saa)
        soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
        terrain_sun = _temperature(config, "terrain_sunlit", leaf_sun - 273.15)
        terrain_shade = _temperature(config, "terrain_shaded", leaf_shade - 273.15)
        if vegetation == "bare":
            terrain_emissivity = _soil_emissivity(config, terrain_vza, spectra["terrain_e"])
            model.set_spectral_input(terrain_emissivity, spectra["roof_e"])
            model.set_thermal_input(soil_sun, soil_shade, terrain_sun, terrain_shade)
            component_temperatures = (soil_sun, soil_shade, terrain_sun, terrain_shade)
        else:
            model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
            model.set_spectral_input(spectra["terrain_e"], spectra["leaf_e"])
            model.set_thermal_input(soil_sun, soil_shade, leaf_sun, leaf_shade)
            component_temperatures = (soil_sun, soil_shade, leaf_sun, leaf_shade)
        emissivity = np.asarray(model.calculate_effective_component_emissivity(2)).reshape(-1)
        radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures))))
        return radiance, float(inv_planck(wavelength, radiance))

    from base.urban import Urban
    from base.urban_veg import Urban_Veg

    model = Urban() if vegetation == "bare" else Urban_Veg()
    model.set_structural_input(_urban_shapes(config))
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    roof_sun = _temperature(config, "roof_sunlit", 45.0)
    roof_shade = _temperature(config, "roof_shaded", 30.0)
    wall_sun = _temperature(config, "wall_sunlit", 45.0)
    wall_shade = _temperature(config, "wall_shaded", 30.0)
    street_sun = _temperature(config, "street_sunlit", 45.0)
    street_shade = _temperature(config, "street_shaded", 30.0)
    model_geometry = _urban_model_geometry(geometry)
    model_vza, model_vaa, model_sza, model_saa = model_geometry
    model.set_angular_input(np.asarray([model_vza]), np.asarray([model_vaa]), model_sza, model_saa)
    model.set_thermal_input(roof_sun, roof_shade, wall_sun, wall_shade, street_sun, street_shade)
    if vegetation == "bare":
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade)
    else:
        model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", "crown_width", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"], spectra["leaf_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade, 0.5 * (leaf_sun + leaf_shade))
    emissivity = np.asarray(model.calculate_effective_component_emissivity(2)).reshape(-1)
    radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures[:len(emissivity)]))))
    return radiance, float(inv_planck(wavelength, radiance))


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
            radiance, brightness_temperature = _surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)
            if (surface != "slope" and
                    (not np.isfinite(radiance) or not np.isfinite(brightness_temperature))):
                radiance, brightness_temperature = _vegetation_result(config, vegetation, float(wavelength), geometry, band_spectra)
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


def run(input_path=DEFAULT_INPUT, output_path=None, model=None, **overrides):
    input_path = Path(input_path).resolve()
    config = read_input(input_path)
    config.update({key: str(item) for key, item in overrides.items() if item is not None})
    surface, vegetation = _canonical(config)
    if model is not None:
        vegetation = {"homogeneous": "hom", "raw": "row", "forest": "crown"}.get(str(model).lower(), str(model).lower())
        if vegetation not in {"bare", "hom", "row", "crown"}:
            raise ValueError("model must be bare, hom, row, or crown")
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
    if snapshot.resolve() != input_path:
        shutil.copy2(input_path, snapshot)

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
    for output_index, row in enumerate(rows):
        row["index"] = output_index
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
    print(run(args.input, args.output, args.model))


if __name__ == "__main__":
    main()
