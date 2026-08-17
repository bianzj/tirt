import argparse
import csv
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
    strict_band_values,
    value,
)
from base.physicsF import inv_planck, planck
from base.row import Row


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input.csv"
OUTPUT_FIELDS = (
    "index", "surface_model", "vegetation_model", "model", "wavelength_um",
    "vza", "vaa", "sza", "saa", "raa", "emissivity_leaf", "emissivity_soil",
    "emissivity_roof", "emissivity_wall", "emissivity_street",
    "radiance", "brightness_temperature_K", "brightness_temperature_C",
)


def _float(config, *keys, default=0.0):
    return float(value(config, *keys, default=default))


def _temperature(config, name, default_c):
    kelvin = value(config, f"{name}_k")
    if kelvin not in (None, ""):
        return float(kelvin)
    return _float(config, f"{name}_c", default=default_c) + 273.15


def _canonical(config):
    surface = str(value(config, "surface_model", "surface", default="plane")).strip().lower()
    surface = {"flat": "plane", "bare": "plane", "building": "urban", "city": "urban"}.get(surface, surface)
    vegetation = str(value(config, "vegetation_model", "vegetation", default="hom")).strip().lower()
    vegetation = {"homogeneous": "hom", "raw": "row", "row_crop": "row", "forest": "crown", "bare_soil": "bare"}.get(vegetation, vegetation)
    if surface not in {"plane", "terrain", "urban"}:
        raise ValueError("surface_model 必须是 plane、terrain 或 urban")
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
            raise ValueError("geometry_mode 必须是 0、1 或 2") from exc
        if mode not in {0, 1, 2}:
            raise ValueError("geometry_mode 必须是 0、1 或 2")
    else:
        source = str(value(config, "geometry_source", "direction_source", default="direct")).strip().lower()
        mode = 2 if source in {"file", "txt", "read"} else 0

    if mode == 2:
        path = resolve_path(input_path, value(config, "geometry_file", "directions_file", default="geometry.txt"))
        vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
        return list(zip(vza, vaa, sza, saa))

    if mode == 1:
        raw_vza = value(config, "principal_vza", "principal_plane_vza")
        if raw_vza not in (None, ""):
            principal_vza = numbers(raw_vza)
        else:
            step = _float(config, "principal_vza_step", default=10.0)
            maximum = _float(config, "principal_vza_max", default=80.0)
            if step <= 0.0 or maximum < 0.0:
                raise ValueError("principal_vza_step 必须大于 0，principal_vza_max 不能小于 0")
            principal_vza = tuple(np.arange(0.0, maximum + 0.5 * step, step))
        if not principal_vza or any(angle < 0.0 or angle > 90.0 for angle in principal_vza):
            raise ValueError("principal_vza 必须位于 0 到 90 度之间")
        return [
            (view_zenith, (solar_azimuth + relative_azimuth) % 360.0, solar_zenith, solar_azimuth)
            for solar_zenith in sza_values
            for solar_azimuth in saa_values
            for relative_azimuth in (0.0, 90.0, 180.0, 270.0)
            for view_zenith in principal_vza
        ]

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


def _spectra(config, input_path, wavelengths):
    size = len(wavelengths)
    file_mode = {"file", "txt", "read"}
    leaf_from_file = str(value(config, "vegetation_spectrum_source", default="direct")).lower() in file_mode
    soil_from_file = str(value(config, "soil_spectrum_source", default="direct")).lower() in file_mode
    urban_from_file = str(value(config, "urban_spectrum_source", default="direct")).lower() in file_mode
    leaf_e = np.full(size, 0.985) if leaf_from_file else strict_band_values(config, ("emissivity_leaf", "leaf_emissivity"), size, 0.985, "emissivity_leaf")
    soil_e = np.full(size, 0.955) if soil_from_file else strict_band_values(config, ("emissivity_soil", "soil_emissivity"), size, 0.955, "emissivity_soil")
    roof_e = np.full(size, 0.95) if urban_from_file else strict_band_values(config, ("emissivity_roof", "urban_emissivity_roof"), size, 0.95, "emissivity_roof")
    wall_e = np.full(size, 0.92) if urban_from_file else strict_band_values(config, ("emissivity_wall", "urban_emissivity_wall"), size, 0.92, "emissivity_wall")
    street_e = np.full(size, 0.955) if urban_from_file else strict_band_values(config, ("emissivity_street", "urban_emissivity_street"), size, 0.955, "emissivity_street")

    if leaf_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "vegetation_spectrum_file")))
        leaf_e = interpolate(table, wavelengths, 1)

    if soil_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "soil_spectrum_file")))
        soil_e = interpolate(table, wavelengths, 1)

    if urban_from_file:
        table = read_spectrum(resolve_path(input_path, value(config, "urban_spectrum_file")))
        if table.shape[1] >= 7:
            roof_e, wall_e, street_e = (interpolate(table, wavelengths, col) for col in (4, 5, 6))
        elif table.shape[1] >= 4:
            roof_e, wall_e, street_e = (interpolate(table, wavelengths, col) for col in (1, 2, 3))
        else:
            raise ValueError("urban_spectrum_file 至少需要 wavelength roof wall street 四列")

    return {
        "leaf_e": np.asarray(leaf_e), "soil_e": np.asarray(soil_e),
        "roof_e": np.asarray(roof_e), "wall_e": np.asarray(wall_e), "street_e": np.asarray(street_e),
    }


def _vegetation_temperatures(config):
    return (
        _temperature(config, "soil_sunlit", 45.0), _temperature(config, "soil_shaded", 30.0),
        _temperature(config, "leaf_sunlit", 33.0), _temperature(config, "leaf_shaded", 30.0),
    )


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


def _vegetation_result(config, vegetation, wavelength, geometry, spectra):
    vza, vaa, sza, saa = geometry
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    if vegetation == "bare":
        radiance = 0.5 * spectra["soil_e"] * (planck(wavelength, soil_sun) + planck(wavelength, soil_shade))
        return float(radiance), float(inv_planck(wavelength, radiance))
    lai = _float(config, "lai", "vegetation_lai", default=1.5)
    hspot = _float(config, "hspot", "vegetation_hspot", default=0.15)
    if vegetation == "hom":
        model = Hom()
        model.set_structure(lai, hspot)
        model.set_angle(np.asarray([vza]), np.asarray([sza]), np.asarray([abs(vaa - saa) % 360.0]))
    elif vegetation == "row":
        model = Row()
        model.set_structure(lai, hspot, _float(config, "row_width", default=1.0), _float(config, "row_blank", default=0.1), _float(config, "row_height", default=1.0))
        model.set_angle(np.asarray([vza]), np.asarray([sza]), np.asarray([vaa]), np.asarray([saa]), np.asarray([_float(config, "row_azimuth", default=0.0)]))
    else:
        model = Crown()
        model.set_structure(lai, hspot, _float(config, "crown_density", default=0.02), _float(config, "crown_height", "crown_width", default=3.0), _float(config, "crown_radius", "crown_width", default=3.0))
        model.set_angle(np.asarray([vza]), np.asarray([sza]), np.asarray([abs(vaa - saa) % 360.0]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance = float(np.asarray(model.run(ifradiance=1)).reshape(-1)[0])
    return radiance, float(inv_planck(wavelength, radiance))


def _surface_result(config, surface, vegetation, wavelength, geometry, spectra):
    vza, vaa, sza, saa = geometry
    if surface == "plane":
        return _vegetation_result(config, vegetation, wavelength, geometry, spectra)
    if surface == "terrain":
        from base.terrain import Terrain
        from base.terrain_veg import Terrain_Veg

        model = Terrain() if vegetation == "bare" else Terrain_Veg()
        model.set_structural_input(_terrain_shapes(config))
        model.set_angular_input(np.asarray([max(abs(vza), 1.0e-3)]), np.asarray([vaa]), sza, saa)
        soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
        terrain_sun = _temperature(config, "terrain_sunlit", leaf_sun - 273.15)
        terrain_shade = _temperature(config, "terrain_shaded", leaf_shade - 273.15)
        if vegetation == "bare":
            model.set_spectral_input(spectra["soil_e"], spectra["roof_e"])
            model.set_thermal_input(soil_sun, soil_shade, terrain_sun, terrain_shade)
            component_temperatures = (soil_sun, soil_shade, terrain_sun, terrain_shade)
        else:
            model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), _float(config, "crown_radius", "crown_width", default=3.0), _float(config, "crown_height", default=3.0), _float(config, "hspot", default=0.15)])
            model.set_spectral_input(spectra["soil_e"], spectra["leaf_e"])
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
    model.set_angular_input(np.asarray([vza]), np.asarray([vaa]), sza, saa)
    model.set_thermal_input(roof_sun, roof_shade, wall_sun, wall_shade, street_sun, street_shade)
    if vegetation == "bare":
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade)
    else:
        model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), _float(config, "crown_radius", "crown_width", default=3.0), _float(config, "crown_height", default=3.0), _float(config, "hspot", default=0.15)])
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"], spectra["leaf_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade, 0.5 * (leaf_sun + leaf_shade))
    emissivity = np.asarray(model.calculate_effective_component_emissivity(2)).reshape(-1)
    radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures[:len(emissivity)]))))
    return radiance, float(inv_planck(wavelength, radiance))


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
    geometries = _geometry(config, input_path)
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
    index = 0
    for geometry in geometries:
        vza, vaa, sza, saa = geometry
        raa = abs(vaa - saa) % 360.0
        raa = min(raa, 360.0 - raa)
        for band_index, wavelength in enumerate(wavelengths):
            band_spectra = {key: values[band_index] for key, values in spectra.items()}
            radiance, brightness_temperature = _surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)
            rows.append({
                "index": index, "surface_model": surface, "vegetation_model": vegetation,
                "model": f"{surface}+{vegetation}", "wavelength_um": float(wavelength),
                "vza": float(vza), "vaa": float(vaa), "sza": float(sza), "saa": float(saa), "raa": float(raa),
                "emissivity_leaf": float(spectra["leaf_e"][band_index]), "emissivity_soil": float(spectra["soil_e"][band_index]),
                "emissivity_roof": float(spectra["roof_e"][band_index]), "emissivity_wall": float(spectra["wall_e"][band_index]), "emissivity_street": float(spectra["street_e"][band_index]),
                "radiance": radiance, "brightness_temperature_K": brightness_temperature, "brightness_temperature_C": brightness_temperature - 273.15,
            })
            index += 1
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
