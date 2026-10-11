import argparse
import csv
import itertools
import shutil
import warnings
from pathlib import Path

import numpy as np

from base.crown import Crown
from base.hom import Hom
from base.inputio import (
    canonical_models,
    geometry as _geometry,
    numbers,
    paired_geometry as _paired_geometry,
    read_input,
    repeat_to,
    spectra as _spectra,
    terrain_shapes as _terrain_shapes,
    value,
    wavelengths as _wavelengths,
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
    ("lai", "LAI"), ("hspot", "Hotspot"), ("sza", "SZA"), ("saa", "SAA"),
    ("slope_angle", "Slope angle"), ("slope_aspect", "Slope aspect"),
    ("row_width", "Row width"), ("row_blank", "Row gap"), ("row_height", "Row height"), ("row_azimuth", "Row azimuth"),
    ("crown_rad_a", "Crown rad-a"), ("crown_rad_b", "Crown rad-b"), ("crown_density", "Crown density"),
)


def _float(config, *keys, default=0.0):
    return float(value(config, *keys, default=default))


def _temperature(config, name, default_c):
    kelvin = value(config, f"{name}_k")
    if kelvin not in (None, ""):
        return float(kelvin)
    return _float(config, f"{name}_c", default=default_c) + 273.15


def _vegetation_temperatures(config):
    return (
        _temperature(config, "soil_sunlit", 45.0), _temperature(config, "soil_shaded", 30.0),
        _temperature(config, "leaf_sunlit", 33.0), _temperature(config, "leaf_shaded", 30.0),
    )


def _soil_emissivity(config, view_zenith, base_emissivity):
    enabled = str(value(config, "hapke_enable", default="0")).strip().lower()
    if enabled in {"0", "false", "no", "off"}:
        return float(base_emissivity)
    albedo = _float(config, "soil_reflectance", default=1.0 - float(base_emissivity))
    K = _float(config, "hapke_k", default=1.0)
    result = hapke_thermal_emissivity(view_zenith, K, albedo)
    return float(np.asarray(result).reshape(-1)[0])


def _shape_values(config, key, default):
    return numbers(value(config, key), (default,))


def _forest_shape(config):
    return np.asarray([
        _float(config, "lai", default=1.5),
        _float(config, "crown_density", default=0.02),
        max(_float(config, "crown_rad_a", default=3.0), 1.0e-6),
        max(_float(config, "crown_rad_b", default=3.0), 1.0e-6),
        _float(config, "hspot", default=0.15),
        max(_float(config, "crown_height", default=3.0), 0.0),
    ])


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
    lai = _float(config, "lai", default=1.5)
    if vegetation == "bare" or lai <= 0.0:
        # A flat bare surface casts no shadows: fully sunlit while the sun is above the horizon.
        soil_emissivity = _soil_emissivity(config, model_vza, spectra["soil_e"])
        soil_temperature = soil_sun if abs(float(sza)) < 90.0 else soil_shade
        radiance = soil_emissivity * planck(wavelength, soil_temperature)
        return float(radiance), float(inv_planck(wavelength, radiance))
    hspot = _float(config, "hspot", default=0.15)
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
        rad_a = max(_float(config, "crown_rad_a", default=3.0), 1.0e-6)
        rad_b = max(_float(config, "crown_rad_b", default=3.0), 1.0e-6)
        model.set_structure(lai, hspot, _float(config, "crown_density", default=0.02), rad_b, rad_a, max(_float(config, "crown_height", default=3.0), 0.0))
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance = float(np.asarray(model.run(ifradiance=1)).reshape(-1)[0])
    return radiance, float(inv_planck(wavelength, radiance))


def _slope_local_geometry(config, geometry):
    """Transform global view and solar directions into slope-local angles."""
    vza, vaa, sza, saa = geometry
    slope_angle = _float(config, "slope_angle", default=0.0)
    slope_aspect = _float(config, "slope_aspect", default=0.0)

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
    if vegetation in {"bare", "hom", "row"}:
        # A bare surface or a turbid/row canopy parallel to the slope is evaluated in slope-local geometry;
        # slope1 folds the local sun zenith into the upper hemisphere, so test for a sun behind the slope here.
        if _slope_cosine(config, geometry[2], geometry[3]) <= 0.0:
            soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
            config = dict(config, soil_sunlit_k=str(soil_shade), leaf_sunlit_k=str(leaf_shade))
        if vegetation == "row":
            slope_aspect = _float(config, "slope_aspect", default=0.0)
            config = dict(config, row_azimuth=str(_float(config, "row_azimuth", default=0.0) - slope_aspect))
        return _vegetation_result(config, vegetation, wavelength, _slope_local_geometry(config, geometry), spectra)

    from base.slope_veg import Slope_Veg

    vza, vaa, sza, saa = geometry
    soil_sun, soil_shade, leaf_sun, leaf_shade = _vegetation_temperatures(config)
    lai = _float(config, "lai", default=1.5)
    hspot = _float(config, "hspot", default=0.15)
    stand_density = _float(config, "crown_density", default=0.02)
    crown_rad_a = max(_float(config, "crown_rad_a", default=3.0), 1.0e-6)
    crown_rad_b = max(_float(config, "crown_rad_b", default=3.0), 1.0e-6)
    slope_angle = _float(config, "slope_angle", default=0.0)
    slope_aspect = _float(config, "slope_aspect", default=0.0)

    model = Slope_Veg()
    model.set_angle(np.asarray([abs(vza)]), np.asarray([sza]), np.asarray([vaa]), np.asarray([saa]))
    model.set_optical(wavelength, spectra["soil_e"], spectra["leaf_e"])
    model.set_structure(lai, hspot, stand_density, crown_rad_b, crown_rad_a, max(_float(config, "crown_height", default=3.0), 0.0))
    model.set_slope(slope_angle, slope_aspect)
    model.set_thermal(soil_sun, soil_shade, leaf_sun, leaf_shade)
    radiance_values = np.asarray(model.run(ifradiance=1)).reshape(-1)
    if not np.all(np.isfinite(radiance_values)):
        return float("nan"), float("nan")
    radiance = float(radiance_values[0])
    return radiance, float(inv_planck(wavelength, radiance))


def _slope_cosine(config, zenith, azimuth):
    """Cosine between a direction and the slope normal; negative means behind the slope."""
    slope_angle = _float(config, "slope_angle", default=0.0)
    slope_aspect = _float(config, "slope_aspect", default=0.0)
    return float(
        np.cos(np.deg2rad(slope_angle)) * np.cos(np.deg2rad(float(zenith)))
        + np.sin(np.deg2rad(slope_angle)) * np.sin(np.deg2rad(float(zenith)))
        * np.cos(np.deg2rad(float(azimuth) - slope_aspect))
    )


def _slope_backface(config, geometry):
    """Return whether the slope-corrected view zenith exceeds 89 degrees."""
    # The dot product keeps the backside information that slope1 folds into the upper hemisphere.
    cosine = _slope_cosine(config, geometry[0], geometry[1])
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
            soil_emissivity = _soil_emissivity(config, terrain_vza, spectra["soil_e"])
            model.set_spectral_input(soil_emissivity, spectra["terrain_e"])
            model.set_thermal_input(soil_sun, soil_shade, terrain_sun, terrain_shade)
            component_temperatures = (soil_sun, soil_shade, terrain_sun, terrain_shade)
        else:
            model.forestshape = _forest_shape(config)
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
        model.forestshape = _forest_shape(config)
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"], spectra["leaf_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade, 0.5 * (leaf_sun + leaf_shade))
    emissivity = np.asarray(model.calculate_effective_component_emissivity(2)).reshape(-1)
    radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures[:len(emissivity)]))))
    return radiance, float(inv_planck(wavelength, radiance))


def _scenario_configs(config):
    """Expand semicolon-separated scenario parameters into a Cartesian batch."""
    entries = []
    paired = _paired_geometry(config)
    for key, label in SCENARIO_PARAMETERS:
        if paired and key in {"sza", "saa"}:
            continue
        raw = value(config, key)
        if raw in (None, ""):
            continue
        values = tuple(float(item) for item in numbers(raw))
        if not values:
            continue
        entries.append((key, label, values))

    if not entries:
        return [(1, "default", config.copy())]

    varying = [entry for entry in entries if len(entry[2]) > 1]
    scenario_configs = []
    products = itertools.product(*(entry[2] for entry in entries))
    for scenario_id, combination in enumerate(products, 1):
        scenario_config = config.copy()
        label_parts = []
        for (key, label, values), selected in zip(entries, combination):
            scenario_config[key] = str(selected)
            if len(values) > 1:
                label_parts.append(f"{label}={selected:g}")
        scenario_label = "; ".join(label_parts) if varying else "default"
        scenario_configs.append((scenario_id, scenario_label, scenario_config))
    return scenario_configs


def _simulate_rows(config, input_path, surface, vegetation, wavelengths, spectra, scenario_id, scenario_label):
    geometries = _geometry(config, input_path)
    rows = []
    for geometry in geometries:
        vza, vaa, sza, saa = geometry
        raa = abs(vaa - saa) % 360.0
        raa = min(raa, 360.0 - raa)
        for band_index, wavelength in enumerate(wavelengths):
            band_spectra = {key: values[band_index] for key, values in spectra.items()}
            radiance, brightness_temperature = _surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)
            rows.append({
                "index": 0, "scenario_id": scenario_id, "scenario_label": scenario_label,
                "surface_model": surface, "vegetation_model": vegetation,
                "model": f"{surface}+{vegetation}", "wavelength_um": float(wavelength),
                "vza": float(vza), "vaa": float(vaa), "sza": float(sza), "saa": float(saa), "raa": float(raa),
                "emissivity_leaf": float(spectra["leaf_e"][band_index]), "emissivity_soil": float(spectra["soil_e"][band_index]),
                "emissivity_roof": float(spectra["roof_e"][band_index]), "emissivity_wall": float(spectra["wall_e"][band_index]), "emissivity_street": float(spectra["street_e"][band_index]),
                "radiance": radiance, "brightness_temperature_K": brightness_temperature, "brightness_temperature_C": brightness_temperature - 273.15,
            })
    return rows


def run(input_path=DEFAULT_INPUT, output_path=None, model=None, **overrides):
    input_path = Path(input_path).resolve()
    config = read_input(input_path)
    config.update({key.lower(): str(item) for key, item in overrides.items() if item is not None})
    if model is not None:
        config["vegetation_model"] = str(model)
    surface, vegetation = canonical_models(config)
    leaf_projection = value(config, "g")
    if leaf_projection not in (None, "") and abs(float(leaf_projection) - 0.5) > 1.0e-9:
        raise ValueError("方向性解析模型基于球形叶倾角分布，叶片投影系数 G 固定为 0.5")
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
    snapshot = (output_path.parent / "input.csv").resolve()
    # Never overwrite the project's own input templates with a case snapshot.
    if snapshot != input_path and snapshot.parent != ROOT:
        shutil.copy2(input_path, snapshot)

    rows = []
    for scenario_id, scenario_label, scenario_config in _scenario_configs(config):
        rows.extend(_simulate_rows(
            scenario_config, input_path, surface, vegetation, wavelengths, spectra,
            scenario_id, scenario_label,
        ))
    # Do not expose invalid radiometric records to the GUI or plot readers.
    valid = [
        row for row in rows
        if np.isfinite(float(row["radiance"])) and np.isfinite(float(row["brightness_temperature_K"]))
    ]
    if len(valid) < len(rows):
        warnings.warn(f"{len(rows) - len(valid)} 条非有限结果已剔除（例如坡面背面不可见的方向）", RuntimeWarning, stacklevel=2)
    rows = valid
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
