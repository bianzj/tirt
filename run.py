import argparse
import csv
import itertools
import warnings
from pathlib import Path

import numpy as np

from base import canopy_closure
from base.crown import Crown
from base.hom import Hom
from base.inputio import (
    BAND_PROPERTIES,
    is_input_snapshot,
    numbers,
    read_band_table,
    read_geometry,
    read_input,
    repeat_to,
    resolve_path,
    value,
    write_input_snapshot,
)
from base.optical import InvisibleDirection
from base.physicsF import inv_planck, planck, slope1
from base.row import Row
from base.soilopt import hapke_directional_emissivity


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input.csv"


class BackfaceOmittedWarning(UserWarning):
    """Some requested directions look at the back of the slope and were left out."""

OUTPUT_FIELDS = (
    "index", "scenario_id", "scenario_label", "surface_model", "vegetation_model", "model", "wavelength_um",
    "vza", "vaa", "sza", "saa", "raa", "emissivity_leaf", "emissivity_soil",
    "emissivity_roof", "emissivity_wall", "emissivity_street",
    "radiance", "brightness_temperature_K", "brightness_temperature_C", "reflectance_brf",
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
    ("crown_rad_b", ("crown_rad_b",), "Crown rad-b"),
    ("crown_height", ("crown_height",), "Crown height"),
    ("crown_density", ("crown_density",), "Crown density"),
)


def _float(config, *keys, default=0.0):
    raw = value(config, *keys, default=default)
    try:
        return float(raw)
    except ValueError as exc:
        raise ValueError(f"{keys[0]} must be a number, got {raw!r}") from exc


def _temperature(config, name, default_c):
    kelvin = value(config, f"{name}_k")
    if kelvin not in (None, ""):
        return _float(config, f"{name}_k")
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
        raise ValueError("surface_model must be plane, slope, terrain or urban")
    if vegetation not in {"bare", "hom", "row", "crown"}:
        raise ValueError("vegetation_model must be bare, hom, row or crown")
    return surface, vegetation


def _geometry(config, input_path):
    geometries = _geometry_unchecked(config, input_path)
    if not geometries:
        raise ValueError("At least one valid viewing direction is required")
    for geometry in geometries:
        if not np.all(np.isfinite(geometry)):
            raise ValueError("Viewing and solar angles must be finite numbers")
        vza, vaa, sza, saa = geometry
        if not (0.0 <= vza <= 90.0 and 0.0 <= sza <= 90.0):
            raise ValueError("vza and sza must be between 0 and 90 degrees")
        if not (0.0 <= vaa <= 360.0 and 0.0 <= saa <= 360.0):
            raise ValueError("vaa and saa must be between 0 and 360 degrees")
    return geometries


def _geometry_unchecked(config, input_path):
    sza_values = numbers(value(config, "sza", "solar_zenith"), (30.0,), key="sza")
    saa_values = numbers(value(config, "saa", "solar_azimuth"), (0.0,), key="saa")
    sza_default, saa_default = sza_values[0], saa_values[0]
    raw_mode = value(config, "geometry_mode", "observation_mode", "direction_mode")
    if raw_mode not in (None, ""):
        try:
            parsed_mode = float(raw_mode)
            if not np.isfinite(parsed_mode) or not parsed_mode.is_integer():
                raise ValueError("geometry_mode must be an integer")
            mode = int(parsed_mode)
        except ValueError as exc:
            raise ValueError("geometry_mode must be 0, 1, 2 or 3") from exc
        if mode not in {0, 1, 2, 3}:
            raise ValueError("geometry_mode must be 0, 1, 2 or 3")
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
            raise ValueError("geometry_view must be parallel, perpendicular or both")
        raw_vza = value(config, "principal_vza", "principal_plane_vza")
        if raw_vza not in (None, ""):
            principal_vza = numbers(raw_vza, key="principal_vza")
        else:
            step = _float(config, "principal_vza_step", default=10.0)
            maximum = _float(config, "principal_vza_max", default=75.0)
            if not np.isfinite(step) or not np.isfinite(maximum) or step <= 0.0 or not 0.0 <= maximum <= 90.0:
                raise ValueError("principal_vza_step must be greater than 0 and principal_vza_max must be between 0 and 90 degrees")
            principal_vza = _stepped_angles(maximum, step)
        if not principal_vza or any(angle < 0.0 or angle > 90.0 for angle in principal_vza):
            raise ValueError("principal_vza must be between 0 and 90 degrees")
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
        if not np.isfinite(raw_step) or not np.isfinite(maximum) or raw_step <= 0.0 or not 0.0 <= maximum <= 90.0:
            raise ValueError("principal_vza_step must be greater than 0 and principal_vza_max must be between 0 and 90 degrees")
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
            path = resolve_path(input_path, value(config, "hemisphere_geometry_file", default="data/directions/hemisphere_directions.txt"))
            vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
            additional = [item for item in zip(vza, vaa, sza, saa) if item[0] <= maximum + 1.0e-9]
        if source in {"file", "txt", "read"}:
            raw_file = value(config, "additional_geometry_file", "observation_geometry_file")
            if raw_file not in (None, ""):
                path = resolve_path(input_path, raw_file)
                vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
                additional.extend(zip(vza, vaa, sza, saa))
        if source in {"direct", "manual", "angles"}:
            vza = numbers(value(config, "observation_vza", "additional_vza"), key="observation_vza")
            vaa = numbers(value(config, "observation_vaa", "additional_vaa"), key="observation_vaa")
            if len(vza) != len(vaa):
                raise ValueError("Manual viewing angles must be paired: observation_vza and observation_vaa have different lengths")
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

    vza = numbers(value(config, "vza", "observation_vza"), (0.0, 30.0, 60.0), key="vza")
    vaa = numbers(value(config, "vaa", "observation_vaa"), (0.0, 90.0, 180.0), key="vaa")
    pair_mode = str(value(config, "geometry_pair_mode", default="cross")).strip().lower()
    if pair_mode not in {"cross", "zip", "pair", "paired"}:
        raise ValueError("geometry_pair_mode must be cross or zip")
    if pair_mode in {"zip", "pair", "paired"}:
        size = max(len(vza), len(vaa), len(sza_values), len(saa_values))
        return list(zip(
            repeat_to(vza, size, "vza"), repeat_to(vaa, size, "vaa"),
            repeat_to(sza_values, size, "sza"), repeat_to(saa_values, size, "saa"),
        ))
    return [(view_zenith, view_azimuth, solar_zenith, solar_azimuth)
            for view_zenith in vza for view_azimuth in vaa
            for solar_zenith in sza_values for solar_azimuth in saa_values]


def _required_properties(surface, vegetation):
    required = set() if vegetation == "bare" else {"leaf_reflectance", "leaf_transmittance"}
    required |= {"street_reflectance", "wall_reflectance", "roof_reflectance"} if surface == "urban" else {"soil_reflectance"}
    return required


def _band_table(config, input_path):
    """Wavelengths and per-band optical properties shared by the shortwave and thermal solvers."""
    source = str(value(config, "spectral_source", default="direct")).strip().lower()
    if source in {"file", "txt", "read"}:
        path = resolve_path(input_path, value(config, "spectrals_file", default="data/spectrals/optical_properties.txt"))
        rows = read_band_table(path, 1 + len(BAND_PROPERTIES))
    else:
        required = _required_properties(config["surface_model"], config["vegetation_model"])
        rows = [numbers(value(config, "wavelengths"), key="wavelengths")]
        for key in BAND_PROPERTIES:
            raw = value(config, key)
            if raw in (None, "") and key in required:
                raise ValueError(f"{key} is required: one value per wavelength")
            # Rows unused by this scene may be omitted; they never enter the result.
            rows.append(numbers(raw, key=key) if raw not in (None, "") else (0.0,) * len(rows[0]))
    wavelengths = np.asarray(rows[0], dtype=float)
    if wavelengths.size == 0 or np.any(wavelengths <= 0.0):
        raise ValueError("wavelengths must contain positive values")
    if np.any(wavelengths < 0.4):
        raise ValueError("Shortwave bands must be at least 0.4 um")
    table = {"wavelengths": wavelengths}
    for key, values in zip(BAND_PROPERTIES, rows[1:]):
        values = np.asarray(values, dtype=float)
        if values.size != wavelengths.size:
            raise ValueError(f"{key} must have one value per wavelength: {values.size} != {wavelengths.size}")
        if np.any(values < 0.0) or np.any(values > 1.0):
            raise ValueError(f"{key} must be between 0 and 1")
        table[key] = values
    if np.any(table["leaf_reflectance"] + table["leaf_transmittance"] > 1.0 + 1.0e-12):
        raise ValueError("leaf_reflectance + leaf_transmittance must not exceed 1")
    return table


def _spectra(table):
    """Kirchhoff emissivities: leaf 1 - rho - tau, opaque surfaces 1 - rho; terrain facets share the soil material."""
    soil_e = 1.0 - table["soil_reflectance"]
    return {
        "leaf_e": np.maximum(1.0 - table["leaf_reflectance"] - table["leaf_transmittance"], 0.0),
        "soil_e": soil_e, "terrain_e": soil_e,
        "roof_e": 1.0 - table["roof_reflectance"], "wall_e": 1.0 - table["wall_reflectance"],
        "street_e": 1.0 - table["street_reflectance"],
    }


def _vegetation_temperatures(config):
    return (
        _temperature(config, "soil_sunlit", 45.0), _temperature(config, "soil_shaded", 30.0),
        _temperature(config, "leaf_sunlit", 33.0), _temperature(config, "leaf_shaded", 30.0),
    )


def _soil_emissivity(config, view_zenith, base_emissivity):
    """Directional soil emissivity; Hapke keeps the hemispherical value equal to the input."""
    enabled = str(value(config, "hapke_enable", "soil_hapke_enable", default="0")).strip().lower()
    if enabled in {"0", "false", "no", "off"}:
        return float(base_emissivity)
    K = _float(config, "hapke_k", "soil_hapke_k", "hapke_K", default=1.0)
    result = hapke_directional_emissivity(min(abs(float(view_zenith)), 89.999), K, float(base_emissivity))
    return float(np.asarray(result).reshape(-1)[0])


def _shape_values(config, key, default):
    return numbers(value(config, key), (default,), key=key)


def _terrain_shapes(config):
    heights = _shape_values(config, "terrain_height", 20.0)
    radii = _shape_values(config, "terrain_radius", 10.0)
    densities = _shape_values(config, "terrain_density", 0.001)
    size = max(len(heights), len(radii), len(densities))
    shapes = np.column_stack([repeat_to(heights, size, "terrain_height"), repeat_to(radii, size, "terrain_radius"), repeat_to(densities, size, "terrain_density")])
    occupied = float(np.sum(shapes[:, 2] * np.pi * shapes[:, 1] ** 2))
    if occupied >= 1.0:
        raise ValueError("Total coverage terrain_density * pi * terrain_radius^2 must be less than 1")
    return shapes


def _urban_shapes(config):
    lengths = _shape_values(config, "building_length", 10.0)
    widths = _shape_values(config, "building_width", 10.0)
    heights = _shape_values(config, "building_height", 20.0)
    densities = _shape_values(config, "building_density", 0.0015)
    azimuths = _shape_values(config, "building_azimuth", 0.0)
    size = max(map(len, (lengths, widths, heights, densities, azimuths)))
    azimuths = repeat_to(azimuths, size, "building_azimuth") % 360.0
    shapes = np.column_stack([
        repeat_to(lengths, size, "building_length"), repeat_to(widths, size, "building_width"),
        repeat_to(heights, size, "building_height"), repeat_to(densities, size, "building_density"),
        azimuths, (azimuths + 90.0) % 360.0,
    ])
    if np.any(shapes[:, :4] < 0.0):
        raise ValueError("building_length, building_width, building_height and building_density must not be negative")
    footprint = float(np.sum(shapes[:, 0] * shapes[:, 1] * shapes[:, 3]))
    if footprint >= 1.0:
        raise ValueError(f"Building footprint fraction (sum of building_density * building_length * building_width) must be less than 1; got {footprint:.4g}")
    return shapes


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
        _float(config, "row_azimuth", default=0.0),
        str(value(config, "thermal_row_geometry", default="parametric")).strip().lower())


def _crown_semiaxes(config):
    """Horizontal radius crown_rad_a and vertical semiaxis crown_rad_b (m)."""
    return (max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6),
            max(_float(config, "crown_rad_b", default=3.0), 1.0e-6))


def _forest_options(config):
    """Crown geometry, inner/outer hotspot widths and terrain transport, defaulting to the validated branch."""
    mode = str(value(config, "thermal_crown_geometry", default="parametric")).strip().lower()
    inner = value(config, "thermal_hotspot_inner")
    outer = value(config, "thermal_hotspot_outer")
    transport = str(value(config, "forest_transport", "thermal_forest_transport", default="legacy")).strip().lower()
    return mode, inner, outer, transport


def _canopy_structure(config, vegetation):
    """(kind, shape, row azimuth, crown base) for canopy_closure."""
    if vegetation == "row":
        shape = (_float(config, "row_width", default=0.5), _float(config, "row_blank", default=0.5),
                 _float(config, "row_height", default=1.0))
        return "row", shape, _float(config, "row_azimuth", default=0.0), None
    if vegetation == "crown":
        rad_a, rad_b = _crown_semiaxes(config)
        return "crown", (_float(config, "crown_density", default=0.02), rad_b, rad_a), 0.0, \
            _float(config, "crown_height", default=3.0)
    return "hom", None, 0.0, None


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
    level = canopy_closure.closure(config)
    if level != "legacy" and not (vegetation == "row" and _row_closure(config)):
        kind, shape, row_azimuth, base = _canopy_structure(config, vegetation)
        weights = canopy_closure.thermal_weights(
            kind, lai, hspot, shape, spectra["leaf_e"], spectra["soil_e"], model_vza, vaa, model_sza, saa,
            row_azimuth, base, _soil_emissivity(config, model_vza, spectra["soil_e"]), level)
        radiance = float(np.sum(weights * planck(wavelength, np.asarray([soil_sun, soil_shade, leaf_sun, leaf_shade]))))
        return _checked_radiance(radiance, wavelength, f"plane+{vegetation}", geometry)
    if vegetation == "hom":
        model = Hom()
        model.set_structure(lai, hspot)
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    elif vegetation == "row":
        model = Row()
        model.set_geometry(str(value(config, "thermal_row_geometry", default="parametric")).strip().lower())
        model.set_closure(_row_closure(config))
        model.set_structure(lai, hspot, _float(config, "row_width", default=0.5), _float(config, "row_blank", default=0.5), _float(config, "row_height", default=1.0))
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([vaa]), np.asarray([saa]), np.asarray([_float(config, "row_azimuth", default=0.0)]))
    else:
        model = Crown()
        model.set_hotspot_geometry(*_forest_options(config)[:3])
        model.crown_base = _float(config, "crown_height", default=3.0)
        rad_a = max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6)
        rad_b = max(_float(config, "crown_rad_b", default=3.0), 1.0e-6)
        model.set_structure(lai, hspot, _float(config, "crown_density", default=0.02), rad_b, rad_a)
        model.set_angle(np.asarray([model_vza]), np.asarray([model_sza]), np.asarray([abs(vaa - saa) % 360.0]))
    model.set_optical(wavelength, _soil_emissivity(config, model_vza, spectra["soil_e"]), spectra["leaf_e"])
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
    standard = vegetation != "bare" and canopy_closure.closure(config) != "legacy"
    if standard or (vegetation == "crown" and _forest_options(config)[0] != "legacy"):
        from base.solar import slope_direction
        slope = _float(config, "slope_angle", "terrain_slope", default=0.0)
        aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)
        local_view, local_vaa, _ = slope_direction(geometry[0], geometry[1], slope, aspect)
        local_sun, local_saa, lit = slope_direction(geometry[2], geometry[3], slope, aspect)
        local_config = dict(config)
        if not lit:
            local_config["soil_sunlit_k"] = _temperature(config, "soil_shaded", 30.0)
            local_config["leaf_sunlit_k"] = _temperature(config, "leaf_shaded", 30.0)
        if vegetation == "row":
            local_config["row_azimuth"] = slope_direction(90.0, _float(config, "row_azimuth", default=0.0),
                                                          slope, aspect)[1]
        return _vegetation_result(local_config, vegetation, wavelength,
                                  (local_view, local_vaa, local_sun, local_saa), spectra)
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
    crown_rad_b = max(_float(config, "slope_crown_rad_b", "crown_rad_b", default=3.0), 1.0e-6)
    slope_angle = _float(config, "slope_angle", "terrain_slope", default=0.0)
    slope_aspect = _float(config, "slope_aspect", "terrain_aspect", default=0.0)

    model = Slope_Veg()
    model.set_angle(np.asarray([abs(vza)]), np.asarray([sza]), np.asarray([vaa]), np.asarray([saa]))
    local_vza = _slope_local_geometry(config, geometry)[0]
    model.set_optical(wavelength, _soil_emissivity(config, local_vza, spectra["soil_e"]), spectra["leaf_e"])
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


def _row_closure(config):
    """Unified plane-row closure for thermal and shortwave bands, or None for the standard model."""
    raw = value(config, "row_closure")
    closure = None if raw in (None, "") else str(raw).strip().lower()
    if closure not in (None, "legacy", "conditioned", "reference"):
        raise ValueError("row_closure must be empty, legacy, conditioned or reference")
    return closure


def _check_row_closure(config, surface, vegetation):
    closure = _row_closure(config)
    if vegetation == "row" and surface != "plane" and closure is not None:
        raise ValueError(f"row_closure={closure} currently supports plane rows only.")


def _surface_result(config, surface, vegetation, wavelength, geometry, spectra):
    if (vegetation == "crown" and surface == "urban" and
            _forest_options(config)[0] not in ("legacy", "clumping_dual", "parametric")):
        raise ValueError("Urban thermal crown geometry supports legacy, clumping_dual or parametric.")
    _check_row_closure(config, surface, vegetation)
    vza, vaa, sza, saa = geometry
    if surface == "plane":
        return _vegetation_result(config, vegetation, wavelength, geometry, spectra)
    if surface == "slope":
        if _slope_backface(config, geometry):
            return float("nan"), float("nan")
        return _slope_result(config, vegetation, wavelength, geometry, spectra)
    if surface == "terrain":
        from base.terrain import Terrain
        from base.terrain_forest import TerrainForest

        shapes = _terrain_shapes(config)
        if np.all(shapes[:, 2] == 0.0):
            return _vegetation_result(config, vegetation, wavelength, geometry, spectra)
        if vegetation == "bare":
            model = Terrain()
        else:
            model = TerrainForest()
            if canopy_closure.closure(config) != "legacy":
                if vegetation == "crown" and _forest_options(config)[3] == "reciprocal_exchange":
                    raise ValueError("forest_transport=reciprocal_exchange requires canopy_closure=legacy.")
                canopy_closure.install(model, canopy_closure.closure(config))
            if vegetation == "crown":
                mode, inner, outer, transport = _forest_options(config)
                if transport == "reciprocal_exchange" and outer is None:
                    # Same R/hcr outer width as the plane and urban whole-crown geometry.
                    outer = _crown_semiaxes(config)[0] / _crown_semiaxes(config)[1]
                model.set_hotspot_geometry(mode, inner, outer)
                model.set_radiative_transport(transport)
        model.set_structural_input(shapes)
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
            model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
            crown_base = _float(config, "crown_height", default=3.0) if vegetation == "crown" else 0.0
            model.forestshape = np.append(model.forestshape, crown_base)
            _set_canopy_geometry(model, config, vegetation)
            model.set_spectral_input(_soil_emissivity(config, terrain_vza, spectra["soil_e"]), spectra["leaf_e"])
            model.set_thermal_input(soil_sun, soil_shade, leaf_sun, leaf_shade)
            component_temperatures = (soil_sun, soil_shade, leaf_sun, leaf_shade)
        emissivity = _checked_components(
            model.calculate_effective_component_emissivity(2), f"terrain+{vegetation}", geometry, wavelength)
        radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures))))
        return _checked_radiance(radiance, wavelength, f"terrain+{vegetation}", geometry)

    from base.urban import Urban
    from base.urban_veg import Urban_Veg

    shapes = _urban_shapes(config)
    if np.all(shapes[:, 3] == 0.0):
        # Empty buildings reduce exactly to a plane canopy over street material.
        plane_config = config.copy()
        plane_config["soil_sunlit_k"] = str(_temperature(config, "street_sunlit", 45.0))
        plane_config["soil_shaded_k"] = str(_temperature(config, "street_shaded", 30.0))
        # Streets are not soil: the Hapke soil model never applies to urban surfaces.
        plane_config["hapke_enable"] = "0"
        plane_spectra = dict(spectra, soil_e=spectra["street_e"])
        return _vegetation_result(plane_config, vegetation, wavelength, geometry, plane_spectra)
    model = Urban() if vegetation == "bare" else Urban_Veg()
    model.set_structural_input(shapes)
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
        model.forestshape = np.asarray([_float(config, "lai", default=1.5), _float(config, "crown_density", default=0.02), max(_float(config, "crown_rad_a", "crown_radius", default=3.0), 1.0e-6), max(_float(config, "crown_rad_b", default=3.0), 1.0e-6), _float(config, "hspot", default=0.15)])
        _set_canopy_geometry(model, config, vegetation)
        if vegetation == "crown":
            model.set_hotspot_geometry(*_forest_options(config)[:3])
            model.crown_base = _float(config, "crown_height", default=3.0)
        if canopy_closure.closure(config) != "legacy":
            canopy_closure.install(model, canopy_closure.closure(config))
        model.set_spectral_input(spectra["street_e"], spectra["wall_e"], spectra["roof_e"], spectra["leaf_e"])
        component_temperatures = (wall_sun, wall_shade, street_sun, street_shade, roof_sun, roof_shade, leaf_sun, leaf_shade)
    # The standard and upgrade closures share the street + four-wall radiosity of the optical urban model.
    default_transport = "legacy" if canopy_closure.closure(config) == "legacy" else "radiosity"
    transport = str(value(config, "urban_transport", default=default_transport) or default_transport).strip().lower()
    if transport == "radiosity":
        from base.urban_radiosity import effective_component_emissivity
        weights = effective_component_emissivity(model)
    elif transport == "legacy":
        weights = model.calculate_effective_component_emissivity(2)
    else:
        raise ValueError("urban_transport must be legacy or radiosity")
    emissivity = _checked_components(weights, f"urban+{vegetation}", geometry, wavelength)
    radiance = float(np.sum(emissivity * planck(wavelength, np.asarray(component_temperatures[:len(emissivity)]))))
    return _checked_radiance(radiance, wavelength, f"urban+{vegetation}", geometry)


def _scenario_configs(config):
    """Expand semicolon-separated scenario parameters into a Cartesian batch."""
    entries = []
    for key, aliases, label in SCENARIO_PARAMETERS:
        raw = value(config, *aliases)
        if raw in (None, ""):
            continue
        values = numbers(raw, key=key)
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


SKY_REFERENCE_TEMPERATURE = 300.0
COMPONENT_TEMPERATURE_KEYS = tuple(f"{material}_{lighting}_k"
                                   for material in ("soil", "leaf", "terrain", "roof", "wall", "street")
                                   for lighting in ("sunlit", "shaded"))


def _sky_temperature(config):
    """Equivalent sky brightness temperature (K) of the downwelling radiance, or None for a black sky."""
    kelvin = value(config, "sky_temperature_k")
    celsius = value(config, "sky_temperature_c")
    if kelvin in (None, "") and celsius in (None, ""):
        return None
    result = _float(config, "sky_temperature_k") if kelvin not in (None, "") else _float(config, "sky_temperature_c") + 273.15
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError("sky_temperature must be a positive finite temperature")
    return result


def _with_sky_reflection(config, surface, vegetation, wavelength, geometry, spectra, radiance, sky_temperature):
    """Add (1 - directional effective emissivity) x B(T_sky) for an isotropic sky (Kirchhoff reciprocity)."""
    isothermal = dict(config, **{key: str(SKY_REFERENCE_TEMPERATURE) for key in COMPONENT_TEMPERATURE_KEYS})
    reference, _ = _surface_result(isothermal, surface, vegetation, wavelength, geometry, spectra)
    emissivity = float(reference) / float(planck(wavelength, SKY_REFERENCE_TEMPERATURE))
    total = float(radiance) + max(0.0, 1.0 - emissivity) * float(planck(wavelength, sky_temperature))
    return _checked_radiance(total, wavelength, f"{surface}+{vegetation}", geometry)


def _optical_band(config, surface, vegetation, wavelength, geometry, band):
    """Shortwave BRF from the explicit VNIR interface, using this band's optical properties."""
    from base.optical import reflectance

    settings = dict(config)
    diffuse_fraction = _float(config, "shortwave_diffuse_fraction", default=0.2)
    if diffuse_fraction < 0.0:
        diffuse_fraction = 0.2
    if not np.isfinite(diffuse_fraction) or not 0.0 <= diffuse_fraction <= 1.0:
        raise ValueError("shortwave_diffuse_fraction must be between 0 and 1 for the directional optical model")
    for key in ("lai", "hspot", "surface_type", "distribution", "row_width", "row_blank",
                "row_height", "row_azimuth", "crown_rad_a", "crown_rad_b", "crown_height",
                "crown_density", "stand_density", "height", "width", "slope_angle",
                "slope_aspect", "terrain_type"):
        if settings.get(key) not in (None, ""):
            settings[key] = _float(settings, key)
    if settings.get("shortwave_sunlit_surface") not in (None, ""):
        settings["shortwave_sunlit_surface"] = str(settings["shortwave_sunlit_surface"]).lower() not in ("0", "false", "no")
    settings.setdefault("surface_type", 0.0)
    if surface == "terrain" and settings.get("terrain_type") in (None, ""):
        settings["terrain_type"] = 2.0
    settings.setdefault("distribution", 2.0 if vegetation in ("row", "raw") else 0.0)
    settings.setdefault("stand_density", float(settings.get("crown_density", 0.02) or 0.02))
    settings.setdefault("height", 2.0 * float(settings.get("crown_rad_b", 3.0) or 3.0))
    settings.setdefault("width", 2.0 * float(settings.get("crown_rad_a", 3.0) or 3.0))
    if canopy_closure.closure(config) != "legacy" and vegetation in ("crown", "forest"):
        # Same crown as the thermal path: semiaxes crown_rad_b / crown_rad_a, density and base height.
        rad_a, rad_b = _crown_semiaxes(config)
        settings.update(height=rad_b, width=rad_a, stand_density=_float(config, "crown_density", default=0.02),
                        crown_bottom=_float(config, "crown_height", default=3.0))
    if vegetation == "bare":
        settings["lai"] = 0.0
    else:
        settings["lai"] = _float(config, "lai", "vegetation_lai", default=1.5)
    if surface == "urban":
        shapes = _urban_shapes(config)
        if len(shapes) > 1:
            raise ValueError("The shortwave urban model supports a single building class; give one value per building_* parameter")
        settings.update(zip(("urban_building_length", "urban_building_width", "urban_building_height",
                             "urban_density", "urban_building_azimuth"), map(float, shapes[0, :5])))
    settings.update(
        surface_model=surface,
        vegetation_model=vegetation,
        wavelengths=[float(wavelength) * 1000.0],
        leaf_reflectance=band["leaf_reflectance"],
        leaf_transmittance=band["leaf_transmittance"],
        soil_reflectance=band["soil_reflectance"],
        leaf_reflectance_spectrum=[band["leaf_reflectance"]],
        leaf_transmittance_spectrum=[band["leaf_transmittance"]],
        soil_reflectance_spectrum=[band["soil_reflectance"]],
        urban_reflectance_roof_spectrum=[band["roof_reflectance"]],
        urban_reflectance_wall_spectrum=[band["wall_reflectance"]],
        urban_reflectance_street_spectrum=[band["street_reflectance"]],
        shortwave_diffuse_fraction=diffuse_fraction,
        hapke_enable=0,
    )
    if surface == "terrain" and vegetation in ("crown", "forest"):
        method = str(value(config, "optical_terrain_method",
                           default="reciprocal_exchange")).strip().lower()
        if method not in ("reciprocal_exchange", "facet_exchange",
                          "finite_crown_ports", "effective_slab"):
            raise ValueError("Unknown optical_terrain_method: " + method)
        settings["optical_terrain_method"] = method
        settings["forest_transport"] = (
            "reciprocal_exchange" if method == "reciprocal_exchange" else "legacy")
        if method in ("reciprocal_exchange", "finite_crown_ports"):
            settings.setdefault("optical_crown_geometry", "clumping_dual")
        if method == "reciprocal_exchange":
            settings.setdefault("optical_terrain_source_depth", "effective_forest")
            settings.setdefault("optical_hotspot_outer", 1.0)
            settings.setdefault("crown_bottom",
                                float(settings.get("crown_height", 3.0) or 3.0))
    vza, vaa, sza, saa = geometry
    if surface == "terrain" and vegetation not in ("crown", "forest"):
        # The public VNIR wrapper restricts compound terrain to forests, but the
        # underlying directional model supports the other covers through the
        # compound-terrain component fractions.
        from base.radiation import directional_reflectance
        return float(directional_reflectance(
            settings, sza, saa, vza, vaa, band["leaf_reflectance"], band["leaf_transmittance"],
            band["soil_reflectance"], use_hapke=False, diffuse_fraction=diffuse_fraction))
    # Single slopes (every cover) use base.optical: local-plane BRF re-referred to horizontal irradiance.
    return float(reflectance(settings, float(wavelength) * 1000.0, vza, vaa, sza, saa))


def _simulate_rows(config, input_path, surface, vegetation, table, spectra, scenario_id, scenario_label):
    _check_row_closure(config, surface, vegetation)
    geometries = _geometry(config, input_path)
    rows = []
    index = 0
    for geometry in geometries:
        vza, vaa, sza, saa = geometry
        raa = abs(vaa - saa) % 360.0
        raa = min(raa, 360.0 - raa)
        for band_index, wavelength in enumerate(table["wavelengths"]):
            band_spectra = {key: values[band_index] for key, values in spectra.items()}
            band = {key: float(table[key][band_index]) for key in BAND_PROPERTIES}
            radiance, brightness_temperature = float("nan"), float("nan")
            reflectance_brf = float("nan")
            try:
                if float(wavelength) < 2.5:
                    try:
                        reflectance_brf = _optical_band(config, surface, vegetation, float(wavelength), geometry, band)
                    except InvisibleDirection:
                        # Directions outside the optical adapter's validity domain stay missing.
                        reflectance_brf = float("nan")
                else:
                    radiance, brightness_temperature = _surface_result(config, surface, vegetation, float(wavelength), geometry, band_spectra)
                    sky_temperature = _sky_temperature(config)
                    if sky_temperature is not None and np.isfinite(radiance):
                        radiance, brightness_temperature = _with_sky_reflection(
                            config, surface, vegetation, float(wavelength), geometry, band_spectra, radiance, sky_temperature)
            except FloatingPointError as exc:
                raise FloatingPointError(
                    f"{surface}+{vegetation}, wavelength={wavelength}, geometry={geometry}: {exc}") from exc
            rows.append({
                "index": index, "scenario_id": scenario_id, "scenario_label": scenario_label,
                "surface_model": surface, "vegetation_model": vegetation,
                "model": f"{surface}+{vegetation}", "wavelength_um": float(wavelength),
                "vza": float(vza), "vaa": float(vaa), "sza": float(sza), "saa": float(saa), "raa": float(raa),
                "emissivity_leaf": round(float(spectra["leaf_e"][band_index]), 12), "emissivity_soil": round(float(spectra["soil_e"][band_index]), 12),
                "emissivity_roof": round(float(spectra["roof_e"][band_index]), 12), "emissivity_wall": round(float(spectra["wall_e"][band_index]), 12), "emissivity_street": round(float(spectra["street_e"][band_index]), 12),
                "radiance": radiance, "brightness_temperature_K": brightness_temperature, "brightness_temperature_C": brightness_temperature - 273.15,
                "reflectance_brf": reflectance_brf,
            })
            index += 1
    return rows


def run(input_path=DEFAULT_INPUT, output_path=None, model=None, report=None, **overrides):
    input_path = Path(input_path).resolve()
    config = read_input(input_path)
    config.update({key.lower(): str(item) for key, item in overrides.items() if item is not None})
    surface, vegetation = _canonical(config)
    if model is not None:
        vegetation = {"homogeneous": "hom", "raw": "row", "forest": "crown"}.get(str(model).lower(), str(model).lower())
        if vegetation not in {"bare", "hom", "row", "crown"}:
            raise ValueError("model must be bare, hom, row, or crown")
    config["surface_model"], config["vegetation_model"] = surface, vegetation
    table = _band_table(config, input_path)
    spectra = _spectra(table)
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
    if snapshot.resolve() == input_path or (snapshot.is_file() and not is_input_snapshot(snapshot)):
        snapshot = output_path.parent / "input_used.csv"

    rows = []
    for scenario_id, scenario_label, scenario_config in _scenario_configs(config):
        rows.extend(_simulate_rows(
            scenario_config, input_path, surface, vegetation, table, spectra,
            scenario_id, scenario_label,
        ))
    # Do not expose invalid radiometric records to the GUI or plot readers.
    requested = len(rows)
    rows = [
        row for row in rows
        if np.isfinite(float(row["brightness_temperature_K"])) or np.isfinite(float(row["reflectance_brf"]))
    ]
    if not rows:
        raise ValueError("No visible viewing directions; check the slope orientation and viewing angles")
    if len(rows) < requested:
        warnings.warn(
            f"{requested - len(rows)} observation records have no valid values (slope backface or optical direction not visible) and were omitted from the output",
            BackfaceOmittedWarning)
    if report is not None:
        report["omitted_backface"] = requested - len(rows)
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
    print(run(args.input, args.output, args.model))


if __name__ == "__main__":
    main()
