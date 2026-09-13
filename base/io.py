import re
import csv
from pathlib import Path

import numpy as np

from .leafopt import fluspect
from .soilopt import bsm, soil_tir_emissivity
from .structs import DefinedStructure


LAND_COVER_TYPES = {
    "water": 0,
    "evergreen_needleleaf_forest": 1,
    "evergreen_broadleaf_forest": 2,
    "deciduous_needleleaf_forest": 3,
    "deciduous_broadleaf_forest": 4,
    "mixed_forest": 5,
    "closed_shrublands": 6,
    "open_shrublands": 7,
    "woody_savannas": 8,
    "savannas": 9,
    "grasslands": 10,
    "permanent_wetlands": 11,
    "permanent_watlands": 11,
    "croplands": 12,
    "urban_and_builtup": 13,
    "urban": 13,
    "builtup": 13,
    "cropland_vegetation_mosaic": 14,
    "snow_and_ice": 15,
    "barren_sparsely_vegetated": 16,
    "unclassified": 0,
}

DISTRIBUTIONS = {
    "homogeneous": 1,
    "hom": 1,
    "row": 2,
    "crown": 3,
    "urban": 4,
}

SURFACE_MODELS = {
    "plane": "plane",
    "bare_soil": "bare_soil",
    "soil": "bare_soil",
    "homogeneous": "homogeneous",
    "hom": "homogeneous",
    "forest": "forest",
    "crown": "forest",
    "row": "row",
    "ridge_row": "row",
    "slope": "terrain",
    "terrain": "terrain",
    "urban": "urban",
    "building": "urban",
    "city": "urban",
}

VEGETATION_MODELS = {
    "bare": "bare",
    "bare_soil": "bare",
    "soil": "bare",
    "barren": "bare",
    "hom": "hom",
    "homogeneous": "hom",
    "row": "row",
    "raw": "row",
    "crown": "crown",
    "forest": "crown",
}


def read_input(path):
    """Read grouped input.csv and return a flat key/value configuration."""
    values = {}
    sections = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as file:
        for row in csv.reader(file):
            fields = [item.strip() for item in row]
            if not fields or not fields[0] or fields[0].startswith("#"):
                continue
            if fields[0].startswith("[") and fields[0].endswith("]"):
                continue
            if len(fields) >= 3:
                section, key, value = fields[0].lower(), fields[1].lower(), fields[2]
                sections.setdefault(section, {})[key] = value
                values[key] = value
            elif len(fields) >= 2 and fields[1]:
                values[fields[0].lower()] = fields[1]
    values["_sections"] = sections
    return values


def _input_number(overrides, key, default):
    value = overrides.get(key, default)
    return float(value) if value not in (None, "") else default


def _surface_model(value):
    text = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if text not in SURFACE_MODELS:
        raise ValueError("surface_model must be plane, terrain or urban")
    return SURFACE_MODELS[text]


def _vegetation_model(value):
    text = str(value or "hom").strip().lower().replace(" ", "_").replace("-", "_")
    if text not in VEGETATION_MODELS:
        raise ValueError("vegetation_model must be bare, hom, row or crown")
    return VEGETATION_MODELS[text]


def _input_value(overrides, *keys):
    for key in keys:
        if key in overrides and overrides[key] != "":
            return overrides[key]
    return None


def _number_list(value):
    if value in (None, ""):
        return []
    values = []
    for item in re.split(r"[;,\s]+", str(value)):
        if item.strip():
            values.append(float(item))
    return values


def _mean_override(value, default):
    values = _number_list(value)
    return float(np.mean(values)) if values else float(default)


def _thermal_band_values(overrides, key, wavelengths, default):
    values = _number_list(_input_value(overrides, key))
    if not values:
        return [float(default)] * len(wavelengths)
    if len(values) == 1:
        return values * len(wavelengths)
    if len(values) != len(wavelengths):
        raise ValueError(f"{key} must contain one value or one value per thermal wavelength")
    return values


def _land_cover_code(value):
    if value is None:
        return None
    text = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if text in LAND_COVER_TYPES:
        return LAND_COVER_TYPES[text]
    return int(round(float(value)))


def _distribution_code(value):
    if value is None:
        return None
    text = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if text in DISTRIBUTIONS:
        return DISTRIBUTIONS[text]
    return int(round(float(value)))


def _apply_input_overrides(canopy, leaf, leafbio, soil_optical, soil_thermal, aero, meta, overrides):
    model = _surface_model(_input_value(overrides, "surface_model", "surface_option") or "plane")
    vegetation_model = _vegetation_model(
        _input_value(overrides, "vegetation_model", "vegetation_type") or "hom"
    )
    if model == "urban":
        # Urban geometry controls the aerodynamic branch. Keep the selected
        # vegetation distribution so the street canopy remains active.
        canopy[9] = 13.0
        canopy[10] = {
            "bare": 0.0,
            "hom": 1.0,
            "row": 2.0,
            "crown": 3.0,
        }[vegetation_model]
    elif vegetation_model == "bare":
        canopy[9], canopy[10] = 0.0, 0.0
    elif vegetation_model == "crown":
        canopy[9], canopy[10] = 1.0, 3.0
    elif vegetation_model == "row":
        canopy[9], canopy[10] = 0.0, 2.0
    else:
        canopy[9], canopy[10] = 0.0, 1.0

    canopy_keys = {
        0: ("lai",), 1: ("stand_density", "density"), 2: ("height",),
        3: ("width", "canopy_width"), 4: ("gleaf", "g"), 5: ("lidfa",),
        6: ("lidfb", "lidb"), 7: ("hspot",), 8: ("leaf_width",),
    }
    for index, keys in canopy_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            canopy[index] = float(value)
    if vegetation_model == "bare":
        canopy[0] = 0.0

    leaf_keys = {
        0: ("cab",), 1: ("cw",), 2: ("cdm",), 3: ("cs",), 4: ("n",),
        5: ("leaf_reflectance", "leaf_reflectance_tir"),
        6: ("leaf_transmittance", "leaf_transmittance_tir"),
    }
    for index, keys in leaf_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            if index in (5, 6):
                value = _mean_override(value, leaf[index])
            leaf[index] = float(value)

    leafbio_keys = {
        0: ("vcmax",), 1: ("ballberry_m", "m"), 2: ("leaf_type",),
        8: ("rdparam",), 9: ("tyear",), 10: ("beta",), 11: ("knpqs",),
        12: ("qls",), 13: ("kv",), 14: ("stressfactor",),
    }
    for index, keys in leafbio_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            leafbio[index] = float(value)
    for index in range(5):
        value = _input_value(overrides, f"tparam{index + 1}")
        if value is not None:
            leafbio[3 + index] = float(value)

    soil_optical_keys = {
        0: ("soil_brightness", "bsm_brightness"),
        1: ("soil_latitude", "bsm_latitude"),
        2: ("soil_longitude", "bsm_longitude"),
        3: ("soil_reflectance", "soil_reflectance_tir"),
    }
    for index, keys in soil_optical_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            if index == 3:
                value = _mean_override(value, soil_optical[index])
            soil_optical[index] = float(value)

    soil_keys = {
        0: ("soil_rss", "rss"), 1: ("soil_cs", "cs_soil"),
        2: ("soil_rhos", "rhos"), 3: ("soil_lambda", "soil_lambdas"),
        4: ("soil_rbs", "rbs"), 5: ("soil_moisture", "smc"),
        6: ("soil_temperature_c", "soil_temperature"),
        7: ("satwater", "soil_satwater"),
    }
    for index, keys in soil_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            soil_thermal[index] = float(value)

    aero_keys = {
        0: ("aero_zo",), 1: ("aero_d",), 2: ("aero_rbc",),
        3: ("aero_cr",), 4: ("aero_cd",), 5: ("aero_cd1",),
        6: ("aero_psicor",), 7: ("aero_cssoil",), 8: ("aero_rwc",),
        9: ("aero_rbs",),
    }
    for index, keys in aero_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None:
            aero[index] = float(value)

    meta_keys = {
        0: ("z",), 2: ("ea",), 3: ("ca",), 4: ("oa",),
        5: ("soil_temperature_c",), 8: ("satwater",),
    }
    for index, keys in meta_keys.items():
        value = _input_value(overrides, *keys)
        if value is not None and index < len(meta):
            meta[index] = float(value)

    return model



def _numeric_rows(path):
    rows = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        values = []
        for token in re.split(r"\s+", line.strip()):
            try:
                values.append(float(token))
            except ValueError:
                pass
        if values:
            rows.append(values)
    return rows


def _read_optical_coefficients(path):
    rows = []
    with Path(path).open(encoding="utf-8-sig", newline="") as file:
        for row in csv.reader(file):
            try:
                values = [float(value) for value in row]
            except ValueError:
                continue
            if len(values) >= 18:
                rows.append(values)
    return rows


def _read_measured_spectrum(path, wavelengths, value_column=1):
    """Read wavelength/value text data and interpolate it to the model grid."""
    samples = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        fields = re.split(r"[,;\s]+", text)
        try:
            wavelength = float(fields[0])
            value = float(fields[value_column])
        except (IndexError, ValueError):
            continue
        if wavelength <= 20.0:
            wavelength *= 1000.0
        samples.append((wavelength, value))
    if len(samples) < 2:
        raise ValueError(f"测量光谱文件至少需要两行波长和数值: {path}")
    samples.sort()
    x = np.asarray([item[0] for item in samples], dtype=float)
    y = np.asarray([item[1] for item in samples], dtype=float)
    unique, indices = np.unique(x, return_index=True)
    return np.interp(np.asarray(wavelengths[:2001], dtype=float), unique, y[indices]).tolist()


def _read_measured_columns(path, wavelengths, columns):
    """Read one wavelength column plus several spectral value columns."""
    samples = {column: [] for column in columns}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        fields = re.split(r"[,;\s]+", text)
        try:
            wavelength = float(fields[0])
        except (IndexError, ValueError):
            continue
        if wavelength <= 20.0:
            wavelength *= 1000.0
        for column in columns:
            try:
                samples[column].append((wavelength, float(fields[column])))
            except (IndexError, ValueError):
                pass
    result = {}
    grid = np.asarray(wavelengths, dtype=float)
    for column, values in samples.items():
        if len(values) < 2:
            raise ValueError(f"urban spectrum column {column} needs at least two samples: {path}")
        values.sort()
        x = np.asarray([item[0] for item in values], dtype=float)
        y = np.asarray([item[1] for item in values], dtype=float)
        unique, indices = np.unique(x, return_index=True)
        result[column] = np.interp(grid, unique, y[indices]).tolist()
    return result


def read_defined(path=None, overrides=None):
    """Build model settings from input.csv; the legacy defined file is not required."""
    overrides = dict(overrides or {})
    input_base = Path(str(overrides.get("_input_base", ""))) if overrides.get("_input_base") else None
    if input_base is None:
        input_base = Path(path).resolve().parent if path else Path.cwd()

    canopy = [
        _input_number(overrides, "lai", 2.0),
        _input_number(overrides, "stand_density", 1.0),
        _input_number(overrides, "height", 3.0),
        _input_number(overrides, "width", 3.0),
        _input_number(overrides, "gleaf", 0.5),
        _input_number(overrides, "lidfa", -0.35),
        _input_number(overrides, "lidfb", -0.15),
        _input_number(overrides, "hspot", 0.2),
        _input_number(overrides, "leaf_width", 0.2), 0.0, 1.0,
    ]
    leaf = [
        _input_number(overrides, "cab", 58.0),
        _input_number(overrides, "cw", 0.013),
        _input_number(overrides, "cdm", 0.0036),
        _input_number(overrides, "cs", 0.0),
        _input_number(overrides, "n", 1.86),
        _mean_override(_input_value(overrides, "leaf_reflectance", "leaf_reflectance_tir"), 0.015),
        _mean_override(_input_value(overrides, "leaf_transmittance", "leaf_transmittance_tir"), 0.01),
    ]
    leafbio = [
        _input_number(overrides, "vcmax", 80.0), _input_number(overrides, "ballberry_m", 9.0),
        _input_number(overrides, "leaf_type", 0.0), _input_number(overrides, "tparam1", 0.2),
        _input_number(overrides, "tparam2", 0.3), _input_number(overrides, "tparam3", 281.0),
        _input_number(overrides, "tparam4", 313.0), _input_number(overrides, "tparam5", 328.0),
        _input_number(overrides, "rdparam", 0.015), _input_number(overrides, "tyear", 20.0),
        _input_number(overrides, "beta", 0.4), _input_number(overrides, "knpqs", 0.0),
        _input_number(overrides, "qls", 1.0), _input_number(overrides, "kv", 0.6396),
        _input_number(overrides, "stressfactor", 1.0),
    ]
    soil_optical = [
        _input_number(overrides, "bsm_brightness", 0.25),
        _input_number(overrides, "bsm_latitude", 25.0),
        _input_number(overrides, "bsm_longitude", 45.0),
        _mean_override(_input_value(overrides, "soil_reflectance", "soil_reflectance_tir"), 0.045),
    ]
    soil_thermal = [
        _input_number(overrides, "soil_rss", 200.0), _input_number(overrides, "soil_cs", 1180.0),
        _input_number(overrides, "soil_rhos", 1100.0), _input_number(overrides, "soil_lambda", 1.55),
        _input_number(overrides, "soil_rbs", 10.0), _input_number(overrides, "soil_moisture", 0.25),
        _input_number(overrides, "soil_temperature_c", 20.0), _input_number(overrides, "satwater", 0.45),
    ]
    aero = [
        _input_number(overrides, "aero_zo", 10.0), _input_number(overrides, "aero_d", 1.34),
        _input_number(overrides, "aero_rbc", 10.0), _input_number(overrides, "aero_cr", 0.35),
        _input_number(overrides, "aero_cd", 0.2), _input_number(overrides, "aero_cd1", 20.6),
        _input_number(overrides, "aero_psicor", 0.2), _input_number(overrides, "aero_cssoil", 0.01),
        _input_number(overrides, "aero_rwc", 0.0), _input_number(overrides, "aero_rbs", 10.0),
    ]
    meta = [
        _input_number(overrides, "z", 15.0), _input_number(overrides, "soil_moisture", 0.25),
        _input_number(overrides, "ea", 15.0), _input_number(overrides, "ca", 380.0),
        _input_number(overrides, "oa", 209.0), _input_number(overrides, "soil_temperature_c", 20.0),
        _input_number(overrides, "satwater", 0.45), _input_number(overrides, "dtime", 1800.0),
    ]
    surface_model = _apply_input_overrides(
        canopy, leaf, leafbio, soil_optical, soil_thermal, aero, meta, overrides,
    )
    vegetation_model = _vegetation_model(
        _input_value(overrides, "vegetation_model", "vegetation_type") or "hom"
    )
    if vegetation_model == "crown":
        canopy[1] = _input_number(overrides, "crown_density", canopy[1])
        canopy[2] = _input_number(overrides, "crown_height", canopy[2])
        canopy[3] = _input_number(overrides, "crown_width", canopy[3])
    if surface_model == "urban" and vegetation_model == "bare":
        canopy[0] = 0.0

    def input_path(*keys, default):
        value = _input_value(overrides, *keys) or default
        file_path = Path(str(value))
        return file_path if file_path.is_absolute() else input_base / file_path

    wavelengths = [float(400 + i) for i in range(2001)]
    wavelengths += [float(2500 + 100 * i) for i in range(126)]
    wavelengths += [float(16000 + 1000 * i) for i in range(35)]

    def spectrum(file_path):
        values = []
        for line in Path(file_path).read_text(encoding="utf-8-sig").splitlines():
            try:
                values.append(float(line.split()[0]))
            except (IndexError, ValueError):
                continue
        return values

    esun = spectrum(input_path("esun_file", "solar_direct_file", default="data/defined/Esun_.dat"))
    esky = spectrum(input_path("esky_file", "solar_sky_file", default="data/defined/Esky_.dat"))
    if len(esun) == len(wavelengths) and len(esky) == len(wavelengths):
        def normalize(values, start, end):
            total = sum(
                0.5 * (values[i] + values[i + 1]) * (wavelengths[i + 1] - wavelengths[i])
                for i in range(start, end - 1)
            )
            return [v / (total * 0.001) for v in values[start:end]]

        fesun = normalize(esun, 0, 2001) + normalize(esun, 2001, len(wavelengths))
        fesky = normalize(esky, 0, 2001) + normalize(esky, 2001, len(wavelengths))
    else:
        fesun = fesky = []

    coeff = _read_optical_coefficients(
        input_path("optical_coefficients_file", "optipar_file", default="data/defined/optipar_new.csv")
    )
    if len(coeff) >= 2001:
        leaf_refl, leaf_trans = fluspect(coeff, leaf)
        soil_refl = bsm(coeff, soil_optical, soil_thermal[5])
    else:
        leaf_refl = leaf_trans = soil_refl = []

    vegetation_source = str(
        _input_value(overrides or {}, "vegetation_spectrum_source", "vegetation_spectral_source") or "model"
    ).strip().lower()
    soil_source = str(
        _input_value(overrides or {}, "soil_spectrum_source", "soil_spectral_source") or "model"
    ).strip().lower()
    if vegetation_source in {"measured", "measurement", "file", "txt"}:
        path = _input_value(overrides or {}, "vegetation_spectrum_file", "vegetation_spectral_file", "leaf_spectrum_file")
        if not path:
            raise ValueError("vegetation_spectrum_source=measured 时必须设置 vegetation_spectrum_file")
        spectrum_path = Path(str(path))
        if not spectrum_path.is_absolute():
            spectrum_path = input_base / spectrum_path
        leaf_refl = _read_measured_spectrum(spectrum_path, wavelengths, 1)
        leaf_trans = _read_measured_spectrum(spectrum_path, wavelengths, 2) if len(
            re.split(r"[,;\s]+", next(
                line for line in spectrum_path.read_text(encoding="utf-8-sig").splitlines()
                if line.strip() and not line.strip().startswith("#")
            ).strip())
        ) > 2 else [0.0] * 2001
    if soil_source in {"measured", "measurement", "file", "txt"}:
        path = _input_value(overrides or {}, "soil_spectrum_file", "soil_spectral_file")
        if not path:
            raise ValueError("soil_spectrum_source=measured 时必须设置 soil_spectrum_file")
        spectrum_path = Path(str(path))
        if not spectrum_path.is_absolute():
            spectrum_path = input_base / spectrum_path
        soil_refl = _read_measured_spectrum(spectrum_path, wavelengths, 1)

    urban_source = str(
        _input_value(overrides or {}, "urban_spectrum_source", "urban_spectral_source") or "model"
    ).strip().lower()
    urban_refl_roof_value = _input_number(overrides or {}, "urban_reflectance_roof", 0.05)
    urban_refl_wall_value = _input_number(overrides or {}, "urban_reflectance_wall", 0.08)
    urban_refl_street_value = _input_number(overrides or {}, "urban_reflectance_street", soil_optical[3])
    urban_refl_roof = [urban_refl_roof_value] * len(wavelengths)
    urban_refl_wall = [urban_refl_wall_value] * len(wavelengths)
    urban_refl_street = [urban_refl_street_value] * len(wavelengths)
    urban_abs_roof_value = _input_number(overrides or {}, "urban_absorption_roof", 1.0 - urban_refl_roof_value)
    urban_abs_wall_value = _input_number(overrides or {}, "urban_absorption_wall", 1.0 - urban_refl_wall_value)
    urban_abs_street_value = _input_number(overrides or {}, "urban_absorption_street", 1.0 - urban_refl_street_value)
    urban_eps_roof = [_input_number(overrides or {}, "urban_emissivity_roof", 0.95)] * len(wavelengths)
    urban_eps_wall = [_input_number(overrides or {}, "urban_emissivity_wall", 0.92)] * len(wavelengths)
    urban_eps_street = [_input_number(overrides or {}, "urban_emissivity_street", 0.955)] * len(wavelengths)
    roof_file = _input_value(overrides or {}, "urban_roof_spectrum_file", "roof_spectrum_file")
    if roof_file:
        roof_path = Path(str(roof_file))
        if not roof_path.is_absolute():
            roof_path = input_base / roof_path
        urban_refl_roof = _read_measured_spectrum(roof_path, wavelengths, 1) + [urban_refl_roof_value] * (len(wavelengths) - 2001)
    wall_file = _input_value(overrides or {}, "urban_wall_spectrum_file", "wall_spectrum_file")
    if wall_file:
        wall_path = Path(str(wall_file))
        if not wall_path.is_absolute():
            wall_path = input_base / wall_path
        urban_refl_wall = _read_measured_spectrum(wall_path, wavelengths, 1) + [urban_refl_wall_value] * (len(wavelengths) - 2001)
    street_file = _input_value(overrides or {}, "urban_street_spectrum_file", "street_spectrum_file")
    if street_file:
        street_path = Path(str(street_file))
        if not street_path.is_absolute():
            street_path = input_base / street_path
        urban_refl_street = _read_measured_spectrum(street_path, wavelengths, 1) + [urban_refl_street_value] * (len(wavelengths) - 2001)

    if urban_source in {"measured", "measurement", "file", "txt"}:
        urban_file = _input_value(overrides or {}, "urban_spectrum_file", "urban_spectral_file")
        if not urban_file:
            raise ValueError("urban_spectrum_source=measured requires urban_spectrum_file")
        urban_path = Path(str(urban_file))
        if not urban_path.is_absolute():
            urban_path = input_base / urban_path
        urban_columns = _read_measured_columns(urban_path, wavelengths, (1, 2, 3, 4, 5, 6))
        urban_refl_roof = urban_columns[1]
        urban_refl_wall = urban_columns[2]
        urban_refl_street = urban_columns[3]
        urban_eps_roof = urban_columns[4]
        urban_eps_wall = urban_columns[5]
        urban_eps_street = urban_columns[6]

    observed_wavelengths = _number_list(_input_value(overrides or {}, "observation_wavelengths", "wavelengths"))
    tir_wavelengths = [wavelength for wavelength in observed_wavelengths if wavelength >= 2.5]
    if not tir_wavelengths:
        tir_wavelengths = [10.5]
    leaf_tir_reflectance = _thermal_band_values(overrides or {}, "leaf_reflectance", tir_wavelengths, leaf[5])
    leaf_tir_transmittance = _thermal_band_values(overrides or {}, "leaf_transmittance", tir_wavelengths, leaf[6])
    soil_tir_reflectance = _thermal_band_values(overrides or {}, "soil_reflectance", tir_wavelengths, soil_optical[3])

    raw_surface = str(_input_value(overrides or {}, "surface_model", "surface_option") or "plane").strip().lower().replace("-", "_")
    terrain_type = 1 if raw_surface == "slope" else int(round(_input_number(overrides or {}, "terrain_type", 1.0)))
    terrain_slope = _input_number(overrides or {}, "terrain_slope", 0.0)
    terrain_aspect = _input_number(overrides or {}, "terrain_aspect", 0.0)
    terrain_radius = _input_number(overrides or {}, "terrain_radius", 100.0)
    terrain_height = _input_number(overrides or {}, "terrain_height", 20.0)
    terrain_density = _input_number(overrides or {}, "terrain_density", 1.0)
    soil_water_enabled = int(round(_input_number(overrides or {}, "soil_tir_water_effect", 0.0)))
    soil_emissivity_value = soil_tir_emissivity(
        soil_thermal[5],
        enabled=soil_water_enabled,
        dry_moisture=_input_number(overrides or {}, "soil_tir_moisture_dry", 0.05),
        wet_moisture=_input_number(overrides or {}, "soil_tir_moisture_wet", soil_thermal[7]),
        dry_emissivity=_input_number(overrides or {}, "soil_tir_emissivity_dry", 1.0 - soil_optical[3]),
        wet_emissivity=_input_number(overrides or {}, "soil_tir_emissivity_wet", 0.99),
    )
    legacy_hapke = max(
        _input_number(overrides or {}, "hapke_optical_enable", 0.0),
        _input_number(overrides or {}, "hapke_thermal_enable", 0.0),
    )

    return DefinedStructure.from_mapping({
        "surface_type": int(round(canopy[9])),
        "distribution": int(round(canopy[10])),
        "surface_model": surface_model or "plane",
        "surface_option": surface_model or "plane",
        "vegetation_model": vegetation_model,
        "lai": canopy[0], "stand_density": canopy[1],
        "height": canopy[2], "width": canopy[3],
        "gleaf": canopy[4], "lidfa": canopy[5],
        "lidfb": canopy[6], "hspot": canopy[7],
        "leaf_width": canopy[8],
        "leaf_reflectance": leaf[5],
        "leaf_transmittance": leaf[6],
        "leaf_emissivity": 1.0 - leaf[5] - leaf[6],
        "vcmax": leafbio[0], "ballberry_m": leafbio[1],
        "leaf_type": leafbio[2], "kv": leafbio[4],
        "rdparam": leafbio[8], "tparam": leafbio[3:8],
        "tyear": leafbio[9], "beta": leafbio[10],
        "knpqs": leafbio[11], "qls": leafbio[12],
        "stressfactor": leafbio[14],
        "stomatal_method": str(_input_value(overrides or {}, "stomatal_method") or "farquhar").strip().lower(),
        "stomatal_g0": _input_number(overrides or {}, "stomatal_g0", 0.01),
        "soil_reflectance": 1.0 - soil_emissivity_value if soil_water_enabled == 1 else soil_optical[3],
        "soil_emissivity": soil_emissivity_value,
        "hapke_enable": int(round(_input_number(overrides or {}, "hapke_enable", legacy_hapke))),
        "hapke_K": _input_number(overrides or {}, "hapke_k", 1.0),
        "urban_spectrum_source": urban_source,
        "urban_spectrum_file": str(_input_value(overrides or {}, "urban_spectrum_file", "urban_spectral_file") or ""),
        "urban_reflectance_roof_spectrum": urban_refl_roof,
        "urban_reflectance_wall_spectrum": urban_refl_wall,
        "urban_reflectance_street_spectrum": urban_refl_street,
        "urban_absorption_roof": urban_abs_roof_value,
        "urban_absorption_wall": urban_abs_wall_value,
        "urban_absorption_street": urban_abs_street_value,
        "urban_emissivity_roof_spectrum": urban_eps_roof,
        "urban_emissivity_wall_spectrum": urban_eps_wall,
        "urban_emissivity_street_spectrum": urban_eps_street,
        "soil_tir_water_effect": soil_water_enabled,
        "soil_tir_emissivity_dry": _input_number(overrides or {}, "soil_tir_emissivity_dry", 1.0 - soil_optical[3]),
        "soil_tir_emissivity_wet": _input_number(overrides or {}, "soil_tir_emissivity_wet", 0.99),
        "soil_tir_moisture_dry": _input_number(overrides or {}, "soil_tir_moisture_dry", 0.05),
        "soil_tir_moisture_wet": _input_number(overrides or {}, "soil_tir_moisture_wet", soil_thermal[7]),
        "urban_rss": _input_number(overrides or {}, "urban_rss", 3662.8562),
        "urban_cs": _input_number(overrides or {}, "urban_cs", 1000000.0),
        "urban_rhos": _input_number(overrides or {}, "urban_rhos", 1800.0),
        "urban_lambda": _input_number(overrides or {}, "urban_lambda", 1.5),
        "urban_rbs": _input_number(overrides or {}, "urban_rbs", 10.0),
        "urban_temperature_c": _input_number(overrides or {}, "urban_temperature_c", 20.0),
        "urban_depth_step": _input_number(overrides or {}, "urban_depth_step", 0.05),
        "urban_wall_thickness": max(_input_number(overrides or {}, "urban_wall_thickness", 0.5), 1.0e-3),
        "urban_wall_layers": max(1, int(round(_input_number(overrides or {}, "urban_wall_layers", 3)))),
        "urban_internal_temperature_c": _input_number(overrides or {}, "urban_internal_temperature_c", 24.0),
        "urban_wall_diffuse_factor": min(max(_input_number(overrides or {}, "urban_wall_diffuse_factor", 0.5), 0.0), 1.0),
        "urban_emissivity_roof": _input_number(overrides or {}, "urban_emissivity_roof", 0.95),
        "urban_emissivity_wall": _input_number(overrides or {}, "urban_emissivity_wall", 0.92),
        "urban_emissivity_street": _input_number(overrides or {}, "urban_emissivity_street", 0.955),
        "urban_street_radiation_factor": _input_number(overrides or {}, "urban_street_radiation_factor", 1.0),
        "soil_rss": _input_number(overrides or {}, "soil_rss", soil_thermal[0]), "soil_cs": soil_thermal[1],
        "soil_rbs": soil_thermal[4],
        "soil_rhos": soil_thermal[2], "soil_lambda": soil_thermal[3],
        "soil_heat_method": str(_input_value(overrides or {}, "soil_heat_method") or "surface").strip().lower(),
        "urban_heat_method": str(_input_value(overrides or {}, "urban_heat_method") or "depth").strip().lower(),
        "temperature_method": str(_input_value(overrides or {}, "temperature_method") or "robust").strip().lower(),
        "soil_heat_fraction": _input_number(overrides or {}, "soil_heat_fraction", 0.25),
        "urban_heat_fraction": _input_number(overrides or {}, "urban_heat_fraction", 0.4),
        "soil_depth_step": _input_number(overrides or {}, "soil_depth_step", 0.05),
        "soil_temperature_method": str(_input_value(overrides or {}, "temperature_method") or "robust").strip().lower(),
        "urban_temperature_method": str(_input_value(overrides or {}, "temperature_method") or "robust").strip().lower(),
        "soil_moisture": min(max(float(soil_thermal[5]), 0.0), 0.45), "soil_temperature_c": meta[5],
        "satwater": 0.45,
        "z": meta[0], "ea": meta[2], "ca": meta[3], "oa": meta[4],
        "shortwave_diffuse_fraction": min(max(_input_number(overrides or {}, "shortwave_diffuse_fraction", -1.0), -1.0), 1.0),
        "aero_zo": aero[0], "aero_d": aero[1], "aero_rbc": aero[2],
        "aero_cr": aero[3], "aero_cd": aero[4], "aero_cd1": aero[5],
        "aero_psicor": aero[6], "aero_cssoil": aero[7],
        "aero_rwc": aero[8], "aero_rbs": aero[9],
        "slope": terrain_slope,
        "aspect": terrain_aspect,
        "terrain_type": terrain_type,
        "terrain_slope": terrain_slope,
        "terrain_aspect": terrain_aspect,
        "terrain_radius": terrain_radius,
        "terrain_height": terrain_height,
        "terrain_density": terrain_density,
        "row_width": _input_number(overrides or {}, "row_width", _input_number(overrides or {}, "ridge_width", canopy[3])),
        "row_blank": _input_number(overrides or {}, "row_blank", _input_number(overrides or {}, "ridge_blank", canopy[3])),
        "row_height": _input_number(overrides or {}, "row_height", _input_number(overrides or {}, "ridge_height", canopy[2])),
        "row_azimuth": _input_number(overrides or {}, "row_azimuth", _input_number(overrides or {}, "ridge_azimuth", 0.0)),
        "crown_radius": _input_number(overrides or {}, "crown_radius", 1.5),
        "crown_height": _input_number(overrides or {}, "crown_height", canopy[2]),
        "crown_width": _input_number(overrides or {}, "crown_width", canopy[3]),
        "crown_density": _input_number(overrides or {}, "crown_density", canopy[1]),
        "urban_density": _input_number(overrides or {}, "building_density", 0.0015),
        "urban_building_length": _input_number(overrides or {}, "building_length", 10.0),
        "urban_building_width": _input_number(overrides or {}, "building_width", 10.0),
        "urban_building_height": _input_number(overrides or {}, "building_height", 20.0),
        "urban_vegetation_fraction": min(max(_input_number(overrides or {}, "urban_vegetation_fraction", 0.30), 0.0), 1.0),
        "urban_roughness_length": max(_input_number(overrides or {}, "urban_roughness_length", 0.15), 0.01),
        "urban_canyon_roughness_ratio": max(_input_number(overrides or {}, "urban_canyon_roughness_ratio", 0.05), 0.005),
        "urban_canyon_min_ratio": max(_input_number(overrides or {}, "urban_canyon_min_ratio", 1.25), 1.0),
        "urban_leaf_wind_decay": max(_input_number(overrides or {}, "urban_leaf_wind_decay", 0.5), 0.0),
        "urban_leaf_boundary_factor": max(_input_number(overrides or {}, "urban_leaf_boundary_factor", 70.0), 0.0),
        "wavelengths": wavelengths, "fesun": fesun, "fesky": fesky,
        "leaf_reflectance_spectrum": leaf_refl,
        "leaf_transmittance_spectrum": leaf_trans,
        "soil_reflectance_spectrum": soil_refl,
        "tir_wavelengths": tir_wavelengths,
        "leaf_reflectance_tir_spectrum": leaf_tir_reflectance,
        "leaf_transmittance_tir_spectrum": leaf_tir_transmittance,
        "soil_reflectance_tir_spectrum": soil_tir_reflectance,
    })


def read_meteo(path):
    rows = []
    for line_number, line in enumerate(
        Path(path).read_text(encoding="utf-8-sig").splitlines(), 1
    ):
        values = []
        for token in re.split(r"\s+", line.strip()):
            if token:
                try:
                    values.append(float(token))
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number} 存在非数值字段") from exc
        if values:
            rows.append(values)
    if len(rows) < 2:
        raise ValueError("meteo.txt 数据不足")

    expected = int(round(rows[0][0]))
    data = []
    for row in rows[1:]:
        if len(row) < 7:
            raise ValueError("meteo.txt 每行至少需要 7 列")
        data.append({
            "time": row[0], "air_temperature_c": row[1],
            "vapor_pressure": row[2], "pressure": row[3],
            "wind_speed": max(row[4], 0.2),
            "shortwave": row[5], "longwave": row[6],
        })
    if expected != len(data):
        raise ValueError(f"meteo.txt 声明 {expected} 行，实际读取 {len(data)} 行")
    return data


def read_directions(path):
    """Read observation view directions as vza vaa pairs in degrees."""
    directions = []
    for line_number, line in enumerate(
        Path(path).read_text(encoding="utf-8-sig").splitlines(), 1
    ):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        fields = re.split(r"[\s,;]+", text)
        try:
            vza, vaa = float(fields[0]), float(fields[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number} 观测方向必须包含 vza vaa 两个数值") from exc
        if not 0.0 <= vza <= 180.0:
            raise ValueError(f"{path}:{line_number} vza 必须在 0 到 180 度之间")
        if not 0.0 <= vaa <= 360.0:
            raise ValueError(f"{path}:{line_number} vaa 必须在 0 到 360 度之间")
        directions.append((vza, vaa))
    if not directions:
        raise ValueError(f"{path} 没有有效的观测方向")
    return tuple(directions)


def read_spectrals(path):
    """Read observation wavelengths in micrometres, one value per line."""
    wavelengths = []
    for line_number, line in enumerate(
        Path(path).read_text(encoding="utf-8-sig").splitlines(), 1
    ):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        try:
            wavelength = float(re.split(r"[\s,;]+", text)[0])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number} 波段必须是微米数值") from exc
        if wavelength <= 0.0:
            raise ValueError(f"{path}:{line_number} 波段必须大于 0")
        wavelengths.append(wavelength)
    if not wavelengths:
        raise ValueError(f"{path} 没有有效的观测波段")
    return tuple(wavelengths)
