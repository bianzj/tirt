import csv
import re
from functools import lru_cache
from pathlib import Path

import numpy as np

from .inputio import FILE_MODES, canonical_models, number, read_spectrum, spectra, terrain_shapes, value, wavelengths
from .leafopt import fluspect
from .soilopt import bsm, soil_tir_emissivity
from .structs import DefinedStructure


# Model grid of the shortwave spectra: 400-2400 nm by 1 nm, then 2.5-15 um and 16-50 um.
SPECTRAL_GRID = [float(400 + i) for i in range(2001)] + [float(2500 + 100 * i) for i in range(126)] + [float(16000 + 1000 * i) for i in range(35)]


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


def _fluspect(coeff, leaf):
    reflectance, transmittance = _fluspect_cached(tuple(map(tuple, coeff[:2001])), tuple(leaf[:5]))
    return list(reflectance), list(transmittance)


@lru_cache(maxsize=32)
def _fluspect_cached(coeff, leaf):
    reflectance, transmittance = fluspect(coeff, leaf)
    return tuple(reflectance), tuple(transmittance)


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


def _column_count(path):
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        text = line.strip()
        if text and not text.startswith("#"):
            return len(re.split(r"[,;\s]+", text))
    return 0


def _aerodynamic_height(config, vegetation):
    """Vegetation top used by the aerodynamic resistances."""
    if vegetation == "row":
        return number(config, "row_height", default=1.0)
    if vegetation == "crown":
        return number(config, "crown_height", default=3.0) + 2.0 * number(config, "crown_rad_b", default=3.0)
    return number(config, "canopy_height", default=3.0)


def read_defined(config):
    """Energy-balance settings of a TIRTEB input written in the TiRT parameter system."""
    input_file = Path(str(config["_input_path"]))
    input_base = input_file.parent
    surface, vegetation = canonical_models(config)

    def num(key, default):
        return number(config, key, default=default)

    def path(key, default=None):
        raw = value(config, key, default=default)
        if raw in (None, ""):
            return None
        file_path = Path(str(raw))
        return file_path if file_path.is_absolute() else input_base / file_path

    leaf = [num("cab", 58.0), num("cw", 0.013), num("cdm", 0.0036), num("cs", 0.0), num("n", 1.86)]
    soil_optical = [num("bsm_brightness", 0.25), num("bsm_latitude", 25.0), num("bsm_longitude", 45.0)]
    soil_moisture = num("soil_moisture", 0.25)
    satwater = num("satwater", 0.45)

    grid = SPECTRAL_GRID
    esun = _read_numbers(path("esun_file", "data/defined/Esun_.dat"))
    esky = _read_numbers(path("esky_file", "data/defined/Esky_.dat"))
    if len(esun) == len(grid) and len(esky) == len(grid):
        def normalize(values, start, end):
            total = sum(
                0.5 * (values[i] + values[i + 1]) * (grid[i + 1] - grid[i])
                for i in range(start, end - 1)
            )
            return [v / (total * 0.001) for v in values[start:end]]

        fesun = normalize(esun, 0, 2001) + normalize(esun, 2001, len(grid))
        fesky = normalize(esky, 0, 2001) + normalize(esky, 2001, len(grid))
    else:
        fesun = fesky = []

    coeff = _read_optical_coefficients(path("optical_coefficients_file", "data/defined/optipar_new.csv"))
    leaf_refl = leaf_trans = soil_refl = []
    if str(value(config, "leaf_optics_source", default="prospect")).strip().lower() in FILE_MODES:
        leaf_file = path("leaf_optics_file")
        if leaf_file is None:
            raise ValueError("leaf_optics_source=file 时必须设置 leaf_optics_file")
        leaf_refl = _read_measured_spectrum(leaf_file, grid, 1)
        leaf_trans = _read_measured_spectrum(leaf_file, grid, 2) if _column_count(leaf_file) > 2 else [0.0] * 2001
    elif len(coeff) >= 2001:
        leaf_refl, leaf_trans = _fluspect(coeff, leaf)
    if str(value(config, "soil_optics_source", default="bsm")).strip().lower() in FILE_MODES:
        soil_file = path("soil_optics_file")
        if soil_file is None:
            raise ValueError("soil_optics_source=file 时必须设置 soil_optics_file")
        soil_refl = _read_measured_spectrum(soil_file, grid, 1)
    elif len(coeff) >= 2001:
        soil_refl = bsm(coeff, soil_optical, soil_moisture)

    urban_reflectance = {
        "roof": num("reflectance_roof", 0.05), "wall": num("reflectance_wall", 0.08), "street": num("reflectance_street", 0.045),
    }
    urban_absorption = {component: num(f"absorption_{component}", 1.0 - reflectance) for component, reflectance in urban_reflectance.items()}
    urban_reflectance_spectra = {component: [reflectance] * len(grid) for component, reflectance in urban_reflectance.items()}
    for component in urban_reflectance:
        optics_file = path(f"{component}_optics_file")
        if optics_file is not None:
            urban_reflectance_spectra[component] = _read_measured_spectrum(optics_file, grid, 1) + [urban_reflectance[component]] * (len(grid) - 2001)
    urban_file = path("urban_spectrum_file") if str(value(config, "urban_spectrum_source", default="direct")).strip().lower() in FILE_MODES else None
    if urban_file is not None and read_spectrum(urban_file).shape[1] >= 7:
        columns = _read_measured_columns(urban_file, grid, (1, 2, 3))
        for component, column in zip(("roof", "wall", "street"), (1, 2, 3)):
            urban_reflectance_spectra[component] = columns[column]

    # Thermal-infrared emissivities at the observed thermal bands; the broadband value used by the
    # energy balance is their mean, and the TIR reflectance is 1 - emissivity (no leaf transmission).
    bands = [float(band) for band in wavelengths(config, input_file) if band >= 2.5] or [10.5]
    emissivity = spectra(config, input_file, bands)
    water_effect = int(round(num("soil_tir_water_effect", 0.0)))
    emissivity["soil_e"] = np.asarray([
        soil_tir_emissivity(
            soil_moisture, enabled=water_effect,
            dry_moisture=num("soil_tir_moisture_dry", 0.05), wet_moisture=num("soil_tir_moisture_wet", satwater),
            dry_emissivity=dry, wet_emissivity=num("soil_tir_emissivity_wet", 0.99),
        )
        for dry in emissivity["soil_e"]
    ])
    tir_emissivity = {name: [float(item) for item in values] for name, values in emissivity.items()}
    leaf_emissivity = float(np.mean(tir_emissivity["leaf_e"]))
    soil_emissivity = float(np.mean(tir_emissivity["soil_e"]))

    return DefinedStructure.from_mapping({
        "surface_model": surface,
        "vegetation_model": vegetation,
        "lai": 0.0 if vegetation == "bare" else num("lai", 2.0),
        "canopy_height": _aerodynamic_height(config, vegetation),
        "gleaf": num("g", 0.5), "lidfa": num("lidfa", -0.35), "lidfb": num("lidfb", -0.15),
        "hspot": num("hspot", 0.2),
        "leaf_width": num("leaf_width", 0.2),
        "leaf_reflectance": 1.0 - leaf_emissivity,
        "leaf_transmittance": 0.0,
        "leaf_emissivity": leaf_emissivity,
        "vcmax": num("vcmax", 80.0), "ballberry_m": num("ballberry_m", 9.0),
        "leaf_type": num("leaf_type", 0.0), "kv": num("kv", 0.6396),
        "rdparam": num("rdparam", 0.015),
        "tparam": [num("tparam1", 0.2), num("tparam2", 0.3), num("tparam3", 281.0), num("tparam4", 313.0), num("tparam5", 328.0)],
        "tyear": num("tyear", 20.0), "beta": num("beta", 0.4),
        "knpqs": num("knpqs", 0.0), "qls": num("qls", 1.0),
        "stressfactor": num("stressfactor", 1.0),
        "stomatal_method": str(value(config, "stomatal_method", default="farquhar")).strip().lower(),
        "stomatal_g0": num("stomatal_g0", 0.01),
        "soil_reflectance": 1.0 - soil_emissivity,
        "soil_emissivity": soil_emissivity,
        "hapke_enable": int(round(num("hapke_enable", 0.0))),
        "hapke_K": num("hapke_k", 1.0),
        "urban_reflectance_roof_spectrum": urban_reflectance_spectra["roof"],
        "urban_reflectance_wall_spectrum": urban_reflectance_spectra["wall"],
        "urban_reflectance_street_spectrum": urban_reflectance_spectra["street"],
        "urban_absorption_roof": urban_absorption["roof"],
        "urban_absorption_wall": urban_absorption["wall"],
        "urban_absorption_street": urban_absorption["street"],
        "urban_emissivity_roof": float(np.mean(tir_emissivity["roof_e"])),
        "urban_emissivity_wall": float(np.mean(tir_emissivity["wall_e"])),
        "urban_emissivity_street": float(np.mean(tir_emissivity["street_e"])),
        "urban_rss": num("urban_rss", 3662.8562),
        "urban_cs": num("urban_cs", 1000000.0),
        "urban_rhos": num("urban_rhos", 1800.0),
        "urban_lambda": num("urban_lambda", 1.5),
        "urban_rbs": num("urban_rbs", 10.0),
        "urban_temperature_c": num("urban_temperature_c", 20.0),
        "urban_depth_step": num("urban_depth_step", 0.05),
        "urban_wall_thickness": max(num("urban_wall_thickness", 0.5), 1.0e-3),
        "urban_wall_layers": max(1, int(round(num("urban_wall_layers", 3)))),
        "urban_internal_temperature_c": num("urban_internal_temperature_c", 24.0),
        "urban_wall_diffuse_factor": min(max(num("urban_wall_diffuse_factor", 0.5), 0.0), 1.0),
        "urban_street_radiation_factor": num("urban_street_radiation_factor", 1.0),
        "soil_rss": num("soil_rss", 200.0), "soil_cs": num("soil_cs", 1180.0),
        "soil_rss_method": str(value(config, "soil_rss_method", default="sellers")).strip().lower(),
        "soil_rbs": num("soil_rbs", 10.0),
        "soil_rhos": num("soil_rhos", 1100.0), "soil_lambda": num("soil_lambda", 1.55),
        "soil_heat_method": str(value(config, "soil_heat_method", default="surface")).strip().lower(),
        "urban_heat_method": str(value(config, "urban_heat_method", default="depth")).strip().lower(),
        "temperature_method": str(value(config, "temperature_method", default="robust")).strip().lower(),
        "soil_heat_fraction": num("soil_heat_fraction", 0.25),
        "urban_heat_fraction": num("urban_heat_fraction", 0.4),
        "soil_depth_step": num("soil_depth_step", 0.05),
        "soil_temperature_method": str(value(config, "soil_temperature_method", "temperature_method", default="robust")).strip().lower(),
        "urban_temperature_method": str(value(config, "temperature_method", default="robust")).strip().lower(),
        "soil_moisture": min(max(soil_moisture, 0.0), 0.45), "soil_temperature_c": num("soil_temperature_c", 20.0),
        "satwater": max(satwater, 1.0e-6),
        "dtime": max(num("dtime", 1800.0), 1.0),
        "z": num("z", 15.0), "ca": num("ca", 380.0), "oa": num("oa", 209.0),
        "shortwave_diffuse_fraction": min(max(num("shortwave_diffuse_fraction", -1.0), -1.0), 1.0),
        "aero_zo": num("aero_zo", 10.0), "aero_d": num("aero_d", 1.34), "aero_rbc": num("aero_rbc", 10.0),
        "aero_cr": num("aero_cr", 0.35), "aero_cd": num("aero_cd", 0.2), "aero_cd1": num("aero_cd1", 20.6),
        "aero_psicor": num("aero_psicor", 0.2), "aero_cssoil": num("aero_cssoil", 0.01),
        "aero_rwc": num("aero_rwc", 0.0), "aero_rbs": num("aero_rbs", 10.0),
        "slope_angle": num("slope_angle", 0.0) if surface == "slope" else 0.0,
        "slope_aspect": num("slope_aspect", 0.0) if surface == "slope" else 0.0,
        "terrain_shapes": tuple(map(tuple, terrain_shapes(config))) if surface == "terrain" else (),
        "row_width": num("row_width", 0.5),
        "row_blank": num("row_blank", 0.5),
        "row_height": num("row_height", 1.0),
        "row_azimuth": num("row_azimuth", 0.0),
        "crown_rad_a": max(num("crown_rad_a", 3.0), 1.0e-6),
        "crown_rad_b": max(num("crown_rad_b", 3.0), 1.0e-6),
        "crown_height": max(num("crown_height", 3.0), 0.0),
        "crown_density": max(num("crown_density", 0.02), 1.0e-6),
        "urban_density": num("building_density", 0.0015),
        "urban_building_length": num("building_length", 10.0),
        "urban_building_width": num("building_width", 10.0),
        "urban_building_height": num("building_height", 20.0),
        "urban_building_azimuth": num("building_azimuth", 0.0),
        "urban_roughness_length": max(num("urban_roughness_length", 0.15), 0.01),
        "urban_canyon_roughness_ratio": max(num("urban_canyon_roughness_ratio", 0.05), 0.005),
        "urban_canyon_min_ratio": max(num("urban_canyon_min_ratio", 1.25), 1.0),
        "urban_leaf_wind_decay": max(num("urban_leaf_wind_decay", 0.5), 0.0),
        "urban_leaf_boundary_factor": max(num("urban_leaf_boundary_factor", 70.0), 0.0),
        "wavelengths": grid, "fesun": fesun, "fesky": fesky,
        "leaf_reflectance_spectrum": leaf_refl,
        "leaf_transmittance_spectrum": leaf_trans,
        "soil_reflectance_spectrum": soil_refl,
        "tir_wavelengths": bands,
        "tir_emissivity": tir_emissivity,
    })


def _read_numbers(file_path):
    values = []
    for line in Path(file_path).read_text(encoding="utf-8-sig").splitlines():
        try:
            values.append(float(line.split()[0]))
        except (IndexError, ValueError):
            continue
    return values


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
