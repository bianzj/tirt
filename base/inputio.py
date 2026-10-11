"""Input readers and the parameter system shared by TiRT (run.py) and TIRTEB (tirteb_run.py)."""

import csv
import re
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FILE_MODES = {"file", "txt", "read"}
SURFACE_ALIASES = {"flat": "plane", "bare": "plane", "single_slope": "slope", "sloped": "slope", "building": "urban", "city": "urban"}
VEGETATION_ALIASES = {"homogeneous": "hom", "raw": "row", "row_crop": "row", "forest": "crown", "bare_soil": "bare"}
PAIRED_GEOMETRY_MODES = {"zip", "pair", "paired"}
# Typical band emissivities at 8, 10.5 and 12 um.
EMISSIVITY_DEFAULTS = {
    "leaf_e": (0.975, 0.985, 0.980), "soil_e": (0.920, 0.955, 0.965), "terrain_e": (0.920, 0.950, 0.960),
    "roof_e": (0.920, 0.950, 0.940), "wall_e": (0.900, 0.930, 0.920), "street_e": (0.940, 0.955, 0.950),
}


def read_input(input_path):
    """Read grouped input.csv into flat and scoped keys."""
    config = {"_input_path": str(Path(input_path).resolve())}
    sections = {}
    with Path(input_path).open(newline="", encoding="utf-8-sig") as file:
        for raw in csv.reader(file):
            if not raw:
                continue
            cells = [cell.strip() for cell in raw]
            if not cells or not cells[0] or cells[0].startswith("#"):
                continue
            if cells[0].startswith("["):
                continue
            if len(cells) < 3 or not cells[1]:
                continue
            section, key, value = cells[0].lower(), cells[1].lower(), cells[2]
            config[key] = value
            config[f"{section}_{key}"] = value
            sections.setdefault(section, {})[key] = value
    config["_sections"] = sections
    return config


def numbers(value, default=()):
    if value in (None, ""):
        return tuple(float(item) for item in default)
    text = str(value).replace(";", ",")
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def value(config, *keys, default=None):
    for key in keys:
        item = config.get(key)
        if item not in (None, ""):
            return item
    return default


def resolve_path(input_path, raw_path, default=None):
    raw_path = raw_path if raw_path not in (None, "") else default
    if raw_path in (None, ""):
        return None
    path = Path(str(raw_path))
    if path.is_absolute():
        return path
    return Path(input_path).resolve().parent / path


def _numeric_rows(path, minimum_columns=2):
    rows = []
    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        fields = re.split(r"[\s,;]+", text)
        try:
            row = [float(item) for item in fields]
        except ValueError:
            continue
        if len(row) < minimum_columns:
            raise ValueError(f"{path}:{line_number} 至少需要 {minimum_columns} 个数值列")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path} 没有有效的数值数据")
    return np.asarray(rows, dtype=float)


def read_geometry(path, default_sza=30.0, default_saa=0.0):
    """Read vza vaa, or vza vaa sza saa, one observation per line."""
    rows = _numeric_rows(path, minimum_columns=2)
    if rows.shape[1] not in (2, 4):
        raise ValueError("geometry file 每行必须是 vza vaa 或 vza vaa sza saa")
    if np.any(rows[:, 0] < 0.0) or np.any(rows[:, 0] > 180.0):
        raise ValueError("geometry file 中 vza 必须在 0 到 180 度之间")
    if np.any(rows[:, 1] < 0.0) or np.any(rows[:, 1] > 360.0):
        raise ValueError("geometry file 中 vaa 必须在 0 到 360 度之间")
    sza = rows[:, 2] if rows.shape[1] == 4 else np.full(len(rows), default_sza)
    saa = rows[:, 3] if rows.shape[1] == 4 else np.full(len(rows), default_saa)
    return rows[:, 0], rows[:, 1], sza, saa


def read_wavelengths(path):
    rows = _numeric_rows(path, minimum_columns=1)
    wavelengths = rows[:, 0]
    if np.any(wavelengths <= 0.0):
        raise ValueError("spectral file 中波长必须大于 0")
    return wavelengths


def read_spectrum(path):
    """Read a wavelength-first spectrum table, allowing a text header."""
    rows = _numeric_rows(path, minimum_columns=2)
    order = np.argsort(rows[:, 0])
    return rows[order]


def interpolate(table, wavelengths, column):
    if table.shape[1] <= column:
        raise ValueError(f"spectrum file 缺少第 {column + 1} 列")
    return np.interp(wavelengths, table[:, 0], table[:, column])


def strict_band_values(config, keys, size, default, label):
    raw = value(config, *keys)
    if raw in (None, ""):
        return np.full(size, float(default), dtype=float)
    result = np.asarray(numbers(raw), dtype=float)
    if result.size != size:
        raise ValueError(f"{label} 的数目必须与 wavelengths 的数目一致: {result.size} != {size}")
    return result


def repeat_to(values, size, label):
    values = np.asarray(values, dtype=float)
    if values.size == size:
        return values
    if values.size == 1:
        return np.repeat(values.item(), size)
    raise ValueError(f"{label} 的数目必须与观测记录数一致: {values.size} != {size}")


def number(config, *keys, default=0.0):
    return float(value(config, *keys, default=default))


def canonical_models(config):
    """Return (surface_model, vegetation_model); terrain and urban vegetation is always discrete crowns."""
    surface = str(value(config, "surface_model", default="plane")).strip().lower()
    surface = SURFACE_ALIASES.get(surface, surface)
    vegetation = str(value(config, "vegetation_model", default="hom")).strip().lower()
    vegetation = VEGETATION_ALIASES.get(vegetation, vegetation)
    if surface not in {"plane", "slope", "terrain", "urban"}:
        raise ValueError("surface_model 必须是 plane、slope、terrain 或 urban")
    if vegetation not in {"bare", "hom", "row", "crown"}:
        raise ValueError("vegetation_model 必须是 bare、hom、row 或 crown")
    if surface in {"terrain", "urban"} and vegetation in {"hom", "row"}:
        warnings.warn(f"{surface} 地表的植被只支持离散树冠模型，{vegetation} 按 crown 计算", RuntimeWarning, stacklevel=2)
        vegetation = "crown"
    return surface, vegetation


def stepped_angles(maximum, step):
    angles = list(np.arange(0.0, maximum + 1.0e-12, step))
    if not angles or not np.isclose(angles[-1], maximum):
        angles.append(maximum)
    return tuple(angles)


def paired_geometry(config):
    return str(value(config, "geometry_pair_mode", default="cross")).strip().lower() in PAIRED_GEOMETRY_MODES


def geometry(config, input_path):
    """Observation geometry as (vza, vaa, sza, saa) tuples, selected by geometry_mode 0-3."""
    sza_values = numbers(value(config, "sza"), (30.0,))
    saa_values = numbers(value(config, "saa"), (0.0,))
    sza_default, saa_default = sza_values[0], saa_values[0]
    try:
        mode = int(float(value(config, "geometry_mode", default=0)))
    except ValueError as exc:
        raise ValueError("geometry_mode 必须是 0、1、2 或 3") from exc
    if mode not in {0, 1, 2, 3}:
        raise ValueError("geometry_mode 必须是 0、1、2 或 3")

    if mode == 2:
        path = resolve_path(input_path, value(config, "geometry_file", default="geometry.txt"))
        vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
        return list(zip(vza, vaa, sza, saa))

    if mode == 1:
        plane = str(value(config, "geometry_view", default="both")).strip().lower()
        relative_azimuths = {
            "parallel": (0.0, 180.0),
            "perpendicular": (90.0, 270.0),
            "both": (0.0, 90.0, 180.0, 270.0),
        }.get(plane)
        if relative_azimuths is None:
            raise ValueError("geometry_view 必须是 parallel、perpendicular 或 both")
        raw_vza = value(config, "principal_vza")
        if raw_vza not in (None, ""):
            principal_vza = numbers(raw_vza)
        else:
            step = number(config, "principal_vza_step", default=10.0)
            maximum = number(config, "principal_vza_max", default=75.0)
            if step <= 0.0 or maximum < 0.0:
                raise ValueError("principal_vza_step 必须大于 0，principal_vza_max 不能小于 0")
            principal_vza = stepped_angles(maximum, step)
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
        raw_step = number(config, "principal_vza_step", default=10.0)
        maximum = number(config, "principal_vza_max", default=75.0)
        if raw_step <= 0.0 or maximum < 0.0:
            raise ValueError("principal_vza_step 必须大于 0，principal_vza_max 不能小于 0")
        principal_vza = stepped_angles(maximum, raw_step)
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
        source = str(value(config, "additional_geometry_source", default="none")).strip().lower()
        additional = []
        if include_hemisphere:
            raw_path = value(config, "additional_geometry_file")
            path = resolve_path(input_path, raw_path) if raw_path not in (None, "") else ROOT / "data" / "directions" / "hemisphere_directions.txt"
            vza, vaa, sza, saa = read_geometry(path, sza_default, saa_default)
            additional = [item for item in zip(vza, vaa, sza, saa) if item[0] <= maximum + 1.0e-9]
        if not include_hemisphere and source in FILE_MODES:
            raw_file = value(config, "additional_geometry_file")
            if raw_file not in (None, ""):
                vza, vaa, sza, saa = read_geometry(resolve_path(input_path, raw_file), sza_default, saa_default)
                additional.extend(zip(vza, vaa, sza, saa))
        if source in {"direct", "manual", "angles"}:
            vza = numbers(value(config, "observation_vza"))
            vaa = numbers(value(config, "observation_vaa"))
            if len(vza) != len(vaa):
                raise ValueError("手动观测角度必须一一对应：observation_vza 和 observation_vaa 数量不同")
            additional.extend((view_zenith, view_azimuth, sza_default, saa_default)
                              for view_zenith, view_azimuth in zip(vza, vaa))
        combined = []
        seen = set()
        for item in principal + additional:
            key = tuple(round(float(item_value), 8) for item_value in item)
            if key not in seen:
                seen.add(key)
                combined.append(item)
        return combined

    vza = numbers(value(config, "vza"), (0.0, 30.0, 60.0))
    vaa = numbers(value(config, "vaa"), (0.0, 90.0, 180.0))
    if paired_geometry(config):
        size = max(len(vza), len(vaa), len(sza_values), len(saa_values))
        return list(zip(
            repeat_to(vza, size, "vza"), repeat_to(vaa, size, "vaa"),
            repeat_to(sza_values, size, "sza"), repeat_to(saa_values, size, "saa"),
        ))
    return [(view_zenith, view_azimuth, solar_zenith, solar_azimuth)
            for view_zenith in vza for view_azimuth in vaa
            for solar_zenith in sza_values for solar_azimuth in saa_values]


def wavelengths(config, input_path):
    if str(value(config, "spectral_source", default="direct")).strip().lower() in FILE_MODES:
        path = resolve_path(input_path, value(config, "spectrals_file", default="wavelengths.txt"))
        return np.asarray(read_wavelengths(path), dtype=float)
    result = np.asarray(numbers(value(config, "wavelengths", default=10.5)), dtype=float)
    if result.size == 0 or np.any(result <= 0.0):
        raise ValueError("wavelengths 必须包含正数波长")
    return result


def _band_values(config, key, bands, defaults):
    raw = value(config, key)
    if raw in (None, ""):
        return np.interp(bands, [8.0, 10.5, 12.0], defaults)
    result = np.asarray(numbers(raw), dtype=float)
    if result.size == 1:
        return np.repeat(result.item(), len(bands))
    if result.size != len(bands):
        raise ValueError(f"{key} 的数目必须是 1 或与热红外波段数一致: {result.size} != {len(bands)}")
    return result


def spectra(config, input_path, bands):
    """Thermal-infrared emissivities of every component at the given bands (um)."""
    bands = np.asarray(bands, dtype=float)
    keys = {"leaf_e": "emissivity_leaf", "soil_e": "emissivity_soil", "terrain_e": "emissivity_terrain",
            "roof_e": "emissivity_roof", "wall_e": "emissivity_wall", "street_e": "emissivity_street"}
    result = {name: _band_values(config, key, bands, EMISSIVITY_DEFAULTS[name]) for name, key in keys.items()}

    def from_file(source_key, file_key):
        if str(value(config, source_key, default="direct")).strip().lower() not in FILE_MODES:
            return None
        return read_spectrum(resolve_path(input_path, value(config, file_key)))

    for name, component in (("leaf_e", "vegetation"), ("soil_e", "soil"), ("terrain_e", "terrain"),
                            ("roof_e", "roof"), ("wall_e", "wall"), ("street_e", "street")):
        table = from_file(f"{component}_spectrum_source", f"{component}_spectrum_file")
        if table is not None:
            result[name] = interpolate(table, bands, 1)
    table = from_file("urban_spectrum_source", "urban_spectrum_file")
    if table is not None:
        # 7 columns: wavelength, roof/wall/street VNIR reflectance, roof/wall/street TIR emissivity.
        columns = (4, 5, 6) if table.shape[1] >= 7 else (1, 2, 3)
        if table.shape[1] < 4:
            raise ValueError("urban_spectrum_file 至少需要 wavelength roof wall street 四列")
        for name, column in zip(("roof_e", "wall_e", "street_e"), columns):
            result[name] = interpolate(table, bands, column)
    return result


def terrain_shapes(config):
    """Compound-terrain cones as rows of (height, radius, density m-2); the cones must cover less than the scene."""
    heights = numbers(value(config, "terrain_height"), (20.0,))
    radii = numbers(value(config, "terrain_radius"), (10.0,))
    densities = numbers(value(config, "terrain_density"), (0.001,))
    size = max(len(heights), len(radii), len(densities))
    shapes = np.column_stack([repeat_to(heights, size, "terrain_height"), repeat_to(radii, size, "terrain_radius"),
                              repeat_to(densities, size, "terrain_density")])
    if float(np.sum(shapes[:, 2] * np.pi * shapes[:, 1] ** 2)) >= 1.0:
        raise ValueError("terrain_density * pi * terrain_radius^2 的总占据率必须小于 1")
    return shapes
