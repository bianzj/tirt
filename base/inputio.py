"""Input and file readers for the unified TiRT runner."""

import csv
import re
import shutil
import warnings
from pathlib import Path

import numpy as np


SNAPSHOT_MARKER = "# Effective TiRT input"
BAND_PROPERTIES = (
    "leaf_reflectance", "leaf_transmittance", "soil_reflectance",
    "street_reflectance", "wall_reflectance", "roof_reflectance",
)
KNOWN_KEYS = frozenset((
    "case_name", "cases_dir", "result_dir",
    "surface_model", "surface", "vegetation_model", "vegetation",
    "slope_angle", "terrain_slope", "slope_aspect", "terrain_aspect",
    "terrain_radius", "terrain_height", "terrain_density", "terrain_type", "surface_type", "distribution",
    "building_density", "building_length", "building_width", "building_height", "building_azimuth",
    "urban_transport",
    "lai", "vegetation_lai", "hspot", "vegetation_hspot", "canopy_closure",
    "height", "width", "leaf_width",
    "row_width", "row_blank", "row_height", "row_azimuth", "thermal_row_geometry", "row_closure",
    "crown_rad_a", "crown_radius", "crown_rad_b", "crown_height", "crown_density", "stand_density",
    "crown_bottom", "crown_base_height",
    "slope_density", "slope_crown_rad_a", "slope_crown_rad_b", "slope_crown_radius",
    "thermal_crown_geometry", "thermal_hotspot_inner", "thermal_hotspot_outer",
    "forest_transport", "thermal_forest_transport",
    "optical_terrain_method", "optical_terrain_source_depth", "optical_terrain_primary",
    "optical_crown_geometry", "optical_crown_method", "optical_crown_quadrature", "optical_crown_source",
    "optical_hotspot_inner", "optical_hotspot_outer", "optical_row_method",
    "shortwave_diffuse_fraction", "shortwave_sunlit_surface", "shortwave_enforce_solar_limit",
    "solar_day_of_year",
    "spectral_source", "spectrals_file", "wavelengths", *BAND_PROPERTIES,
    "hapke_enable", "soil_hapke_enable", "hapke_k", "soil_hapke_k", "hapke_optical_enable",
    *(f"{material}_{lighting}_{unit}" for material in ("soil", "leaf", "terrain", "roof", "wall", "street")
      for lighting in ("sunlit", "shaded") for unit in ("c", "k")),
    "sky_temperature_k", "sky_temperature_c",
    "geometry_mode", "observation_mode", "direction_mode", "geometry_source", "direction_source",
    "geometry_file", "directions_file", "vza", "vaa", "sza", "solar_zenith", "saa", "solar_azimuth",
    "observation_vza", "observation_vaa", "additional_vza", "additional_vaa", "geometry_pair_mode",
    "geometry_view", "principal_plane", "principal_vza", "principal_plane_vza", "principal_vza_step",
    "principal_vza_max", "include_principal", "include_hemisphere", "hemisphere_geometry_file",
    "additional_geometry_source", "observation_source", "additional_geometry_file",
    "observation_geometry_file",
))


def read_input(input_path):
    """Read grouped input.csv into flat and scoped keys."""
    config = {"_input_path": str(Path(input_path).resolve())}
    sections = {}
    try:
        with Path(input_path).open(newline="", encoding="utf-8-sig") as file:
            reader = csv.reader(file)
            for raw in reader:
                cells = [cell.strip() for cell in raw]
                if not cells or not cells[0] or cells[0].startswith("#"):
                    continue
                if cells[0].startswith("["):
                    continue
                if len(cells) < 3:
                    warnings.warn(f"{input_path}:{reader.line_num}: row needs section,key,value; ignored: {raw}")
                    continue
                if not cells[1]:
                    continue
                section, key, value = cells[0].lower(), cells[1].lower(), cells[2]
                if section != "effective" and key not in KNOWN_KEYS:
                    warnings.warn(f"{input_path}:{reader.line_num}: unknown parameter '{cells[1]}' is ignored")
                config[key] = value
                # Snapshots contain exact lookup keys, including scoped aliases.
                if section != "effective":
                    config[f"{section}_{key}"] = value
                sections.setdefault(section, {})[key] = value
    except UnicodeDecodeError as exc:
        raise ValueError(f"{input_path} is not UTF-8 encoded; please save it as UTF-8 (请保存为 UTF-8 编码)") from exc
    config["_sections"] = sections
    return config


def is_input_snapshot(path):
    with Path(path).open(encoding="utf-8-sig", errors="replace") as file:
        return file.readline().startswith(SNAPSHOT_MARKER)


def numbers(value, default=(), key=None):
    if value in (None, ""):
        return tuple(float(item) for item in default)
    text = str(value).replace(";", ",")
    try:
        result = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    except ValueError as exc:
        raise ValueError(f"{key or 'Numeric parameter'} must be numbers separated by ';', got {value!r}") from exc
    if not all(np.isfinite(item) for item in result):
        raise ValueError(f"{key or 'Numeric parameters'} must be finite; NaN and Inf are not allowed")
    return result


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
        text = line.split("#", 1)[0].strip()
        if not text or text.startswith("#"):
            continue
        fields = re.split(r"[\s,;]+", text)
        try:
            row = [float(item) for item in fields]
        except ValueError as exc:
            numeric_token = False
            for field in fields:
                try:
                    float(field)
                    numeric_token = True
                    break
                except ValueError:
                    pass
            if rows or numeric_token:
                raise ValueError(f"{path}:{line_number} numeric data row contains an invalid field") from exc
            continue  # Permit text headers before the first data row.
        if len(row) < minimum_columns:
            raise ValueError(f"{path}:{line_number} needs at least {minimum_columns} numeric columns")
        if not all(np.isfinite(item) for item in row):
            raise ValueError(f"{path}:{line_number} must not contain NaN or Inf")
        if rows and len(row) != len(rows[0]):
            raise ValueError(f"{path}:{line_number} has a different number of numeric columns than the previous row")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path} contains no valid numeric data")
    return np.asarray(rows, dtype=float)


def read_geometry(path, default_sza=30.0, default_saa=0.0):
    """Read vza vaa, or vza vaa sza saa, one observation per line."""
    rows = _numeric_rows(path, minimum_columns=2)
    if rows.shape[1] not in (2, 4):
        raise ValueError("Each geometry file row must be vza vaa or vza vaa sza saa")
    if np.any(rows[:, 0] < 0.0) or np.any(rows[:, 0] > 90.0):
        raise ValueError("vza in the geometry file must be between 0 and 90 degrees")
    if np.any(rows[:, 1] < 0.0) or np.any(rows[:, 1] > 360.0):
        raise ValueError("vaa in the geometry file must be between 0 and 360 degrees")
    sza = rows[:, 2] if rows.shape[1] == 4 else np.full(len(rows), default_sza)
    saa = rows[:, 3] if rows.shape[1] == 4 else np.full(len(rows), default_saa)
    if np.any(sza < 0.0) or np.any(sza > 90.0) or not np.all(np.isfinite(sza)):
        raise ValueError("sza in the geometry file must be between 0 and 90 degrees")
    if np.any(saa < 0.0) or np.any(saa > 360.0) or not np.all(np.isfinite(saa)):
        raise ValueError("saa in the geometry file must be between 0 and 360 degrees")
    return rows[:, 0], rows[:, 1], sza, saa


_SEPARATORS = r"[;,，；]"


def band_cells(text):
    """One band-table row: a float per band, None for a blank or "-" cell."""
    text = str(text if text is not None else "").strip()
    if not text:
        return []
    parts = re.split(_SEPARATORS, text) if re.search(_SEPARATORS, text) else text.split()
    cells = []
    for part in parts:
        part = part.strip()
        if part in ("", "-"):
            cells.append(None)
            continue
        number = float(part)
        if not np.isfinite(number):
            raise ValueError("Band table values must be finite")
        cells.append(number)
    return cells


def _is_label(token):
    if token == "-" or re.search(_SEPARATORS, token):
        return False
    try:
        float(token)
        return False
    except ValueError:
        return True


def read_band_cells(path, row_count):
    """Read one row per property (first row = wavelengths); each row may start with a text label."""
    rows = []
    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        text = line.split("#", 1)[0].strip()
        if re.search(r"[:：]", text):
            text = re.split(r"[:：]", text, maxsplit=1)[1]
        tokens = text.split()
        while tokens and _is_label(tokens[0]):
            tokens.pop(0)
        if not tokens:
            continue
        try:
            rows.append(band_cells(" ".join(tokens)))
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number} contains a non-numeric value") from exc
    if len(rows) != row_count:
        raise ValueError(f"{path} must contain {row_count} rows, found {len(rows)}")
    return rows


def read_band_table(path, row_count):
    """Band table in which every cell must hold a number."""
    rows = read_band_cells(path, row_count)
    if any(cell is None for row in rows for cell in row):
        raise ValueError(f"{path}: every band table cell needs a value")
    return [tuple(row) for row in rows]


def strict_band_values(config, keys, size, default, label):
    raw = value(config, *keys)
    if raw in (None, ""):
        return np.full(size, float(default), dtype=float)
    result = np.asarray(numbers(raw), dtype=float)
    if result.size != size:
        raise ValueError(f"{label} must have one value per wavelength: {result.size} != {size}")
    return result


def repeat_to(values, size, label):
    values = np.asarray(values, dtype=float)
    if values.size == size:
        return values
    if values.size == 1:
        return np.repeat(values.item(), size)
    raise ValueError(f"{label} must have one value per observation record: {values.size} != {size}")


def _active_input_files(config):
    """Yield the files actually used by the TiRT directional runner."""
    file_modes = {"file", "txt", "read"}
    if str(value(config, "spectral_source", default="direct")).strip().lower() in file_modes:
        yield "spectrals_file", value(config, "spectrals_file", default="data/spectrals/optical_properties.txt")
    raw_mode = value(config, "geometry_mode", "observation_mode", "direction_mode")
    mode = int(float(raw_mode)) if raw_mode not in (None, "") else (
        2 if str(value(config, "geometry_source", "direction_source", default="direct")).strip().lower() in file_modes else 0
    )
    if mode == 2:
        yield "geometry_file", value(config, "geometry_file", "directions_file", default="geometry.txt")
    elif mode == 3:
        if str(value(config, "include_hemisphere", default="0")).strip().lower() not in {"0", "false", "no", "off"}:
            yield "hemisphere_geometry_file", value(config, "hemisphere_geometry_file", default="data/directions/hemisphere_directions.txt")
        if str(value(config, "additional_geometry_source", "observation_source", default="none")).strip().lower() in file_modes:
            path = value(config, "additional_geometry_file", "observation_geometry_file")
            if path not in (None, ""):
                yield "additional_geometry_file", path


def write_input_snapshot(config, input_path, snapshot_path, output_dir=None):
    """Save effective overrides and local copies of referenced inputs for replay."""
    input_path, snapshot_path = Path(input_path).resolve(), Path(snapshot_path).resolve()
    if input_path == snapshot_path:
        raise ValueError("The input snapshot must not overwrite the source input.csv")
    flat = {key: str(item) for key, item in config.items()
            if not key.startswith("_") and item is not None}
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    for index, (key, raw_path) in enumerate(_active_input_files(config), 1):
        source = resolve_path(input_path, raw_path)
        if source is None or not source.is_file():
            raise ValueError(f"File referenced by {key} does not exist: {raw_path}")
        relative = Path("inputs") / f"{index:02d}_{key}{source.suffix or '.txt'}"
        destination = snapshot_path.parent / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)
        flat[key] = relative.as_posix()
    if output_dir is not None:
        # A moved case must never replay into the original result directory.
        flat["result_dir"] = "."
    with snapshot_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([f"{SNAPSHOT_MARKER} snapshot written by TiRT; overwritten on rerun", "", "", ""])
        for key, item in flat.items():
            writer.writerow(["effective", key, item, ""])
    return snapshot_path
