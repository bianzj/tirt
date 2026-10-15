"""Input and file readers for the unified TiRT runner."""

import csv
import re
import shutil
from pathlib import Path

import numpy as np


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
            # Snapshots contain exact lookup keys, including scoped aliases.
            if section != "effective":
                config[f"{section}_{key}"] = value
            sections.setdefault(section, {})[key] = value
    config["_sections"] = sections
    return config


def numbers(value, default=()):
    if value in (None, ""):
        return tuple(float(item) for item in default)
    text = str(value).replace(";", ",")
    result = tuple(float(item.strip()) for item in text.split(",") if item.strip())
    if not all(np.isfinite(item) for item in result):
        raise ValueError("数值参数必须为有限数值，不能使用 NaN 或 Inf")
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
                raise ValueError(f"{path}:{line_number} 数值数据行包含无效字段") from exc
            continue  # Permit text headers before the first data row.
        if len(row) < minimum_columns:
            raise ValueError(f"{path}:{line_number} 至少需要 {minimum_columns} 个数值列")
        if not all(np.isfinite(item) for item in row):
            raise ValueError(f"{path}:{line_number} 不能包含 NaN 或 Inf")
        if rows and len(row) != len(rows[0]):
            raise ValueError(f"{path}:{line_number} 数值列数与前一行不一致")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path} 没有有效的数值数据")
    return np.asarray(rows, dtype=float)


def read_geometry(path, default_sza=30.0, default_saa=0.0):
    """Read vza vaa, or vza vaa sza saa, one observation per line."""
    rows = _numeric_rows(path, minimum_columns=2)
    if rows.shape[1] not in (2, 4):
        raise ValueError("geometry file 每行必须是 vza vaa 或 vza vaa sza saa")
    if np.any(rows[:, 0] < 0.0) or np.any(rows[:, 0] > 90.0):
        raise ValueError("geometry file 中 vza 必须在 0 到 90 度之间")
    if np.any(rows[:, 1] < 0.0) or np.any(rows[:, 1] > 360.0):
        raise ValueError("geometry file 中 vaa 必须在 0 到 360 度之间")
    sza = rows[:, 2] if rows.shape[1] == 4 else np.full(len(rows), default_sza)
    saa = rows[:, 3] if rows.shape[1] == 4 else np.full(len(rows), default_saa)
    if np.any(sza < 0.0) or np.any(sza > 90.0) or not np.all(np.isfinite(sza)):
        raise ValueError("geometry file 中 sza 必须在 0 到 90 度之间")
    if np.any(saa < 0.0) or np.any(saa > 360.0) or not np.all(np.isfinite(saa)):
        raise ValueError("geometry file 中 saa 必须在 0 到 360 度之间")
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


def _active_input_files(config):
    """Yield the files actually used by the TiRT directional runner."""
    file_modes = {"file", "txt", "read"}
    if str(value(config, "spectral_source", "wavelength_source", default="direct")).strip().lower() in file_modes:
        yield "spectrals_file", value(config, "spectrals_file", "wavelengths_file", default="wavelengths.txt")
    for component in ("vegetation", "soil", "terrain", "urban", "roof", "wall", "street"):
        if str(value(config, f"{component}_spectrum_source", default="direct")).strip().lower() in file_modes:
            yield f"{component}_spectrum_file", value(config, f"{component}_spectrum_file")
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
        raise ValueError("输入快照不能覆盖源 input.csv")
    flat = {key: str(item) for key, item in config.items()
            if not key.startswith("_") and item is not None}
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    for index, (key, raw_path) in enumerate(_active_input_files(config), 1):
        source = resolve_path(input_path, raw_path)
        if source is None or not source.is_file():
            raise ValueError(f"{key} 引用的文件不存在: {raw_path}")
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
        writer.writerow(["# Effective TiRT input", "", "", ""])
        for key, item in flat.items():
            writer.writerow(["effective", key, item, ""])
    return snapshot_path
