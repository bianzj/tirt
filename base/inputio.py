"""Input and file readers for the unified TiRT runner."""

import csv
import re
from pathlib import Path

import numpy as np


def read_input(input_path):
    """Read tirteb-style grouped input.csv into flat and scoped keys."""
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
