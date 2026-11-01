"""DEM facet ensemble for the TiRT directional models.

A DEM is decomposed into slope/aspect classes plus a cast-shadow mask.  Each
class is evaluated by the existing single-slope model, so the DEM mode uses the
same canopy, optical and thermal closures as ``surface_model=slope`` without a
new radiative-transfer equation.  It is intended for directional TiRT runs
(``run.py`` / ``thermal_run.py`` / ``optical_run.py``); TIRTEB's point energy
balance keeps its four shared temperatures and does not consume this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from .solar import slope_direction


@dataclass(frozen=True)
class Grid:
    z: np.ndarray
    cell: float
    path: Path


def _resolve_path(config):
    raw = config.get("dem_file")
    if raw in (None, ""):
        raise ValueError("surface_model=dem requires dem_file (an ESRI ASCII grid .asc)")
    path = Path(str(raw))
    if not path.is_absolute():
        base = Path(str(config.get("_input_path") or ".")).resolve().parent
        path = base / path
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f"dem_file does not exist: {path}")
    return path


@lru_cache(maxsize=8)
def _grid_cached(path_text, mtime):
    lines = Path(path_text).read_text(encoding="utf-8", errors="replace").splitlines()
    header = {}
    first_data = None
    for index, line in enumerate(lines):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        match = text.split()
        key = match[0].upper()
        if key in {"NCOLS", "NROWS", "XLLCORNER", "YLLCORNER", "XLLCENTER", "YLLCENTER",
                   "CELLSIZE", "NODATA_VALUE"} and len(match) >= 2:
            header[key] = float(match[1])
            continue
        first_data = index
        break
    if first_data is None or not {"NCOLS", "NROWS", "CELLSIZE"} <= set(header):
        raise ValueError(f"{path_text}: not an ESRI ASCII grid with NCOLS/NROWS/CELLSIZE")
    ncols, nrows = int(header["NCOLS"]), int(header["NROWS"])
    values = []
    for line in lines[first_data:first_data + nrows]:
        values.extend(float(item) for item in line.split())
    if len(values) < ncols * nrows:
        raise ValueError(f"{path_text}: expected {ncols * nrows} values, found {len(values)}")
    z = np.asarray(values[:ncols * nrows], dtype=float).reshape(nrows, ncols)
    nodata = header.get("NODATA_VALUE")
    if nodata is not None:
        z = np.where(z == nodata, np.nan, z)
    return Grid(z=z, cell=float(header["CELLSIZE"]), path=Path(path_text))


def read_grid(config):
    path = _resolve_path(config)
    return _grid_cached(str(path), path.stat().st_mtime_ns)


def slope_aspect(grid):
    z = np.nan_to_num(grid.z, nan=0.0)
    dzdx = np.gradient(z, grid.cell, axis=1)          # east
    dzdy_north = -np.gradient(z, grid.cell, axis=0)   # rows north -> south
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy_north)))
    aspect = np.degrees(np.arctan2(-dzdx, -dzdy_north)) % 360.0
    slope[0, :] = slope[-1, :] = slope[:, 0] = slope[:, -1] = 0.0
    return slope, aspect


def cast_shadow(grid, sza, saa):
    """Cast-shadow mask by a horizon ray march along the solar azimuth."""
    z = np.nan_to_num(grid.z, nan=0.0)
    nrows, ncols = z.shape
    if sza >= 90.0:
        return np.ones_like(z, dtype=bool)
    sun_elevation = math.radians(90.0 - sza)
    dx, dy = math.sin(math.radians(saa)), math.cos(math.radians(saa))
    rows, cols = np.mgrid[0:nrows, 0:ncols]
    max_tan = np.zeros_like(z)
    for step in range(1, int(math.hypot(nrows, ncols) * 1.1) + 2):
        r, c = rows - step * dy, cols + step * dx
        valid = (r >= 0.0) & (r <= nrows - 1.001) & (c >= 0.0) & (c <= ncols - 1.001)
        if not np.any(valid):
            break
        r0 = np.clip(np.floor(r).astype(int), 0, nrows - 1)
        c0 = np.clip(np.floor(c).astype(int), 0, ncols - 1)
        fr, fc = np.clip(r - r0, 0.0, 1.0), np.clip(c - c0, 0.0, 1.0)
        r1, c1 = np.clip(r0 + 1, 0, nrows - 1), np.clip(c0 + 1, 0, ncols - 1)
        zs = (z[r0, c0] * (1 - fr) * (1 - fc) + z[r0, c1] * (1 - fr) * fc +
              z[r1, c0] * fr * (1 - fc) + z[r1, c1] * fr * fc)
        max_tan = np.maximum(max_tan, np.where(valid, (zs - z) / (step * grid.cell), 0.0))
    return max_tan > math.tan(sun_elevation)


@lru_cache(maxsize=32)
def _classes_cached(path_text, mtime, sza, saa, slope_bin, aspect_bin, shadow_enabled):
    grid = _grid_cached(path_text, mtime)
    slope, aspect = slope_aspect(grid)
    shadow = cast_shadow(grid, float(sza), float(saa)) if shadow_enabled else np.zeros_like(slope, dtype=bool)
    inner = (slice(1, -1), slice(1, -1))
    slope_values = np.floor(slope[inner] / slope_bin) * slope_bin
    aspect_values = np.floor(aspect[inner] / aspect_bin) * aspect_bin
    code = (np.rint(slope_values).astype(np.int64) * 1000 + np.rint(aspect_values).astype(np.int64))
    sun_code, sun_count = np.unique(code[~shadow[inner]], return_counts=True)
    shade_code, shade_count = np.unique(code[shadow[inner]], return_counts=True)
    counts = {}
    for value, count in zip(sun_code, sun_count):
        counts[int(value)] = [int(count), 0]
    for value, count in zip(shade_code, shade_count):
        entry = counts.setdefault(int(value), [0, 0])
        entry[1] = int(count)
    classes = []
    for value in sorted(counts):
        slope_value = (value // 1000) / 1.0
        aspect_value = (value % 1000) / 1.0
        sun, shade = counts[value]
        if sun or shade:
            classes.append((slope_value, aspect_value % 360.0, sun, shade))
    return tuple(classes)


def facet_classes(config, sza, saa):
    grid = read_grid(config)
    slope_bin = max(float(config.get("dem_slope_bin", 2.0)), 0.25)
    aspect_bin = max(float(config.get("dem_aspect_bin", 10.0)), 0.25)
    shadow_enabled = str(config.get("dem_shadow", "1")).strip().lower() not in {"0", "false", "no", "off"}
    return _classes_cached(str(grid.path), grid.path.stat().st_mtime_ns, float(sza), float(saa),
                           slope_bin, aspect_bin, shadow_enabled)


def _facet_config(config, slope, aspect, shaded, *, optical=False):
    local = dict(config, surface_model="slope", terrain_type="1",
                 slope_angle=f"{slope:.8f}", slope_aspect=f"{aspect:.8f}")
    if optical:
        local["shortwave_sunlit_surface"] = "0" if shaded else "1"
    if shaded:
        from run import _temperature
        for component, default in (("soil", 30.0), ("leaf", 30.0), ("terrain", 30.0)):
            local[f"{component}_sunlit_k"] = str(_temperature(config, f"{component}_shaded", default))
    return local


def _facet_weights(config, classes, vza, vaa):
    total = []
    for slope, aspect, sun, shade in classes:
        local_vza, _, _ = slope_direction(vza, vaa, slope, aspect)
        if local_vza >= 89.5:
            continue
        projected = math.cos(math.radians(local_vza))
        if projected <= 0.0:
            continue
        total.append((slope, aspect, sun, shade, projected))
    return total


def radiance(config, vegetation, wavelength, geometry, spectra):
    """Projected-area weighted radiance of the DEM facet ensemble."""
    from run import _surface_result, _checked_radiance

    vza, vaa, sza, saa = geometry
    classes = _facet_weights(config, facet_classes(config, sza, saa), vza, vaa)
    if not classes:
        return float("nan"), float("nan")
    total = 0.0
    weight = 0.0
    for slope, aspect, sun, shade, projected in classes:
        for shaded, count in ((False, sun), (True, shade)):
            if count == 0:
                continue
            local = _facet_config(config, slope, aspect, shaded)
            value, _ = _surface_result(local, "slope", vegetation, wavelength, geometry, spectra)
            if not np.isfinite(value):
                continue
            contribution = count * projected
            total += contribution * float(value)
            weight += contribution
    if weight <= 0.0:
        return float("nan"), float("nan")
    mean = total / weight
    return _checked_radiance(mean, wavelength, f"dem+{vegetation}", geometry)


def reflectance(config, wavelength, geometry, band):
    """Projected-area weighted BRF of the DEM facet ensemble."""
    from run import _optical_band

    vza, vaa, sza, saa = geometry
    classes = _facet_weights(config, facet_classes(config, sza, saa), vza, vaa)
    if not classes:
        return float("nan")
    total = 0.0
    weight = 0.0
    for slope, aspect, sun, shade, projected in classes:
        for shaded, count in ((False, sun), (True, shade)):
            if count == 0:
                continue
            local = _facet_config(config, slope, aspect, shaded, optical=True)
            value = _optical_band(local, "slope", str(config.get("vegetation_model", "bare")),
                                  float(wavelength), geometry, band)
            if not np.isfinite(value):
                continue
            contribution = count * projected
            total += contribution * float(value)
            weight += contribution
    return total / weight if weight > 0.0 else float("nan")


__all__ = ["Grid", "read_grid", "slope_aspect", "cast_shadow", "facet_classes", "radiance", "reflectance"]
