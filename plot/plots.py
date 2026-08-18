"""Polar and solar-principal-plane plots for TiRT output.csv."""

import csv
import math
from pathlib import Path

import numpy as np


def _matplotlib(show=False):
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def load_output(input_path, wavelength=None, quantity="brightness_temperature_C"):
    """Load one wavelength from the unified runner output.

    If wavelength is omitted, the first wavelength in the file is selected.
    The returned value is a list of numeric row dictionaries.
    """
    input_path = Path(input_path)
    with input_path.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise ValueError(f"{input_path} 没有可绘图的数据")
    fields = set(rows[0])
    required = {"vza", "vaa", "sza", "saa", "wavelength_um", quantity}
    missing = sorted(required - fields)
    if missing:
        raise ValueError(f"output.csv 缺少绘图字段: {', '.join(missing)}")
    wavelengths = np.asarray([float(row["wavelength_um"]) for row in rows])
    selected = float(np.unique(wavelengths)[0] if wavelength is None else wavelength)
    mask = np.isclose(wavelengths, selected, rtol=0.0, atol=1.0e-7)
    if not np.any(mask):
        available = ", ".join(f"{item:g}" for item in sorted(set(wavelengths)))
        raise ValueError(f"波段 {selected:g} 不存在，可选波段: {available}")
    result = []
    for row, keep in zip(rows, mask):
        if not keep:
            continue
        result.append({
            **row,
            "wavelength_um": float(row["wavelength_um"]),
            "vza": float(row["vza"]),
            "vaa": float(row["vaa"]),
            "sza": float(row["sza"]),
            "saa": float(row["saa"]),
            "value": float(row[quantity]),
        })
    return result, selected, quantity


def _save(fig, output_path, show):
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=240, bbox_inches="tight")
    if show:
        fig.show()
    return output_path


def _label(quantity):
    labels = {
        "brightness_temperature_C": "Brightness temperature (C)",
        "brightness_temperature_K": "Brightness temperature (K)",
        "radiance": "Radiance",
    }
    return labels.get(quantity, quantity)


def plot_polar(rows, output_path=None, quantity="brightness_temperature_C", title=None, show=False):
    """Draw a dense, periodically interpolated view-hemisphere contour map."""
    plt = _matplotlib(show)
    zeniths = np.sort(np.unique([float(row["vza"]) for row in rows]))
    azimuths = np.sort(np.unique([float(row["vaa"]) % 360.0 for row in rows]))
    azimuths = azimuths[azimuths < 360.0]
    if zeniths.size < 2 or azimuths.size < 2:
        raise ValueError("极坐标填充至少需要两个天顶角和两个方位角")

    lookup = {
        (float(row["vza"]), float(row["vaa"]) % 360.0): float(row["value"])
        for row in rows
    }
    values = np.asarray(
        [[lookup.get((zenith, azimuth), np.nan) for azimuth in azimuths] for zenith in zeniths],
        dtype=float,
    )
    if not np.any(np.isfinite(values)):
        raise ValueError("极坐标数据没有有限值")

    # Interpolate around the periodic azimuth axis first, then along VZA.
    fine_azimuths = np.linspace(0.0, 360.0, 361)
    fine_zeniths = np.linspace(float(zeniths.min()), float(zeniths.max()), 121)
    azimuth_grid = np.deg2rad(fine_azimuths)
    radial_values = np.empty((zeniths.size, fine_azimuths.size), dtype=float)
    for index, row in enumerate(values):
        valid = np.isfinite(row)
        if not np.any(valid):
            radial_values[index] = np.nan
            continue
        x = azimuths[valid]
        y = row[valid]
        x_extended = np.r_[x, x[0] + 360.0]
        y_extended = np.r_[y, y[0]]
        radial_values[index] = np.interp(fine_azimuths, x_extended, y_extended)

    polar_values = np.empty((fine_zeniths.size, fine_azimuths.size), dtype=float)
    for index in range(fine_azimuths.size):
        valid = np.isfinite(radial_values[:, index])
        if not np.any(valid):
            polar_values[:, index] = np.nan
            continue
        polar_values[:, index] = np.interp(
            fine_zeniths, zeniths[valid], radial_values[valid, index]
        )
    polar_values[0, :] = polar_values[0, 0]

    fig, ax = plt.subplots(figsize=(7.2, 6.4), subplot_kw={"projection": "polar"})
    finite = polar_values[np.isfinite(polar_values)]
    value_min, value_max = float(finite.min()), float(finite.max())
    levels = (
        np.linspace(value_min - 1e-6, value_max + 1e-6, 20)
        if math.isclose(value_min, value_max, rel_tol=0.0, abs_tol=1e-12)
        else np.linspace(value_min, value_max, 60)
    )
    contour = ax.contourf(
        np.tile(azimuth_grid, (fine_zeniths.size, 1)),
        np.tile(fine_zeniths[:, None], (1, fine_azimuths.size)),
        polar_values,
        levels=levels,
        cmap="turbo",
        extend="both",
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rmax(max(float(zeniths.max()), 1.0))
    ax.set_rticks(zeniths.tolist())
    ax.set_rlabel_position(22.5)
    ax.set_title(title or "TiRT directional response")
    fig.colorbar(contour, ax=ax, pad=0.12, label=_label(quantity))
    fig.tight_layout()
    result = _save(fig, output_path, show)
    plt.close(fig)
    return result


def _principal_plane(rows, plane, tolerance):
    delta = np.asarray([(row["vaa"] - row["saa"]) % 360.0 for row in rows])
    if plane == "parallel":
        mask = (np.minimum(delta, 360.0 - delta) <= tolerance) | (np.abs(delta - 180.0) <= tolerance)
        sign = np.where(np.abs(delta - 180.0) <= tolerance, -1.0, 1.0)
    else:
        mask = (np.abs(delta - 90.0) <= tolerance) | (np.abs(delta - 270.0) <= tolerance)
        sign = np.where(np.abs(delta - 270.0) <= tolerance, -1.0, 1.0)
    if not np.any(mask):
        raise ValueError(f"没有找到太阳主平面{plane}方向，请检查 vaa、saa 或 angle_tolerance")
    selected = [row for row, keep in zip(rows, mask) if keep]
    x = np.asarray([row["vza"] for row in selected]) * sign[mask]
    y = np.asarray([row["value"] for row in selected])
    order = np.argsort(x)
    return x[order], y[order]


def _plot_plane(rows, plane, output_path=None, quantity="brightness_temperature_C", title=None, tolerance=1.0e-6, show=False):
    plt = _matplotlib(show)
    x, y = _principal_plane(rows, plane, tolerance)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(x, y, "o-", linewidth=1.4, markersize=4)
    ax.axvline(0.0, color="0.5", linewidth=0.8)
    ax.set_xlabel("Signed view zenith angle (degree)")
    ax.set_ylabel(_label(quantity))
    ax.grid(True, alpha=0.3)
    ax.set_title(title or f"TiRT solar principal-plane {plane}")
    fig.tight_layout()
    result = _save(fig, output_path, show)
    plt.close(fig)
    return result


def plot_parallel(rows, output_path=None, quantity="brightness_temperature_C", title=None, tolerance=1.0e-6, show=False):
    return _plot_plane(rows, "parallel", output_path, quantity, title, tolerance, show)


def plot_perpendicular(rows, output_path=None, quantity="brightness_temperature_C", title=None, tolerance=1.0e-6, show=False):
    return _plot_plane(rows, "perpendicular", output_path, quantity, title, tolerance, show)


def plot_all(input_path, output_dir=None, wavelength=None, quantity="brightness_temperature_C", tolerance=1.0e-6, show=False):
    rows, selected, quantity = load_output(input_path, wavelength, quantity)
    output_dir = Path(output_dir) if output_dir is not None else Path(input_path).parent / "plots"
    tag = f"{selected:g}".replace(".", "p")
    result = {
        "polar": plot_polar(rows, output_dir / f"polar_{tag}.png", quantity, show=show),
    }
    for name, function in (("parallel", plot_parallel), ("perpendicular", plot_perpendicular)):
        try:
            result[name] = function(rows, output_dir / f"{name}_{tag}.png", quantity, tolerance=tolerance, show=show)
        except ValueError:
            result[name] = None
    return result
