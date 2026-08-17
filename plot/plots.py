"""Polar and solar-principal-plane plots for TiRT output.csv."""

import csv
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
    """Draw a view-hemisphere polar scatter plot."""
    plt = _matplotlib(show)
    theta = np.radians([row["vaa"] for row in rows])
    radius = np.asarray([row["vza"] for row in rows])
    values = np.asarray([row["value"] for row in rows])
    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    points = ax.scatter(theta, radius, c=values, cmap="turbo", s=34, edgecolors="none")
    ax.set_ylim(0.0, max(5.0, float(np.max(radius)) * 1.05))
    ax.set_title(title or "TiRT directional response")
    fig.colorbar(points, ax=ax, pad=0.1, label=_label(quantity))
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
