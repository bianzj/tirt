"""Plot directional observations as polar maps and line plots.

Examples:
    python analysis/plot_observations.py
    python analysis/plot_observations.py --time 172.0 --bands 0p55 0p86 10p5
    python analysis/plot_observations.py --result-dir cases/homogeneous_multispectral_rt
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULT_DIR = ROOT / "cases" / "hom_urban"
ANALYSIS_DIR = Path(__file__).resolve().parent


def _number(value: str) -> float:
    return float(value.strip())


def read_observations(path: Path) -> list[dict[str, float | str]]:
    """Read observations.csv and convert numeric columns to floats."""
    rows: list[dict[str, float | str]] = []
    with path.open("r", newline="", encoding="utf-8-sig") as file:
        for row in csv.DictReader(file):
            converted: dict[str, float | str] = {}
            for key, value in row.items():
                if value is None:
                    continue
                try:
                    converted[key] = _number(value)
                except ValueError:
                    converted[key] = value
            rows.append(converted)
    if not rows:
        raise ValueError(f"没有读取到观测数据: {path}")
    return rows


def observation_fields(rows: list[dict[str, float | str]]) -> list[str]:
    """Return observation value fields in their CSV order."""
    fields = list(rows[0])
    metadata = {"node", "time", "veg", "surface", "vza", "vaa", "sza", "saa"}
    return [field for field in fields if field not in metadata]


def field_label(field: str) -> str:
    return re.sub("_+", " ", field).replace(" temperature K", " temperature (K)")


def field_stem(field: str) -> str:
    if field.endswith("_temperature_K"):
        return field[: -len("_temperature_K")]
    return field.replace("_reflectance", "")


def time_tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def select_rows(rows: list[dict[str, float | str]], selected_time: float | None) -> list[float]:
    times = sorted({float(row["time"]) for row in rows})
    if selected_time is None:
        return times
    nearest = min(times, key=lambda value: abs(value - selected_time))
    if not math.isclose(nearest, selected_time, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError(f"时间 {selected_time} 不在 observations.csv 中，可选时间: {times}")
    return [nearest]


def rows_at_time(rows: list[dict[str, float | str]], time_value: float) -> list[dict[str, float | str]]:
    return [
        row for row in rows
        if math.isclose(float(row["time"]), time_value, rel_tol=0.0, abs_tol=1e-8)
    ]


def plot_polar(rows: list[dict[str, float | str]], field: str, time_value: float, output: Path) -> None:
    zeniths = np.sort(np.unique([float(row["vza"]) for row in rows]))
    azimuths = np.sort(np.unique([float(row["vaa"]) % 360.0 for row in rows]))
    azimuths = azimuths[azimuths < 360.0]
    if zeniths.size < 2 or azimuths.size < 2:
        return

    lookup = {
        (float(row["vza"]), float(row["vaa"]) % 360.0): float(row[field])
        for row in rows
    }
    values = np.asarray(
        [[lookup.get((zenith, azimuth), np.nan) for azimuth in azimuths] for zenith in zeniths],
        dtype=float,
    )
    if not np.any(np.isfinite(values)):
        return

    # Interpolate onto a dense rectangular (azimuth, zenith) grid before
    # projecting it into polar coordinates. The azimuth is periodic.
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
        polar_values[:, index] = np.interp(
            fine_zeniths, zeniths, radial_values[:, index]
        )
    polar_values[0, :] = polar_values[0, 0]

    fig, ax = plt.subplots(figsize=(7.2, 6.4), subplot_kw={"projection": "polar"})
    finite = polar_values[np.isfinite(polar_values)]
    value_min, value_max = float(finite.min()), float(finite.max())
    if math.isclose(value_min, value_max, rel_tol=0.0, abs_tol=1e-12):
        levels = np.linspace(value_min - 1e-6, value_max + 1e-6, 20)
    else:
        levels = np.linspace(value_min, value_max, 60)
    contour = ax.contourf(
        np.tile(azimuth_grid, (fine_zeniths.size, 1)),
        np.tile(fine_zeniths[:, None], (1, fine_azimuths.size)),
        polar_values,
        levels=levels,
        cmap="viridis",
        extend="both",
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rmax(max(float(zeniths.max()), 1.0))
    ax.set_rticks(zeniths.tolist())
    ax.set_rlabel_position(22.5)
    ax.set_title(f"{field_label(field)}\ntime = {time_value:g}")
    colorbar = fig.colorbar(contour, ax=ax, pad=0.12)
    colorbar.set_label(field_label(field))
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_lines(rows: list[dict[str, float | str]], field: str, time_value: float, output: Path) -> None:
    grouped: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["vza"])].append(row)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for vza in sorted(grouped):
        values = sorted(grouped[vza], key=lambda row: float(row["vaa"]))
        x = np.asarray([float(row["vaa"]) for row in values])
        y = np.asarray([float(row[field]) for row in values])
        valid = np.isfinite(y)
        if np.any(valid):
            ax.plot(x[valid], y[valid], marker="o", ms=2.8, lw=1.2, label=f"vza={vza:g} deg")
    ax.set_xlabel("Viewing azimuth angle (deg)")
    ax.set_ylabel(field_label(field))
    ax.set_xlim(0.0, 360.0)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{field_label(field)} vs viewing azimuth\ntime = {time_value:g}")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_zenith_lines(rows: list[dict[str, float | str]], field: str, time_value: float, output: Path) -> None:
    grouped: dict[float, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[float(row["vaa"])].append(row)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for vaa in sorted(grouped):
        values = sorted(grouped[vaa], key=lambda row: float(row["vza"]))
        x = np.asarray([float(row["vza"]) for row in values])
        y = np.asarray([float(row[field]) for row in values])
        valid = np.isfinite(y)
        if np.any(valid):
            ax.plot(x[valid], y[valid], marker="o", ms=2.8, lw=1.0, label=f"vaa={vaa:g} deg")
    ax.set_xlabel("Viewing zenith angle (deg)")
    ax.set_ylabel(field_label(field))
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{field_label(field)} vs viewing zenith\ntime = {time_value:g}")
    ax.legend(ncol=3, fontsize=7)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot directional observations.csv results")
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--output-dir", type=Path, help="direct output directory; default is analysis/<case_name>")
    parser.add_argument("--time", type=float, help="only plot the nearest available time")
    parser.add_argument(
        "--bands",
        nargs="+",
        help="wavelength tags, for example 0p55 0p86 10p5",
    )
    args = parser.parse_args()

    observations_path = args.result_dir / "observations.csv"
    rows = read_observations(observations_path)
    fields = observation_fields(rows)
    if args.bands:
        requested = set(args.bands)
        fields = [field for field in fields if field_stem(field) in requested]
        if not fields:
            raise ValueError(f"没有匹配的波段: {sorted(requested)}")

    case_name = safe_name(args.result_dir.resolve().name)
    output_root = args.output_dir.resolve() if args.output_dir else ANALYSIS_DIR / case_name
    output_root.mkdir(parents=True, exist_ok=True)

    for time_value in select_rows(rows, args.time):
        time_rows = rows_at_time(rows, time_value)
        tag = time_tag(time_value)
        for field in fields:
            stem = safe_name(field)
            plot_polar(time_rows, field, time_value, output_root / f"time_{tag}_{stem}.png")
            # Line plots are kept available for later use but are disabled for now.

    print(f"读取: {observations_path}")
    print(f"时间数: {len(select_rows(rows, args.time))}, 波段结果数: {len(fields)}")
    print(f"输出: {output_root}")


if __name__ == "__main__":
    main()
