"""Forest (T3) component fractions: TiRT crown model vs HiSTREAM eRaytracing decomposition, both principal planes.

Usage: python crown_fractions_plot.py [--old previous_TiRT_fractions.csv]
Writes runs/tirt/T3_fractions_TiRT.csv and runs/compare/crown_fractions.png.
"""
import argparse
import csv
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = HERE / "runs"
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

from base.proportion import proportion_bidirectional_crown_one  # noqa: E402

LAI, DENSITY, RADIUS, TRUNK, SZA, CROWN_DEPTH = 1.5, 0.02, 3.0, 3.0, 30.0, 6.0
VZA = np.arange(0, 76, 5)
COMPONENTS = ["soil_sunlit", "soil_shaded", "leaf_sunlit", "leaf_shaded"]
PANELS = [("土壤总比例", ("soil_sunlit", "soil_shaded")), ("光照土壤", ("soil_sunlit",)), ("阴影土壤", ("soil_shaded",)),
          ("光照植被", ("leaf_sunlit",)), ("阴影植被", ("leaf_shaded",))]
PLANES = {"太阳主平面（+：vaa 0，太阳一侧）": (0.0, 180.0), "垂直主平面（+：vaa 90）": (90.0, 270.0)}
LEAVES = {0.9: ("", "tab:green"), 0.3: ("_leaf0.3", "tab:blue"), 0.15: ("_leaf0.15", "tab:purple")}


def read_fractions(path):
    with path.open(encoding="utf-8") as file:
        return {(float(r["vza"]), float(r["vaa"]) % 360): {c: float(r[f"f_{c}"]) for c in COMPONENTS} for r in csv.DictReader(file)}


def tirt_fractions(hspot):
    table = {}
    for vaa in (0.0, 90.0, 180.0, 270.0):
        parts = proportion_bidirectional_crown_one(LAI, DENSITY, hspot, RADIUS, RADIUS, np.maximum(VZA, 1.0e-3), np.full(VZA.size, SZA),
                                                   np.full(VZA.size, vaa), hcb=TRUNK)
        for index, vza in enumerate(VZA):
            table[(float(vza), vaa)] = {c: float(parts[k][index]) for k, c in enumerate(COMPONENTS)}
    return table


def curve(table, keys, positive, negative):
    points = {}
    for (vza, vaa), values in table.items():
        value = sum(values[k] for k in keys)
        if np.isclose(vaa, positive) or (vza == 0 and np.isclose(vaa, negative)):
            points[vza] = value
        elif np.isclose(vaa, negative):
            points[-vza] = value
    x = np.array(sorted(points))
    return x, np.array([points[v] for v in x])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, help="previous TiRT fractions CSV to overlay")
    cli = parser.parse_args()
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    series = []
    for leaf, (suffix, color) in LEAVES.items():
        path = RUNS / f"T3_fractions_eRaytracing{suffix}.csv"
        if path.exists():
            series.append((f"HiSTREAM 叶片 {leaf:g} m", read_fractions(path), color, "o-"))
            series.append((f"TiRT hspot {leaf / CROWN_DEPTH:.3g}", tirt_fractions(leaf / CROWN_DEPTH), color, "-"))
    if cli.old:
        series.append(("TiRT 旧模型 hspot 0.15", read_fractions(cli.old), "0.4", "--"))

    target = RUNS / "tirt" / "T3_fractions_TiRT.csv"
    current = tirt_fractions(0.15)
    with target.open("w", newline="", encoding="utf-8") as sink:
        writer = csv.writer(sink)
        writer.writerow(["vza", "vaa", *[f"f_{c}" for c in COMPONENTS], "sum"])
        for (vza, vaa), values in sorted(current.items(), key=lambda item: (item[0][0], item[0][1])):
            writer.writerow([f"{vza:g}", f"{vaa:g}", *[f"{values[c]:.5f}" for c in COMPONENTS], f"{sum(values.values()):.5f}"])

    fig, axes = plt.subplots(len(PLANES), len(PANELS), figsize=(22, 8.5), squeeze=False)
    for i, (plane, (positive, negative)) in enumerate(PLANES.items()):
        for j, (title, keys) in enumerate(PANELS):
            ax = axes[i][j]
            for label, table, color, style in series:
                x, y = curve(table, keys, positive, negative)
                ax.plot(x, y, style, color=color, ms=3, lw=1.3, label=label)
            ax.set_title(f"{title} — {plane}", fontsize=9)
            ax.set_xlabel("带符号观测天顶角 (°)")
            ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=8)
    fig.suptitle(f"离散树冠组分比例（LAI {LAI}，树冠半径 {RADIUS:g} m，密度 {DENSITY}，树干 {TRUNK:g} m，SZA {SZA:g}°）", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = RUNS / "compare" / "crown_fractions.png"
    fig.savefig(out, dpi=110)
    print(target)
    print(out)


if __name__ == "__main__":
    main()
