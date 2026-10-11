"""Compare TiRT with StreamSim/HiSTREAM outputs: per-test principal/perpendicular-plane plots and statistics."""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
OUT = RUNS / "compare"
PLANES = {"principal (+: vaa 0, sun side)": (0.0, 180.0), "perpendicular (+: vaa 90)": (90.0, 270.0)}


def load():
    """Series label = file-name suffix after the test name (variants such as _seed2 stay separate)."""
    data = defaultdict(lambda: defaultdict(dict))
    for path in RUNS.rglob("*.csv"):
        if OUT in path.parents or "_frac_" in path.name or "_fractions_" in path.name or path.name.endswith("_streamsim.csv"):
            continue
        with path.open(encoding="utf-8-sig") as file:
            reader = csv.DictReader(file)
            if not reader.fieldnames or "bt_K" not in reader.fieldnames:
                continue
            for row in reader:
                test = row["test"]
                label = "TiRT" if row["mode"] == "TiRT" else path.stem[len(test) + 1:] or row["mode"]
                key = (round(float(row["vza"]), 3), round(float(row["vaa"]) % 360.0, 3))
                data[test][label][key] = float(row["bt_K"])
    return data


def signed_curve(values, positive, negative):
    points = []
    for (vza, vaa), bt in values.items():
        if np.isclose(vaa, positive) or (vza == 0 and np.isclose(vaa, negative)):
            points.append((vza, bt))
        elif np.isclose(vaa, negative):
            points.append((-vza, bt))
    points = sorted(dict(points).items())
    return np.array([p[0] for p in points]), np.array([p[1] for p in points])


def main():
    data = load()
    OUT.mkdir(parents=True, exist_ok=True)
    tests = sorted(data)
    fig, axes = plt.subplots(len(tests), 2, figsize=(11, 3.2 * len(tests)), squeeze=False)
    summary = []
    for i, test in enumerate(tests):
        modes = data[test]
        reference = next((m for m in ("eRaytracing", "eFacetRT", "eVoxelRT") if m in modes), None)
        for j, (title, (positive, negative)) in enumerate(PLANES.items()):
            ax = axes[i][j]
            for mode, values in sorted(modes.items()):
                if "_" in mode:
                    continue
                x, y = signed_curve(values, positive, negative)
                ax.plot(x, y - 273.15, marker="o" if mode != "TiRT" else None, ms=3, lw=1.6 if mode == "TiRT" else 1.0, label=mode)
            ax.set_title(f"{test} — {title}", fontsize=9)
            ax.set_xlabel("signed VZA (deg)")
            ax.set_ylabel("BT (°C)")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        if reference is None or "TiRT" not in modes:
            continue
        for mode in modes:
            if mode == reference:
                continue
            common = sorted(set(modes[mode]) & set(modes[reference]))
            d = np.array([modes[mode][k] - modes[reference][k] for k in common])
            ref = np.array([modes[reference][k] for k in common])
            own = np.array([modes[mode][k] for k in common])
            summary.append({
                "test": test, "mode": mode, "reference": reference, "n": len(common),
                "bias_K": float(np.mean(d)), "rmse_K": float(np.sqrt(np.mean(d * d))), "max_abs_K": float(np.max(np.abs(d))),
                "range_mode_K": float(own.max() - own.min()), "range_ref_K": float(ref.max() - ref.min()),
                "nadir_mode_K": modes[mode].get((0.0, 0.0), np.nan), "nadir_ref_K": modes[reference].get((0.0, 0.0), np.nan),
            })
    fig.tight_layout()
    fig.savefig(OUT / "comparison.png", dpi=130)
    with (OUT / "summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(summary[0]) if summary else ["test"])
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(f"{row['test']:16s} {row['mode']:10s} vs {row['reference']:11s} n={row['n']:3d} "
              f"bias {row['bias_K']:+6.2f} K  rmse {row['rmse_K']:5.2f}  max {row['max_abs_K']:5.2f}  "
              f"angular range {row['range_mode_K']:5.2f} vs {row['range_ref_K']:5.2f}  "
              f"nadir {row['nadir_mode_K'] - 273.15:6.2f} vs {row['nadir_ref_K'] - 273.15:6.2f} °C")
    print(f"plot: {OUT / 'comparison.png'}")


if __name__ == "__main__":
    sys.exit(main())
