"""Summarise StreamSim benchmark CSVs: noise estimates, mode spread and TiRT differences."""

import csv
import math
import sys
from pathlib import Path

RUNS = Path(__file__).resolve().parent / "runs"
TESTS = ["T1_plane_bare", "T2_plane_hom", "T3_plane_crown", "T4_urban_bare", "T5_row_az0", "T5_row_az90",
         "S1_slope_hom", "S1_slope_crown", "TV_terrain_crown", "UV_urban_crown"]


def load(path):
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as source:
        rows = {(float(r["vza"]), float(r["vaa"])): float(r["bt_K"]) for r in csv.DictReader(source)}
    # TiRT reports NaN where a slope faces away from the sensor.
    return {key: value for key, value in rows.items() if math.isfinite(value)}


def diff_stats(a, b, keys=None):
    keys = [k for k in (keys or a) if k in a and k in b]
    values = [a[k] - b[k] for k in keys]
    if not values:
        return "n/a"
    rms = (sum(v * v for v in values) / len(values)) ** 0.5
    return f"mean {sum(values) / len(values):+.3f}  rms {rms:.3f}  max|d| {max(abs(v) for v in values):.3f} K"


def symmetry(a):
    pairs = [((z, 90.0), (z, 270.0)) for z in range(5, 80, 5)]
    values = [abs(a[p] - a[q]) for p, q in pairs if p in a and q in a]
    return f"mean {sum(values) / len(values):.3f}  max {max(values):.3f} K" if values else "n/a"


def combine(test):
    """runs/<test>_streamsim.csv: all main-mode rows of one test in TiRT's column layout."""
    rows = []
    for mode in ("eRaytracing", "eVoxelRT"):
        path = RUNS / f"{test}_{mode}.csv"
        if path.exists():
            with path.open(encoding="utf-8") as source:
                rows.extend(csv.DictReader(source))
    if rows:
        with (RUNS / f"{test}_streamsim.csv").open("w", newline="", encoding="utf-8") as sink:
            writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    out = []
    for test in TESTS:
        combine(test)
        tirt = load(RUNS / "tirt" / f"{test}.csv")
        out.append(f"## {test}")
        for mode in ("eRaytracing", "eVoxelRT"):
            base = load(RUNS / f"{test}_{mode}.csv")
            if base is None:
                continue
            out.append(f"- {mode}: nadir {base[(0.0, 0.0)]:.3f} K, hotspot(30,0) {base[(30.0, 0.0)]:.3f} K, "
                       f"(75,180) {base[(75.0, 180.0)]:.3f} K")
            out.append(f"  - symmetry |BT(vaa=90)-BT(vaa=270)|: {symmetry(base)}")
            for tag, label in (("_seed2", "scene seed 2 - seed 1"), ("_spp128", "128 spp - 32 spp"),
                               ("_depth1", "depth 1 (no scattering) - depth 4"), ("_leaf0.3", "leaf 0.3 m - default")):
                other = load(RUNS / f"{test}_{mode}{tag}.csv")
                if other:
                    out.append(f"  - {label}: {diff_stats(other, base)}")
            if tirt:
                out.append(f"  - {mode} - TiRT: {diff_stats(base, tirt)}")
        rt, vx = load(RUNS / f"{test}_eRaytracing.csv"), load(RUNS / f"{test}_eVoxelRT.csv")
        if rt and vx:
            out.append(f"- eVoxelRT - eRaytracing: {diff_stats(vx, rt)}")
        out.append("")
    text = "\n".join(out)
    (RUNS / "summary.md").write_text(text, encoding="utf-8")
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
