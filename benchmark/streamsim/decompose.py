"""Component-fraction decomposition (T3 crowns, T5 rows) with eRaytracing (blackbody surfaces, emission only).

For each leaf size, four indicator runs put one component (sunlit soil, shaded soil,
sunlit leaf, shaded leaf) at T_HOT and the rest at T_COLD; with eps=1 (no reflected
paths at any depth) and a 3 K sky, L = sum f_i B(T_i), so f_i = (L_i - B_cold) / (B_hot - B_cold).

Usage: python decompose.py [--test T3|T5az0|T5az90] [--leaf ...] [--samples 32]
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_benchmark as rb  # noqa: E402

T_HOT, T_COLD = 320.0, 280.0
COMPONENTS = ["soil_sunlit", "soil_shaded", "leaf_sunlit", "leaf_shaded"]


def planck(temperature):
    return rb.C1 / (rb.WAVELENGTH_UM ** 5 * (math.exp(rb.C2 / (rb.WAVELENGTH_UM * temperature)) - 1.0))


def indicator_temps(component):
    temps = [T_COLD] * 4
    temps[COMPONENTS.index(component)] = T_HOT
    return ",".join(f"{t:g}" for t in temps)


def crown_cover(run_dir, length, radius=3.0, cells=2000):
    """Fraction of the domain covered by >=1 crown in vertical projection (periodic copies included)."""
    centres = []
    for path in (run_dir / "input" / "positions").glob("crown*_position.txt"):
        centres += [tuple(map(float, line.split()[:2])) for line in path.read_text().splitlines() if line.strip()]
    axis = (np.arange(cells) + 0.5) * length / cells
    covered = np.zeros((cells, cells), dtype=bool)
    for x, z in centres:
        if -radius <= x <= length + radius and -radius <= z <= length + radius:
            ix = (axis >= x - radius) & (axis <= x + radius)
            iz = (axis >= z - radius) & (axis <= z + radius)
            dx, dz = np.meshgrid(axis[ix] - x, axis[iz] - z, indexing="ij")
            covered[np.ix_(ix, iz)] |= dx * dx + dz * dz <= radius * radius
    return float(covered.mean())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", default="T3", choices=["T3", "T5az0", "T5az90"])
    parser.add_argument("--leaf", nargs="*", type=float, default=None, help="default: 0.9 0.3 for T3, 0.15 for rows")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600)
    # eps=1 zeroes every reflected path, so depth>1 stays emission-only; it lets
    # camera rays spend a path step on the periodic wrap (depth 1 leaves ~0.5% zero samples).
    parser.add_argument("--depth", type=int, default=4)
    cli = parser.parse_args()
    args = argparse.Namespace(samples=cli.samples, seed=cli.seed, res=512, depth=cli.depth, timeout=cli.timeout)
    default_leaf = 0.9 if cli.test == "T3" else 0.15
    leaves = cli.leaf or ([0.9, 0.3] if cli.test == "T3" else [0.15])
    prefix = "T3" if cli.test == "T3" else rb.NAMES[cli.test]
    b_hot, b_cold = planck(T_HOT), planck(T_COLD)
    geometry = {}
    for leaf in leaves:
        suffix = "" if abs(leaf - default_leaf) < 1e-9 else f"_leaf{leaf:g}"
        radiance = {}
        for component in COMPONENTS:
            extra = ["--leaf", f"{leaf:g}", "--blackbody", "--temps", indicator_temps(component)]
            tag = f"_frac_{component}{suffix}"
            result = rb.run_case(cli.test, "eRaytracing", args, extra, tag)
            radiance[component] = {(r["vza"], r["vaa"]): r["radiance"] for r in result["rows"]}
            meta = result["meta"]
        keys = sorted(radiance[COMPONENTS[0]], key=lambda k: (k[1], k[0]))
        target = rb.RUNS / f"{prefix}_fractions_eRaytracing{suffix}.csv"
        with target.open("w", newline="", encoding="utf-8") as sink:
            writer = csv.writer(sink)
            writer.writerow(["vza", "vaa", *[f"f_{c}" for c in COMPONENTS], "sum"])
            for key in keys:
                fractions = [(radiance[c][key] - b_cold) / (b_hot - b_cold) for c in COMPONENTS]
                writer.writerow([f"{key[0]:g}", f"{key[1]:g}", *[f"{f:.5f}" for f in fractions], f"{sum(fractions):.5f}"])
        run_dir = rb.RUNS / f"{rb.NAMES[cli.test]}_eRaytracing_frac_soil_sunlit{suffix}"
        nadir = radiance["soil_sunlit"][(0.0, 0.0)], radiance["soil_shaded"][(0.0, 0.0)]
        info = {"realised_lai": meta["leafAreaIndex"], "leaf_size": meta["leafSize"],
                "nadir_gap_fraction": sum((value - b_cold) / (b_hot - b_cold) for value in nadir)}
        if cli.test == "T3":
            info.update({"leaves_per_crown": meta["leavesPerCrown"], "crowns": meta["crowns"],
                         "nominal_crown_cover_1_minus_exp": 1 - math.exp(-meta["crowns"] / meta["L"] ** 2 * math.pi * 9),
                         "crown_cover_fraction": crown_cover(run_dir, meta["L"])})
        else:
            info.update({"row_azimuth": meta["rowAzimuth"], "rows": meta["rows"], "in_row_lad": meta["inRowLad"],
                         "row_cover_fraction": meta["rowCoverEnvelope"], "leaf_pieces": meta["leafPieces"]})
        geometry[f"leaf{leaf:g}"] = info
        print(f"leaf {leaf:g}: {target.name}", json.dumps(info), flush=True)
    (rb.RUNS / f"{prefix}_fractions_geometry.json").write_text(json.dumps(geometry, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
