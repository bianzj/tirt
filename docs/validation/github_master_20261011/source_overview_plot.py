"""Validation figures: for every benchmarked scene, HiSTREAM eRaytracing (3D path tracing), a 1D homogeneous-layer run
and TiRT, drawn as scene sketch | solar principal plane | perpendicular plane (both planes on one brightness-temperature
axis), plus one statistics table.

The 1D homogeneous run is TiRT with the structure removed: flat ground (no slope, terrain or buildings) and, if the scene
has vegetation, a horizontally homogeneous layer with the same LAI, emissivities and temperatures.

Usage: python overview_plot.py   (after tirt_side.py and the eRaytracing runs)
Writes docs/validation/<nn>_<scene>.png and docs/validation/overview.csv.
"""
import csv
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "docs" / "validation"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
warnings.simplefilter("ignore")

import compare  # noqa: E402
import scene_sketch  # noqa: E402
import tirt_side  # noqa: E402
from run import run  # noqa: E402

UV_LEAF = {"leaf_sunlit_c": 31.5, "leaf_shaded_c": 31.5}


def tirt(test, **extra):
    path = run(ROOT / "input.csv", Path(tempfile.mkdtemp(prefix="tirt-ov-")) / "o.csv", **{**tirt_side.COMMON, **tirt_side.TESTS[test], **extra})
    with path.open(encoding="utf-8") as file:
        return {(round(float(r["vza"]), 3), round(float(r["vaa"]) % 360, 3)): float(r["brightness_temperature_K"]) for r in csv.DictReader(file)}


def homogeneous(test, **extra):
    vegetation = "bare" if tirt_side.TESTS[test]["vegetation_model"] == "bare" else "hom"
    return tirt(test, surface_model="plane", vegetation_model=vegetation, hspot=tirt_side.COMMON["hspot"], **extra)


def errors(model, ref):
    keys = sorted(k for k in set(ref) & set(model) if np.isfinite(model[k]))
    d = np.array([model[k] - ref[k] for k in keys])
    return keys, d


def main():
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    data = compare.load()
    # (file stem, scene label, test, HiSTREAM reference, TiRT, 1D homogeneous)
    scenes = [
        ("01_T1_plane_bare", "平地裸土", "T1_plane_bare", data["T1_plane_bare"]["eRaytracing"], data["T1_plane_bare"]["TiRT"], homogeneous("T1_plane_bare")),
        ("02_T2_plane_hom", "平地均质冠层（LAI 1.5）", "T2_plane_hom", data["T2_plane_hom"]["eRaytracing"], data["T2_plane_hom"]["TiRT"], homogeneous("T2_plane_hom")),
        ("03_T3_plane_crown_leaf0.9", "平地离散森林（叶片 0.9 m）", "T3_plane_crown", data["T3_plane_crown"]["eRaytracing"], data["T3_plane_crown"]["TiRT"], homogeneous("T3_plane_crown")),
        ("04_T3_plane_crown_leaf0.15", "平地离散森林（叶片 0.15 m）", "T3_plane_crown", data["T3_plane_crown"]["eRaytracing_leaf0.15"], tirt("T3_plane_crown", hspot=0.025), homogeneous("T3_plane_crown")),
        ("05_T5_row_az0", "垄行作物（垄向 0°）", "T5_row_az0", data["T5_row_az0"]["eRaytracing"], data["T5_row_az0"]["TiRT"], homogeneous("T5_row_az0")),
        ("06_T5_row_az90", "垄行作物（垄向 90°）", "T5_row_az90", data["T5_row_az90"]["eRaytracing"], data["T5_row_az90"]["TiRT"], homogeneous("T5_row_az90")),
        ("07_S1_slope_hom", "单坡 30° 朝东 + 均质冠层", "S1_slope_hom", data["S1_slope_hom"]["eRaytracing"], data["S1_slope_hom"]["TiRT"], homogeneous("S1_slope_hom")),
        ("08_S1_slope_crown", "单坡 30° 朝东 + 森林（叶片 0.15 m）", "S1_slope_crown", data["S1_slope_crown"]["eRaytracing"], data["S1_slope_crown"]["TiRT"], homogeneous("S1_slope_crown")),
        ("09_TV_terrain_crown", "复合地形 + 森林（叶片 0.15 m）", "TV_terrain_crown", data["TV_terrain_crown"]["eRaytracing"], data["TV_terrain_crown"]["TiRT"], homogeneous("TV_terrain_crown")),
        ("10_T4_urban_bare", "城市（无植被）", "T4_urban_bare", data["T4_urban_bare"]["eRaytracing"], data["T4_urban_bare"]["TiRT"], homogeneous("T4_urban_bare")),
        ("11_UV_urban_crown", "城市 + 森林（叶片 0.15 m，31.5 °C）", "UV_urban_crown", data["UV_urban_crown"]["eRaytracing"], data["UV_urban_crown"]["TiRT"], homogeneous("UV_urban_crown", **UV_LEAF)),
    ]
    rows = []
    OUT.mkdir(parents=True, exist_ok=True)
    for stem, name, test, ref, model, flat in scenes:
        keys, d = errors(model, ref)
        _, h = errors(flat, ref)
        r = np.array([ref[k] for k in keys])
        rows.append({"scene": name, "test": test, "figure": f"{stem}.png", "n": len(keys),
                     "tirt_bias_K": d.mean(), "tirt_rmse_K": np.sqrt((d * d).mean()), "tirt_max_abs_K": np.abs(d).max(),
                     "hom1d_bias_K": h.mean(), "hom1d_rmse_K": np.sqrt((h * h).mean()), "hom1d_max_abs_K": np.abs(h).max(),
                     "range_ref_K": r.max() - r.min(), "nadir_ref_C": ref[(0.0, 0.0)] - 273.15, "nadir_tirt_C": model[(0.0, 0.0)] - 273.15,
                     "nadir_hom1d_C": flat[(0.0, 0.0)] - 273.15})
        row = rows[-1]
        series = [
            (ref, "o-", {"color": "tab:orange", "ms": 3, "lw": 1.0}, "HiSTREAM 三维路径追踪"),
            (flat, "--", {"color": "0.4", "lw": 1.4}, "一维均质假设"),
            (model, "-", {"color": "tab:blue", "lw": 1.8}, "TiRT"),
        ]
        note = (f"相对 HiSTREAM 的 RMSE：TiRT {row['tirt_rmse_K']:.2f} K（偏差 {row['tirt_bias_K']:+.2f} K），"
                f"一维均质 {row['hom1d_rmse_K']:.2f} K（偏差 {row['hom1d_bias_K']:+.2f} K）")
        if all(np.isclose(flat[k], model[k], atol=1e-6) for k in keys):
            note += "；本场景无三维结构，一维均质与 TiRT 重合"
        scene_sketch.scene_figure(test, name, series, OUT / f"{stem}.png", note=note)
    with (OUT / "overview.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in row.items()} for row in rows)
    for row in rows:
        print(f"{row['test']:18s} TiRT rmse {row['tirt_rmse_K']:.2f} bias {row['tirt_bias_K']:+.2f}   "
              f"1D hom rmse {row['hom1d_rmse_K']:.2f} bias {row['hom1d_bias_K']:+.2f}")
    print(OUT)


if __name__ == "__main__":
    main()
