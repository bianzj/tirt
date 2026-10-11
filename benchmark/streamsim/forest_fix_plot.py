"""Forest + compound terrain (TV) or urban (UV): TiRT before/after the forest-background scattering step vs HiSTREAM eRaytracing.

Usage: python forest_fix_plot.py [TV|UV ...]   (default: both)
Writes runs/compare/terrain_forest_fix.png and runs/compare/urban_forest_fix.png.
"Before" is the current model with the forest-background term switched off.
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
warnings.simplefilter("ignore")

import base.terrain_veg as terrain_veg  # noqa: E402
import base.urban_veg as urban_veg  # noqa: E402
import compare  # noqa: E402
import scene_sketch  # noqa: E402
import tirt_side  # noqa: E402
from run import run  # noqa: E402

SCENES = {
    "TV": ("TV_terrain_crown", terrain_veg, "复合地形 + 森林：森林–背景散射修正前后", "terrain_forest_fix.png"),
    "UV": ("UV_urban_crown", urban_veg, "城市 + 森林（叶片 31.5 °C）：森林–背景散射修正前后", "urban_forest_fix.png"),
}


def tirt(test, module, forest_background):
    original = module.emissivity_scatter_analytical
    if not forest_background:
        module.emissivity_scatter_analytical = lambda lai, el, es, vza, sza: (np.zeros((4, np.size(vza))), np.zeros(np.size(vza)))
    try:
        path = run(ROOT / "input.csv", Path(tempfile.mkdtemp(prefix="tirt-fix-")) / "o.csv", **{**tirt_side.COMMON, **tirt_side.TESTS[test]})
    finally:
        module.emissivity_scatter_analytical = original
    with path.open(encoding="utf-8") as file:
        return {(round(float(r["vza"]), 3), round(float(r["vaa"]) % 360, 3)): float(r["brightness_temperature_K"]) for r in csv.DictReader(file)}


def stats(model, ref):
    d = np.array([model[k] - ref[k] for k in ref if k in model])
    return f"偏差 {d.mean():+.2f} K，RMSE {np.sqrt((d * d).mean()):.2f} K，最大 {np.abs(d).max():.2f} K"


def plot(key):
    test, module, title, filename = SCENES[key]
    ref = compare.load()[test]["eRaytracing"]
    before, after = tirt(test, module, False), tirt(test, module, True)
    out = HERE / "runs" / "compare" / filename
    scene_sketch.scene_figure(test, title, [
        (ref, "o-", {"color": "tab:orange", "ms": 3, "lw": 1.0}, "HiSTREAM 路径追踪"),
        (before, "--", {"color": "0.45", "lw": 1.2}, f"TiRT 修正前（{stats(before, ref)}）"),
        (after, "-", {"color": "tab:blue", "lw": 1.8}, f"TiRT 修正后（{stats(after, ref)}）"),
    ], out)
    print(f"{key}: before {stats(before, ref)}; after {stats(after, ref)}")
    print(out)


if __name__ == "__main__":
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    for key in sys.argv[1:] or SCENES:
        plot(key)
