"""Three levels on the HiSTREAM terrain + crown (TV) and urban + crown (UV) scenes, 10.5 um.

python validation/compare_complex.py

The macro exchange stays the model's own: terrain interreflection (forest_transport=legacy) and the
street + four-wall radiosity (urban_transport=radiosity). Only the crown closures change:
standard   : dual-scale soil hotspot + covered-phase two-layer leaves + node radiosity (crown2)
parametric : parametric crown four components + node radiosity (crown2)
current    : the current default model, for reference.
Writes validation/results/complex/{directions,summary}.csv and one figure per scene.
"""
import functools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation/results/complex"
SCENES = {"09_TV_terrain_crown": "复合地形 + 树冠", "11_UV_urban_crown": "城市 + 树冠"}
SZA, SAA = 30.0, 0.0
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation"))

import compare_all_scenes as cas  # noqa: E402
import run  # noqa: E402
from base.terrain_forest import TerrainForest  # noqa: E402
from base.urban_veg import Urban_Veg  # noqa: E402

LEVELS = {"current": "现有标准模型", "standard": cas.LEVELS["standard"], "parametric": cas.LEVELS["parametric"]}


def use_level(level):
    fraction = None if level == "current" else functools.partial(cas.crown_level, level)
    scatter = None if level == "current" else staticmethod(cas.crown_scatter)
    for model in (TerrainForest, Urban_Veg):
        model.canopy_fraction_model = fraction
        model.canopy_scatter_model = scatter


def thermal(scene, level, vza, vaa):
    config = dict(scene["thermal_config"])
    if level != "current" and config["surface_model"] == "urban":
        config["urban_transport"] = "radiosity"
    use_level(level)
    try:
        return run._surface_result(config, config["surface_model"], "crown", 10.5, (vza, vaa, SZA, SAA),
                                   scene["emissivity"])[1]
    finally:
        use_level("current")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(ROOT / "validation/results/reference_comparison.csv")
    rows, summary = [], []
    for stem, label in SCENES.items():
        scene = json.loads((ROOT / "validation/reference/inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        frame = reference[(reference.scene == stem) & (reference.wavelength_nm == 10500) &
                          (reference.status == "ok")].dropna(subset=["reference"])
        values = {level: np.array([thermal(scene, level, r.vza, r.vaa) for r in frame.itertuples()]) for level in LEVELS}
        for level, v in values.items():
            e = v - frame.reference.to_numpy()
            summary.append(dict(scene=stem, label=label, level=level, rmse=float(np.sqrt(np.mean(e ** 2))),
                                bias=float(e.mean()), max_abs=float(np.abs(e).max())))
        for i, r in enumerate(frame.itertuples()):
            rows.append(dict(scene=stem, vza=r.vza, vaa=r.vaa, histream=r.reference, **{k: v[i] for k, v in values.items()}))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False, encoding="utf-8-sig")
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 200):
        print(summary.drop(columns="label").to_string(index=False))
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    methods = (("histream", "HiSTREAM 模拟", dict(ls="none", marker="o", ms=4.5)),
               ("standard", LEVELS["standard"], dict(ls="-", lw=1.8)),
               ("parametric", LEVELS["parametric"], dict(ls="--", lw=1.5)),
               ("current", LEVELS["current"], dict(ls=":", lw=1.5)))
    planes = (("太阳主平面", (0.0, 180.0), "tab:red"), ("垂直主平面", (90.0, 270.0), "tab:blue"))
    for stem, label in SCENES.items():
        frame = data[data.scene == stem]
        fig, panel = plt.subplots(figsize=(9, 6))
        for plane, (positive, negative), color in planes:
            ordered = pd.concat([
                frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                frame[frame.vza == 0].assign(sign=1),
                frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
            for key, name, style in methods:
                panel.plot(ordered.vza * ordered.sign, ordered[key] - 273.15, color=color, label=f"{plane} | {name}", **style)
        panel.set(xlim=(-80, 80), ylim=(31.0, 43.0), ylabel="亮温 (°C)", xlabel="带符号观测天顶角 (°)",
                  title=f"{label} | 热红外 10.5 µm | 太阳天顶角 30°，太阳方位角 0°")
        panel.grid(alpha=0.3)
        panel.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(OUTPUT / f"complex_{stem}.png", dpi=120)
        plt.close(fig)


if __name__ == "__main__":
    main()
