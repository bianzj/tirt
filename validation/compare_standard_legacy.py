"""canopy_closure standard vs legacy against HiSTREAM on all validation scenes (670 nm, 860 nm, 10.5 um).

python validation/compare_standard_legacy.py

Writes validation/results/standard_vs_legacy/{directions,summary}.csv and scene_<stem>.png.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/standard_vs_legacy"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation"))

import run  # noqa: E402
from validate_current_model import optical  # noqa: E402

LABELS = {
    "01_T1_plane_bare": "平面裸土", "02_T2_plane_hom": "平面均质", "04_T3_plane_crown_leaf0.15": "平面离散树冠",
    "05_T5_row_az0": "垄行（垄向 0°）", "06_T5_row_az90": "垄行（垄向 90°）", "07_S1_slope_hom": "坡面均质",
    "08_S1_slope_crown": "坡面离散树冠", "09_TV_terrain_crown": "复合地形 + 树冠", "10_T4_urban_bare": "城市（无植被）",
    "11_UV_urban_crown": "城市 + 树冠",
}
BANDS = ((670, "红光 670 nm", "BRF", (0.0, 0.20)), (860, "近红外 860 nm", "BRF", (0.10, 0.60)),
         (10500, "热红外 10.5 µm", "亮温 (°C)", (31.0, 43.0)))
CLOSURES = {"standard": "标准版（双尺度热点 + 双层 + 节点辐射度）", "legacy": "legacy（原有模型）"}


def evaluate():
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    inputs = {stem: json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8")) for stem in LABELS}
    rows = []
    for r in reference.itertuples():
        inp = inputs[r.scene]
        values = {}
        for closure in CLOSURES:
            config = dict(inp["thermal_config"], canopy_closure=closure)
            settings = dict(inp["optical_settings"], canopy_closure=closure)
            if r.wavelength_nm == 10500:
                values[closure] = run._surface_result(config, config["surface_model"], config["vegetation_model"], 10.5,
                                                      (r.vza, r.vaa, 30.0, 0.0), inp["emissivity"])[1]
            else:
                values[closure] = optical(settings, config, r.wavelength_nm, r.vza, r.vaa)
        rows.append(dict(scene=r.scene, wavelength_nm=r.wavelength_nm, vza=r.vza, vaa=r.vaa, histream=r.reference, **values))
    return pd.DataFrame(rows)


def summarize(data):
    rows = []
    for (scene, wave), group in data.groupby(["scene", "wavelength_nm"], sort=False):
        for closure in CLOSURES:
            e = group[closure] - group.histream
            rows.append(dict(scene=scene, wavelength_nm=wave, closure=closure,
                             rmse=float(np.sqrt(np.mean(e ** 2))), bias=float(e.mean())))
    return pd.DataFrame(rows)


def plot(data, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    methods = (("histream", "HiSTREAM 模拟", dict(ls="none", marker="o", ms=4.5)),
               ("standard", CLOSURES["standard"], dict(ls="-", lw=1.8)),
               ("legacy", CLOSURES["legacy"], dict(ls="--", lw=1.4)))
    planes = (("太阳主平面", (0.0, 180.0), "tab:red"), ("垂直主平面", (90.0, 270.0), "tab:blue"))
    for stem, label in LABELS.items():
        fig, axes = plt.subplots(1, len(BANDS), figsize=(19, 5.8))
        for panel, (wave, band, unit, limits) in zip(axes, BANDS):
            frame = data[(data.scene == stem) & (data.wavelength_nm == wave)]
            shift = 273.15 if wave == 10500 else 0.0
            for plane, (positive, negative), color in planes:
                ordered = pd.concat([
                    frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                    frame[frame.vza == 0].assign(sign=1),
                    frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
                for key, name, style in methods:
                    panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, color=color,
                               label=f"{plane} | {name}", **style)
            score = summary[(summary.scene == stem) & (summary.wavelength_nm == wave)].set_index("closure").rmse
            fmt = "{:.3f} K" if wave == 10500 else "{:.4f}"
            panel.set(xlim=(-80, 80), ylim=limits, ylabel=unit, xlabel="带符号观测天顶角 (°)",
                      title=f"{band}\nRMSE 标准版 {fmt.format(score['standard'])} | legacy {fmt.format(score['legacy'])}")
            panel.grid(alpha=0.3)
        axes[0].legend(fontsize=7.5, ncol=1, loc="upper left")
        fig.suptitle(f"{label} | 太阳天顶角 30°，太阳方位角 0°", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(OUTPUT / f"scene_{stem}.png", dpi=110)
        plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = evaluate()
    summary = summarize(data)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    plot(data, summary)


if __name__ == "__main__":
    main()
