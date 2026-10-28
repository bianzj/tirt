"""Emission-only (depth 1) thermal term: HiSTREAM square/disk leaves vs the model at two hotspot widths.

python validation/compare_hotspot_scale.py

Model first-order term: K-integral (homogeneous) or parametric rows, hspot = leaf/height (0.15) and
0.424 x leaf/height (finite-leaf overlap integral). Uses the depth-1 runs in validation/streamsim/runs.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation"))

import compare_routes as cr  # noqa: E402
import run  # noqa: E402
from base import row_closure  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402

RUNS = ROOT / "validation/streamsim/runs"
OUTPUT = ROOT / "validation/results/routes"
HSPOTS = (0.15, 0.424 * 0.15)
Y_LIMITS = (31.0, 40.0)


def first_order(scene, kind, hspot, vza, vaa):
    lai, _, shape, row_az = cr.structure(scene, kind)
    emis = scene["emissivity"]
    E = np.array([emis["soil_e"], emis["soil_e"], emis["leaf_e"], emis["leaf_e"]])
    parts = (cr.k_fractions(lai, hspot, shape, row_az, vza, vaa)[0] if kind == "hom" else
             row_closure.fractions("reference", lai, hspot, *shape, vza, vaa, 30.0, 0.0, row_az))
    temperatures = np.asarray(run._vegetation_temperatures(scene["thermal_config"]))
    return float(inv_planck(10.5, np.sum(parts * E * planck(10.5, temperatures))))


def main():
    rows = []
    for stem, label, kind, _ in cr.SCENES:
        scene = json.loads((ROOT / "validation/reference/inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        square = pd.read_csv(RUNS / f"{cr.RUN_NAMES[stem]}_eRaytracing_ss_d1.csv").set_index(["vza", "vaa"]).bt_K
        disk = pd.read_csv(RUNS / f"{cr.RUN_NAMES[stem]}_eRaytracing_disk_d1.csv").set_index(["vza", "vaa"]).bt_K
        for (vza, vaa) in square.index:
            rows.append(dict(scene=stem, vza=vza, vaa=vaa, square=square[(vza, vaa)], disk=disk[(vza, vaa)],
                             **{f"model_{h:.3f}": first_order(scene, kind, h, vza, vaa) for h in HSPOTS}))
    data = pd.DataFrame(rows)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUTPUT / "hotspot_scale_direct.csv", index=False)
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    styles = (("square", "HiSTREAM 方形叶", "o-", "tab:orange"), ("disk", "HiSTREAM 圆形叶", "s--", "tab:brown"),
              (f"model_{HSPOTS[0]:.3f}", "模型 hspot = 叶片/冠层高 = 0.150", "-", "tab:red"),
              (f"model_{HSPOTS[1]:.3f}", "模型 hspot = 0.424 × 叶片/冠层高 = 0.064", "-", "tab:blue"))
    fig, axes = plt.subplots(1, 2 * len(cr.SCENES), figsize=(26, 5.2), sharex=True, sharey=True)
    for block, (stem, label, kind, _) in enumerate(cr.SCENES):
        frame = data[data.scene == stem]
        for offset, ((positive, negative), title) in enumerate(zip(((0.0, 180.0), (90.0, 270.0)),
                                                                    ("平行主平面", "垂直主平面"))):
            ordered = pd.concat([
                frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                frame[frame.vza == 0].assign(sign=1),
                frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
            panel = axes[2 * block + offset]
            for key, name, style, color in styles:
                panel.plot(ordered.vza * ordered.sign, ordered[key] - 273.15, style, color=color, label=name, lw=1.3, ms=3)
            first = "K积分" if kind == "hom" else "parametric"
            panel.set(title=f"{label}（{first}）| {title}", xlim=(-80, 80), ylim=Y_LIMITS)
            panel.set_xlabel("带符号观测天顶角 (°)")
            panel.grid(alpha=0.3)
    axes[0].set_ylabel("热红外 10.5 µm 一次项亮温 (°C)")
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("一次项（只发射，HiSTREAM depth 1）| 热点宽度的影响 | 太阳天顶角 30°，太阳方位角 0°")
    fig.tight_layout()
    fig.savefig(OUTPUT / "hotspot_scale_direct.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
