"""Figures for the test report (current standard model).

1. TiRT vs HiSTREAM: one figure per scene (3 bands x principal / perpendicular plane), from validation/results.
2. TIRTEB vs SCOPE 1.73: component temperature time series and energy fluxes, from validation/scope/runs.
Outputs: docs/report/figures/*.png (all < 3000 px).
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures"
plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"], "axes.unicode_minus": False})
BANDS = (("670", "BRF 0.67 µm"), ("860", "BRF 0.86 µm"), ("10500", "亮温 10.5 µm (°C)"))


def histream():
    manifest = json.loads((ROOT / "validation/reference/manifest.json").read_text(encoding="utf-8"))
    labels = {s["stem"]: s["label"] for s in manifest["scenes"]}
    data = defaultdict(list)
    for r in csv.DictReader(open(ROOT / "validation/results/reference_comparison.csv", encoding="utf-8")):
        data[(r["scene"], r["wavelength_nm"])].append(r)
    acc = {(r["scene"], r["wavelength_nm"]): r for r in csv.DictReader(open(ROOT / "validation/results/accuracy.csv"))}
    for stem, label in labels.items():
        fig, axes = plt.subplots(3, 2, figsize=(11, 7.4), sharex=True)
        for row, (band, name) in enumerate(BANDS):
            rows = data[(stem, band)]
            offset = 273.15 if band == "10500" else 0.0
            for col, (pos, neg, title) in enumerate(((0., 180., "太阳主平面"), (90., 270., "垂直太阳主平面"))):
                pts = []
                for r in rows:
                    vza, rel = float(r["vza"]), (float(r["vaa"]) - float(r["saa"])) % 360
                    sign = 1 if vza == 0 or abs(rel - pos) < 1e-6 else (-1 if abs(rel - neg) < 1e-6 else 0)
                    if sign:
                        pts.append((sign * vza, float(r["value"]) - offset, float(r["reference"]) - offset))
                pts.sort()
                ax = axes[row, col]
                if pts:
                    x, model, ref = map(np.array, zip(*pts))
                    ax.plot(x, ref, "o-", color="tab:orange", ms=3, lw=1.2, label="HiSTREAM")
                    ax.plot(x, model, "-", color="tab:blue", lw=1.6, label="TiRT")
                ax.grid(alpha=0.3)
                ax.set_xlim(-80, 80)
                if row == 0:
                    ax.set_title(title, fontsize=11)
                if col == 0:
                    ax.set_ylabel(name)
                if row == 2:
                    ax.set_xlabel("带符号观测天顶角 (°)，正值为太阳一侧")
            a = acc[(stem, band)]
            unit = " K" if band == "10500" else ""
            axes[row, 0].text(0.02, 0.96, f"RMSE {float(a['rmse']):.4g}{unit}  偏差 {float(a['bias']):+.4g}{unit}",
                              transform=axes[row, 0].transAxes, va="top", fontsize=8.5,
                              bbox=dict(fc="white", ec="none", alpha=0.85))
            for ax in axes[row]:
                ax.ticklabel_format(axis="y", useOffset=False)
        axes[0, 1].legend(fontsize=9, loc="best")
        fig.suptitle(f"{label}：TiRT 与 HiSTREAM（太阳天顶角 30°，方位 0°）", fontsize=13)
        fig.tight_layout()
        fig.savefig(OUT / f"histream_{stem}.png", dpi=130)
        plt.close(fig)
    # summary bar chart of RMSE
    stems = list(labels)
    fig, axes = plt.subplots(3, 1, figsize=(9, 11))
    for ax, (band, name) in zip(axes, BANDS):
        values = [float(acc[(s, band)]["rmse"]) for s in stems]
        ax.barh(range(len(stems)), values, color="tab:blue")
        ax.set_yticks(range(len(stems)), [labels[s] for s in stems], fontsize=9)
        ax.invert_yaxis()
        ax.set_title(f"RMSE：{name.replace(' (°C)', '')}" + (" (K)" if band == "10500" else ""), fontsize=11)
        ax.grid(axis="x", alpha=0.3)
        for i, v in enumerate(values):
            ax.text(v, i, f" {v:.3g}", va="center", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "histream_rmse.png", dpi=120)
    plt.close(fig)


def scope():
    rows = list(csv.DictReader(open(ROOT / "validation/scope/runs/scope_comparison.csv", encoding="utf-8")))
    t = np.array([float(r["time"]) for r in rows])
    hour = (t - np.floor(t[0])) * 24.0

    def series(key, source):
        return np.array([float(r[f"{key}_{source}"]) for r in rows])

    stats = {r["variable"] + "|" + r["period"]: r for r in
             csv.DictReader(open(ROOT / "validation/scope/runs/scope_summary.csv", encoding="utf-8"))}

    def panel(ax, key, title):
        s, m = series(key, "scope"), series(key, "tirteb")
        ax.plot(hour, s, "o-", color="tab:orange", ms=3, lw=1.2, label="SCOPE 1.73")
        ax.plot(hour, m, "-", color="tab:blue", lw=1.6, label="TIRTEB")
        st = stats[key + "|all"]
        ax.set_title(f"{title}\n偏差 {float(st['bias']):+.2f}，RMSE {float(st['rmse']):.2f}", fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_xticks(range(0, 25, 3))

    temps = (("Tsoil_sunlit", "光照土壤温度 (°C)"), ("Tsoil_shaded", "阴影土壤温度 (°C)"),
             ("Tsoil", "土壤平均温度 (°C)"), ("Tcanopy", "冠层温度 (°C)"), ("BT_nadir", "天底亮温 10.5 µm (°C)"))
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    for ax, (key, title) in zip(axes.flat, temps):
        panel(ax, key, title)
    ax = axes.flat[5]
    ax.plot(hour, series("Tcanopy", "tirteb") - series("Tcanopy", "scope"), lw=1.5, label="冠层")
    ax.plot(hour, series("Tsoil", "tirteb") - series("Tsoil", "scope"), lw=1.5, label="土壤平均")
    ax.plot(hour, series("BT_nadir", "tirteb") - series("BT_nadir", "scope"), lw=1.5, label="天底亮温")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_title("温度差 TIRTEB − SCOPE (K)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_xticks(range(0, 25, 3))
    axes[0, 0].legend(fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel("当地时间 (h)，DOY 169")
    fig.suptitle("TIRTEB 与 SCOPE 1.73 验证算例：温度日变化（均质冠层，LAI 3）", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "scope_temperature.png", dpi=130)
    plt.close(fig)

    fluxes = (("Rn", "净辐射 Rn (W/m²)"), ("LE", "潜热 LE (W/m²)"), ("H", "感热 H (W/m²)"),
              ("G", "土壤热通量 G (W/m²)"), ("Rn_canopy", "冠层净辐射 (W/m²)"), ("Rn_soil", "土壤净辐射 (W/m²)"),
              ("LE_canopy", "冠层潜热 (W/m²)"), ("LE_soil", "土壤潜热 (W/m²)"))
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), sharex=True)
    for ax, (key, title) in zip(axes.flat, fluxes):
        panel(ax, key, title)
    axes[0, 0].legend(fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel("当地时间 (h)")
    fig.suptitle("TIRTEB 与 SCOPE 1.73：能量通量日变化", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "scope_fluxes.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    histream()
    scope()
    print("ok")
