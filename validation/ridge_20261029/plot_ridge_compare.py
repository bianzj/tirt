"""Plot the HiSTREAM vs TiRT multi-ridge comparison."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "ridge_compare.csv").open(encoding="utf-8")))
plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
colors = {90.0: "tab:blue", 270.0: "tab:red"}
for ax, sza in zip(axes, (30.0, 70.0)):
    for vaa in (90.0, 270.0):
        data = [r for r in rows if float(r["sza"]) == sza and float(r["vaa"]) == vaa]
        data.sort(key=lambda r: float(r["vza"]))
        x = [float(r["vza"]) for r in data]
        tirt = [float(r["tirt_bt"]) for r in data]
        hist = [float(r["histream_bt"]) for r in data]
        ax.plot(x, tirt, "-", color=colors[vaa], lw=1.8, label=f"TiRT 面元集合, VAA {vaa:g}°")
        ax.plot(x, hist, "o", color=colors[vaa], ms=5, mfc="white", label=f"HiSTREAM, VAA {vaa:g}°")
    d = np.asarray([float(r["difference"]) for r in rows if float(r["sza"]) == sza])
    ax.set_title(f"SZA {sza:.0f}°，SAA 270°；RMSE {np.sqrt(np.mean(d ** 2)):.3f} K，偏差 {d.mean():+.3f} K")
    ax.set_xlabel("观测天顶角 VZA (°)")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("方向亮温 (K)")
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle("多条横向沙丘（100.5 m × 100.5 m，201×201，起伏 6.1 m）：TiRT surface_model=dem vs HiSTREAM eRaytracing", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(HERE / "ridge_vs_histream.png", dpi=110)
plt.close(fig)
print("saved", HERE / "ridge_vs_histream.png")
