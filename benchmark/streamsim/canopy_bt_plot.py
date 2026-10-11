"""Canopy brightness temperature only: TiRT with hotspot_analytical vs exact leaf integral vs HiSTREAM (T2 hom, T3 crown)."""
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

import base.proportion as proportion  # noqa: E402
from base.hotspot import hotspot_analytical as old_leaf  # noqa: E402
import compare  # noqa: E402
import hotspot_experiment as H  # noqa: E402
import tirt_side  # noqa: E402
from run import run  # noqa: E402

TESTS = {"T2_plane_hom": "均质冠层", "T3_plane_crown": "离散树冠"}


def tirt_bt(test, patched, **extra):
    proportion.hotspot_analytical = H.exact_leaf if patched else old_leaf
    try:
        work = Path(tempfile.mkdtemp(prefix="tirt-can-"))
        path = run(ROOT / "input.csv", work / "o.csv", **{**tirt_side.COMMON, **tirt_side.TESTS[test], **extra})
        with path.open(encoding="utf-8") as f:
            return {(round(float(r["vza"]), 3), round(float(r["vaa"]) % 360, 3)): float(r["brightness_temperature_K"]) for r in csv.DictReader(f)}
    finally:
        proportion.hotspot_analytical = old_leaf


def main():
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    data = compare.load()
    planes = {"太阳主平面": (0.0, 180.0), "垂直主平面": (90.0, 270.0)}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for i, (test, name) in enumerate(TESTS.items()):
        ref = data[test]["eRaytracing"]
        old, new = tirt_bt(test, False), tirt_bt(test, True)
        for j, (title, (pos, neg)) in enumerate(planes.items()):
            ax = axes[i][j]
            for label, vals, style in (("HiSTREAM 路径追踪", ref, "g-o"), ("TiRT hotspot_analytical", old, "r--"), ("TiRT 精确叶片积分", new, "b-")):
                x, y = compare.signed_curve(vals, pos, neg)
                ax.plot(x, y - 273.15, style, ms=3, lw=1.4, label=label)
            ax.set_title(f"{name} — {title}", fontsize=10)
            ax.set_xlabel("带符号观测天顶角 (°)"); ax.set_ylabel("亮温 (°C)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        keys = sorted(set(ref) & set(old) & set(new))
        text = []
        for label, vals in (("analytical", old), ("精确", new)):
            d = np.array([vals[k] - ref[k] for k in keys])
            text.append(f"{label}: 偏差 {d.mean():+.2f} K, RMSE {np.sqrt((d ** 2).mean()):.2f} K")
        axes[i][0].text(0.02, 0.03, "\n".join(text), transform=axes[i][0].transAxes, fontsize=8, va="bottom",
                        bbox=dict(boxstyle="round", fc="white", alpha=0.85))

        ax = axes[i][2]
        for lai, color in ((1.5, "tab:blue"), (3.0, "tab:green"), (5.0, "tab:orange")):
            extra = {"lai": lai}
            a, b = tirt_bt(test, False, **extra), tirt_bt(test, True, **extra)
            x, ya = compare.signed_curve(a, 90.0, 270.0)
            _, yb = compare.signed_curve(b, 90.0, 270.0)
            ax.plot(x, yb - ya, color=color, label=f"LAI {lai:g}")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{name} — 精确 − analytical 的亮温差（垂直主平面）", fontsize=10)
        ax.set_xlabel("带符号观测天顶角 (°)"); ax.set_ylabel("ΔBT (K)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("冠层亮温：hotspot_analytical 与精确叶片积分（10.5 µm，SZA 30°，叶片 33/30 °C，土壤 45/30 °C）", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = HERE / "runs" / "compare" / "canopy_bt.png"
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
