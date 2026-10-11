"""Experiment: replace hotspot_analytical (two-layer leaf approximation) by the exact bidirectional-gap integral.

Leaf sunlit fraction = integral of k_o * P_so(x) dx over canopy depth, with the same P_so as hotspot_layer
(so soil and leaves share one hotspot model). Nothing in base/ is modified; the function is patched in memory.
"""
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

import base.proportion as proportion  # noqa: E402
from base.hotspot import hotspot_analytical as old_leaf  # noqa: E402
import compare  # noqa: E402
import tirt_side  # noqa: E402
from run import run  # noqa: E402

G = 0.5


def exact_leaf_scalar(lai, hspot, vza, sza, raa, ci_i, ci_v, n=400):
    rd = np.pi / 180.0
    uv, us = np.cos(vza * rd), np.cos(sza * rd)
    ko, ks = ci_v * G / uv, ci_i * G / us
    dso2 = np.tan(sza * rd) ** 2 + np.tan(vza * rd) ** 2 - 2 * np.tan(sza * rd) * np.tan(vza * rd) * np.cos(raa * rd)
    alf = np.sqrt(max(dso2, 0.0)) / hspot * 2.0 / (G / uv + G / us)
    x = np.linspace(0.0, lai, n + 1)
    shape = x * (1 - alf * x / (2 * lai)) if alf < 1e-9 else lai / alf * (-np.expm1(-alf * x / lai))
    pso = np.exp(-(ks + ko) * x + np.sqrt(ko * ks) * shape)
    return float(np.trapezoid(ko * pso, x))


def exact_leaf(lai, hspot, vza, sza, vsa, CIi=1, CIv=1, G=0.5):
    """Drop-in for hotspot_analytical (same signature, returns the sunlit-leaf fraction of the scene)."""
    vza, sza, vsa = np.broadcast_arrays(np.atleast_1d(vza), np.atleast_1d(sza), np.atleast_1d(vsa))
    ci_i, ci_v = np.broadcast_to(np.asarray(CIi, float), vza.shape), np.broadcast_to(np.asarray(CIv, float), vza.shape)
    out = np.array([exact_leaf_scalar(lai, hspot, float(a), float(b), float(c), float(i), float(v))
                    for a, b, c, i, v in zip(vza.ravel(), sza.ravel(), vsa.ravel(), ci_i.ravel(), ci_v.ravel())])
    return out.reshape(vza.shape)


def tirt_curves(test):
    warnings.simplefilter("ignore", RuntimeWarning)
    work = Path(tempfile.mkdtemp(prefix="tirt-hs-"))
    path = run(ROOT / "input.csv", work / test / "output.csv", **{**tirt_side.COMMON, **tirt_side.TESTS[test]})
    import csv
    values = {}
    with path.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            values[(round(float(r["vza"]), 3), round(float(r["vaa"]) % 360, 3))] = float(r["brightness_temperature_K"])
    return values


TESTS = ("T2_plane_hom", "T3_plane_crown", "T5_row_az0", "T5_row_az90")
PLANES = {"principal (+: vaa 0, sun side)": (0.0, 180.0), "perpendicular (+: vaa 90)": (90.0, 270.0)}


def main():
    data = compare.load()
    fig, axes = plt.subplots(len(TESTS), 2, figsize=(11, 3.2 * len(TESTS)), squeeze=False)
    stats = []
    for i, test in enumerate(TESTS):
        reference, old = data[test]["eRaytracing"], data[test]["TiRT"]
        proportion.hotspot_analytical = exact_leaf
        try:
            new = tirt_curves(test)
        finally:
            proportion.hotspot_analytical = old_leaf
        for j, (title, (pos, neg)) in enumerate(PLANES.items()):
            ax = axes[i][j]
            for label, vals, style in (("HiSTREAM eRaytracing", reference, "g-o"), ("TiRT current", old, "r--"), ("TiRT exact leaf integral", new, "b-")):
                x, y = compare.signed_curve(vals, pos, neg)
                ax.plot(x, y - 273.15, style, ms=3, lw=1.4, label=label)
            ax.set_title(f"{test} — {title}", fontsize=9)
            ax.set_xlabel("signed VZA (deg)"); ax.set_ylabel("BT (°C)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        keys = sorted(set(reference) & set(old) & set(new))
        for name, vals in (("TiRT current", old), ("TiRT exact leaf integral", new)):
            d = np.array([vals[k] - reference[k] for k in keys])
            stats.append((test, name, d.mean(), np.sqrt((d ** 2).mean()), np.abs(d).max()))
    fig.tight_layout()
    out = HERE / "runs" / "compare" / "hotspot_experiment.png"
    fig.savefig(out, dpi=130)
    for s in stats:
        print(f"{s[0]:15s} {s[1]:26s} bias {s[2]:+.2f} K  rmse {s[3]:.2f}  max {s[4]:.2f}")
    print(out)


if __name__ == "__main__":
    main()
