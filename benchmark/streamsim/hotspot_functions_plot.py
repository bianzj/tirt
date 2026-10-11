"""Plot hotspot_layer (sunlit soil) vs hotspot_analytical (sunlit leaf) against the exact bidirectional-gap integral."""
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
warnings.simplefilter("ignore")

import base.hotspot as hs  # noqa: E402
import hotspot_experiment as H  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402

SZA, HSPOT, G = 30.0, 0.15, 0.5
arr = lambda v: np.array([float(v)])


def layer(lai, vza, raa):
    return float(np.asarray(hs.hotspot_layer(lai, HSPOT, arr(vza), arr(SZA), arr(raa))).reshape(-1)[0])


def analytical(lai, vza, raa):
    return float(np.asarray(hs.hotspot_analytical(lai, HSPOT, arr(vza), arr(SZA), arr(raa))).reshape(-1)[0])


def exact_leaf(lai, vza, raa):
    return H.exact_leaf_scalar(lai, HSPOT, vza, SZA, raa, 1.0, 1.0)


def main():
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    vza = np.arange(0, 76, 2.5)
    vza_plot = np.maximum(vza, 1e-3)

    # (0,0): sunlit leaf, principal plane signed VZA, LAI 1.5
    ax = axes[0][0]
    xs = np.concatenate([-vza[::-1], vza])
    for label, fn, style in (("exact integral", exact_leaf, "k-"), ("hotspot_analytical (leaf)", analytical, "r--")):
        ys = [fn(1.5, max(abs(v), 1e-3), 180 if v < 0 else 0) for v in xs]
        ax.plot(xs, ys, style, label=label)
    ys = [layer(1.5, max(abs(v), 1e-3), 180 if v < 0 else 0) for v in xs]
    ax.plot(xs, ys, "b-", label="hotspot_layer (soil, = exact)")
    ax.set_title("LAI 1.5, solar principal plane: sunlit fractions of the scene", fontsize=10)
    ax.set_xlabel("signed VZA (deg)"); ax.set_ylabel("fraction of scene"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (0,1): leaf error vs VZA for several LAI (raa 90)
    ax = axes[0][1]
    for lai in (0.5, 1.5, 3.0, 5.0):
        ys = [analytical(lai, v, 90) - exact_leaf(lai, v, 90) for v in vza_plot]
        ax.plot(vza, ys, label=f"LAI {lai}")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_title("hotspot_analytical − exact (sunlit leaf), raa 90°", fontsize=10)
    ax.set_xlabel("VZA (deg)"); ax.set_ylabel("difference in scene fraction"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (0,2): max error vs LAI for both functions at vza 0 and 60
    ax = axes[0][2]
    lais = np.linspace(0.25, 6.0, 24)
    for v, style in ((0.001, "-"), (60, "--")):
        ax.plot(lais, [analytical(l, v, 90) - exact_leaf(l, v, 90) for l in lais], "r" + style, label=f"analytical − exact, VZA {int(v)}°")
        ax.plot(lais, [layer(l, v, 90) - layer(l, v, 90) for l in lais], "b" + style, label=f"layer − exact, VZA {int(v)}°" if v > 1 else "layer − exact")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_title("error vs LAI (raa 90°)", fontsize=10)
    ax.set_xlabel("LAI"); ax.set_ylabel("difference in scene fraction"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (1,0): sunlit share of the visible component (soil vs leaf), LAI 1.5 and 5
    ax = axes[1][0]
    for lai, color in ((1.5, "tab:blue"), (5.0, "tab:orange")):
        vis_soil = np.exp(-G * lai / np.cos(np.deg2rad(vza_plot)))
        ax.plot(vza, [layer(lai, v, 90) for v in vza_plot] / vis_soil, color=color, ls="-", label=f"soil share (layer), LAI {lai}")
        ax.plot(vza, [analytical(lai, v, 90) for v in vza_plot] / (1 - vis_soil), color=color, ls="--", label=f"leaf share (analytical), LAI {lai}")
        ax.plot(vza, [exact_leaf(lai, v, 90) for v in vza_plot] / (1 - vis_soil), color=color, ls=":", label=f"leaf share (exact), LAI {lai}")
    ax.set_title("sunlit share of the visible component, raa 90°", fontsize=10)
    ax.set_xlabel("VZA (deg)"); ax.set_ylabel("sunlit / visible"); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(fontsize=7)

    # (1,1): BT effect for leaf 33/30 C and for a stronger contrast 45/30 C
    ax = axes[1][1]
    for (ts, th), style in (((33.0, 30.0), "-"), ((45.0, 30.0), "--")):
        Bs, Bh = planck(10.5, ts + 273.15), planck(10.5, th + 273.15)
        for lai, color in ((1.5, "tab:blue"), (5.0, "tab:orange")):
            # sunlit-leaf error moves the same amount of leaf area from sunlit to shaded; evaluate on a half sunlit / half shaded leaf mix
            d = np.array([analytical(lai, v, 90) - exact_leaf(lai, v, 90) for v in vza_plot])
            vis_leaf = 1 - np.exp(-G * lai / np.cos(np.deg2rad(vza_plot)))
            exact_share = np.array([exact_leaf(lai, v, 90) for v in vza_plot]) / vis_leaf
            share_new = exact_share + d / vis_leaf
            # leaf contributes with weight vis_leaf; BT effect of the leaf term alone, soil assumed at the mean leaf temperature
            radiance = lambda s: 0.985 * (s * Bs + (1 - s) * Bh)
            ax.plot(vza, [float(inv_planck(10.5, radiance(a)) - inv_planck(10.5, radiance(b))) * w
                          for a, b, w in zip(share_new, exact_share, vis_leaf)],
                    color=color, ls=style, label=f"LAI {lai}, leaf {int(ts)}/{int(th)} °C")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_title("leaf-term brightness-temperature error (leaf only)", fontsize=10)
    ax.set_xlabel("VZA (deg)"); ax.set_ylabel("ΔBT (K) from the sunlit-leaf error"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (1,2): heat map of the error over LAI and VZA
    ax = axes[1][2]
    lai_grid = np.linspace(0.25, 6.0, 24)
    err = np.array([[analytical(l, v, 90) - exact_leaf(l, v, 90) for v in vza_plot[::2]] for l in lai_grid])
    im = ax.pcolormesh(vza[::2], lai_grid, err, shading="auto", cmap="RdBu", vmin=-0.3, vmax=0.3)
    ax.set_title("hotspot_analytical − exact over LAI × VZA (raa 90°)", fontsize=10)
    ax.set_xlabel("VZA (deg)"); ax.set_ylabel("LAI"); fig.colorbar(im, ax=ax, label="difference in scene fraction")

    fig.suptitle(f"hotspot_layer (sunlit soil) vs hotspot_analytical (sunlit leaf) — homogeneous canopy, hspot {HSPOT}, SZA {SZA:.0f}°", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = HERE / "runs" / "compare" / "hotspot_layer_vs_analytical.png"
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
