"""Independent check of the row-crop fractions: exact 2-D periodic ray integration (turbid rows) vs eRaytracing vs TiRT."""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
warnings.simplefilter("ignore")

from base.proportion import proportion_bidirectional_row_one  # noqa: E402

K = 1.5          # G * in-row LAD (m-1)
WIDTH, GAP, HEIGHT = 0.5, 0.5, 1.0
PERIOD = WIDTH + GAP
SZA = 30.0


def row_optical_depth(x0, tan_theta, cos_theta, direction, across):
    """Optical depth of the ray from soil point x0 (cross-row coordinate) up to the canopy top."""
    z = (np.arange(200) + 0.5) / 200 * HEIGHT
    if not across:
        return np.where((x0 % PERIOD) < WIDTH, K * HEIGHT / cos_theta, 0.0)
    x = x0[:, None] + direction * tan_theta * z[None, :]
    return K * HEIGHT / cos_theta * ((x % PERIOD) < WIDTH).mean(axis=1)


def exact_fractions(row_az, vza, vaa):
    """Visible soil and sunlit visible soil for rows along compass azimuth row_az (sun at saa 0, sza 30)."""
    x0 = (np.arange(400) + 0.5) / 400 * PERIOD
    rd = np.deg2rad
    # sign of the cross-row component of the horizontal direction
    def components(az):
        rel = rd(az - row_az)
        return np.cos(rel), np.sin(rel)  # along-row, across-row
    _, v_across = components(vaa)
    _, s_across = components(0.0)
    tv = row_optical_depth(x0, np.tan(rd(vza)) * abs(v_across), np.cos(rd(vza)), np.sign(v_across), abs(v_across) > 1e-9)
    ts = row_optical_depth(x0, np.tan(rd(SZA)) * abs(s_across), np.cos(rd(SZA)), np.sign(s_across), abs(s_across) > 1e-9)
    return float(np.mean(np.exp(-tv))), float(np.mean(np.exp(-tv - ts)))


def tirt_fractions(row_az, vza, vaa):
    one = lambda v: np.array([float(v)])
    p = proportion_bidirectional_row_one(1.5, 0.15, 0.5, 0.5, 1.0, one(max(vza, 1e-3)), one(vaa), one(SZA), one(0.0), one(row_az))
    p = [float(np.asarray(x).reshape(-1)[0]) for x in p]
    return p[0] + p[1], p[0]


def rt_table(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[(float(r["vza"]), float(r["vaa"]))] = (float(r["f_soil_sunlit"]) + float(r["f_soil_shaded"]), float(r["f_soil_sunlit"]))
    return out


def main():
    for row_az in (0, 90):
        rt = rt_table(HERE / "runs" / f"T5_row_az{row_az}_fractions_eRaytracing.csv")
        print(f"\n=== rows along azimuth {row_az}  (visible soil / sunlit soil; exact | eRaytracing | TiRT) ===")
        err = {"rt_vis": [], "rt_sun": [], "tirt_vis": [], "tirt_sun": []}
        for vaa in (0, 90, 180):
            for vza in (0, 20, 30, 40, 60, 75):
                ev, es = exact_fractions(row_az, vza, vaa)
                rv, rs = rt[(float(vza), float(vaa))]
                tv, ts = tirt_fractions(row_az, vza, vaa)
                tag = "  (hotspot region: exact ignores leaf overlap)" if abs(vza - SZA) < 12 and vaa == 0 else ""
                print(f"vaa {vaa:3d} vza {vza:2d}: vis {ev:.3f} | {rv:.3f} | {tv:.3f}    sunlit {es:.3f} | {rs:.3f} | {ts:.3f}{tag}")
                if not tag:
                    for key, a, b in (("rt_vis", rv, ev), ("rt_sun", rs, es), ("tirt_vis", tv, ev), ("tirt_sun", ts, es)):
                        err[key].append(a - b)
        print("deviation from exact (mean / rms, hotspot region excluded): " + "   ".join(
            f"{k} {np.mean(v):+.3f}/{np.sqrt(np.mean(np.square(v))):.3f}" for k, v in err.items()))


if __name__ == "__main__":
    main()
