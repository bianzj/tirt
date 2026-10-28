"""Three-level comparison on the HiSTREAM validation scenes (670 nm, 860 nm, 10.5 um).

python validation/compare_all_scenes.py

histream   : reference
standard   : dual-scale soil hotspot + two-layer sunlit leaves + node radiosity
             (soil hotspot: hom = layer hotspot, rows = two-phase strip overlap, crowns = lens overlap)
parametric : exact geometric-optics four components (hom K-integral, rows cross-section integral,
             crowns Poisson spheroids) + the same node radiosity
Nodes: hom3 / row3 / crown2, with bi-Lambertian (back/forward) leaf scattering. Slopes are evaluated in the slope-local frame HiSTREAM renders.
Terrain and urban scenes are not yet in this framework.
Writes validation/results/all_scenes/{directions,summary}.csv and scene_<stem>.png.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/all_scenes"
SZA, SAA, G = 30.0, 0.0, 0.5
NODES = {"hom": "hom3", "row": "row3", "crown": "crown2"}
LEVELS = {"standard": "双尺度热点 + 双层 + 节点辐射度", "parametric": "parametric + 节点辐射度"}
SCENES = {
    "01_T1_plane_bare": "平面裸土", "02_T2_plane_hom": "平面均质", "04_T3_plane_crown_leaf0.15": "平面离散树冠",
    "05_T5_row_az0": "垄行（垄向 0°）", "06_T5_row_az90": "垄行（垄向 90°）", "07_S1_slope_hom": "坡面均质",
    "08_S1_slope_crown": "坡面离散树冠",
}
BANDS = ((670, "红光 670 nm", "BRF", (0.0, 0.20)), (860, "近红外 860 nm", "BRF", (0.10, 0.60)),
         (10500, "热红外 10.5 µm", "亮温 (°C)", (31.0, 43.0)))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation/streamsim"))

import run  # noqa: E402
import run_benchmark  # noqa: E402
from base import node_radiosity, row_closure  # noqa: E402
from base.crown_parametric import crown_fractions  # noqa: E402
from base.hotspot import two_layer_sunlit_leaf_covered  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402
from base.radiation import relative_azimuth  # noqa: E402
from base.vnir import joint_transmission, spherical_coefficients  # noqa: E402


def structure(config):
    kind = config["vegetation_model"]
    if kind == "hom":
        return kind, (0.5, 0.0, 1.0), 0.0
    if kind == "row":
        return kind, tuple(float(config[k]) for k in ("row_width", "row_blank", "row_height")), float(config["row_azimuth"])
    return kind, tuple(float(config[k]) for k in ("crown_density", "crown_rad_b", "crown_rad_a")), 0.0


def crown_level(level, lai, density, hspot, vertical, horizontal, base, vza, vaa, sza, saa):
    """Crown four components of a level for one direction (degrees)."""
    zv, zs = min(abs(float(vza)), 89.9), min(abs(float(sza)), 89.9)
    raa = float(relative_azimuth(saa, vaa))
    args = (lai, density, hspot, vertical, horizontal, base)
    parts = crown_fractions(*args, zv, zs, raa)[:, 0]
    if level == "parametric":
        return parts
    pv = parts[0] + parts[1]
    sun = crown_fractions(*args, zs, zs, 0.0)[:, 0]
    leaf = two_layer_sunlit_leaf_covered(lai, pv, sun[0] + sun[1], math.cos(math.radians(zv)),
                                         math.cos(math.radians(zs)), row_closure.leaf_correlation(hspot, zv, vaa, zs, saa))
    return np.array([parts[0], parts[1], leaf, 1.0 - pv - leaf])


def crown_scatter(kind, lai, emis_leaf, emis_soil, vza, sza, density, vertical, horizontal):
    """Node-radiosity canopy scattering (4, N) with the canopy_scatter signature (crown schedule is azimuth-free)."""
    el, es = float(np.ravel(emis_leaf)[0]), float(np.ravel(emis_soil)[0])
    vza, sza = np.broadcast_arrays(np.atleast_1d(vza).astype(float), np.atleast_1d(sza).astype(float))
    return np.asarray([node_radiosity.thermal_scatter(NODES["crown"], float(lai), (float(density), float(vertical),
                                                      float(horizontal)), el, es, v, 0.0, s, 0.0)
                       for v, s in zip(vza.ravel(), sza.ravel())]).T


def fractions(level, config, vza, vaa, sza, saa):
    """Sunlit soil, shaded soil, sunlit leaf, shaded leaf of one level."""
    lai, hspot = float(config["lai"]), float(config["hspot"])
    kind, shape, az = structure(config)
    raa = float(relative_azimuth(saa, vaa))
    if kind == "crown":
        density, vertical, horizontal = shape
        return crown_level(level, lai, density, hspot, vertical, horizontal, float(config["crown_height"]),
                           vza, vaa, sza, saa)
    if level == "parametric":
        if kind == "row":
            return row_closure.fractions("reference", lai, hspot, *shape, vza, vaa, sza, saa, az)
        zv, zs = min(abs(vza), 89.999), min(sza, 89.999)
        ks, kv = G / math.cos(math.radians(zs)), G / math.cos(math.radians(zv))
        joint, path = joint_transmission(lai, ks, kv, hspot, zs, zv, raa)
        gap = math.exp(-kv * lai)
        return np.array([joint, gap - joint, kv * path, 1.0 - gap - kv * path])
    soil = row_closure.fractions("conditioned", lai, hspot, *shape, vza, vaa, sza, saa, az)
    pv = soil[0] + soil[1]
    ps = row_closure.sun_gap(lai, *shape, sza, saa, az)
    zv, zs = min(abs(vza), 89.999), min(sza, 89.999)
    leaf = two_layer_sunlit_leaf_covered(lai, pv, ps, math.cos(math.radians(zv)), math.cos(math.radians(zs)),
                                         row_closure.leaf_correlation(hspot, vza, vaa, sza, saa))
    return np.array([soil[0], soil[1], leaf, 1.0 - pv - leaf])


def local(config, vza, vaa):
    if config["surface_model"] != "slope":
        return vza, vaa, SZA, SAA
    frame = run_benchmark.slope_frame(float(config["slope_angle"]), float(config["slope_aspect"]))
    (v, a), (s, sa) = run_benchmark.to_local(vza, vaa, frame), run_benchmark.to_local(SZA, SAA, frame)
    return v, a, s, sa


def thermal(scene, level, vza, vaa):
    config, emis = scene["thermal_config"], scene["emissivity"]
    if config["vegetation_model"] == "bare":
        return run._surface_result(config, "plane", "bare", 10.5, (vza, vaa, SZA, SAA), emis)[1]
    kind, shape, az = structure(config)
    v, a, s, sa = local(config, vza, vaa)
    el, es = emis["leaf_e"], emis["soil_e"]
    weights = fractions(level, config, v, a, s, sa) * np.array([es, es, el, el]) + \
        node_radiosity.thermal_scatter(NODES[kind], float(config["lai"]), shape, el, es, v, a, s, sa, az)
    radiance = float(np.sum(weights * planck(10.5, np.asarray(run._vegetation_temperatures(config)))))
    return float(inv_planck(10.5, radiance))


def optical(scene, level, wavelength, vza, vaa):
    config, settings = scene["thermal_config"], scene["optical_settings"]
    band = settings["wavelengths"].index(wavelength)
    rho, tau, soil = (settings[key][band] for key in
                      ("leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"))
    if config["vegetation_model"] == "bare":
        return float(soil)
    kind, shape, az = structure(config)
    v, a, s, sa = local(config, vza, vaa)
    parts = fractions(level, config, v, a, s, sa)
    _, ko, sob, sof = spherical_coefficients(float(s), float(v), float(relative_azimuth(sa, a)))
    single = (sob * rho + sof * tau) / ko * parts[2] + soil * parts[0]
    return single + node_radiosity.shortwave_multiple(NODES[kind], float(config["lai"]), shape, rho, tau, soil,
                                                      v, a, s, sa, az)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(ROOT / "validation/results/reference_comparison.csv")
    rows, summary = [], []
    for stem, label in SCENES.items():
        scene = json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        for wavelength, *_ in BANDS:
            frame = reference[(reference.scene == stem) & (reference.wavelength_nm == wavelength) &
                              (reference.status == "ok")].dropna(subset=["reference"])
            values = {level: np.array([thermal(scene, level, r.vza, r.vaa) if wavelength == 10500 else
                                       optical(scene, level, wavelength, r.vza, r.vaa) for r in frame.itertuples()])
                      for level in LEVELS}
            for level, v in values.items():
                e = v - frame.reference.to_numpy()
                summary.append(dict(scene=stem, label=label, wavelength_nm=wavelength, level=level,
                                    rmse=float(np.sqrt(np.mean(e ** 2))), bias=float(e.mean()),
                                    max_abs=float(np.abs(e).max())))
            for i, r in enumerate(frame.itertuples()):
                rows.append(dict(scene=stem, wavelength_nm=wavelength, vza=r.vza, vaa=r.vaa, histream=r.reference,
                                 **{level: v[i] for level, v in values.items()}))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False, encoding="utf-8-sig")
    table = summary.pivot_table(index=["scene", "wavelength_nm"], columns="level", values=["rmse", "bias"])
    with pd.option_context("display.float_format", "{:.4f}".format, "display.width", 200):
        print(table.to_string())
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    methods = (("histream", "HiSTREAM 模拟", dict(ls="none", marker="o", ms=4.5)),
               ("standard", LEVELS["standard"], dict(ls="-", lw=1.8)),
               ("parametric", LEVELS["parametric"], dict(ls="--", lw=1.5)))
    planes = (("太阳主平面", (0.0, 180.0), "tab:red"), ("垂直主平面", (90.0, 270.0), "tab:blue"))
    for stem, label in SCENES.items():
        fig, axes = plt.subplots(1, len(BANDS), figsize=(19, 5.6))
        for panel, (wavelength, band, unit, limits) in zip(axes, BANDS):
            frame = data[(data.scene == stem) & (data.wavelength_nm == wavelength)]
            shift = 273.15 if wavelength == 10500 else 0.0
            for plane, (positive, negative), color in planes:
                ordered = pd.concat([
                    frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                    frame[frame.vza == 0].assign(sign=1),
                    frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
                for key, name, style in methods:
                    panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, color=color,
                               label=f"{plane} | {name}", **style)
            panel.set(xlim=(-80, 80), ylim=limits, ylabel=unit, xlabel="带符号观测天顶角 (°)", title=band)
            panel.grid(alpha=0.3)
        axes[0].legend(fontsize=8, ncol=2)
        fig.suptitle(f"{label} | 太阳天顶角 30°，太阳方位角 0°", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(OUTPUT / f"scene_{stem}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
