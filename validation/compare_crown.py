"""Discrete-forest scenes (HiSTREAM T3 plane / S1 slope crowns): hotspot x multiple-scattering combinations.

python validation/compare_crown.py

hotspot   : parametric (closed-form Poisson spheroid, legacy-style calibrated) or clumping_dual
scattering: crown3 / crown2 node radiosity, or the equivalent-homogeneous closure (thermal parametric
            closure with budget limit; optics four-stream slab on the crown solar gap)
standard  : current default model. Slopes are evaluated in the slope-local frame HiSTREAM renders.
Writes validation/results/crown/{directions,summary}.csv and one figure per scene.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation/results/crown"
SZA, SAA = 30.0, 0.0
SCENES = {"04_T3_plane_crown_leaf0.15": "平面离散树冠", "08_S1_slope_crown": "坡面离散树冠"}
METHODS = {"parametric+等效均质": ("parametric", "legacy"), "parametric+crown2": ("parametric", "crown2"),
           "双尺度双层+等效均质": ("dual_two_layer", "legacy"), "双尺度双层+crown2": ("dual_two_layer", "crown2"),
           "双尺度双层+crown3": ("dual_two_layer", "crown3"), "legacy+等效均质": ("legacy", "legacy")}
BANDS = ((670, "红光 670 nm", "BRF", (0.0, 0.20)), (860, "近红外 860 nm", "BRF", (0.10, 0.60)),
         (10500, "热红外 10.5 µm", "亮温 (°C)", (31.0, 43.0)))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation/streamsim"))

import run  # noqa: E402
import run_benchmark  # noqa: E402
from base import node_radiosity  # noqa: E402
from base.crown_parametric import crown_fractions_two_layer  # noqa: E402
from base.emissivity import canopy_scatter  # noqa: E402
from base.gap import gap_probability_crown_analytical  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402
from base.radiation import relative_azimuth  # noqa: E402
from base.thermal_crown import component_fractions  # noqa: E402
from base.vegetation import vegetation_proportions  # noqa: E402
from base.vnir import canopy_reflectance, spherical_coefficients  # noqa: E402


def crown(config):
    return dict(lai=float(config["lai"]), hspot=float(config["hspot"]), density=float(config["crown_density"]),
                vertical=float(config["crown_rad_b"]), horizontal=float(config["crown_rad_a"]),
                base=float(config["crown_height"]))


def fractions(c, geometry, vza, vaa, sza, saa):
    raa = float(relative_azimuth(saa, vaa))
    if geometry == "dual_two_layer":
        return crown_fractions_two_layer(c["lai"], c["density"], c["hspot"], c["vertical"], c["horizontal"], c["base"],
                                         min(abs(vza), 89.9), min(sza, 89.9), raa)[:, 0]
    if geometry == "legacy":
        return vegetation_proportions("crown", c["lai"], c["hspot"], vza, vaa, sza, saa, density=c["density"],
                                      crown_height=c["vertical"], crown_radius=c["horizontal"])[:, 0]
    settings = dict(lai=c["lai"], stand_density=c["density"], height=c["vertical"], width=c["horizontal"],
                    hspot=c["hspot"], crown_base=c["base"], thermal_crown_geometry=geometry)
    return component_fractions(settings, sza, vza, raa)[:, 0]


def local(config, vza, vaa):
    if config["surface_model"] != "slope":
        return vza, vaa, SZA, SAA
    frame = run_benchmark.slope_frame(float(config["slope_angle"]), float(config["slope_aspect"]))
    (v, a), (s, sa) = run_benchmark.to_local(vza, vaa, frame), run_benchmark.to_local(SZA, SAA, frame)
    return v, a, s, sa


def thermal(scene, method, vza, vaa):
    config, emis = scene["thermal_config"], scene["emissivity"]
    geometry, scatter = METHODS[method]
    c = crown(config)
    v, a, s, sa = local(config, vza, vaa)
    el, es = emis["leaf_e"], emis["soil_e"]
    direct = fractions(c, geometry, v, a, s, sa) * np.array([es, es, el, el])
    shape = (c["density"], c["vertical"], c["horizontal"])
    if scatter == "legacy":
        extra = np.asarray(canopy_scatter("crown", c["lai"], np.array([el]), np.array([es]), np.array([min(v, 89.9)]),
                                          np.array([min(s, 89.9)]), *shape), float).reshape(4)
        budget = max(0.0, 1.0 - float(direct.sum()))
        extra = extra * min(1.0, budget / extra.sum()) if extra.sum() > 0 else extra
    else:
        extra = node_radiosity.thermal_scatter(scatter, c["lai"], shape, el, es, v, a, s, sa)
    radiance = float(np.sum((direct + extra) * planck(10.5, np.asarray(run._vegetation_temperatures(config)))))
    return float(inv_planck(10.5, radiance))


def optical(scene, method, wavelength, vza, vaa):
    config, settings = scene["thermal_config"], scene["optical_settings"]
    geometry, scatter = METHODS[method]
    band = settings["wavelengths"].index(wavelength)
    rho, tau, soil = (settings[key][band] for key in
                      ("leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"))
    c = crown(config)
    v, a, s, sa = local(config, vza, vaa)
    parts = fractions(c, geometry, v, a, s, sa)
    raa = float(relative_azimuth(sa, a))
    shape = (c["density"], c["vertical"], c["horizontal"])
    if scatter == "legacy":
        solar_gap = float(gap_probability_crown_analytical(c["lai"], *shape, np.array([s]))[0])
        return canopy_reflectance(c["lai"], c["hspot"], rho, tau, soil, float(s), float(v), raa, parts, solar_gap,
                                  False, 0.0, True, soil)
    _, ko, sob, sof = spherical_coefficients(float(s), float(v), raa)
    single = (sob * rho + sof * tau) / ko * parts[2] + soil * parts[0]
    return single + node_radiosity.shortwave_multiple(scatter, c["lai"], shape, rho, tau, soil, v, a, s, sa)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    standard = pd.read_csv(ROOT / "validation/results/reference_comparison.csv")
    rows, summary = [], []
    for stem, label in SCENES.items():
        scene = json.loads((ROOT / "validation/reference/inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        for wavelength, *_ in BANDS:
            frame = standard[(standard.scene == stem) & (standard.wavelength_nm == wavelength) &
                             (standard.status == "ok")].dropna(subset=["reference"])
            values = {"standard": frame.value.to_numpy()}
            for method in METHODS:
                values[method] = np.array([thermal(scene, method, r.vza, r.vaa) if wavelength == 10500 else
                                           optical(scene, method, wavelength, r.vza, r.vaa) for r in frame.itertuples()])
            for method, v in values.items():
                e = v - frame.reference.to_numpy()
                summary.append(dict(scene=stem, wavelength_nm=wavelength, method=method,
                                    rmse=float(np.sqrt(np.mean(e ** 2))), bias=float(e.mean())))
            for i, r in enumerate(frame.itertuples()):
                rows.append(dict(scene=stem, wavelength_nm=wavelength, vza=r.vza, vaa=r.vaa, histream=r.reference,
                                 **{m: v[i] for m, v in values.items()}))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    order = ["standard"] + list(METHODS)
    with pd.option_context("display.float_format", "{:.4f}".format, "display.width", 220):
        for column in ("rmse", "bias"):
            print(f"--- {column} ---")
            print(summary.pivot_table(index=["scene", "wavelength_nm"], columns="method", values=column)[order].to_string())
    plot(data)


def plot(data, shown=("standard", "双尺度双层+等效均质", "双尺度双层+crown2")):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    names = {"histream": "HiSTREAM 模拟", "standard": "现有标准模型（parametric）",
             "parametric+等效均质": "parametric + 等效均质", "parametric+crown2": "parametric + 节点辐射度(2)",
             "双尺度双层+等效均质": "双尺度热点 + 双层 + 等效均质", "双尺度双层+crown2": "双尺度热点 + 双层 + 节点辐射度(2)",
             "双尺度双层+crown3": "双尺度热点 + 双层 + 节点辐射度(3)", "legacy+等效均质": "legacy（单尺度 + 双层）"}
    styles = [dict(ls="none", marker="o", ms=4.5), dict(ls="--", lw=1.4), dict(ls="-", lw=1.8), dict(ls=":", lw=1.8)]
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
                for key, style in zip(("histream",) + tuple(shown), styles):
                    panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, color=color,
                               label=f"{plane} | {names[key]}", **style)
            panel.set(xlim=(-80, 80), ylim=limits, ylabel=unit, xlabel="带符号观测天顶角 (°)", title=band)
            panel.grid(alpha=0.3)
        axes[0].legend(fontsize=7.5, ncol=2)
        fig.suptitle(f"{label} | 太阳天顶角 30°，太阳方位角 0°", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(OUTPUT / f"crown_{stem}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
