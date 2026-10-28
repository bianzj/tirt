"""Homogeneous canopy (HiSTREAM T2): standard model vs the new route vs the analytical reference.

python validation/compare_hom_route.py

standard  : two-layer CI hotspot + parametric thermal closure; optics K-integral + four-stream
route     : Beer-Lambert + K-integral hotspot + four-stream (thermal four-stream by reciprocity)
reference : K-integral + slab cross-section radiosity
Writes validation/results/hom_route/{directions,summary}.csv and a figure.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/hom_route"
STEM = "02_T2_plane_hom"
BANDS = ((670, "红光 670 nm", "BRF"), (860, "近红外 860 nm", "BRF"), (10500, "热红外 10.5 µm", "K"))
Y_LIMITS = {670: (0.0, 0.14), 860: (0.24, 0.50), 10500: (31.0, 40.0)}
G, SZA, SAA = 0.5, 30.0, 0.0
sys.path.insert(0, str(ROOT))

import run  # noqa: E402
from base import row_radiosity  # noqa: E402
from base.optical import reflectance  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402
from base.radiation import relative_azimuth  # noqa: E402
from base.vnir import diffuse_field, joint_transmission, spherical_coefficients  # noqa: E402


def k_fractions(lai, hspot, vza, raa):
    ks, ko = G / math.cos(math.radians(SZA)), G / math.cos(math.radians(min(vza, 89.999)))
    joint, path = joint_transmission(lai, ks, ko, hspot, SZA, min(vza, 89.999), raa)
    gap = math.exp(-ko * lai)
    return np.array([joint, gap - joint, ko * path, 1.0 - gap - ko * path])


def four_stream_scatter(lai, el, es, vza):
    """Effective emissivities of multiple scattering: absorptance of a beam from the view (reciprocity)."""
    ko = G / math.cos(math.radians(min(vza, 89.999)))
    gap = math.exp(-ko * lai)
    _, _, _, ground, top = diffuse_field(lai, ko, 1.0 - el, 0.0, 1.0 - es, 1.0, 0.0)
    soil, leaf = es * ground, (1.0 - top - es * (ground + gap)) - el * (1.0 - gap)
    sun_gap = math.exp(-G * lai / math.cos(math.radians(SZA)))
    sunlit = (1.0 - sun_gap) * math.cos(math.radians(SZA)) / (G * lai)
    return np.array([soil * sun_gap, soil * (1.0 - sun_gap), leaf * sunlit, leaf * (1.0 - sunlit)])


def thermal(scene, method, vza, vaa):
    config, emis = scene["thermal_config"], scene["emissivity"]
    if method == "standard":
        return run._surface_result(config, "plane", "hom", 10.5, (vza, vaa, SZA, SAA), emis)[1]
    lai, hspot = float(config["lai"]), float(config["hspot"])
    el, es = emis["leaf_e"], emis["soil_e"]
    parts = k_fractions(lai, hspot, vza, float(relative_azimuth(SAA, vaa)))
    direct = parts * np.array([es, es, el, el])
    if method == "route":
        weights = direct + four_stream_scatter(lai, el, es, vza)
    else:
        weights = direct + row_radiosity.scatter_emissivity(lai, 0.5, 0.0, 1.0, el, es, vza, vaa, SZA, SAA, 90.0,
                                                            soil_view=parts[0] + parts[1])
    radiance = float(np.sum(weights * planck(10.5, np.asarray(run._vegetation_temperatures(config)))))
    return float(inv_planck(10.5, radiance))


def optical(scene, method, wavelength, vza, vaa):
    settings = scene["optical_settings"]
    if method in ("standard", "route"):
        return reflectance(settings, wavelength, vza, vaa, SZA, SAA)
    band = settings["wavelengths"].index(wavelength)
    rho, tau, soil = (settings[key][band] for key in
                      ("leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"))
    lai, hspot, raa = float(settings["lai"]), float(settings["hspot"]), float(relative_azimuth(SAA, vaa))
    ks, ko, sob, sof = spherical_coefficients(SZA, float(vza), raa)
    joint, path = joint_transmission(lai, ks, ko, hspot, SZA, float(vza), raa)
    return (sob * rho + sof * tau) * path + soil * joint + \
        row_radiosity.shortwave_multiple(lai, 0.5, 0.0, 1.0, rho, tau, soil, vza, vaa, SZA, SAA, 90.0)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    scene = json.loads((BUNDLE / "inputs" / f"{STEM}.json").read_text(encoding="utf-8"))
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    rows, summary = [], []
    for wavelength, _, unit in BANDS:
        targets = reference[(reference.scene == STEM) & (reference.wavelength_nm == wavelength)]
        values = {}
        for method in ("reference", "standard", "route"):
            values[method] = np.array([thermal(scene, method, r.vza, r.vaa) if wavelength == 10500 else
                                       optical(scene, method, wavelength, r.vza, r.vaa) for r in targets.itertuples()])
            error = values[method] - targets.reference.to_numpy()
            summary.append(dict(wavelength_nm=wavelength, unit=unit, method=method,
                                rmse_histream=float(np.sqrt(np.mean(error ** 2))), bias_histream=float(error.mean()),
                                rmse_reference=float(np.sqrt(np.mean((values[method] - values["reference"]) ** 2)))))
        for i, r in enumerate(targets.itertuples()):
            rows.append(dict(wavelength_nm=wavelength, vza=r.vza, vaa=r.vaa, histream=r.reference,
                             **{method: values[method][i] for method in values}))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    with pd.option_context("display.float_format", "{:.4f}".format):
        print(summary.to_string(index=False))
    for wavelength, _, unit in BANDS:
        frame = data[data.wavelength_nm == wavelength][["histream", "standard", "route", "reference"]]
        shift = 273.15 if unit == "K" else 0.0
        print(f"{wavelength} range: {frame.min().min() - shift:.4f} .. {frame.max().max() - shift:.4f}")
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    styles = (("histream", "HiSTREAM 模拟", "o-", "tab:orange"), ("standard", "现有标准模型", "--", "tab:gray"),
              ("route", "新路线（K积分 + 四流）", "-", "tab:red"), ("reference", "解析基准", ":", "tab:blue"))
    fig, axes = plt.subplots(len(BANDS), 2, figsize=(12, 11.5), sharex=True, sharey="row")
    for row, (wavelength, band, unit) in enumerate(BANDS):
        frame = data[data.wavelength_nm == wavelength]
        for col, ((positive, negative), title) in enumerate(zip(((0.0, 180.0), (90.0, 270.0)),
                                                                 ("平行主平面", "垂直主平面"))):
            ordered = pd.concat([
                frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                frame[frame.vza == 0].assign(sign=1),
                frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
            panel = axes[row, col]
            shift = 273.15 if unit == "K" else 0.0
            for key, name, style, color in styles:
                panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, style, color=color, label=name, lw=1.4, ms=3)
            panel.set(title=f"均质冠层 | {title}", xlim=(-80, 80), ylim=Y_LIMITS[wavelength])
            panel.grid(alpha=0.3)
        axes[row, 0].set_ylabel(f"{band}\n" + ("亮温 (°C)" if unit == "K" else "BRF"))
    for panel in axes[-1]:
        panel.set_xlabel("带符号观测天顶角 (°)")
    axes[0, 0].legend()
    fig.suptitle("均质 T2 | LAI 1.5，太阳天顶角 30°，太阳方位角 0°")
    fig.tight_layout()
    fig.savefig(OUTPUT / "hom_route_comparison.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
