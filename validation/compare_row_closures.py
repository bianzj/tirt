"""Unified row closures against the analytical reference and HiSTREAM (670 nm, 860 nm, 10.5 um).

python validation/compare_row_closures.py

Each closure is a (hotspot, scattering) pair from base/row_closure.py, identical in all bands:
legacy        : CI transmittance + CI layer hotspot + equivalent-homogeneous scattering
conditioned   : legacy + hotspot correlation compensation (scattering as legacy)
reference     : parametric geometry + row cross-section radiosity
exact_hotspot : reference hotspot + equivalent-homogeneous scattering (isolates scattering error)
exact_scatter : conditioned hotspot + radiosity scattering (isolates hotspot error)
Writes validation/results/row_closures/{directions,summary}.csv and a figure.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/row_closures"
SCENES = (("05_T5_row_az0", "垄行（垄向 0°）"), ("06_T5_row_az90", "垄行（垄向 90°）"))
BANDS = ((670, "红光 670 nm", "BRF"), (860, "近红外 860 nm", "BRF"), (10500, "热红外 10.5 µm", "K"))
Y_LIMITS = {670: (0.0, 0.14), 860: (0.24, 0.50), 10500: (31.0, 40.0)}
METHODS = {"reference": ("reference", "reference"), "legacy": ("legacy", "legacy"),
           "conditioned": ("conditioned", "conditioned"), "exact_hotspot": ("reference", "legacy"),
           "exact_scatter": ("conditioned", "reference")}
sys.path.insert(0, str(ROOT))

import run  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402
from base.row_closure import optical_reflectance, thermal_weights  # noqa: E402


def predict(scene, method, wavelength, vza, vaa):
    hotspot, scatter = METHODS[method]
    config, settings = scene["thermal_config"], scene["optical_settings"]
    shape = [float(settings[key]) for key in ("row_width", "row_blank", "row_height", "row_azimuth")]
    if wavelength == 10500:
        weights = thermal_weights(hotspot, scatter, float(settings["lai"]), float(settings["hspot"]), *shape[:3],
                                  scene["emissivity"]["leaf_e"], scene["emissivity"]["soil_e"],
                                  vza, vaa, 30.0, 0.0, shape[3])
        radiance = float(np.sum(weights * planck(10.5, np.asarray(run._vegetation_temperatures(config)))))
        return float(inv_planck(10.5, radiance))
    band = settings["wavelengths"].index(wavelength)
    rho, tau, soil = (settings[key][band] for key in
                      ("leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"))
    return optical_reflectance(hotspot, scatter, settings, 30.0, 0.0, vza, vaa, rho, tau, soil,
                               float(settings["shortwave_diffuse_fraction"]))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    rows, summary = [], []
    for stem, label in SCENES:
        scene = json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        for wavelength, band, unit in BANDS:
            targets = reference[(reference.scene == stem) & (reference.wavelength_nm == wavelength)]
            values = {}
            for method in METHODS:
                start = time.perf_counter()
                values[method] = np.array([predict(scene, method, wavelength, r.vza, r.vaa) for r in targets.itertuples()])
                seconds = time.perf_counter() - start
                to_histream = values[method] - targets.reference.to_numpy()
                summary.append(dict(scene=stem, wavelength_nm=wavelength, unit=unit, method=method,
                                    rmse_histream=float(np.sqrt(np.mean(to_histream ** 2))),
                                    bias_histream=float(to_histream.mean()),
                                    rmse_reference=float(np.sqrt(np.mean((values[method] - values["reference"]) ** 2))),
                                    seconds=seconds))
            for i, r in enumerate(targets.itertuples()):
                rows.append(dict(scene=stem, wavelength_nm=wavelength, vza=r.vza, vaa=r.vaa, histream=r.reference,
                                 **{method: values[method][i] for method in METHODS}))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    with pd.option_context("display.float_format", "{:.4f}".format):
        print(summary.to_string(index=False))
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    styles = (("histream", "HiSTREAM 模拟", "o-", "tab:orange"), ("legacy", "一般闭合", "--", "tab:gray"),
              ("conditioned", "条件化闭合", "-", "tab:red"), ("reference", "解析基准", ":", "tab:blue"))
    fig, axes = plt.subplots(len(BANDS), 4, figsize=(20, 11.5), sharex=True, sharey="row")
    for row, (wavelength, band, unit) in enumerate(BANDS):
        for block, (stem, label) in enumerate(SCENES):
            frame = data[(data.scene == stem) & (data.wavelength_nm == wavelength)]
            for offset, ((positive, negative), title) in enumerate(zip(((0.0, 180.0), (90.0, 270.0)),
                                                                        ("平行主平面", "垂直主平面"))):
                ordered = pd.concat([
                    frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                    frame[frame.vza == 0].assign(sign=1),
                    frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
                panel = axes[row, 2 * block + offset]
                shift = 273.15 if unit == "K" else 0.0
                for key, name, style, color in styles:
                    panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, style, color=color, label=name, lw=1.4, ms=3)
                panel.set(title=f"{label} | {title}", xlim=(-80, 80), ylim=Y_LIMITS[wavelength])
                panel.grid(alpha=0.3)
        axes[row, 0].set_ylabel(f"{band}\n" + ("亮温 (°C)" if unit == "K" else "BRF"))
    for panel in axes[-1]:
        panel.set_xlabel("带符号观测天顶角 (°)")
    axes[0, 0].legend()
    fig.suptitle("垄行 T5 | 太阳天顶角 30°，太阳方位角 0°")
    fig.tight_layout()
    fig.savefig(OUTPUT / "row_closures_comparison.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
