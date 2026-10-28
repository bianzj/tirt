"""Compare the urban thermal scattering closures against HiSTREAM.

python validation/compare_urban_radiosity.py

legacy    : single-bounce urban exchange (Urban / Urban_Veg)
radiosity : street + four-wall radiosity, all orders (base/urban_radiosity.py)
Both share the same direct (first-order) geometric-optical term. Writes per-direction
values, summary statistics, an isothermal closure check and principal-plane figures to
validation/results/urban_radiosity/.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/urban_radiosity"
SCENES = (("10_T4_urban_bare", "城市（无植被）"), ("11_UV_urban_crown", "城市＋街道森林"))
METHODS = ("legacy", "radiosity")
TEMPERATURE_KEYS = tuple(f"{m}_{l}_k" for m in ("soil", "leaf", "terrain", "roof", "wall", "street")
                         for l in ("sunlit", "shaded"))
sys.path.insert(0, str(ROOT))

import run  # noqa: E402
from base.physicsF import planck  # noqa: E402


def brightness(config, method, geometry, spectra):
    local = dict(config, urban_transport=method)
    return run._surface_result(local, "urban", local["vegetation_model"], 10.5, geometry, spectra)[1]


def isothermal_emissivity(config, method, geometry, spectra, temperature=300.0):
    local = dict(config, urban_transport=method, **{key: temperature for key in TEMPERATURE_KEYS})
    radiance, _ = run._surface_result(local, "urban", local["vegetation_model"], 10.5, geometry, spectra)
    return radiance / float(planck(10.5, temperature))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    rows, summary, closure = [], [], []
    for stem, label in SCENES:
        scene = json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        config, spectra = scene["thermal_config"], scene["emissivity"]
        sza, saa = float(config["sza"]), float(config["saa"])
        targets = reference[(reference.scene == stem) & (reference.wavelength_nm == 10500)]
        values = {}
        for method in METHODS:
            start = time.perf_counter()
            values[method] = np.array([brightness(config, method, (r.vza, r.vaa, sza, saa), spectra)
                                       for r in targets.itertuples()])
            seconds = time.perf_counter() - start
            difference = values[method] - targets.reference.to_numpy()
            summary.append(dict(scene=stem, label=label, method=method, n=len(difference),
                                rmse_K=float(np.sqrt(np.mean(difference ** 2))), bias_K=float(difference.mean()),
                                max_abs_K=float(np.abs(difference).max()), seconds=seconds))
            for view in (0.0, 30.0, 60.0):
                eps = isothermal_emissivity(config, method, (view, 0.0, sza, saa), spectra)
                closure.append(dict(scene=stem, method=method, vza=view, emissivity=eps, closed=bool(eps <= 1.0 + 1e-9)))
        for i, r in enumerate(targets.itertuples()):
            rows.append(dict(scene=stem, vza=r.vza, vaa=r.vaa, histream_K=r.reference,
                             legacy_K=values["legacy"][i], radiosity_K=values["radiosity"][i]))
    data, summary, closure = pd.DataFrame(rows), pd.DataFrame(summary), pd.DataFrame(closure)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    closure.to_csv(OUTPUT / "isothermal_closure.csv", index=False)
    print(summary.to_string(index=False))
    print(closure.to_string(index=False))
    plot(data)


def plot(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    styles = (("histream_K", "HiSTREAM", "o-", "tab:orange"), ("legacy_K", "当前（单次交换）", "--", "tab:gray"),
              ("radiosity_K", "5×5 radiosity", "-", "tab:blue"))
    for row, (stem, label) in enumerate(SCENES):
        frame = data[data.scene == stem]
        for col, ((positive, negative), title) in enumerate(zip(((0.0, 180.0), (90.0, 270.0)),
                                                                 ("太阳主平面", "垂直太阳主平面"))):
            ordered = pd.concat([
                frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
                frame[frame.vza == 0].assign(sign=1),
                frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])
            panel = axes[row, col]
            for key, name, style, color in styles:
                panel.plot(ordered.vza * ordered.sign, ordered[key] - 273.15, style, color=color, label=name, lw=1.5, ms=3)
            panel.set(title=f"{label} | {title}", xlim=(-80, 80))
            panel.grid(alpha=0.3)
        axes[row, 0].set_ylabel("亮温 (°C)")
    for panel in axes[1]:
        panel.set_xlabel("带符号观测天顶角 (°)")
    axes[0, 0].legend()
    fig.suptitle("10.5 µm | SZA 30°, SAA 0°")
    fig.tight_layout()
    fig.savefig(OUTPUT / "urban_radiosity_comparison.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
