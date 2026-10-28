"""Basic route vs new route on HiSTREAM T2 (homogeneous) and T5 (rows), 670 nm, 860 nm, 10.5 um.

python validation/compare_routes.py

standard : current default model
basic    : Beer-Lambert / CI transmittance + K-integral hotspot (on the equivalent layer) + four-stream
new      : dual-scale hotspot (layer + correlation compensation) + node radiosity (scheduled nodes)
reference: exact geometry + cross-section radiosity
Writes validation/results/routes/{directions,summary}.csv and a figure.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/routes"
SCENES = (("02_T2_plane_hom", "均质", "hom", ("hom2", "hom3")),
          ("05_T5_row_az0", "垄行（垄向 0°）", "row", ("row2", "row3")),
          ("06_T5_row_az90", "垄行（垄向 90°）", "row", ("row2", "row3")))
NEW = {"hom": "hom3", "row": "row3"}
BANDS = ((670, "红光 670 nm", "BRF"), (860, "近红外 860 nm", "BRF"), (10500, "热红外 10.5 µm", "K"))
Y_LIMITS = {670: (0.0, 0.14), 860: (0.24, 0.50), 10500: (31.0, 40.0)}
G, SZA, SAA = 0.5, 30.0, 0.0
sys.path.insert(0, str(ROOT))

import run  # noqa: E402
from base import node_radiosity, row_closure, row_radiosity  # noqa: E402
from base.optical import reflectance  # noqa: E402
from base.physicsF import inv_planck, planck  # noqa: E402
from base.radiation import relative_azimuth  # noqa: E402
from base.vnir import canopy_reflectance, diffuse_field, joint_transmission, spherical_coefficients  # noqa: E402


def structure(scene, kind):
    s = scene["optical_settings"]
    lai, hspot = float(s["lai"]), float(s["hspot"])
    if kind == "hom":
        return lai, hspot, (0.5, 0.0, 1.0), 0.0
    return lai, hspot, (float(s["row_width"]), float(s["row_blank"]), float(s["row_height"])), float(s["row_azimuth"])


def gaps(lai, shape, zenith, azimuth, row_azimuth):
    return row_closure._gap(lai, *shape, row_closure._zenith(zenith), float(azimuth) - row_azimuth)


def k_fractions(lai, hspot, shape, row_azimuth, vza, vaa):
    """K-integral joint transmission on the CI-equivalent layer (Beer-Lambert when homogeneous)."""
    zv, zs = row_closure._zenith(vza), SZA
    pv, ps = gaps(lai, shape, zv, vaa, row_azimuth), gaps(lai, shape, zs, SAA, row_azimuth)
    kv, ks = -math.log(max(pv, 1e-12)) / lai, -math.log(max(ps, 1e-12)) / lai
    joint, path = joint_transmission(lai, ks, kv, hspot, zs, zv, float(relative_azimuth(SAA, vaa)))
    soil_sun, leaf_sun = min(joint, pv), min(kv * path, 1.0 - pv)
    return np.array([soil_sun, pv - soil_sun, leaf_sun, 1.0 - pv - leaf_sun]), pv, ps


def four_stream_thermal(lai, el, es, pv, ps, vza):
    """Four-stream slab scattering by reciprocity: absorptance of a beam from the view direction."""
    ko = G / math.cos(math.radians(row_closure._zenith(vza)))
    depth = max(-math.log(max(pv, 1e-12)) / ko, 1e-8)
    _, _, _, ground, top = diffuse_field(depth, ko, 1.0 - el, 0.0, 1.0 - es, 1.0, 0.0)
    soil, leaf = es * ground, (1.0 - top - es * (ground + pv)) - el * (1.0 - pv)
    sunlit = (1.0 - ps) * math.cos(math.radians(SZA)) / (G * lai)
    return np.array([soil * ps, soil * (1.0 - ps), leaf * sunlit, leaf * (1.0 - sunlit)])


def thermal(scene, kind, method, vza, vaa):
    config, emis = scene["thermal_config"], scene["emissivity"]
    if method == "standard":
        local = dict(config, vegetation_model=kind)
        return run._surface_result(local, "plane", kind, 10.5, (vza, vaa, SZA, SAA), emis)[1]
    lai, hspot, shape, row_az = structure(scene, kind)
    el, es = emis["leaf_e"], emis["soil_e"]
    E = np.array([es, es, el, el])
    if method == "basic":
        parts, pv, ps = k_fractions(lai, hspot, shape, row_az, vza, vaa)
        weights = parts * E + four_stream_thermal(lai, el, es, pv, ps, vza)
    elif method == "reference" and kind == "row":
        weights = row_closure.thermal_weights("reference", "reference", lai, hspot, *shape, el, es, vza, vaa, SZA, SAA, row_az)
    elif method == "reference":
        parts, _, _ = k_fractions(lai, hspot, shape, row_az, vza, vaa)
        weights = parts * E + row_radiosity.scatter_emissivity(lai, *shape, el, es, vza, vaa, SZA, SAA, 90.0,
                                                               soil_view=parts[0] + parts[1])
    else:
        parts = row_closure.fractions("conditioned", lai, hspot, *shape, vza, vaa, SZA, SAA, row_az)
        weights = parts * E + node_radiosity.thermal_scatter(method, lai, shape, el, es, vza, vaa, SZA, SAA, row_az)
    radiance = float(np.sum(weights * planck(10.5, np.asarray(run._vegetation_temperatures(config)))))
    return float(inv_planck(10.5, radiance))


def optical(scene, kind, method, wavelength, vza, vaa):
    settings = scene["optical_settings"]
    if method == "standard":
        return reflectance(dict(settings, vegetation_model=kind), wavelength, vza, vaa, SZA, SAA)
    band = settings["wavelengths"].index(wavelength)
    rho, tau, soil = (settings[key][band] for key in
                      ("leaf_reflectance_spectrum", "leaf_transmittance_spectrum", "soil_reflectance_spectrum"))
    lai, hspot, shape, row_az = structure(scene, kind)
    raa = float(relative_azimuth(SAA, vaa))
    if method == "basic":
        parts, _, ps = k_fractions(lai, hspot, shape, row_az, vza, vaa)
        return canopy_reflectance(lai, hspot, rho, tau, soil, SZA, float(vza), raa, parts, ps, False, 0.0, True, soil)
    _, ko, sob, sof = spherical_coefficients(SZA, float(vza), raa)
    if method == "reference":
        if kind == "row":
            return row_closure.optical_reflectance("reference", "reference", settings, SZA, SAA, vza, vaa, rho, tau, soil, 0.0)
        parts, _, _ = k_fractions(lai, hspot, shape, row_az, vza, vaa)
        return (sob * rho + sof * tau) / ko * parts[2] + soil * parts[0] + \
            row_radiosity.shortwave_multiple(lai, *shape, rho, tau, soil, vza, vaa, SZA, SAA, 90.0)
    parts = row_closure.fractions("conditioned", lai, hspot, *shape, vza, vaa, SZA, SAA, row_az)
    single = (sob * rho + sof * tau) / ko * parts[2] + soil * parts[0]
    return single + node_radiosity.shortwave_multiple(method, lai, shape, rho, tau, soil, vza, vaa, SZA, SAA, row_az)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    rows, summary = [], []
    for stem, label, kind, schedules in SCENES:
        scene = json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8"))
        for wavelength, _, unit in BANDS:
            targets = reference[(reference.scene == stem) & (reference.wavelength_nm == wavelength)]
            values = {}
            for method in ("reference", "standard", "basic") + schedules:
                values[method] = np.array([thermal(scene, kind, method, r.vza, r.vaa) if wavelength == 10500 else
                                           optical(scene, kind, method, wavelength, r.vza, r.vaa)
                                           for r in targets.itertuples()])
                error = values[method] - targets.reference.to_numpy()
                summary.append(dict(scene=stem, wavelength_nm=wavelength, method=method,
                                    rmse_histream=float(np.sqrt(np.mean(error ** 2))), bias_histream=float(error.mean()),
                                    rmse_reference=float(np.sqrt(np.mean((values[method] - values["reference"]) ** 2)))))
            for i, r in enumerate(targets.itertuples()):
                rows.append(dict(scene=stem, wavelength_nm=wavelength, vza=r.vza, vaa=r.vaa, histream=r.reference,
                                 standard=values["standard"][i], basic=values["basic"][i],
                                 new=values[NEW[kind]][i], reference=values["reference"][i]))
    data, summary = pd.DataFrame(rows), pd.DataFrame(summary)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False)
    with pd.option_context("display.float_format", "{:.4f}".format):
        for column in ("rmse_histream", "bias_histream", "rmse_reference"):
            print(f"--- {column} ---")
            table = summary.pivot_table(index=["wavelength_nm", "scene"], columns="method", values=column)
            print(table[[c for c in ("standard", "basic", "hom2", "hom3", "row2", "row3", "reference") if c in table]].to_string())
    plot(data)


RUNS = ROOT / "validation/streamsim/runs"
RUN_NAMES = {"02_T2_plane_hom": "T2_plane_hom", "05_T5_row_az0": "T5_row_az0", "06_T5_row_az90": "T5_row_az90"}


def with_thermal_runs(data, tag):
    """Replace the 10.5 um HiSTREAM column by fresh engine runs (validation/streamsim/runs/*<tag>.csv)."""
    data = data.copy()
    for stem, name in RUN_NAMES.items():
        runs = pd.read_csv(RUNS / f"{name}_eRaytracing{tag}.csv").set_index(["vza", "vaa"]).bt_K
        rows = (data.scene == stem) & (data.wavelength_nm == 10500)
        data.loc[rows, "histream"] = [runs[(r.vza, r.vaa)] for r in data[rows].itertuples()]
    return data


def plot(data, filename="routes_comparison.png", heading="HiSTREAM 与 双尺度热点 + 节点辐射度"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    styles = (("histream", "HiSTREAM 模拟", "o-", "tab:orange"),
              ("new", "双尺度热点 + 节点辐射度", "-", "tab:red"))
    fig, axes = plt.subplots(len(BANDS), 2 * len(SCENES), figsize=(26, 11.5), sharex=True, sharey="row")
    for row, (wavelength, band, unit) in enumerate(BANDS):
        for block, (stem, label, _, _) in enumerate(SCENES):
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
                    panel.plot(ordered.vza * ordered.sign, ordered[key] - shift, style, color=color, label=name, lw=1.3, ms=3)
                panel.set(title=f"{label} | {title}", xlim=(-80, 80), ylim=Y_LIMITS[wavelength])
                panel.grid(alpha=0.3)
        axes[row, 0].set_ylabel(f"{band}\n" + ("亮温 (°C)" if unit == "K" else "BRF"))
    for panel in axes[-1]:
        panel.set_xlabel("带符号观测天顶角 (°)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"{heading} | 太阳天顶角 30°，太阳方位角 0°")
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
