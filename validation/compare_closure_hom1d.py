"""canopy_closure standard / upgrade against HiSTREAM and a 1D homogeneous baseline (670 nm, 860 nm, 10.5 um).

python validation/compare_closure_hom1d.py

The 1D homogeneous baseline is the scene with the structure removed: the T2 plane homogeneous scene (T1 plane bare
for scenes without vegetation), which shares LAI, spectra and temperatures with every scene, with the scene's own
emissivities and the same closure as the model.
Writes validation/results/closure_vs_hom1d/{directions,summary}.csv and <closure>/scene_<stem>.png.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "validation/reference"
OUTPUT = ROOT / "validation/results/closure_vs_hom1d"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "validation"))

import run  # noqa: E402
from validate_current_model import optical  # noqa: E402

LABELS = {
    "01_T1_plane_bare": "平面裸土", "02_T2_plane_hom": "平面均质", "04_T3_plane_crown_leaf0.15": "平面离散树冠",
    "05_T5_row_az0": "垄行（垄向 0°）", "06_T5_row_az90": "垄行（垄向 90°）", "07_S1_slope_hom": "坡面均质",
    "08_S1_slope_crown": "坡面离散树冠", "09_TV_terrain_crown": "复合地形 + 树冠", "10_T4_urban_bare": "城市（无植被）",
    "11_UV_urban_crown": "城市 + 树冠",
}
BANDS = ((670, "红光 670 nm", "BRF"), (860, "近红外 860 nm", "BRF"), (10500, "热红外 10.5 µm", "亮温 (°C)"))
CLOSURES = {"standard": "标准版（双尺度热点 + 双层 + 节点辐射度）", "upgrade": "升级版（精确几何光学 + 节点辐射度）"}
FLAT = {False: "02_T2_plane_hom", True: "01_T1_plane_bare"}
TEMPERATURES = ("soil_sunlit_c", "soil_shaded_c", "leaf_sunlit_c", "leaf_shaded_c", "lai")


def value(inp, closure, wave, vza, vaa):
    config, settings = dict(inp["thermal_config"]), dict(inp["optical_settings"])
    config["canopy_closure"] = settings["canopy_closure"] = closure
    if wave == 10500:
        return run._surface_result(config, config["surface_model"], config["vegetation_model"], 10.5,
                                   (vza, vaa, 30.0, 0.0), inp["emissivity"])[1]
    return optical(settings, config, wave, vza, vaa)


def evaluate():
    reference = pd.read_csv(BUNDLE / "histream_reference.csv")
    comparison = pd.read_csv(ROOT / "validation/results/reference_comparison.csv")
    valid = comparison[comparison.status == "ok"][["scene", "wavelength_nm", "vza", "vaa"]]
    reference = reference.merge(valid, on=["scene", "wavelength_nm", "vza", "vaa"], validate="one_to_one")
    inputs = {stem: json.loads((BUNDLE / "inputs" / f"{stem}.json").read_text(encoding="utf-8")) for stem in LABELS}
    flats = {}
    for stem, inp in inputs.items():
        flat = inputs[FLAT[inp["thermal_config"]["vegetation_model"] == "bare"]]
        spectra = {k: v for k, v in flat["optical_settings"].items() if "spectrum" in k}
        assert all(inp["optical_settings"][k] == v for k, v in spectra.items()), stem
        assert all(float(inp["thermal_config"][k]) == float(flat["thermal_config"][k]) for k in TEMPERATURES), stem
        flats[stem] = dict(flat, emissivity=inp["emissivity"])
    rows = []
    for r in reference.itertuples():
        inp = inputs[r.scene]
        row = dict(scene=r.scene, wavelength_nm=r.wavelength_nm, vza=r.vza, vaa=r.vaa, histream=r.reference)
        for closure in CLOSURES:
            row[closure] = value(inp, closure, r.wavelength_nm, r.vza, r.vaa)
            row[f"hom1d_{closure}"] = value(flats[r.scene], closure, r.wavelength_nm, r.vza, r.vaa)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(data):
    rows = []
    for (scene, wave), group in data.groupby(["scene", "wavelength_nm"], sort=False):
        for closure in CLOSURES:
            for key, series in ((closure, "model"), (f"hom1d_{closure}", "hom1d")):
                e = group[key] - group.histream
                rows.append(dict(scene=scene, wavelength_nm=wave, closure=closure, series=series, n=len(e),
                                 rmse=float(np.sqrt(np.mean(e ** 2))), bias=float(e.mean()),
                                 max_abs=float(np.abs(e).max())))
    return pd.DataFrame(rows)


def ordered(frame, positive, negative):
    return pd.concat([
        frame[(frame.vaa == negative) & (frame.vza > 0)].sort_values("vza", ascending=False).assign(sign=-1),
        frame[frame.vza == 0].assign(sign=1),
        frame[(frame.vaa == positive) & (frame.vza > 0)].sort_values("vza").assign(sign=1)])


def plot(data, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"], "axes.unicode_minus": False})
    planes = (("太阳主平面", (0.0, 180.0), "tab:red"), ("垂直主平面", (90.0, 270.0), "tab:blue"))
    keys = ["histream"] + [k for c in CLOSURES for k in (c, f"hom1d_{c}")]
    for closure, closure_name in CLOSURES.items():
        (OUTPUT / closure).mkdir(parents=True, exist_ok=True)
        methods = (("histream", "HiSTREAM 模拟", dict(ls="none", marker="o", ms=4.5)),
                   (closure, closure_name.split("（")[0], dict(ls="-", lw=1.8)),
                   (f"hom1d_{closure}", "一维均质", dict(ls="--", lw=1.3, alpha=0.75)))
        for stem, label in LABELS.items():
            fig, axes = plt.subplots(1, len(BANDS), figsize=(19, 5.8))
            for panel, (wave, band, unit) in zip(axes, BANDS):
                frame = data[(data.scene == stem) & (data.wavelength_nm == wave)]
                shift = 273.15 if wave == 10500 else 0.0
                for plane, (positive, negative), color in planes:
                    line = ordered(frame, positive, negative)
                    for key, name, style in methods:
                        panel.plot(line.vza * line.sign, line[key] - shift, color=color,
                                   label=f"{plane} | {name}", **style)
                # Same y range for both closures so the standard and upgrade figures compare directly.
                values = frame[keys].to_numpy() - shift
                low, high = np.nanmin(values), np.nanmax(values)
                pad = max(0.08 * (high - low), 0.3 if wave == 10500 else 0.005)
                score = summary[(summary.scene == stem) & (summary.wavelength_nm == wave) &
                                (summary.closure == closure)].set_index("series").rmse
                fmt = "{:.3f} K" if wave == 10500 else "{:.4f}"
                panel.set(xlim=(-80, 80), ylim=(low - pad, high + pad), ylabel=unit, xlabel="带符号观测天顶角 (°)",
                          title=f"{band}\nRMSE {closure_name.split('（')[0]} {fmt.format(score['model'])}"
                                f" | 一维均质 {fmt.format(score['hom1d'])}")
                panel.grid(alpha=0.3)
            axes[0].legend(fontsize=7.5, ncol=1, loc="best")
            fig.suptitle(f"{label} | {closure_name} vs 一维均质 | 太阳天顶角 30°，太阳方位角 0°", fontsize=14)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(OUTPUT / closure / f"scene_{stem}.png", dpi=110)
            plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = evaluate()
    summary = summarize(data)
    data.to_csv(OUTPUT / "directions.csv", index=False)
    summary.to_csv(OUTPUT / "summary.csv", index=False, encoding="utf-8-sig")
    table = summary.pivot_table(index=["scene", "wavelength_nm"], columns=["closure", "series"], values="rmse", sort=False)
    with pd.option_context("display.float_format", "{:.4f}".format, "display.width", 200):
        print(table.to_string())
    plot(data, summary)


if __name__ == "__main__":
    main()
