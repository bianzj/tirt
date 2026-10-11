"""TIRTEB (plane + homogeneous canopy) against the SCOPE 2 verification run.

SCOPE 2's verification run (output/verificationdata) is 30 independent scenarios from a Latin hypercube
(input/dataset for_verification): LAI, leaf pigments, Vcmax, Ball-Berry slope, leaf angles, Rin, Ta, ea, sun zenith,
view zenith and relative azimuth vary; every other input is input/input_data_default.csv. Each scenario is run
through TIRTEB as one time step whose clock time reproduces the scenario's sun zenith (equator, equinox).
A second check feeds SCOPE's own canopy and soil temperatures to TiRT to separate the radiative transfer.

Usage: python compare_scope2.py [--scope <SCOPE2 folder containing SCOPE.m>]
Writes runs/scope2_tirteb/<nn>/..., runs/scope2_comparison.csv, runs/scope2_summary.csv, runs/scope2_comparison.png.
"""
import argparse
import csv
import math
import sys
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
warnings.simplefilter("ignore")

import run as tirt  # noqa: E402
from base.physicsF import inv_planck  # noqa: E402
from base.solar import solar_angles  # noqa: E402
from compare_scope import tirteb_series  # noqa: E402
from tirteb_run import run  # noqa: E402

DEFAULT_SCOPE = Path(r"C:\BaiduSyncdisk\app\modelling\SCOPE2-master\SCOPE2-master")
YEAR, DOY, LAT, LON = 2021, 80, 0.0, 0.0
KAPPA = 0.4


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig") as file:
        rows = [row for row in csv.reader(file) if row and not row[0].startswith("#")]
    header, body = rows[0], rows[1:]
    return {name: np.asarray([float(row[index]) for row in body]) for index, name in enumerate(header) if name}


def defaults(scope):
    values = {}
    for row in csv.reader((scope / "input" / "input_data_default.csv").open(encoding="utf-8-sig")):
        if len(row) > 1 and row[0] and row[1]:
            try:
                values[row[0]] = [float(item) for item in row[1:] if item != ""]
            except ValueError:
                continue
    return {key: (items[0] if len(items) == 1 else items) for key, items in values.items()}


def zo_and_d(lai, hc, par):
    """SCOPE 2 zo_and_d.m (Verhoef, McNaughton & Jacobs 1997)."""
    sq = math.sqrt(par["CD1"] * lai / 2.0)
    g1 = max(3.3, (par["CSSOIL"] + par["CR"] * lai / 2.0) ** -0.5)
    d = hc * (1.0 - (1.0 - math.exp(-sq)) / sq) if lai > 1e-7 and hc > 1e-7 else 0.0
    return (hc - d) * math.exp(-KAPPA * g1 + par["Psicor"]), d


def time_for_zenith(tts):
    """Morning time (decimal DOY, UTC at lon 0) at which TIRTEB's solar routine gives the zenith `tts`."""
    low, high = DOY + 0.25, DOY + 0.5
    if solar_angles(YEAR, DOY, high, LAT, LON)[0] >= tts:
        return high
    for _ in range(60):
        middle = 0.5 * (low + high)
        if solar_angles(YEAR, DOY, middle, LAT, LON)[0] > tts:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def scope_reference(scope):
    out = scope / "output" / "verificationdata"
    fluxes = read_csv(out / "fluxes.csv")
    radiation = read_csv(out / "radiation.csv")
    vegetation = read_csv(out / "vegetation.csv")
    wl = np.asarray([float(item) for item in (out / "wlS.txt").read_text().split()])
    spectra = np.loadtxt(out / "Lo_spectrum.csv", delimiter=",", comments="#")
    band = int(np.argmin(np.abs(wl - (10500.0 if wl.max() > 100 else 10.5))))
    optical = (wl >= 400.0) & (wl <= 2400.0)
    direct = np.trapezoid(np.loadtxt(out / "Esun.csv", delimiter=",", comments="#")[:, optical], wl[optical], axis=1)
    diffuse = np.trapezoid(np.loadtxt(out / "Esky.csv", delimiter=",", comments="#")[:, optical], wl[optical], axis=1)
    return {
        "diffuse_fraction": diffuse / (direct + diffuse),
        # SCOPE's ebal net radiation does not close with its own radiation output; both are compared.
        "Rn_budget": radiation["ShortIn"] + radiation["LongIn"] - radiation["HemisOutShort"] - radiation["HemisOutLong"],
        "Rn": fluxes["Rntot"], "LE": fluxes["lEtot"], "H": fluxes["Htot"], "G": fluxes["Gtot"],
        "Rn_canopy": fluxes["Rnctot"], "Rn_soil": fluxes["Rnstot"], "LE_canopy": fluxes["lEctot"], "LE_soil": fluxes["lEstot"],
        "Tcanopy": fluxes["Tcave"], "Tsoil": fluxes["Tsave"],
        "BT": np.asarray([float(inv_planck(10.5, value)) for value in spectra[:, band]]) - 273.15,
        "LST": vegetation["LST"] - 273.15,
    }


def soil_spectrum_file(scope, column, target):
    table = np.loadtxt(scope / "input" / "soil_spectra" / "soilnew.txt")
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(target, table[:, [0, column]], fmt="%.6f")
    return target


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", type=Path, default=DEFAULT_SCOPE)
    args = parser.parse_args()
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    par = defaults(args.scope)
    scenarios = read_csv(args.scope / "input" / "dataset for_verification" / "input_data_latin_hypercube_ts.csv")
    ref = scope_reference(args.scope)
    count = ref["Rn"].size
    out_dir = HERE / "runs"
    soil_file = soil_spectrum_file(args.scope, int(par["spectrum"]), out_dir / "scope2_soil_spectrum.txt")
    hc, leafwidth = par["hc"], par["leafwidth"]
    tparam = par["Tparam"]

    model = {key: np.full(count, np.nan) for key in ref if key != "diffuse_fraction"}
    rt_only = np.full(count, np.nan)
    geometry = []
    for index in range(count):
        s = {key: float(values[index]) for key, values in scenarios.items()}
        time = time_for_zenith(s["tts"])
        sza, saa = solar_angles(YEAR, DOY, time, LAT, LON)
        vaa = (saa + s["psi"]) % 360.0
        geometry.append((s["tts"], sza, s["tto"], s["psi"]))
        zo, d = zo_and_d(s["LAI"], hc, par)
        case = out_dir / "scope2_tirteb" / f"{index + 1:02d}"
        meteo = case / "meteo.txt"
        meteo.parent.mkdir(parents=True, exist_ok=True)
        meteo.write_text(f"1\n{time:.8f} {s['Ta']:.4f} {s['ea']:.4f} {par['p']:.2f} {par['u']:.4f} {s['Rin']:.4f} {par['Rli']:.4f}\n", encoding="utf-8")
        overrides = {
            "meteo": str(meteo), "year": YEAR, "doy": DOY, "start_time": time - 1e-6, "end_time": time + 1e-6,
            "latitude": LAT, "longitude": LON, "utc_offset_hours": 0, "dtime": 1800,
            "surface_model": "plane", "vegetation_model": "hom",
            "lai": s["LAI"], "canopy_height": hc, "leaf_width": leafwidth, "hspot": leafwidth / hc,
            "lidfa": s["LIDFa"], "lidfb": s["LIDFb"],
            "cab": s["Cab"], "cdm": s["Cdm"], "cw": s["Cw"], "cs": s["Cs"], "n": s["N"],
            "soil_optics_source": "file", "soil_optics_file": str(soil_file),
            "emissivity_leaf": 1.0 - par["rho_thermal"] - par["tau_thermal"], "emissivity_soil": 1.0 - par["rs_thermal"],
            "spectral_source": "direct", "wavelengths": "0.55;0.86;10.5",
            "geometry_mode": 0, "geometry_pair_mode": "zip", "vza": s["tto"], "vaa": vaa, "observation_time_mode": "all",
            "vcmax": s["Vcmo"], "ballberry_m": s["m"], "leaf_type": par["Type"], "kv": par["kV"], "rdparam": par["Rdparam"],
            "tparam1": tparam[0], "tparam2": tparam[1], "tparam3": tparam[2], "tparam4": tparam[3], "tparam5": tparam[4],
            "tyear": par["Tyear"], "beta": par["beta"], "knpqs": par["kNPQs"], "qls": par["qLs"],
            "stressfactor": par["stressfactor"], "stomatal_method": "farquhar", "stomatal_g0": par["BallBerry0"],
            "soil_rss_method": "constant", "soil_rss": par["rss"], "soil_rbs": par["rbs"], "soil_cs": par["cs"],
            "soil_rhos": par["rhos"], "soil_lambda": par["lambdas"], "soil_moisture": par["SMC"] / 100.0,
            "soil_heat_method": "surface",
            "z": par["z"], "ca": par["Ca"], "oa": par["Oa"],
            "aero_zo": zo, "aero_d": d, "aero_rbc": par["rb"], "aero_cr": par["CR"], "aero_cd": par["Cd"],
            "aero_cd1": par["CD1"], "aero_psicor": par["Psicor"], "aero_cssoil": par["CSSOIL"], "aero_rwc": par["rwc"],
            "aero_rbs": par["rbs"], "shortwave_diffuse_fraction": ref["diffuse_fraction"][index],
            "output_intermediate": "yes",
        }
        output = run(ROOT / "input_tirteb.csv", case / "output.csv", **overrides)
        series, _ = tirteb_series(output)
        for key in ("Rn", "LE", "H", "G", "Rn_canopy", "Rn_soil", "LE_canopy", "LE_soil", "Tcanopy", "Tsoil"):
            model[key][index] = series[key][0]
        model["Rn_budget"][index] = series["Rn"][0]
        with output.with_name("observations.csv").open(encoding="utf-8") as file:
            model["BT"][index] = float(next(csv.DictReader(file))["10p5"]) - 273.15
        model["LST"][index] = model["BT"][index]

        # Radiative transfer only: TiRT's homogeneous canopy at SCOPE's own canopy and soil temperatures.
        temperatures = {"lai": str(s["LAI"]), "hspot": str(leafwidth / hc), "hapke_enable": "0",
                        "soil_sunlit_k": str(ref["Tsoil"][index] + 273.15), "soil_shaded_k": str(ref["Tsoil"][index] + 273.15),
                        "leaf_sunlit_k": str(ref["Tcanopy"][index] + 273.15), "leaf_shaded_k": str(ref["Tcanopy"][index] + 273.15)}
        band = {"leaf_e": overrides["emissivity_leaf"], "soil_e": overrides["emissivity_soil"], "terrain_e": 1, "roof_e": 1, "wall_e": 1, "street_e": 1}
        _, temperature = tirt._surface_result(temperatures, "plane", "hom", 10.5, (s["tto"], s["psi"], s["tts"], 0.0), band)
        rt_only[index] = temperature - 273.15

    keys = ("Rn_budget", "Rn", "LE", "H", "G", "Rn_canopy", "Rn_soil", "LE_canopy", "LE_soil", "Tcanopy", "Tsoil", "BT")
    with (out_dir / "scope2_comparison.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.writer(sink)
        writer.writerow(["scenario", "tts", "tirteb_sza", "tto", "psi"] + [f"{key}_{source}" for key in keys for source in ("scope", "tirteb")]
                        + ["BT_tirt_with_scope_temperatures"])
        for index in range(count):
            writer.writerow([index + 1, *(f"{item:.2f}" for item in geometry[index])]
                            + [f"{series[key][index]:.3f}" for key in keys for series in (ref, model)] + [f"{rt_only[index]:.3f}"])
    # The hypercube samples Ta and ea independently; ea above saturation (condensing air) is not physical.
    saturation = 6.107 * 10.0 ** (7.5 * scenarios["Ta"] / (237.3 + scenarios["Ta"]))
    unsaturated = scenarios["ea"] <= saturation
    rows = []
    comparisons = [(key, model[key]) for key in keys] + [("BT_tirt_with_scope_temperatures", rt_only)]
    for subset, mask in (("all", np.ones(count, bool)), ("unsaturated", unsaturated)):
        for key, values in comparisons:
            reference = ref["BT"] if key.startswith("BT_") else ref[key]
            diff = values - reference
            good = np.isfinite(diff) & mask
            r = np.corrcoef(reference[good], values[good])[0, 1]
            rows.append({"variable": key, "subset": subset, "n": int(good.sum()), "bias": diff[good].mean(),
                         "rmse": np.sqrt((diff[good] ** 2).mean()), "max_abs": np.abs(diff[good]).max(), "r": r,
                         "scope_mean": reference[good].mean(), "tirteb_mean": values[good].mean()})
    with (out_dir / "scope2_summary.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in row.items()} for row in rows)
    for row in rows:
        print(f"{row['variable']:32s} {row['subset']:11s} n {row['n']:2d}  bias {row['bias']:+8.2f}  rmse {row['rmse']:7.2f}  max {row['max_abs']:7.2f}  "
              f"r {row['r']:5.2f}  SCOPE {row['scope_mean']:8.2f}  TIRTEB {row['tirteb_mean']:8.2f}")
    print("sun zenith reproduced to", max(abs(a - b) for a, b, _, _ in geometry), "deg")

    panels = (("Rn_budget", "净辐射，对 SCOPE 辐射收支 (W/m²)", model["Rn_budget"]), ("LE", "潜热 LE (W/m²)", model["LE"]), ("H", "感热 H (W/m²)", model["H"]),
              ("G", "土壤热通量 G (W/m²)", model["G"]), ("Tcanopy", "冠层温度 (°C)", model["Tcanopy"]),
              ("Tsoil", "土壤温度 (°C)", model["Tsoil"]), ("BT", "观测方向亮温 10.5 µm (°C)", model["BT"]),
              ("BT", "TiRT（输入 SCOPE 温度）亮温 (°C)", rt_only))
    fig, axes = plt.subplots(2, 4, figsize=(9.8, 5.6), dpi=100)
    for ax, (key, title, values) in zip(axes.flat, panels):
        reference = ref[key]
        name = key if values is not rt_only else "BT_tirt_with_scope_temperatures"
        stats = [row for row in rows if row["variable"] == name and row["subset"] == "unsaturated"][0]
        color = "tab:blue" if values is not rt_only else "tab:green"
        ax.scatter(reference[unsaturated], values[unsaturated], s=12, color=color)
        ax.scatter(reference[~unsaturated], values[~unsaturated], s=12, facecolors="none", edgecolors=color, linewidths=0.8)
        low = np.nanmin([reference.min(), np.nanmin(values)])
        high = np.nanmax([reference.max(), np.nanmax(values)])
        pad = 0.05 * (high - low)
        ax.plot([low - pad, high + pad], [low - pad, high + pad], color="0.5", lw=0.8, ls="--")
        ax.set_xlim(low - pad, high + pad)
        ax.set_ylim(low - pad, high + pad)
        ax.set_title(f"{title}\n未饱和：偏差 {stats['bias']:+.1f}，RMSE {stats['rmse']:.1f}，r {stats['r']:.2f}", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("SCOPE 2", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("TIRTEB / TiRT", fontsize=8)
    fig.suptitle(f"TIRTEB 与 SCOPE 2 验证算例（30 个拉丁超立方情景，均质冠层；实心 = 未饱和 {int(unsaturated.sum())} 个，空心 = ea 超过饱和）", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    docs = ROOT / "docs" / "validation"
    docs.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "scope2_comparison.png", dpi=100)
    fig.savefig(docs / "tirteb_scope2.png", dpi=100)
    (docs / "tirteb_scope2_summary.csv").write_bytes((out_dir / "scope2_summary.csv").read_bytes())
    print(docs / "tirteb_scope2.png")


if __name__ == "__main__":
    main()
