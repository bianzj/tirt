"""TIRTEB (plane + homogeneous canopy) against the SCOPE 1.73 verification time series.

SCOPE's verification run (output/verificationdata) is a 24 h, half-hourly simulation with the energy balance on.
Its meteorology and parameters (pars_and_input.dat) drive TIRTEB; fluxes, component temperatures and the nadir
10.5 um brightness temperature are compared.

Usage: python compare_scope.py [--scope <SCOPE_v1.73 folder>]
Writes runs/<files>, runs/scope_comparison.csv, runs/scope_summary.csv and runs/scope_comparison.png.
"""
import argparse
import csv
import re
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
warnings.simplefilter("ignore")

from base.physicsF import inv_planck  # noqa: E402
from tirteb_run import run  # noqa: E402

DEFAULT_SCOPE = Path(r"D:\BaiduNetdiskDownload\SCOPE1.73-master\SCOPE-master\SCOPE_v1.73")


def table(path):
    """Header row + numeric rows of a SCOPE .dat file (unit rows are skipped)."""
    lines = [line for line in Path(path).read_text(encoding="latin-1").splitlines() if line.strip()]
    rows = []
    for line in lines:
        try:
            rows.append([float(item) for item in line.split()])
        except ValueError:
            continue
    header = [] if not lines or re.match(r"^\s*[-+0-9.]", lines[0]) else lines[0].split()
    return header, np.asarray(rows)


def column(path, name):
    header, rows = table(path)
    return rows[:, header.index(name)]


def scope_case(scope):
    out = scope / "output" / "verificationdata"
    pars = {name: column(out / "pars_and_input.dat", name) for name in (
        "Rin", "Ta", "Rli", "p", "ea", "u", "tts", "LAI", "hc", "leafwidth", "Cab", "Cdm", "Cw", "Cs", "N",
        "rho_thermal", "tau_thermal", "rs_thermal", "Vcmo", "m", "kV", "Rdparam", "Tyear", "beta", "kNPQs", "qLs",
        "stressfactor", "rss", "rbs", "cs", "rhos", "lambdas", "SMC", "z", "Ca", "Oa", "zo", "d", "rb", "Cd", "CR",
        "CD1", "Psicor", "CSSOIL", "rwc", "LAT", "LON", "timezn", "LIDFa", "LIDFb", "BSMBrightness", "BSMlat", "BSMlon",
        "BallBerry0")}
    time = column(out / "fluxes.dat", "t")
    header, surf = table(out / "surftemp.dat")
    wl = np.concatenate(table(out / "wl.dat")[1]) if table(out / "wl.dat")[1].ndim > 1 else table(out / "wl.dat")[1]
    _, spectrum = table(out / "spectrum_obsdir_BlackBody.dat")
    band = int(np.argmin(np.abs(wl - 10500.0)))
    reference = {
        "time": time,
        "Rn": column(out / "fluxes.dat", "Rntot"), "LE": column(out / "fluxes.dat", "lEtot"),
        "H": column(out / "fluxes.dat", "Htot"), "G": column(out / "fluxes.dat", "Gtot"),
        "Rn_canopy": column(out / "fluxes.dat", "Rnctot"), "Rn_soil": column(out / "fluxes.dat", "Rnstot"),
        "LE_canopy": column(out / "fluxes.dat", "lEctot"), "LE_soil": column(out / "fluxes.dat", "lEstot"),
        # SCOPE stores the soil as Ts = [shaded; sunlit].
        "Tsoil_sunlit": surf[:, 5], "Tsoil_shaded": surf[:, 4], "Tcanopy": surf[:, 6], "Tsoil": surf[:, 7],
        "BT_nadir": np.asarray([float(inv_planck(10.5, value)) for value in spectrum[:, band]]) - 273.15,
        "sza": pars["tts"], "Rin": pars["Rin"],
    }
    return pars, reference


def tirteb_overrides(pars, meteo_path):
    first = {key: float(values[0]) for key, values in pars.items()}
    return {
        "meteo": str(meteo_path), "year": 2006, "doy": 169, "start_time": 169.0, "end_time": 170.0,
        "latitude": first["LAT"], "longitude": first["LON"], "utc_offset_hours": first["timezn"], "dtime": 1800,
        "surface_model": "plane", "vegetation_model": "hom",
        "lai": first["LAI"], "canopy_height": first["hc"], "leaf_width": first["leafwidth"],
        # SCOPE's hotspot parameter is leaf width / canopy height.
        "hspot": first["leafwidth"] / first["hc"], "lidfa": first["LIDFa"], "lidfb": first["LIDFb"],
        "cab": first["Cab"], "cdm": first["Cdm"], "cw": first["Cw"], "cs": first["Cs"], "n": first["N"],
        "spectral_source": "direct", "wavelengths": "0.55;0.86;10.5",
        "leaf_reflectance": f"-;-;{first['rho_thermal']}", "leaf_transmittance": f"-;-;{first['tau_thermal']}",
        "soil_reflectance": f"-;-;{first['rs_thermal']}",
        "direction_source": "direct", "vza": 0, "vaa": 0, "observation_time_mode": "all",
        "vcmax": first["Vcmo"], "ballberry_m": first["m"], "leaf_type": 0, "kv": first["kV"], "rdparam": first["Rdparam"],
        "tyear": first["Tyear"], "beta": first["beta"], "knpqs": first["kNPQs"], "qls": first["qLs"],
        "stressfactor": first["stressfactor"], "stomatal_method": "farquhar", "stomatal_g0": first["BallBerry0"],
        "soil_rss_method": "constant", "soil_rss": first["rss"], "soil_rbs": first["rbs"], "soil_cs": first["cs"],
        "soil_rhos": first["rhos"], "soil_lambda": first["lambdas"], "soil_moisture": first["SMC"],
        "soil_heat_method": "surface",
        "z": first["z"], "ca": first["Ca"], "oa": first["Oa"],
        "aero_zo": first["zo"], "aero_d": first["d"], "aero_rbc": first["rb"], "aero_cr": first["CR"],
        "aero_cd": first["Cd"], "aero_cd1": first["CD1"], "aero_psicor": first["Psicor"],
        "aero_cssoil": first["CSSOIL"], "aero_rwc": first["rwc"], "aero_rbs": first["rbs"],
        "bsm_brightness": first["BSMBrightness"], "bsm_latitude": first["BSMlat"], "bsm_longitude": first["BSMlon"],
        "output_intermediate": "yes",
    }


def write_meteo(pars, time, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(time))]
    for index, moment in enumerate(time):
        lines.append(f"{moment:.4f} {pars['Ta'][index]:.4f} {pars['ea'][index]:.4f} {pars['p'][index]:.2f} "
                     f"{pars['u'][index]:.4f} {pars['Rin'][index]:.3f} {pars['Rli'][index]:.3f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def tirteb_series(output):
    with output.open(encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    with output.with_name("intermediate.csv").open(encoding="utf-8") as file:
        middle = list(csv.DictReader(file))

    def series(rows_, key):
        return np.asarray([float(row[key]) if row.get(key) not in (None, "") else np.nan for row in rows_])

    # Soil fluxes are per unit soil area, weighted by the sunlit/shaded soil fraction; leaf fluxes are per unit
    # ground of the big leaf, weighted by the sunlit/shaded share of the leaf area (TIRTEB's scene aggregation).
    soil_w = np.stack([series(middle, "soil_sunlit_fraction"), series(middle, "soil_shaded_fraction")])
    leaf_w = np.stack([series(middle, "lai_sunlit"), series(middle, "lai_shaded")])
    lai = np.maximum(leaf_w.sum(axis=0), 1e-9)
    result = {
        "time": series(rows, "time"), "Rn": series(rows, "Rn"), "LE": series(rows, "LE"), "H": series(rows, "H"), "G": series(rows, "G"),
        "Tsoil_sunlit": series(rows, "Tss_C"), "Tsoil_shaded": series(rows, "Tsh_C"),
        "Tcanopy": (leaf_w[0] * series(rows, "Tls_C") + leaf_w[1] * series(rows, "Tlh_C")) / lai,
        "Tsoil": soil_w[0] * series(rows, "Tss_C") + soil_w[1] * series(rows, "Tsh_C"),
        "BT_nadir": series(rows, "Tv_C"), "sza": series(middle, "sza_global") if "sza_global" in middle[0] else series(middle, "sza"),
    }
    for flux, prefix in (("Rn", "N_"), ("LE", "LE_")):
        result[f"{flux}_soil"] = soil_w[0] * series(middle, prefix + "soilsunlit") + soil_w[1] * series(middle, prefix + "soilshaded")
        result[f"{flux}_canopy"] = (leaf_w[0] * series(middle, prefix + "leafsunlit") + leaf_w[1] * series(middle, prefix + "leafshaded")) / lai
    return result, middle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", type=Path, default=DEFAULT_SCOPE)
    args = parser.parse_args()
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    pars, ref = scope_case(args.scope)
    out_dir = HERE / "runs"
    meteo = out_dir / "scope_meteo.txt"
    write_meteo(pars, ref["time"], meteo)
    output = run(ROOT / "input_tirteb.csv", out_dir / "tirteb" / "output.csv", **tirteb_overrides(pars, meteo))
    model, middle = tirteb_series(output)

    keys = ("Rn", "LE", "H", "G", "Rn_canopy", "Rn_soil", "LE_canopy", "LE_soil",
            "Tsoil_sunlit", "Tsoil_shaded", "Tsoil", "Tcanopy", "BT_nadir", "sza")
    with (out_dir / "scope_comparison.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.writer(sink)
        writer.writerow(["time"] + [f"{key}_{source}" for key in keys for source in ("scope", "tirteb")])
        for index, moment in enumerate(ref["time"]):
            writer.writerow([f"{moment:.4f}"] + [f"{series[key][index]:.3f}" for key in keys for series in (ref, model)])
    day = ref["Rin"] > 0.0
    summary = []
    for key in keys:
        for label, mask in (("all", np.ones_like(day)), ("day", day), ("night", ~day)):
            d = model[key][mask] - ref[key][mask]
            d = d[np.isfinite(d)]
            summary.append({"variable": key, "period": label, "n": d.size, "bias": d.mean(), "rmse": np.sqrt((d * d).mean()),
                            "max_abs": np.abs(d).max(), "scope_mean": np.nanmean(ref[key][mask]), "tirteb_mean": np.nanmean(model[key][mask])})
    with (out_dir / "scope_summary.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in row.items()} for row in summary)
    for row in summary:
        print(f"{row['variable']:13s} {row['period']:5s} bias {row['bias']:+8.2f}  rmse {row['rmse']:7.2f}  "
              f"max {row['max_abs']:7.2f}  SCOPE {row['scope_mean']:8.2f}  TIRTEB {row['tirteb_mean']:8.2f}")

    hours = (ref["time"] - 169.0) * 24.0
    panels = (("Rn", "净辐射 (W/m²)"), ("LE", "潜热 LE (W/m²)"), ("H", "感热 H (W/m²)"), ("G", "土壤热通量 G (W/m²)"),
              ("Tsoil", "土壤温度 (°C)"), ("Tcanopy", "冠层温度 (°C)"), ("BT_nadir", "天底亮温 10.5 µm (°C)"))
    fig, axes = plt.subplots(2, 4, figsize=(9.8, 5.6), dpi=100, sharex=True)
    for ax, (key, title) in zip(axes.flat, panels):
        ax.plot(hours, ref[key], "o-", color="tab:orange", ms=2.5, lw=1.0, label="SCOPE 1.73")
        ax.plot(hours, model[key], "-", color="tab:blue", lw=1.6, label="TIRTEB")
        rows = [row for row in summary if row["variable"] == key and row["period"] == "all"][0]
        ax.set_title(f"{title}\n偏差 {rows['bias']:+.1f}，RMSE {rows['rmse']:.1f}", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        ax.grid(alpha=0.3)
    axes.flat[-1].axis("off")
    axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(), loc="center", fontsize=9, frameon=False)
    for ax in axes[1]:
        ax.set_xlabel("当地时间 (h)", fontsize=8)
    fig.suptitle("TIRTEB 与 SCOPE 1.73 验证时间序列（DOY 169，LAI 3，均质冠层）", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    docs = ROOT / "docs" / "validation"
    docs.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "scope_comparison.png", dpi=100)
    fig.savefig(docs / "tirteb_scope173.png", dpi=100)
    (docs / "tirteb_scope173_summary.csv").write_bytes((out_dir / "scope_summary.csv").read_bytes())
    print(docs / "tirteb_scope173.png")


if __name__ == "__main__":
    main()
