"""Model-level scenario forward calls, export and optional plotting.

Parameters may come from JSON or a Python dictionary; no GUI is required.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SURFACES = ("plane", "slope", "terrain", "urban")
VEGETATION = ("bare", "hom", "row", "crown")


def signed_geometry(angle, plane, solar_azimuth):
    """Positive angles face SAA (+90 for the perpendicular plane)."""
    if plane not in ("parallel", "perpendicular"):
        raise ValueError(f"Unknown principal plane: {plane}")
    offset = 0 if plane == "parallel" else 90
    return abs(angle), (solar_azimuth + offset + (180 if angle < 0 else 0)) % 360


def validate_case(case):
    domain = case["domain"]
    surface, vegetation = case["surface_model"], case["vegetation_model"]
    if domain not in ("optical", "thermal") or surface not in SURFACES or vegetation not in VEGETATION:
        raise ValueError("Unknown domain/surface/vegetation")
    for name in ("wavelengths", "lai_values", "signed_vza"):
        values = case[name]
        if not values or not all(math.isfinite(float(v)) for v in values) or len(set(values)) != len(values):
            raise ValueError(f"{name} must be nonempty, finite and unique")
    waves = case["wavelengths"]
    if waves != sorted(waves) or min(waves) <= 0:
        raise ValueError("wavelengths must be positive and increasing")
    if domain == "optical" and not all(400 <= v <= 2500 for v in waves):
        raise ValueError("Optical wavelengths must be 400..2500 nm")
    if domain == "thermal" and not all(3 <= v <= 100 for v in waves):
        raise ValueError("Thermal example wavelengths must be 3..100 um (check units)")
    if case["angle_wavelength"] not in waves:
        raise ValueError("angle_wavelength must occur in wavelengths")
    if min(case["lai_values"]) < 0 or (vegetation == "bare" and case["lai_values"] != [0]):
        raise ValueError("LAI must be nonnegative; bare scenes must use lai_values=[0]")
    if any(abs(v) >= 90 for v in case["signed_vza"]):
        raise ValueError("Signed view angles must be inside (-90, 90)")
    for name in ("sza", "saa", "spectrum_vza", "spectrum_vaa"):
        if not math.isfinite(case[name]):
            raise ValueError(f"{name} must be finite")
    if not 0 <= case["sza"] < 90 or not 0 <= case["spectrum_vza"] < 90:
        raise ValueError("Solar and spectral view zenith must be in [0, 90)")
    for name, values in case["materials"].items():
        if len(values) != len(waves) or not all(math.isfinite(v) and 0 <= v <= 1 for v in values):
            raise ValueError(f"{name}: supply one coefficient in [0, 1] per wavelength")
    if domain == "optical":
        material = case["materials"]
        if any(r + t > 1 for r, t in zip(material["leaf_reflectance_spectrum"], material["leaf_transmittance_spectrum"])):
            raise ValueError("Leaf reflectance + transmittance must not exceed 1")
    else:
        if any(not math.isfinite(v) or v <= -273.15 for v in case["temperatures"].values()):
            raise ValueError("Prescribed temperatures must be finite Celsius values above absolute zero")
        material = case["materials"]
        if (surface == "terrain" and "emissivity_terrain" in material
                and material["emissivity_terrain"] != material.get("emissivity_soil")):
            raise ValueError("emissivity_terrain is not supported separately: terrain facets use "
                             "emissivity_soil; remove it or set it equal to emissivity_soil")


def observations(case):
    spectrum = [dict(sweep="spectrum", plane="fixed", signed_vza=None,
                     vza=case["spectrum_vza"], vaa=case["spectrum_vaa"] % 360,
                     wavelength=wave) for wave in case["wavelengths"]]
    angles = []
    for plane in ("parallel", "perpendicular"):
        for angle in sorted(case["signed_vza"]):
            vza, vaa = signed_geometry(angle, plane, case["saa"])
            angles.append(dict(sweep="angles", plane=plane, signed_vza=angle,
                               vza=vza, vaa=vaa, wavelength=case["angle_wavelength"]))
    return spectrum, angles


def thermal_forward(case, lai, requests, directory):
    """Use run.run, preserving its input snapshots and omitted-view semantics."""
    from run import run, BackfaceOmittedWarning, _slope_backface, OUTPUT_FIELDS
    waves = sorted({row["wavelength"] for row in requests})
    views = list(dict.fromkeys((row["vza"], row["vaa"]) for row in requests))
    config = dict(case["structure"], **case["temperatures"], lai=lai,
                  surface_model=case["surface_model"], vegetation_model=case["vegetation_model"],
                  wavelengths=waves, sza=case["sza"], saa=case["saa"],
                  vza=[v[0] for v in views], vaa=[v[1] for v in views],
                  geometry_mode=0, geometry_pair_mode="zip", spectral_source="direct",
                  hapke_enable=0, sky_temperature_k="")
    indices = [case["wavelengths"].index(wave) for wave in waves]
    for key, values in case["materials"].items():
        values = [values[i] for i in indices]
        component = key.removeprefix("emissivity_")
        if component == key:
            config[key] = values
        elif component == "leaf":
            # Opaque-leaf Kirchhoff conversion: rho = 1 - emissivity, tau = 0.
            config["leaf_reflectance"] = [1.0 - item for item in values]
            config["leaf_transmittance"] = [0.0] * len(values)
        elif component != "terrain":
            config[f"{component}_reflectance"] = [1.0 - item for item in values]
    directory.mkdir(parents=True, exist_ok=True)
    source = directory / "input.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        for key, value in config.items():
            serialized = ";".join(map(str, value)) if isinstance(value, (list, tuple)) else value
            writer.writerow(["example", key, serialized, ""])
    if case["surface_model"] == "slope" and all(
            _slope_backface(case["structure"], (vza, vaa, case["sza"], case["saa"]))
            for vza, vaa in views):
        # run.run raises on an empty visible set. Keep this sweep as missing
        # so another visible sweep can still be calculated and plotted.
        with (directory / "output.csv").open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle).writerow(OUTPUT_FIELDS)
        (directory / "input_used.csv").write_bytes(source.read_bytes())
        return {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BackfaceOmittedWarning)
        run(source, directory / "output.csv")
    for warning in caught:
        if not issubclass(warning.category, BackfaceOmittedWarning):
            warnings.warn(warning.message, warning.category)
    with (directory / "output.csv").open(encoding="utf-8", newline="") as handle:
        raw = list(csv.DictReader(handle))
    return {(round(float(r["wavelength_um"]), 9), round(float(r["vza"]), 9), round(float(r["vaa"]), 9)):
            {key: float(r[key]) for key in ("radiance", "brightness_temperature_K", "brightness_temperature_C")}
            for r in raw}


def simulate_case(case, output_dir):
    validate_case(case)
    output_dir = Path(output_dir).resolve()
    requests = observations(case)
    results = ([], [])
    if case["domain"] == "optical":
        from base.optical import reflectance
        from base.solar import slope_direction
    else:
        from run import _slope_backface
    for lai_index, lai in enumerate(case["lai_values"]):
        settings = dict(case["structure"], **case["materials"], lai=lai,
                        surface_model=case["surface_model"], vegetation_model=case["vegetation_model"],
                        wavelengths=case["wavelengths"])
        for sweep_index, sweep in enumerate(requests):
            if case["domain"] == "thermal":
                lookup = thermal_forward(case, lai, sweep, output_dir / "thermal_inputs" /
                                         f"lai_{lai_index}" / ("spectrum" if sweep_index == 0 else "angles"))
            for obs in sweep:
                row = dict(scene=case["name"], domain=case["domain"], lai=lai, **obs,
                           wavelength_unit="nm" if case["domain"] == "optical" else "um",
                           sza=case["sza"], saa=case["saa"], status="ok")
                if case["domain"] == "optical":
                    hidden = case["surface_model"] == "slope" and not slope_direction(
                        obs["vza"], obs["vaa"], settings["slope_angle"], settings["slope_aspect"])[2]
                    row["reflectance_brf"] = None if hidden else float(reflectance(
                        settings, obs["wavelength"], obs["vza"], obs["vaa"], case["sza"], case["saa"]))
                    if hidden:
                        row["status"] = "slope_backface"
                else:
                    key = tuple(round(obs[k], 9) for k in ("wavelength", "vza", "vaa"))
                    if key in lookup:
                        row.update(lookup[key])
                    elif case["surface_model"] == "slope" and _slope_backface(
                            case["structure"], (obs["vza"], obs["vaa"], case["sza"], case["saa"])):
                        row.update(status="slope_backface", radiance=None,
                                   brightness_temperature_K=None, brightness_temperature_C=None)
                    else:
                        raise RuntimeError(f"Thermal backend omitted an expected observation: {key}")
                key = "reflectance_brf" if case["domain"] == "optical" else "brightness_temperature_K"
                if row["status"] == "ok" and (not math.isfinite(row[key]) or row[key] < 0):
                    raise FloatingPointError(f"Invalid model result: {row}")
                if case["domain"] == "thermal" and row["status"] == "ok":
                    if not math.isfinite(row["radiance"]) or row["radiance"] < 0:
                        raise FloatingPointError(f"Invalid radiance: {row}")
                results[sweep_index].append(row)
    return results


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_case(case, spectrum, angles, output, show=False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"],
                         "axes.unicode_minus": False, "font.size": 10})
    optical = case["domain"] == "optical"
    quantity = "reflectance_brf" if optical else "brightness_temperature_C"
    ylabel = "方向反射因子 BRF" if optical else "方向亮温 (°C)"
    unit = "nm" if optical else "μm"
    colors = ("#007c91", "#d56820", "#62499d", "#477a35", "#b44b69")
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    angular_fig, panels = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True, layout="constrained")
    for index, lai in enumerate(case["lai_values"]):
        style = dict(color=colors[index % len(colors)], marker="o", markersize=4,
                     linewidth=1.8, label=f"LAI = {lai:g}" if case["vegetation_model"] != "bare" else "无植被（LAI = 0）")
        spectral_rows = [r for r in spectrum if r["lai"] == lai]
        ax.plot([r["wavelength"] for r in spectral_rows],
                [float("nan") if r[quantity] is None else r[quantity] for r in spectral_rows], **style)
        for panel, plane in zip(panels, ("parallel", "perpendicular")):
            selected = [r for r in angles if r["lai"] == lai and r["plane"] == plane]
            panel.plot([r["signed_vza"] for r in selected],
                       [float("nan") if r[quantity] is None else r[quantity] for r in selected], **style)
    ax.set(xlabel=f"波长 ({unit})", ylabel=ylabel,
           title=f"{case['label']} · {'光学' if optical else '热红外'}\n"
                 f"SZA={case['sza']:g}°，VZA={case['spectrum_vza']:g}°，VAA={case['spectrum_vaa']:g}°")
    if not any(r["status"] == "ok" for r in spectrum):
        ax.text(0.5, 0.5, "该观测方向位于坡面背面，无有效值", transform=ax.transAxes, ha="center")
    for panel, title in zip(panels, ("太阳主平面", "垂直太阳主平面")):
        panel.set(xlabel="带符号观测天顶角 (°)", title=title)
        panel.axvline(0, color="#85929e", linewidth=0.6)
    panels[0].set_ylabel(ylabel)
    angular_fig.suptitle(f"{case['label']} · {case['angle_wavelength']:g} {unit} · SZA={case['sza']:g}°")
    for panel in (ax, *panels):
        panel.grid(alpha=0.22)
        panel.spines[["top", "right"]].set_visible(False)
        panel.legend(frameon=False, fontsize=9)
    fig.savefig(output / "spectrum.png", dpi=180)
    angular_fig.savefig(output / "angles.png", dpi=180)
    if show:
        plt.show()
    plt.close(fig)
    plt.close(angular_fig)


def run_case(case, output_dir=None, show=False, *, plots=True, input_path=None):
    """Run one example, save its two plots, CSVs and reproducibility metadata."""
    validate_case(case)
    if show and not plots:
        raise ValueError("show requires plots=True")
    output = (Path(output_dir) if output_dir is not None else
              ROOT / "cases" / case["domain"] / case["name"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # A failed rerun must not retain a success receipt from an earlier run.
    (output / "verification.json").unlink(missing_ok=True)
    start = time.perf_counter()
    spectrum, angles = simulate_case(case, output)
    write_csv(output / "spectrum.csv", spectrum)
    write_csv(output / "angles.csv", angles)
    (output / "parameters.json").write_text(json.dumps(case, ensure_ascii=False, indent=2), encoding="utf-8")
    if plots:
        plot_case(case, spectrum, angles, output, show)
    else:
        # Avoid retaining figures from a previous run with different parameters.
        for name in ("spectrum.png", "angles.png"):
            (output / name).unlink(missing_ok=True)
    files = [ROOT / ("optical_run.py" if case["domain"] == "optical" else "run.py"), Path(__file__)]
    if case["domain"] == "thermal":
        files.append(ROOT / "thermal_run.py")
    files += sorted((ROOT / "base").glob("*.py"))
    summary = dict(scene=case["name"], domain=case["domain"], energy_balance=False,
                   spectrum_rows=len(spectrum), angle_rows=len(angles),
                   valid_rows=sum(r["status"] == "ok" for r in spectrum + angles),
                   omitted_backface=sum(r["status"] != "ok" for r in spectrum + angles),
                   seconds=round(time.perf_counter() - start, 3), plots=plots,
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    if input_path is not None:
        source = Path(input_path).resolve()
        summary.update(input_file=str(source), input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (output / "verification.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"{case['domain']}/{case['name']}: {summary['valid_rows']} valid rows; {output}", flush=True)
    return summary


def example_main(case):
    parser = argparse.ArgumentParser(description=f"{case['label']} / {case['domain']} — 无能量平衡")
    parser.add_argument("--output-dir", type=Path, help="输出目录；默认 cases/<domain>/<scene>")
    parser.add_argument("--show", action="store_true", help="保存图片后弹出绘图窗口")
    args = parser.parse_args()
    return run_case(case, args.output_dir, args.show)
