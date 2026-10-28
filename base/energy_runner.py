"""JSON adaptation and time-series export; numerical physics stays in tirteb_run."""
import csv
import hashlib
import json
import math
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
TEMPERATURE_FIELDS = ("Ta_C", "Tss_C", "Tsh_C", "Tls_C", "Tlh_C")
URBAN_FIELDS = ("Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C",
                "Twall_shaded_C", "Tstreet_sunlit_C", "Tstreet_shaded_C",
                "Tleaf_sunlit_C", "Tleaf_shaded_C")
CLOSURE_FIELDS = ("closed", "iterations", "energy_balance_residual_W_m2", "L", "ustar",
                  "stability_closed", "stability_iterations", "stability_residual_m_inv")


def _csv_value(value):
    if isinstance(value, list):
        if any(isinstance(item, list) for item in value):
            raise ValueError("Energy parameter arrays must be flat lists")
        return ";".join(_csv_value(item) for item in value)
    if isinstance(value, bool):
        return "1" if value else "0"
    if value is None or isinstance(value, dict):
        raise ValueError("Energy parameters must be strings, finite numbers, booleans or flat lists")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Energy parameters must be finite")
    return str(value)


def _validate_time_range(parameters):
    start, end = (float(parameters[key]) for key in ("start_time", "end_time"))
    if not math.isfinite(start) or not math.isfinite(end) or start > end:
        raise ValueError("start_time and end_time must be finite, with start_time <= end_time")


def load_config(input_path):
    config = json.loads(Path(input_path).read_text(encoding="utf-8-sig"))
    if not isinstance(config, dict):
        raise ValueError("Energy-balance JSON must contain an object")
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise ValueError("Energy-balance JSON requires schema_version=1")
    if config.get("domain") != "energy_balance":
        raise ValueError("tirteb_run.py requires domain='energy_balance'")
    if not isinstance(config.get("parameters"), dict):
        raise ValueError("Energy-balance JSON requires a parameters object")
    parameters = config["parameters"]
    for key in ("name", "label"):
        if not isinstance(config.get(key), str) or not config[key].strip():
            raise ValueError(f"Energy-balance JSON requires {key}")
    for key in ("surface_model", "vegetation_model", "meteo", "esun_file", "esky_file",
                "optical_coefficients_file", "start_time", "end_time"):
        if key not in parameters or parameters[key] == "":
            raise ValueError(f"Missing energy-balance parameter: {key}")
    if parameters["surface_model"] not in ("plane", "slope", "terrain", "urban"):
        raise ValueError("surface_model must be plane, slope, terrain or urban")
    if parameters["vegetation_model"] not in ("bare", "hom", "row", "crown"):
        raise ValueError("vegetation_model must be bare, hom, row or crown")
    for key, value in parameters.items():
        if key.startswith("_") or key != key.strip().lower():
            raise ValueError(f"Use lowercase public parameter names: {key}")
        if key in ("case_result_dir", "cases_dir", "output_intermediate", "plots", "show"):
            raise ValueError(f"Output option {key} belongs in Python, not the JSON parameters")
        _csv_value(value)
    _validate_time_range(parameters)
    return config


def _write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows, temperature_fields, output, label, closed_rows, show):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = [float(row["time"]) for row in rows]
    figures = []
    for stem, fields, ylabel in (
        ("temperatures", temperature_fields, "Temperature (°C)"),
        ("fluxes", ("Rn", "LE", "H", "G"), "Energy flux (W m$^{-2}$)"),
    ):
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
        for field in fields:
            values = [float(row[field]) if row[field] else math.nan for row in rows]
            if any(math.isfinite(value) for value in values):
                ax.plot(x, values, marker=".", linewidth=1.4, label=field)
        ax.set(xlabel="Time (day of year, local standard time)", ylabel=ylabel,
               title=f"{label} | energy balance closed: {closed_rows}/{len(rows)}")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
        fig.savefig(output / f"{stem}.png", dpi=160)
        figures.append(fig)
    if show:
        plt.show()
    for fig in figures:
        plt.close(fig)


def run_job(source, solver, file_keys, *, output_dir=None, plots=True, show=False, overrides=None):
    source = Path(source).resolve()
    config = load_config(source)
    if show and not plots:
        raise ValueError("show requires plots=True")
    parameters = dict(config["parameters"])
    parameters.update({key: value for key, value in (overrides or {}).items() if value is not None})
    if "case_result_dir" in parameters or "cases_dir" in parameters:
        raise ValueError("JSON jobs use output_dir to choose the output location")
    parameters = {key: _csv_value(value) for key, value in parameters.items()}
    # Reference paths are relative to the JSON, never the process working directory.
    for key in file_keys:
        if parameters.get(key):
            parameters[key] = str((source.parent / parameters[key]).resolve())
    for key in ("meteo", "esun_file", "esky_file", "optical_coefficients_file"):
        if not Path(parameters[key]).is_file():
            raise FileNotFoundError(f"{key}: {parameters[key]}")
    _validate_time_range(parameters)
    output = (Path(output_dir) if output_dir is not None else source.parent / "results" / source.stem).resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "verification.json").unlink(missing_ok=True)
    for name in ("temperatures.png", "fluxes.png"):
        (output / name).unlink(missing_ok=True)
    parameters["case_result_dir"] = str(output)
    parameters.setdefault("output_intermediate", "no")
    start = time.perf_counter()
    # Keep an executable CSV snapshot for existing APIs and method comparisons.
    effective = output / "job_input.csv"
    _write_csv(effective, ("section", "parameter", "value"),
               [dict(section="json", parameter=key, value=value) for key, value in parameters.items()])
    result = solver(effective, output / "output.csv")
    with result.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    temperature_fields = TEMPERATURE_FIELDS + (URBAN_FIELDS if parameters["surface_model"] == "urban" else ())
    _write_csv(output / "temperatures.csv", ("node", "time") + temperature_fields + CLOSURE_FIELDS, rows)
    _write_csv(output / "fluxes.csv", ("node", "time", "Rn", "LE", "H", "G") + CLOSURE_FIELDS, rows)
    closed_rows = sum(row["closed"].strip().lower() in ("true", "1") for row in rows)
    if plots:
        _plot(rows, temperature_fields, output, config["name"], closed_rows, show)
    (output / "parameters.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "effective_parameters.json").write_text(json.dumps(parameters, ensure_ascii=False, indent=2), encoding="utf-8")
    files = [ROOT / "tirteb_run.py", *sorted((ROOT / "base").glob("*.py"))]
    receipt = dict(domain="energy_balance", scene=config["name"], energy_balance=True,
                   rows=len(rows), closed_rows=closed_rows, all_closed=closed_rows == len(rows),
                   max_energy_balance_residual_W_m2=max(float(r["energy_balance_residual_W_m2"]) for r in rows),
                   plots=plots, seconds=round(time.perf_counter() - start, 3),
                   input_file=str(source), input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                   data_sha256={key: hashlib.sha256(Path(parameters[key]).read_bytes()).hexdigest()
                                for key in sorted(file_keys) if parameters.get(key) and Path(parameters[key]).is_file()})
    (output / "verification.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(f"energy_balance/{config['name']}: {len(rows)} rows, {closed_rows} closed; {output}", flush=True)
    return result
