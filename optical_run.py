"""Run TIRT optical JSON scenes, with LAI/spectral/angular sweeps and plots.

python optical_run.py --input <scene.json>
python optical_run.py --input <scene.json> --no-plots
python optical_run.py --input <scene.json> --output-dir results/forest

Version 1 jobs produce spectrum.csv, angles.csv and optional PNGs. Previous
settings-based JSON files and their --output CSV option remain supported.
"""
import argparse
import csv
import json
from pathlib import Path

from base.optical import reflectance


def simulate(config):
    settings = config["settings"]
    sza, saa = config.get("sza", 30), config.get("saa", 0)
    maximum, step = config.get("max_vza", 60), config.get("vza_step", 5)
    if not 0 <= maximum < 90 or step <= 0:
        raise ValueError("Require 0 <= max_vza < 90 and vza_step > 0.")
    import numpy as np
    grid = np.unique(np.append(np.arange(0, maximum + 1e-9, step), maximum))
    for plane, offset in (("parallel", 0), ("perpendicular", 90)):
        for angle in np.r_[-grid[:0:-1], grid]:
            vaa = (saa + offset + (180 if angle < 0 else 0)) % 360
            for wave in settings["wavelengths"]:
                yield dict(scene=config.get("name", "optical"), plane=plane,
                           signed_vza=float(angle), vza=abs(float(angle)), vaa=vaa,
                           sza=sza, saa=saa, wavelength_nm=wave,
                           tirt_reflectance=reflectance(settings, wave, abs(float(angle)), vaa, sza, saa),
                           method="TIRT optical directional reflectance")


def load_config(input_path):
    """Read a versioned optical job (or the earlier settings JSON format)."""
    source = Path(input_path).resolve()
    config = json.loads(source.read_text(encoding="utf-8-sig"))
    if not isinstance(config, dict):
        raise ValueError("Optical JSON must contain an object")
    if "schema_version" in config or "domain" in config:
        if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
            raise ValueError("Optical JSON requires schema_version=1")
        if config.get("domain") != "optical":
            raise ValueError("optical_run.py requires domain='optical'")
        if "settings" in config:
            raise ValueError("Do not mix version 1 structure/materials with legacy settings")
        from base.scenario_runner import validate_case
        required = ("name", "label", "surface_model", "vegetation_model", "lai_values",
                    "wavelengths", "angle_wavelength", "sza", "saa", "spectrum_vza",
                    "spectrum_vaa", "signed_vza", "structure", "materials")
        missing = [key for key in required if key not in config]
        if missing:
            raise ValueError("Missing optical JSON fields: " + ", ".join(missing))
        validate_case(config)
    elif not isinstance(config.get("settings"), dict):
        raise ValueError("Use a version 1 optical job or a legacy settings object")
    return config


def run(input_path, output_path=None, *, output_dir=None, plots=True, show=False):
    """Run a JSON file; return spectrum.csv (v1) or the legacy output CSV path.

Version 1 defaults to <JSON directory>/results/<JSON filename without suffix>/.
An explicitly supplied output directory is relative to the caller's cwd.
"""
    source = Path(input_path).resolve()
    config = load_config(source)
    if "schema_version" in config:
        if output_path is not None:
            raise ValueError("Version 1 jobs use output_dir/--output-dir, not --output")
        from base.scenario_runner import run_case
        output = (Path(output_dir) if output_dir is not None else source.parent / "results" / source.stem).resolve()
        run_case(config, output, show, plots=plots, input_path=source)
        return output / "spectrum.csv"
    if output_dir is not None or show:
        raise ValueError("Legacy settings JSON uses --output; use version 1 for plots/output directories")
    output = Path(output_path) if output_path else source.parent / "tirt_optical.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = list(simulate(config))
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output.resolve()


def main(argv=None, *, default_input=None, default_output_dir=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input, required=default_input is None,
                        help="Version 1 optical JSON scene (required)")
    parser.add_argument("--output", type=Path, help="Legacy settings JSON: output CSV file")
    parser.add_argument("--output-dir", type=Path, help="Version 1 JSON: directory for CSVs, figures and metadata")
    parser.add_argument("--no-plots", action="store_true", help="Save numerical results without plotting")
    parser.add_argument("--show", action="store_true", help="Display figures after saving")
    args = parser.parse_args(argv)
    if args.show and args.no_plots:
        parser.error("--show cannot be combined with --no-plots")
    if args.output_dir is None and default_output_dir is not None:
        if "schema_version" in load_config(args.input):
            args.output_dir = default_output_dir
    print(run(args.input, args.output, output_dir=args.output_dir, plots=not args.no_plots, show=args.show))


if __name__ == "__main__":
    # 无内置场景：必须用 --input 指定光学 JSON。
    main()
