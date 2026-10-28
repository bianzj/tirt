"""TIRT thermal JSON scenes: prescribed temperatures, LAI, spectra and angles.

python thermal_run.py
python thermal_run.py --input <scene.json>
Legacy CSV inputs still use the run.py solver and write under cases/<case_name>/.
No energy balance is solved.
"""
import argparse
import json
from pathlib import Path

from run import DEFAULT_INPUT, run as _run_csv


def load_config(input_path):
    """Read and validate a thermal scene without creating output files."""
    config = json.loads(Path(input_path).read_text(encoding="utf-8-sig"))
    if not isinstance(config, dict):
        raise ValueError("Thermal JSON must contain an object")
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise ValueError("Thermal JSON requires schema_version=1")
    if config.get("domain") != "thermal":
        raise ValueError("thermal_run.py requires domain='thermal'")
    required = ("name", "label", "surface_model", "vegetation_model", "lai_values",
                "wavelengths", "angle_wavelength", "sza", "saa", "spectrum_vza",
                "spectrum_vaa", "signed_vza", "structure", "materials", "temperatures")
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError("Missing thermal JSON fields: " + ", ".join(missing))
    from base.scenario_runner import validate_case
    validate_case(config)
    return config


def run(input_path=DEFAULT_INPUT, output_path=None, model=None, report=None, *,
        output_dir=None, plots=True, show=False, **overrides):
    """JSON: save CSVs and optional figures; CSV: preserve the legacy solver API."""
    source = Path(input_path).resolve()
    if source.suffix.lower() != ".json":
        if output_dir is not None or show:
            raise ValueError("CSV inputs use output_path/--output; JSON jobs support output_dir and plots")
        return _run_csv(source, output_path, model, report, **overrides)
    config = load_config(source)
    if output_path is not None or model is not None or overrides:
        raise ValueError("JSON jobs use output_dir and model parameters inside the JSON")
    from base.scenario_runner import run_case
    output = (Path(output_dir) if output_dir is not None else source.parent / "results" / source.stem).resolve()
    summary = run_case(config, output, show, plots=plots, input_path=source)
    if report is not None:
        report.update(summary)
    return output / "spectrum.csv"


def main(argv=None, *, default_input=None, default_output_dir=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input or DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, help="Legacy CSV output file")
    parser.add_argument("--model", choices=("bare", "hom", "row", "crown"), help="Legacy CSV model override")
    parser.add_argument("--output-dir", type=Path, help="JSON output directory")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args(argv)
    if args.show and args.no_plots:
        parser.error("--show cannot be combined with --no-plots")
    if args.output_dir is None and args.input.suffix.lower() == ".json":
        args.output_dir = default_output_dir
    print(run(args.input, args.output, args.model, output_dir=args.output_dir,
              plots=not args.no_plots, show=args.show))


if __name__ == "__main__":
    # 默认读取根目录 input.csv；JSON 场景用 --input 指定。
    main()
