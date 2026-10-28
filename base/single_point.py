"""Single-scene energy balance with one directional thermal observation per time step.

python -m base.single_point --input input_tirteb.csv --vza 0 --vaa 0 --wavelength 10.5

Runs tirteb_run.run for one view direction and one band and writes one row per time step with the
component temperatures, fluxes and the TiRT directional brightness temperature.
"""
import argparse
import csv
from pathlib import Path

TIRTEB_ROOT = Path(__file__).resolve().parents[1]
COLUMNS = ("node", "time", "veg", "surface", "Ta_C", "Tss_C", "Tsh_C", "Tls_C", "Tlh_C",
           "Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
           "LE", "H", "G", "closed", "iterations", "energy_balance_residual_W_m2")


def run_single(input_path, output_path, vza=0.0, vaa=0.0, wavelength=10.5, **overrides):
    import tirteb_run
    from .io import _observation_band_table, read_input

    overrides = {key: str(item) for key, item in overrides.items() if item is not None}
    # The band table also sets the energy-balance emissivities, so it is used unchanged.
    config = dict(read_input(input_path), **overrides)
    bands = _observation_band_table(config, Path(input_path).resolve().parent)["wavelengths"]
    if not any(abs(band - float(wavelength)) < 1.0e-9 for band in bands):
        raise ValueError(f"wavelength {wavelength:g} um is not in the band table ({';'.join(f'{band:g}' for band in bands)}); add it to wavelengths")
    result = Path(output_path)
    full = tirteb_run.run(
        input_path, result.parent / "_tirteb" / "output.csv",
        direction_source="direct", vza=str(vza), vaa=str(vaa), observation_time_mode="all", **overrides)
    with full.open(encoding="utf-8") as file:
        steps = list(csv.DictReader(file))
    with (full.parent / "observations.csv").open(encoding="utf-8") as file:
        observed = {row["time"]: row for row in csv.DictReader(file)}
    tag = tirteb_run._wavelength_tag(wavelength)
    rows = []
    for step in steps:
        row = {key: step.get(key, "") for key in COLUMNS}
        observation = observed.get(step["time"], {})
        value = observation.get(tag, "")
        row.update(vza=vza, vaa=vaa, sza=observation.get("sza", ""), saa=observation.get("saa", ""),
                   wavelength_um=wavelength, brightness_temperature_K=value,
                   brightness_temperature_C=float(value) - 273.15 if value not in ("", None) else "")
        rows.append(row)
    result.parent.mkdir(parents=True, exist_ok=True)
    with result.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=TIRTEB_ROOT / "input_tirteb.csv")
    parser.add_argument("--output", type=Path, default=TIRTEB_ROOT / "cases" / "single_point" / "output.csv")
    parser.add_argument("--meteo", type=Path)
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--latitude", type=float)
    parser.add_argument("--longitude", type=float)
    parser.add_argument("--surface-model")
    parser.add_argument("--vegetation-model")
    parser.add_argument("--vza", type=float, default=0.0)
    parser.add_argument("--vaa", type=float, default=0.0)
    parser.add_argument("--wavelength", type=float, default=10.5)
    parser.add_argument("--max-iter", type=int)
    args = parser.parse_args()
    overrides = dict(meteo=args.meteo.resolve() if args.meteo else None, start_time=args.start_time,
                     end_time=args.end_time, latitude=args.latitude, longitude=args.longitude,
                     surface_model=args.surface_model, vegetation_model=args.vegetation_model, max_iter=args.max_iter)
    print(run_single(args.input, args.output, args.vza, args.vaa, args.wavelength, **overrides))


if __name__ == "__main__":
    main()
