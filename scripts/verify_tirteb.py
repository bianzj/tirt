"""Run a compact TIRTEB regression matrix and validate generated CSV files."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tirteb_run import OUTPUT_FIELDS, run, run_comparison


SURFACES = ("plane", "terrain", "urban")
VEGETATION = ("bare", "hom", "row", "crown")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise AssertionError(f"empty output: {path}")
    missing = [field for field in OUTPUT_FIELDS if field not in rows[0]]
    if missing:
        raise AssertionError(f"missing fields in {path}: {missing}")
    for row in rows:
        for key in ("Tv_C", "LE", "H", "G", "red", "nir"):
            if not math.isfinite(float(row[key])):
                raise AssertionError(f"non-finite {key} in {path}: {row[key]}")
    return rows


def read_intermediate(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise AssertionError(f"empty intermediate: {path}")
    return rows


def assert_urban_resistance(root: Path) -> None:
    rows = read_intermediate(root / "matrix_urban_hom" / "intermediate.csv")
    for row in rows:
        values = {
            key: float(row[key])
            for key in ("raa_roof", "raa_wall", "raa_street", "raa_leaf")
        }
        if not all(math.isfinite(value) and value > 0.0 for value in values.values()):
            raise AssertionError(f"invalid urban aerodynamic resistance: {values}")
        if values["raa_street"] <= values["raa_roof"]:
            raise AssertionError(f"canyon resistance must exceed roof resistance: {values}")
        if values["raa_leaf"] < values["raa_street"]:
            raise AssertionError(f"leaf resistance must include canyon and leaf boundary layers: {values}")
    print("PASS urban resistance ordering: roof < canyon surfaces <= canyon leaves")


def check(name: str, input_path: Path, root: Path, **overrides: object) -> int:
    output = root / name / "output.csv"
    observations = output.with_name("observations.csv")
    if output.exists() and observations.exists():
        try:
            existing = read_rows(output)
            with observations.open(newline="", encoding="utf-8-sig") as file:
                cached_observations = list(csv.DictReader(file))
            if cached_observations:
                print(f"PASS {name} (cached): series={len(existing)} observations={len(cached_observations)}")
                return len(existing)
        except (OSError, ValueError, AssertionError):
            pass
    result = Path(run(input_path, output, **overrides))
    rows = read_rows(result)
    observations = result.with_name("observations.csv")
    with observations.open(newline="", encoding="utf-8-sig") as file:
        observation_rows = list(csv.DictReader(file))
    if not observation_rows:
        raise AssertionError(f"empty observations: {observations}")
    print(f"PASS {name}: series={len(rows)} observations={len(observation_rows)}")
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "input_tirteb.csv")
    parser.add_argument("--output-root", type=Path, default=ROOT / "cases" / "_verify_tirteb")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_root = args.output_root.resolve()
    if output_root.exists() and not args.resume:
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    short = {"start_time": 214.0139, "end_time": 214.0764, "max_iter": 20}
    short_all = {**short, "observation_time_mode": "all"}
    check("default", input_path, output_root, **short_all)
    for surface in SURFACES:
        for vegetation in VEGETATION:
            check(
                f"matrix_{surface}_{vegetation}",
                input_path,
                output_root,
                surface_model=surface,
                vegetation_model=vegetation,
                **short_all,
            )

    check("terrain_single_slope", input_path, output_root, **{**short_all, "surface_model": "terrain", "terrain_type": 1, "terrain_slope": 30, "terrain_aspect": 45})
    check("terrain_compound", input_path, output_root, **{**short_all, "surface_model": "terrain", "terrain_type": 2, "terrain_radius": 60, "terrain_height": 25, "terrain_density": 1.2})
    check("direction_file", input_path, output_root, **{**short_all, "observation_direction_mode": "file", "direction_source": "file", "directions_file": "data/directions/observation_directions.txt"})
    check("spectral_file", input_path, output_root, **{**short_all, "observation_spectral_mode": "file", "spectral_source": "file", "spectrals_file": "data/spectrals/observation_spectrals.txt"})
    check("selected_times", input_path, output_root, **{**short, "observation_time_mode": "time", "observation_times": "214.0139;214.0556"})
    check("selected_hours", input_path, output_root, observation_time_mode="hour", observation_hours="0.0139", start_time=214.0, end_time=215.0, max_iter=20)
    check("intermediate", input_path, output_root, **{**short_all, "output_intermediate": "yes"})

    urban_bare = read_rows(output_root / "matrix_urban_bare" / "output.csv")
    urban_hom = read_rows(output_root / "matrix_urban_hom" / "output.csv")
    if all(abs(float(a["Tv_C"]) - float(b["Tv_C"])) < 1.0e-8 for a, b in zip(urban_bare, urban_hom)):
        raise AssertionError("urban+hom is identical to urban+bare; street canopy is inactive")
    if any(abs(float(row["Tv_C"])) > 100.0 for row in urban_hom):
        raise AssertionError("urban+hom brightness temperature is outside the physical test range")
    print("PASS urban vegetation separation and brightness range")
    assert_urban_resistance(output_root)

    comparison = run_comparison(input_path, output_root / "method_comparison", start_time=214.0139, end_time=214.0764, max_iter=15)
    comparison_rows = read_rows(Path(comparison))
    if len(comparison_rows) != 12 * 4:
        raise AssertionError(f"unexpected comparison rows: {len(comparison_rows)}")
    print(f"PASS method_comparison: rows={len(comparison_rows)}")
    print(f"PASS all TIRTEB checks: {output_root}")


if __name__ == "__main__":
    main()
