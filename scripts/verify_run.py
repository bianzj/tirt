from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from run import OUTPUT_FIELDS, run


SURFACES = ("plane", "slope", "terrain", "urban")
VEGETATIONS = ("bare", "hom", "row", "crown")


def _row_count(path: Path) -> int:
    with path.open(newline="", encoding="utf-8-sig") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise RuntimeError(f"{path} has no data rows")
    missing = [field for field in OUTPUT_FIELDS if field not in rows[0]]
    if missing:
        raise RuntimeError(f"{path} missing output fields: {', '.join(missing)}")
    return len(rows)


def _check(name: str, input_path: Path, output_root: Path, **overrides: object) -> None:
    output_path = output_root / name / "output.csv"
    result = Path(run(input_path, output_path, **overrides))
    rows = _row_count(result)
    print(f"OK {name}: {rows} rows -> {result}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify TiRT run.py and input.csv combinations")
    parser.add_argument("--input", type=Path, default=Path("input.csv"))
    parser.add_argument("--output-root", type=Path, default=Path("cases/_verify_run"))
    args = parser.parse_args()

    input_path = args.input.resolve()
    output_root = args.output_root.resolve()

    _check("default", input_path, output_root)

    for surface in SURFACES:
        for vegetation in VEGETATIONS:
            _check(
                f"matrix_{surface}_{vegetation}",
                input_path,
                output_root,
                surface_model=surface,
                vegetation_model=vegetation,
                geometry_mode=0,
                vza="0;30",
                vaa="0;180",
            )

    _check("geometry_mode_1", input_path, output_root, geometry_mode=1, principal_vza="0;30;60")
    _check(
        "geometry_mode_2",
        input_path,
        output_root,
        geometry_mode=2,
        geometry_file="data/directions/observation_directions.txt",
    )
    _check(
        "geometry_mode_3",
        input_path,
        output_root,
        geometry_mode=3,
        principal_vza_step=15,
        principal_vza_max=30,
        include_hemisphere=1,
    )
    _check(
        "spectral_files",
        input_path,
        output_root,
        spectral_source="file",
        vegetation_spectrum_source="file",
        soil_spectrum_source="file",
        urban_spectrum_source="file",
    )
    _check(
        "zip_angles",
        input_path,
        output_root,
        geometry_mode=0,
        geometry_pair_mode="zip",
        vza="0;30;60",
        vaa="0;90;180",
    )


if __name__ == "__main__":
    main()
