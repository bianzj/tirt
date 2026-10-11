"""TiRT side of the TiRT-vs-StreamSim pilot benchmark: same angles, band, emissivities and temperatures."""

import csv
import sys
import tempfile
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

from run import run  # noqa: E402

VZA = ";".join(str(v) for v in range(0, 76, 5))
COMMON = {
    "geometry_mode": 0, "geometry_pair_mode": "cross", "vza": VZA, "vaa": "0;90;180;270",
    "sza": 30, "saa": 0, "wavelengths": "10.5", "spectral_source": "direct",
    "emissivity_leaf": 0.985, "emissivity_soil": 0.955, "emissivity_terrain": 0.950,
    "emissivity_roof": 0.950, "emissivity_wall": 0.930, "emissivity_street": 0.955,
    "hapke_enable": 0, "lai": 1.5, "hspot": 0.15,
    "crown_rad_a": 3.0, "crown_rad_b": 3.0, "crown_height": 3.0, "crown_density": 0.02,
    "building_density": 0.0015, "building_length": 10, "building_width": 10, "building_height": 20, "building_azimuth": 0,
    "row_width": 0.5, "row_blank": 0.5, "row_height": 1.0, "row_azimuth": 0,
}
TESTS = {
    "T1_plane_bare": {"surface_model": "plane", "vegetation_model": "bare"},
    "T2_plane_hom": {"surface_model": "plane", "vegetation_model": "hom"},
    "T3_plane_crown": {"surface_model": "plane", "vegetation_model": "crown"},
    "T4_urban_bare": {"surface_model": "urban", "vegetation_model": "bare"},
    "T5_row_az0": {"surface_model": "plane", "vegetation_model": "row", "row_azimuth": 0},
    "T5_row_az90": {"surface_model": "plane", "vegetation_model": "row", "row_azimuth": 90},
    "S1_slope_hom": {"surface_model": "slope", "vegetation_model": "hom", "slope_angle": 30, "slope_aspect": 90},
    "S1_slope_crown": {"surface_model": "slope", "vegetation_model": "crown", "slope_angle": 30, "slope_aspect": 90, "hspot": 0.025},
    "TV_terrain_crown": {"surface_model": "terrain", "vegetation_model": "crown", "terrain_height": 20, "terrain_radius": 30,
                         "terrain_density": 0.00015, "hspot": 0.025},
    "UV_urban_crown": {"surface_model": "urban", "vegetation_model": "crown", "hspot": 0.025},
}


def main():
    out_dir = HERE / "runs" / "tirt"
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="tirt-bench-"))
    for name, overrides in TESTS.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            result = run(ROOT / "input.csv", work / name / "output.csv", **{**COMMON, **overrides})
        with result.open(encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        target = out_dir / f"{name}.csv"
        with target.open("w", newline="", encoding="utf-8") as sink:
            writer = csv.writer(sink)
            writer.writerow(["test", "mode", "vza", "vaa", "sza", "saa", "radiance", "bt_K"])
            for row in rows:
                writer.writerow([name, "TiRT", row["vza"], row["vaa"], row["sza"], row["saa"], row["radiance"], row["brightness_temperature_K"]])
        print(f"{name}: {len(rows)} rows -> {target}")


if __name__ == "__main__":
    main()
