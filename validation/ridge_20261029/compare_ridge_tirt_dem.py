"""Compare the integrated TiRT surface_model=dem with the HiSTREAM multi-ridge reference."""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(r"C:\work\tirt")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

import run  # noqa: E402


def main():
    cfg = run.read_input(ROOT / "input.csv")
    cfg.update(surface_model="dem", vegetation_model="bare",
               dem_file=str(Path(r"C:\work\streamstudio\assets\dem-library\dune\transverse_dunes_100m.asc")),
               dem_slope_bin="2", dem_aspect_bin="10", dem_shadow="1",
               sky_temperature_k="none", hapke_enable="0",
               soil_sunlit_k="318.15", soil_shaded_k="303.15")
    spectra = run._spectra(run._band_table(cfg, ROOT / "input.csv"))
    spectra = {key: values[2] for key, values in spectra.items()}
    spectra["soil_e"] = 0.95

    reference = {}
    with (HERE / "histream_ridge.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            reference[(float(row["sza"]), float(row["saa"]), float(row["vza"]), float(row["vaa"]))] = float(row["bt_K"])
    rows = []
    for sza, saa in ((30.0, 270.0), (70.0, 270.0)):
        for vaa in (90.0, 270.0):
            for vza in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0):
                _, tirt = run._surface_result(cfg, "dem", "bare", 10.5, (vza, vaa, sza, saa), spectra)
                hist = reference[(sza, saa, vza, vaa)]
                rows.append(dict(sza=sza, saa=saa, vza=vza, vaa=vaa, tirt_bt=tirt, histream_bt=hist,
                                 difference=tirt - hist))
                print(f"SZA {sza:g} VAA {vaa:g} VZA {vza:g}: TiRT {tirt:8.4f}  HiSTREAM {hist:8.4f}  d {tirt-hist:+7.4f}")
    with (HERE / "ridge_compare.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for sza in (30.0, 70.0):
        for vaa in (90.0, 270.0):
            d = np.asarray([r["difference"] for r in rows if r["sza"] == sza and r["vaa"] == vaa])
            print(f"RMSE SZA {sza:g} VAA {vaa:g}: {np.sqrt(np.mean(d ** 2)):.4f} K  bias {d.mean():+.4f} K")


if __name__ == "__main__":
    main()
