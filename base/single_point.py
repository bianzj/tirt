import argparse
import csv
from pathlib import Path

import numpy as np

TIRTEB_ROOT = Path(__file__).resolve().parents[1]

from .brightness import urban_brightness, vegetation_brightness
from .modelio import read_defined, read_input, read_meteo
from .model import run_node


def model_name(settings):
    if settings.get("surface_model") in {"homogeneous", "row", "forest", "terrain", "urban"}:
        return settings["surface_model"]
    if settings["surface_type"] == 13 or settings["distribution"] == 4:
        return "urban"
    if settings["distribution"] == 2 or settings["surface_type"] in (12, 14):
        return "row"
    if settings["surface_type"] in (1, 2, 3, 4, 5):
        return "crown"
    return "hom"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=TIRTEB_ROOT / "input_tirteb.csv")
    parser.add_argument("--meteo", type=Path)
    parser.add_argument("--output", type=Path, default=TIRTEB_ROOT / "single_point_iterative_output.csv")
    parser.add_argument("--year", type=int)
    parser.add_argument("--doy", type=int)
    parser.add_argument("--lat", type=float)
    parser.add_argument("--lon", type=float)
    parser.add_argument("--dem", type=float, default=0.0)
    parser.add_argument("--slope", type=float, default=0.0)
    parser.add_argument("--vza", type=float, default=0.0)
    parser.add_argument("--vaa", type=float, default=0.0)
    parser.add_argument("--wavelength", type=float, default=10.5)
    parser.add_argument("--max-iter", type=int, default=50)
    args = parser.parse_args()

    config = read_input(args.input)
    config["_input_base"] = str(args.input.resolve().parent)
    if args.meteo is not None:
        config["meteo"] = str(args.meteo.resolve())
    settings = read_defined(None, config)
    meteo_path = Path(str(config["meteo"]))
    if not meteo_path.is_absolute():
        meteo_path = args.input.resolve().parent / meteo_path
    records = read_meteo(meteo_path)
    year = args.year if args.year is not None else int(float(config.get("year", 2024)))
    doy = args.doy if args.doy is not None else int(float(config.get("doy", 172)))
    lat = args.lat if args.lat is not None else float(config.get("latitude", 39.9))
    lon = args.lon if args.lon is not None else float(config.get("longitude", 116.4))
    name = model_name(settings)
    previous = None
    output = []

    for node, meteo in enumerate(records):
        state = run_node(settings, meteo, node, year, doy, lat, lon, args.dem, args.slope, args.vza, args.vaa, previous, args.max_iter)
        if name == "urban":
            radiance, bt = urban_brightness(settings, state, args.wavelength, args.vza, args.vaa)
        else:
            radiance, bt = vegetation_brightness(settings, state, args.wavelength, args.vza, args.vaa)
        row = {
            "node": node, "time": meteo["time"], "model": name,
            "surface_type": settings["surface_type"], "distribution": settings["distribution"],
            "sza": state["sza"], "saa": state["saa"], "vza": args.vza, "vaa": args.vaa,
            "iterations": state["iterations"], "closed": int(state["closed"]),
            "Tss_K": state["Tss"], "Tsh_K": state["Tsh"], "Tls_K": state["Tls"], "Tlh_K": state["Tlh"],
            "ustar": state.get("ustar", np.nan), "raa_leaf": state.get("raa_leaf", 0.0),
            "raa_soil": state.get("raa_soil", 0.0), "raa_roof": state.get("raa_roof", 0.0),
            "raa_wall": state.get("raa_wall", 0.0), "raa_street": state.get("raa_street", 0.0),
            "rssunlit": state.get("rssunlit", np.nan), "rsshaded": state.get("rsshaded", np.nan),
            "shortwave_direct": state.get("shortwave_direct", np.nan),
            "shortwave_diffuse": state.get("shortwave_diffuse", np.nan),
            "longwave_down": state.get("longwave_down", np.nan),
            "tir_scatter_soil": state.get("tir_scatter_soil", np.nan),
            "tir_scatter_leaf": state.get("tir_scatter_leaf", np.nan),
            "gap_view": state.get("gap_view", np.nan), "gap_sun": state.get("gap_sun", np.nan),
            "hotspot_leaf": state.get("hotspot_leaf", np.nan), "hotspot_soil": state.get("hotspot_soil", np.nan),
            "tirt_direct_soil_sun": state.get("tirt_direct_soil_sun", np.nan),
            "tirt_direct_soil_shade": state.get("tirt_direct_soil_shade", np.nan),
            "tirt_direct_leaf_sun": state.get("tirt_direct_leaf_sun", np.nan),
            "tirt_direct_leaf_shade": state.get("tirt_direct_leaf_shade", np.nan),
            "tirt_scatter_soil_sun": state.get("tirt_scatter_soil_sun", np.nan),
            "tirt_scatter_soil_shade": state.get("tirt_scatter_soil_shade", np.nan),
            "tirt_scatter_leaf_sun": state.get("tirt_scatter_leaf_sun", np.nan),
            "tirt_scatter_leaf_shade": state.get("tirt_scatter_leaf_shade", np.nan),
            "LE_leaf_sunlit": state.get("LE_leafsunlit", np.nan), "LE_leaf_shaded": state.get("LE_leafshaded", np.nan),
            "LE_soil_sunlit": state.get("LE_soilsunlit", np.nan), "LE_soil_shaded": state.get("LE_soilshaded", np.nan),
            "H_leaf_sunlit": state.get("H_leafsunlit", np.nan), "H_leaf_shaded": state.get("H_leafshaded", np.nan),
            "H_soil_sunlit": state.get("H_soilsunlit", np.nan), "H_soil_shaded": state.get("H_soilshaded", np.nan),
            "radiance": radiance, "brightness_temperature_K": bt, "brightness_temperature_C": bt - 273.15,
        }
        output.append(row)
        previous = state

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
