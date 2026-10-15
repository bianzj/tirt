"""Run TIRT's exploratory VNIR approximation from a reproducible JSON input."""
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


def run(input_path, output_path=None):
    source = Path(input_path).resolve()
    config = json.loads(source.read_text(encoding="utf-8"))
    output = Path(output_path) if output_path else source.parent / "tirt_optical.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = list(simulate(config))
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output.resolve()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    print(run(args.input, args.output))
