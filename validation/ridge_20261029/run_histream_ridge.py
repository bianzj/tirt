"""Run HiSTREAM eRaytracing on the transverse-dune heightfield for a small view grid.

Builds the case once (Node), then rewrites project.json for each (SZA, SAA, VZA, VAA),
runs the StreamStudio HiSTREAM engine, and reduces each TIF to mean radiance/BT.
"""
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = Path(r"C:\work\tirt")
ENGINE = Path(r"C:\work\streamstudio\models\bin_x64\Release\histream.exe")
NODE = "node"
C1 = 11910.439340652 * 1.0e4
C2 = 14388.291040407
WAVELENGTH_UM = 10.5
VIEWS = [(sza, 270.0, vza, vaa)
         for sza in (30.0, 70.0)
         for vaa in (90.0, 270.0)
         for vza in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0)]


def bt(radiance):
    return C2 / (WAVELENGTH_UM * np.log(C1 / (WAVELENGTH_UM ** 5 * radiance) + 1.0))


def build(out_dir):
    env = dict(os.environ, STREAMSIM_ROOT=r"C:\work\streamstudio")
    subprocess.run([NODE, str(HERE / "build_ridge.mjs"), "--out", str(out_dir), "--res", "256",
                    "--samples", "32", "--sza", "30", "--saa", "270", "--vza", "0", "--vaa", "0"],
                   check=True, cwd=HERE, env=env, capture_output=True, text=True)


def reduce_image(path):
    image = np.asarray(Image.open(path), dtype=np.float64)
    valid = np.isfinite(image) & (image > 0)
    return float(image[valid].mean()) if valid.any() else float("nan")


def main():
    out_dir = HERE / "histream"
    build(out_dir)
    input_dir = out_dir / "input"
    project_path = input_dir / "project.json"
    base = json.loads(project_path.read_text(encoding="utf-8"))
    rows = []
    for sza, saa, vza, vaa in VIEWS:
        project = json.loads(json.dumps(base))
        project["configuration"]["light"]["zenith"] = sza
        project["configuration"]["light"]["azimuth"] = saa
        project["configuration"]["sensor"]["vza"] = vza
        project["configuration"]["sensor"]["vaa"] = vaa
        view_dir = out_dir / "output" / f"sza{sza:g}_saa{saa:g}_vza{vza:g}_vaa{vaa:g}"
        view_dir.mkdir(parents=True, exist_ok=True)
        project["configuration"]["outDir"] = str(view_dir)
        project_path.write_text(json.dumps(project, ensure_ascii=False), encoding="utf-8")
        log_path = out_dir / f"engine_sza{sza:g}_saa{saa:g}_vza{vza:g}_vaa{vaa:g}.log"
        with log_path.open("w", encoding="utf-8", errors="replace") as log:
            result = subprocess.run([str(ENGINE), "eRaytracing", "project.json"], cwd=input_dir,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=900)
        if result.returncode != 0:
            raise RuntimeError(f"engine exit {result.returncode} for {sza}/{saa}/{vza}/{vaa}: {log_path}")
        images = sorted(view_dir.glob("SZA=*.tif"))
        if not images:
            raise RuntimeError(f"no image for {sza}/{saa}/{vza}/{vaa}: {view_dir}")
        radiance = reduce_image(images[-1])
        rows.append(dict(sza=sza, saa=saa, vza=vza, vaa=vaa, radiance=radiance, bt_K=bt(radiance)))
        print(f"SZA {sza:g} SAA {saa:g} VZA {vza:g} VAA {vaa:g}: {radiance:.6f}  {rows[-1]['bt_K']:.4f} K", flush=True)
    with (HERE / "histream_ridge.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
