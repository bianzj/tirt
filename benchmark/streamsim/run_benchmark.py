"""StreamSim/HiSTREAM side of the TiRT pilot benchmark: build, run, and reduce each case.

Usage (from this folder):
    python run_benchmark.py                 # all default cases
    python run_benchmark.py --cases T2:eRaytracing T2:eVoxelRT --samples 32
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
SNAPSHOT = RUNS / "_engine" / "histream.exe"
# Prefer the frozen engine snapshot (exe + DLLs + shader + defined) so a rebuild of
# the StreamSim tree cannot change results mid-benchmark.
ENGINE = Path(os.environ.get("HISTREAM_EXE") or (SNAPSHOT if SNAPSHOT.exists()
              else r"C:\work\streamsim\models\bin_x64\Release\histream.exe"))
NODE = os.environ.get("NODE", "node")
C1 = 11910.439340652 * 1.0e4  # W um^4 m-2 sr-1 (TiRT c1 x 1e4)
C2 = 14388.291040407          # um K
WAVELENGTH_UM = 10.5
DEFAULT_CASES = ["T1:eRaytracing", "T2:eRaytracing", "T2:eVoxelRT", "T3:eRaytracing", "T3:eVoxelRT", "T4:eRaytracing"]
NAMES = {"T1": "T1_plane_bare", "T2": "T2_plane_hom", "T3": "T3_plane_crown", "T4": "T4_urban_bare",
         "T5az0": "T5_row_az0", "T5az90": "T5_row_az90", "S1hom": "S1_slope_hom", "S1crown": "S1_slope_crown",
         "TV": "TV_terrain_crown", "UV": "UV_urban_crown"}
# Slope cases reuse a flat-scene builder; the scene is rendered in slope-local coordinates.
SLOPE_ANGLE, SLOPE_ASPECT = 30.0, 90.0
SLOPE_TESTS = {"S1hom": ("T2", []), "S1crown": ("T3", ["--leaf", "0.15"])}
TIFF_NAME = re.compile(r"SZA=([\d.]+)_SAA=([\d.]+)_VZA=([\d.]+)_VAA=([\d.]+)(?:_[a-z])?\.tif$")


def brightness_temperature(radiance):
    return C2 / (WAVELENGTH_UM * math.log(C1 / (WAVELENGTH_UM ** 5 * radiance) + 1.0))


def run_engine(mode, input_dir, log_path, timeout):
    command = [str(ENGINE), mode, str(input_dir / "project.json")]
    with open(log_path, "w", encoding="utf-8", errors="replace") as log:
        process = subprocess.Popen(command, cwd=input_dir, stdout=log, stderr=subprocess.STDOUT)
        try:
            return process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            subprocess.run(["taskkill", "/pid", str(process.pid), "/T", "/F"], capture_output=True)
            process.wait()
            raise RuntimeError(f"engine timed out after {timeout}s: {command}")


def reduce_images(output_dir, keep_zero=False):
    """Mean radiance per view. Zero pixels are treated as missing unless keep_zero (cold/black decomposition runs,
    where surfaces legitimately render as exactly 0)."""
    rows = []
    for path in sorted(output_dir.glob("SZA=*.tif")):
        match = TIFF_NAME.search(path.name)
        if not match:
            continue
        sza, saa, vza, vaa = (float(value) for value in match.groups())
        image = np.asarray(Image.open(path), dtype=np.float64)
        valid = np.isfinite(image) & ((image >= 0) if keep_zero else (image > 0))
        values = image[valid]
        rows.append({
            "vza": vza, "vaa": vaa, "sza": sza, "saa": saa,
            "radiance": float(values.mean()),
            "pixel_std": float(values.std()), "valid_fraction": float(valid.mean()),
            "half_mean_difference": float(abs(values[: values.size // 2].mean() - values[values.size // 2:].mean())),
        })
    return rows


def write_case_csv(test, mode, rows, target):
    nadir = [row for row in rows if row["vza"] == 0]
    expanded = [row for row in rows if row["vza"] != 0]
    # The engine renders nadir once; TiRT lists it for every azimuth.
    for vaa in (0.0, 90.0, 180.0, 270.0):
        for row in nadir:
            expanded.append({**row, "vaa": vaa})
    expanded.sort(key=lambda row: (row["vaa"], row["vza"]))
    with target.open("w", newline="", encoding="utf-8") as sink:
        writer = csv.writer(sink)
        writer.writerow(["test", "mode", "vza", "vaa", "sza", "saa", "radiance", "bt_K"])
        for row in expanded:
            writer.writerow([NAMES[test], mode, f"{row['vza']:g}", f"{row['vaa']:g}", f"{row['sza']:g}", f"{row['saa']:g}",
                             f"{row['radiance']:.6f}", f"{brightness_temperature(row['radiance']):.4f}"])
    return expanded


def inserted_voxels(log_path):
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    return sum(int(value) for value in re.findall(r"OBJ placement \d+ .*?inserted=(\d+)", text))


def unit(zenith, azimuth):
    """Direction (north, east, up) for a zenith and compass azimuth in degrees."""
    z, a = math.radians(zenith), math.radians(azimuth)
    return np.array([math.sin(z) * math.cos(a), math.sin(z) * math.sin(a), math.cos(z)])


def slope_frame(slope, aspect):
    """Orthonormal (north', east', normal) basis of a plane tilted by `slope` towards compass azimuth `aspect` (TiRT slope1)."""
    normal = unit(slope, aspect)
    north = np.array([1.0, 0.0, 0.0]) - normal[0] * normal
    if np.linalg.norm(north) < 1e-9:
        north = np.array([0.0, 1.0, 0.0]) - normal[1] * normal
    north /= np.linalg.norm(north)
    return north, np.cross(normal, north), normal


def to_local(zenith, azimuth, frame):
    d = unit(zenith, azimuth)
    north, east, normal = frame
    x, y, z = d @ north, d @ east, d @ normal
    return math.degrees(math.acos(max(-1.0, min(1.0, z)))), math.degrees(math.atan2(y, x)) % 360.0


def principal_views(saa=0.0):
    return [(float(abs(a)), (saa + offset + (180.0 if a < 0 else 0.0)) % 360.0) for offset in (0.0, 90.0) for a in range(-75, 76, 5)
            if not (a == 0 and offset == 90.0)]


def run_slope_case(test, mode, args, extra, tag=""):
    """Infinite slope rendered as the equivalent flat scene in slope-local coordinates, one engine run per view."""
    base_test, build_extra = SLOPE_TESTS[test]
    name = f"{NAMES[test]}_{mode}{tag}"
    run_dir = RUNS / name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    build = [NODE, str(HERE / "build_case.mjs"), "--test", base_test, "--mode", mode, "--out", str(run_dir),
             "--samples", str(args.samples), "--seed", str(args.seed), "--res", str(args.res), "--depth", str(args.depth),
             *build_extra, *extra]
    subprocess.run(build, check=True, cwd=HERE, capture_output=True, text=True)
    project_path = run_dir / "input" / "project.json"
    project = json.loads(project_path.read_text(encoding="utf-8"))
    frame = slope_frame(SLOPE_ANGLE, SLOPE_ASPECT)
    light = project["configuration"]["light"]
    sza, saa = light["zenith"], light["azimuth"]
    local_sza, local_saa = to_local(sza, saa, frame)
    light["zenith"], light["azimuth"] = round(local_sza, 6), round(local_saa, 6)
    rows, started = [], time.time()
    for index, (vza, vaa) in enumerate(principal_views(saa)):
        local_vza, local_vaa = to_local(vza, vaa, frame)
        if local_vza > 89.0:
            continue  # slope not visible: TiRT reports NaN
        out_dir = run_dir / "output" / f"view{index:02d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        configuration = project["configuration"]
        configuration["outDir"] = str(out_dir)
        configuration["sensor"].update({"principalPlane": False, "hemisphere": False, "vza": round(local_vza, 6), "vaa": round(local_vaa, 6)})
        project_path.write_text(json.dumps(project, ensure_ascii=False, indent=2), encoding="utf-8")
        code = run_engine(mode, run_dir / "input", run_dir / f"engine_view{index:02d}.log", args.timeout)
        reduced = reduce_images(out_dir, getattr(args, "keep_zero", False))
        if code != 0 or len(reduced) != 1:
            raise RuntimeError(f"{name} view {vza}/{vaa}: engine exit {code}, {len(reduced)} images")
        rows.append({**reduced[0], "vza": vza, "vaa": vaa, "sza": sza, "saa": saa, "local_vza": local_vza, "local_vaa": local_vaa})
    elapsed = time.time() - started
    meta = json.loads((run_dir / "case.json").read_text(encoding="utf-8"))
    meta.update({"slope": SLOPE_ANGLE, "aspect": SLOPE_ASPECT, "local_sza": local_sza, "local_saa": local_saa, "base_test": base_test,
                 "engine": {"path": str(ENGINE), "sha256": hashlib.sha256(ENGINE.read_bytes()).hexdigest()}})
    (run_dir / "case.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    with (run_dir / "diagnostics.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    expanded = write_case_csv(test, mode, rows, RUNS / f"{name}.csv")
    print(f"{name}: {len(rows)} views, {elapsed:.1f}s", flush=True)
    return {"name": name, "meta": meta, "elapsed_s": elapsed, "rows": expanded}


def run_case(test, mode, args, extra, tag=""):
    if test in SLOPE_TESTS:
        return run_slope_case(test, mode, args, extra, tag)
    name = f"{NAMES[test]}_{mode}{tag}"
    run_dir = RUNS / name

    def build_and_run(build_extra):
        if run_dir.exists():
            shutil.rmtree(run_dir)
        build = [NODE, str(HERE / "build_case.mjs"), "--test", test, "--mode", mode, "--out", str(run_dir),
                 "--samples", str(args.samples), "--seed", str(args.seed), "--res", str(args.res),
                 "--depth", str(args.depth), *build_extra]
        subprocess.run(build, check=True, cwd=HERE, capture_output=True, text=True)
        started = time.time()
        code = run_engine(mode, run_dir / "input", run_dir / "engine.log", args.timeout)
        return code, time.time() - started

    code, elapsed = build_and_run(extra)
    if test == "T3" and mode == "eVoxelRT" and "--lad" not in extra:
        # Overlapping crowns share voxels; rescale LAD so scene LAI stays 1.5.
        meta = json.loads((run_dir / "case.json").read_text(encoding="utf-8"))
        lad = 1.5 * meta["L"] ** 2 / (inserted_voxels(run_dir / "engine.log") * meta["voxel"] ** 3)
        code, elapsed = build_and_run([*extra, "--lad", f"{lad:.6f}"])
    rows = reduce_images(run_dir / "output", getattr(args, "keep_zero", False))
    if code != 0 or not rows:
        raise RuntimeError(f"{name}: engine exit {code}, {len(rows)} images; see {run_dir / 'engine.log'}")
    meta = json.loads((run_dir / "case.json").read_text(encoding="utf-8"))
    meta["engine"] = {"path": str(ENGINE), "sha256": hashlib.sha256(ENGINE.read_bytes()).hexdigest()}
    (run_dir / "case.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    with (run_dir / "diagnostics.csv").open("w", newline="", encoding="utf-8") as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    expanded = write_case_csv(test, mode, rows, RUNS / f"{name}.csv")
    print(f"{name}: exit {code}, {len(rows)} views, {elapsed:.1f}s", flush=True)
    return {"name": name, "meta": meta, "elapsed_s": elapsed, "rows": expanded}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=DEFAULT_CASES, help="TEST:MODE entries")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--res", type=int, default=512)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--tag", default="", help="suffix for run folder / CSV name")
    parser.add_argument("--extra", default="", help="extra build_case.mjs args, e.g. '--leaf 0.3'")
    parser.add_argument("--keep-zero", action="store_true", help="count zero-radiance pixels (needed with --cold/--black)")
    args = parser.parse_args()
    if not ENGINE.exists():
        sys.exit(f"HiSTREAM engine not found: {ENGINE}")
    RUNS.mkdir(parents=True, exist_ok=True)
    summary = []
    for entry in args.cases:
        test, mode = entry.split(":")
        result = run_case(test, mode, args, args.extra.split(), args.tag)
        summary.append({"name": result["name"], "elapsed_s": round(result["elapsed_s"], 2), **result["meta"]})
    target = RUNS / f"summary{args.tag}.json"
    previous = json.loads(target.read_text(encoding="utf-8")) if target.exists() else []
    names = {entry["name"] for entry in summary}
    target.write_text(json.dumps([entry for entry in previous if entry["name"] not in names] + summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
