"""Compare the TiRT DEM facet ensemble with the HiSTREAM multi-ridge reference."""
import csv
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(r"C:\work\tirt")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

import run  # noqa: E402
from base.physicsF import inv_planck  # noqa: E402
from base.solar import slope_direction  # noqa: E402

DEM = Path(r"C:\work\streamstudio\assets\dem-library\dune\transverse_dunes_100m.asc")


def read_grid(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = {}
    for line in lines[:6]:
        key, value = line.split()
        header[key.upper()] = float(value)
    ncols, nrows = int(header["NCOLS"]), int(header["NROWS"])
    z = np.array([[float(v) for v in line.split()] for line in lines[6:6 + nrows]])
    return z, header["CELLSIZE"]


def geometry(z, cell):
    dzdx = np.gradient(z, cell, axis=1)
    dzdy_north = -np.gradient(z, cell, axis=0)
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy_north)))
    aspect = np.degrees(np.arctan2(-dzdx, -dzdy_north)) % 360.0
    slope[0, :] = slope[-1, :] = slope[:, 0] = slope[:, -1] = 0.0
    return slope, aspect


def cast_shadow(z, cell, sza, saa):
    nrows, ncols = z.shape
    if sza >= 90.0:
        return np.ones_like(z, dtype=bool)
    sun_elevation = math.radians(90.0 - sza)
    dx, dy = math.sin(math.radians(saa)), math.cos(math.radians(saa))
    rows, cols = np.mgrid[0:nrows, 0:ncols]
    max_tan = np.zeros_like(z)
    for step in range(1, int(math.hypot(nrows, ncols) * 1.1) + 2):
        r, c = rows - step * dy, cols + step * dx
        valid = (r >= 0.0) & (r <= nrows - 1.001) & (c >= 0.0) & (c <= ncols - 1.001)
        if not np.any(valid):
            break
        r0 = np.clip(np.floor(r).astype(int), 0, nrows - 1)
        c0 = np.clip(np.floor(c).astype(int), 0, ncols - 1)
        fr, fc = np.clip(r - r0, 0.0, 1.0), np.clip(c - c0, 0.0, 1.0)
        r1, c1 = np.clip(r0 + 1, 0, nrows - 1), np.clip(c0 + 1, 0, ncols - 1)
        zs = (z[r0, c0] * (1 - fr) * (1 - fc) + z[r0, c1] * (1 - fr) * fc +
              z[r1, c0] * fr * (1 - fc) + z[r1, c1] * fr * fc)
        max_tan = np.maximum(max_tan, np.where(valid, (zs - z) / (step * cell), 0.0))
    return max_tan > math.tan(sun_elevation)


def bins(slope, aspect, shadow, step=2.0, aspect_step=10.0):
    sb = np.floor(slope / step) * step
    ab = np.floor(aspect / aspect_step) * aspect_step
    result = {}
    for r in range(1, slope.shape[0] - 1):
        for c in range(1, slope.shape[1] - 1):
            key = (float(sb[r, c]), float(ab[r, c]))
            entry = result.setdefault(key, [0, 0])
            entry[0 if not shadow[r, c] else 1] += 1
    return [(s, a, n_sun, n_shade) for (s, a), (n_sun, n_shade) in result.items()]


def config_and_spectra():
    cfg = run.read_input(ROOT / "input.csv")
    cfg.update(surface_model="slope", terrain_type="1", vegetation_model="bare", hapke_enable="0",
               sky_temperature_k="none", shortwave_diffuse_fraction="0",
               soil_sunlit_k="318.15", soil_shaded_k="303.15")
    table = run._band_table(cfg, ROOT / "input.csv")
    spectra = run._spectra(table)
    spectra = {key: values[2] for key, values in spectra.items()}
    spectra["soil_e"] = 0.95
    return cfg, spectra


def facet_bt(cfg, spectra, settings, sza, saa, vza, vaa, shaded):
    local = dict(cfg, slope_angle=str(settings[0]), slope_aspect=str(settings[1]))
    if shaded:
        local["soil_sunlit_k"] = "303.15"
    radiance, bt = run._surface_result(local, "slope", "bare", 10.5, (vza, vaa, sza, saa), spectra)
    return bt if np.isfinite(bt) else None


def ensemble(cfg, spectra, classes, sza, saa, vza, vaa):
    total = 0.0
    count = 0
    for slope, aspect, n_sun, n_shade in classes:
        local_vza, _, _ = slope_direction(vza, vaa, slope, aspect)
        if local_vza >= 89.5:
            continue
        # Radiance from a facet contributes in proportion to its projected area.
        projected = math.cos(math.radians(local_vza))
        for shaded, weight in ((False, n_sun * projected), (True, n_shade * projected)):
            if weight == 0:
                continue
            bt = facet_bt(cfg, spectra, (slope, aspect), sza, saa, vza, vaa, shaded)
            if bt is None:
                continue
            total += weight * bt
            count += weight
    return total / count if count else float("nan")


def main():
    z, cell = read_grid(DEM)
    slope, aspect = geometry(z, cell)
    cfg, spectra = config_and_spectra()
    reference = {}
    with (HERE / "histream_ridge.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            reference[(float(row["sza"]), float(row["saa"]), float(row["vza"]), float(row["vaa"]))] = float(row["bt_K"])

    rows = []
    for sza, saa in ((30.0, 270.0), (70.0, 270.0)):
        shadow = cast_shadow(z, cell, sza, saa)
        classes = bins(slope, aspect, shadow)
        for vaa in (90.0, 270.0):
            for vza in (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0):
                tirt = ensemble(cfg, spectra, classes, sza, saa, vza, vaa)
                hist = reference[(sza, saa, vza, vaa)]
                rows.append(dict(sza=sza, saa=saa, vza=vza, vaa=vaa, tirt_bt=tirt, histream_bt=hist,
                                 difference=tirt - hist))
                print(f"SZA {sza:g} VAA {vaa:g} VZA {vza:g}: TiRT {tirt:8.4f}  HiSTREAM {hist:8.4f}  d {tirt-hist:+7.4f}", flush=True)
    with (HERE / "ridge_compare.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for sza, saa in ((30.0, 270.0), (70.0, 270.0)):
        for vaa in (90.0, 270.0):
            selected = [r for r in rows if r["sza"] == sza and r["vaa"] == vaa]
            d = np.asarray([r["difference"] for r in selected])
            print(f"RMSE SZA {sza:g} VAA {vaa:g}: {np.sqrt(np.mean(d ** 2)):.4f} K  bias {d.mean():+.4f} K")


if __name__ == "__main__":
    main()
