"""Check that TIRTEB directional observations equal TiRT for the same scene and temperatures.

TIRTEB computes the energy balance but its directional shortwave and thermal observations use the
same standard closure as ``run.py``. This script drives both interfaces with identical component
temperatures, geometry and band optics and requires them to agree to ``--tolerance``.

Every scene is written with ``input_tirteb.csv`` keys and goes through ``read_input`` /
``read_defined``. The TiRT reference takes its geometry (crown, slope and terrain keys) directly
from the scene keys rather than from the TIRTEB settings, so a wrong TIRTEB parameter mapping
shows up as a difference.

Usage:
    python validation/compare_tirt_tirteb.py
    python validation/compare_tirt_tirteb.py --full

``--full`` adds more sun and view angles; all suns are above the horizon (no sub-horizon case).
A direction behind a single slope is expected to be empty in TIRTEB and rejected by TiRT; both
are treated as equal.
"""
import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run as tirt  # noqa: E402
import tirteb_run  # noqa: E402
from base import tirteb_tirt  # noqa: E402
from base.io import read_defined  # noqa: E402
from base.modelio import read_input  # noqa: E402
from base.solar import slope_direction  # noqa: E402

SCENES = {
    "plane_bare": dict(surface_model="plane", vegetation_model="bare", lai=3.0),
    "plane_hom": dict(surface_model="plane", vegetation_model="hom", lai=3.0),
    "plane_row": dict(surface_model="plane", vegetation_model="row", lai=3.0),
    "plane_crown": dict(surface_model="plane", vegetation_model="crown", lai=3.0,
                        crown_rad_a=1.5, crown_rad_b=3.0, crown_density=0.1, crown_height=1.0),
    "slope_bare": dict(surface_model="terrain", terrain_type=1, terrain_slope=20.0,
                       terrain_aspect=135.0, vegetation_model="bare", lai=3.0),
    "slope_hom": dict(surface_model="terrain", terrain_type=1, slope_angle=20.0,
                      slope_aspect=135.0, vegetation_model="hom", lai=3.0),
    "slope_row": dict(surface_model="terrain", terrain_type=1, terrain_slope=20.0,
                      terrain_aspect=135.0, vegetation_model="row", lai=3.0),
    "slope_crown": dict(surface_model="terrain", terrain_type=1, slope_angle=20.0,
                        slope_aspect=135.0, vegetation_model="crown", lai=3.0,
                        crown_rad_a=1.5, crown_rad_b=3.0, crown_density=0.1, crown_height=1.0),
    "terrain_bare": dict(surface_model="terrain", terrain_type=2, terrain_radius=10.0,
                         terrain_height=20.0, terrain_density=0.001, vegetation_model="bare", lai=3.0),
    "terrain_hom": dict(surface_model="terrain", terrain_type=2, terrain_radius=10.0,
                        terrain_height=20.0, terrain_density=0.001, vegetation_model="hom", lai=3.0),
    "terrain_row": dict(surface_model="terrain", terrain_type=2, terrain_radius=10.0,
                        terrain_height=20.0, terrain_density=0.001, vegetation_model="row", lai=3.0),
    "terrain_crown": dict(surface_model="terrain", terrain_type=2, terrain_radius=10.0,
                          terrain_height=20.0, terrain_density=0.001, vegetation_model="crown", lai=3.0,
                          crown_rad_a=1.5, crown_rad_b=3.0, crown_density=0.1, crown_height=1.0),
    "urban_bare": dict(surface_model="urban", vegetation_model="bare", lai=3.0),
    "urban_hom": dict(surface_model="urban", vegetation_model="hom", lai=3.0),
    "urban_row": dict(surface_model="urban", vegetation_model="row", lai=3.0),
    "urban_crown": dict(surface_model="urban", vegetation_model="crown", lai=3.0,
                        crown_rad_a=1.5, crown_rad_b=3.0, crown_density=0.1, crown_height=1.0),
}

STATE = {
    "Tss": 318.15, "Tsh": 303.15, "Tls": 306.15, "Tlh": 303.15,
    "Troof_sunlit": 318.15, "Troof_shaded": 303.15,
    "Twall_sunlit": 318.15, "Twall_shaded": 303.15,
    "Tstreet_sunlit": 318.15, "Tstreet_shaded": 303.15,
    "sunlit_surface": True, "shortwave_diffuse_fraction": 0.2,
    "sza": 40.0, "saa": 150.0,
}
METEO = {"shortwave": 700.0, "longwave": 350.0, "air_temperature_c": 25.0, "pressure": 101325.0}
SHORTWAVE = (0.55, 0.86)
THERMAL = 10.5
TIRT_GEOMETRY_KEYS = ("crown_rad_a", "crown_rad_b", "crown_height", "crown_density",
                      "terrain_height", "terrain_radius", "terrain_density")


def _slope(scene):
    return (float(scene.get("slope_angle", scene.get("terrain_slope", 0.0))),
            float(scene.get("slope_aspect", scene.get("terrain_aspect", 0.0))))


def _settings(config, scene):
    config = dict(config)
    if "slope_angle" in scene:
        config.pop("terrain_slope", None)
        config.pop("terrain_aspect", None)
    return read_defined(None, dict(config, **{k: str(v) for k, v in scene.items()}))


def _reference(config, scene):
    """TiRT config with the geometry taken from the scene's own keys."""
    config = dict(config, **{key: float(scene[key]) for key in TIRT_GEOMETRY_KEYS if key in scene})
    if int(scene.get("terrain_type", 0)) == 1:
        config["slope_angle"], config["slope_aspect"] = _slope(scene)
    return config


def _state_for(overrides, sza, saa):
    state = dict(STATE, sza=sza, saa=saa, sza_global=sza, saa_global=saa)
    if int(overrides.get("terrain_type", 0)) == 1:
        state["sza"], state["saa"], state["sunlit_surface"] = slope_direction(sza, saa, *_slope(overrides))
    return state


def _ti_band(settings, state):
    band = tirteb_tirt.tir_band_settings(settings, THERMAL)
    surface, vegetation, config = tirteb_tirt.tirt_config(band, state)
    return band, surface, vegetation, config, tirteb_tirt.band_emissivities(band)


def _difference(a, b):
    if math.isfinite(a) and math.isfinite(b):
        return abs(a - b)
    if math.isfinite(a) != math.isfinite(b):
        return math.inf
    return 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=ROOT / "input_tirteb.csv")
    parser.add_argument("--tolerance", type=float, default=1.0e-9)
    parser.add_argument("--full", action="store_true", help="also test sub-horizon and extra angles")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    config = read_input(input_path)
    config["_input_base"] = str(input_path.parent)
    suns = ((20.0, 120.0), (50.0, 240.0)) if args.full else ((40.0, 150.0),)
    views = ((0.0, 0.0), (30.0, 90.0), (60.0, 180.0), (78.0, 270.0)) if args.full else ((0.0, 0.0), (30.0, 60.0))
    failures = 0
    worst = {"thermal": (0.0, None), "shortwave": (0.0, None)}
    for name, overrides in SCENES.items():
        settings = _settings(config, overrides)
        scene_worst = {"thermal": 0.0, "shortwave": 0.0}
        for sza, saa in suns:
            state = _state_for(overrides, sza, saa)
            band, surface, vegetation, _, spectra = _ti_band(settings, state)
            reference_config = _reference(_ti_band(settings, state)[3], overrides)
            for vza, vaa in views:
                lib, bt = tirteb_tirt.tirt_brightness(band, state, THERMAL, vza, vaa)
                lib2, bt2 = tirt._surface_result(dict(reference_config), surface, vegetation,
                                                 THERMAL, (vza, vaa, sza, saa), spectra)
                delta = max(_difference(lib, lib2), _difference(bt, bt2))
                scene_worst["thermal"] = max(scene_worst["thermal"], delta)
                if delta > worst["thermal"][0]:
                    worst["thermal"] = (delta, (name, sza, saa, vza, vaa))

                for wavelength in SHORTWAVE:
                    observed = tirteb_run._band_observation(
                        settings, METEO, state, tirteb_run.model_name(settings), wavelength, vza, vaa)
                    rho, tau, soil = tirteb_run._spectral_values(
                        settings, wavelength, tirteb_run.model_name(settings))
                    optics = tirteb_tirt.band_optics(settings, wavelength)
                    band_optics = dict(optics, leaf_reflectance=rho,
                                       leaf_transmittance=tau, soil_reflectance=soil)
                    scfg = dict(reference_config)
                    try:
                        reference = tirt._optical_band(scfg, surface, vegetation, wavelength,
                                                       (vza, vaa, sza, saa), band_optics)
                    except ValueError:
                        reference = float("nan")
                    delta = _difference(float(observed["reflectance"]), float(reference))
                    scene_worst["shortwave"] = max(scene_worst["shortwave"], delta)
                    if delta > worst["shortwave"][0]:
                        worst["shortwave"] = (delta, (name, wavelength, sza, saa, vza, vaa))
        ok = all(value <= args.tolerance for value in scene_worst.values())
        failures += not ok
        print(f"{name:15s} thermal {scene_worst['thermal']:.2e}  shortwave {scene_worst['shortwave']:.2e}"
              f"  {'ok' if ok else 'FAIL'}")

    print(f"worst thermal   {worst['thermal'][0]:.3e} at {worst['thermal'][1]}")
    print(f"worst shortwave {worst['shortwave'][0]:.3e} at {worst['shortwave'][1]}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
