import argparse
import csv
import json
import shutil
import warnings
from pathlib import Path

import numpy as np

from base import tirteb_radiation
from base.brightness import urban_canopy_fraction, urban_component_fractions, urban_directional_reflectance
from base.tirteb_tirt import band_optics, tir_band_settings, tirt_brightness, tirt_emissivity
from base.modelio import read_defined, read_directions, read_input, read_meteo
from base.model import run_node
from base.solar import slope_direction
from base.structs import InputStructure, ObservationStructure, RTStructure
from base.radiation import clearness_index_warning, directional_reflectance


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input_tirteb.csv"
DEFAULT_OBSERVATION_VZA = (0.0, 30.0, 60.0)
DEFAULT_OBSERVATION_VAA = (0.0, 90.0, 180.0)
OUTPUT_FIELDS = (
    "node", "time", "veg", "surface", "Ta_C",
    "Tss_C", "Tsh_C", "Tls_C", "Tlh_C", "Rn", "LE", "H", "G",
    "Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
    "Tstreet_sunlit_C", "Tstreet_shaded_C", "Tleaf_sunlit_C", "Tleaf_shaded_C",
    "frac_soil_sunlit", "frac_leaf_sunlit", "frac_roof_sunlit", "frac_wall_sunlit", "frac_street_sunlit",
    "Tv_C", "green", "nir",
    "closed", "iterations", "energy_balance_residual_W_m2",
    "L", "ustar", "stability_closed", "stability_iterations", "stability_residual_m_inv",
    "stability_limited", "stability_flux_area_factor", "fallback_solver_used", "temperature_solver",
    "clearness_index", "clearness_index_exceeds_one", "solar_direct_normal_limit_W_m2",
    "shortwave_direct_normal_W_m2", "air_density_kg_m3",
    "aerodynamic_roughness_length", "aerodynamic_displacement", "aerodynamic_canopy_weight",
)
BOOLEAN_DIAGNOSTICS = ("stability_limited", "fallback_solver_used", "clearness_index_exceeds_one")
NUMERIC_DIAGNOSTICS = ("stability_flux_area_factor", "clearness_index", "solar_direct_normal_limit_W_m2",
                       "shortwave_direct_normal_W_m2", "air_density_kg_m3")
AERODYNAMIC_DIAGNOSTICS = ("aerodynamic_roughness_length", "aerodynamic_displacement", "aerodynamic_canopy_weight")
OBSERVATION_BASE_FIELDS = (
    "node", "time", "veg", "surface", "vza", "vaa", "sza", "saa",
)
INPUT_FILE_KEYS = {
    "meteo", "defined", "defined_path", "defined_file",
    "esun_file", "solar_direct_file", "esky_file", "solar_sky_file",
    "optical_coefficients_file", "optipar_file", "directions_file", "spectrals_file",
    "vegetation_spectrum_file", "vegetation_spectral_file", "leaf_spectrum_file",
    "soil_spectrum_file", "soil_spectral_file", "urban_spectrum_file", "urban_spectral_file",
    "urban_roof_spectrum_file", "roof_spectrum_file",
    "urban_wall_spectrum_file", "wall_spectrum_file",
    "urban_street_spectrum_file", "street_spectrum_file",
}


def model_name(settings):
    surface = settings.get("surface_model")
    if surface == "urban" or settings["surface_type"] == 13 or settings["distribution"] == 4:
        return "urban"
    vegetation = settings.get("vegetation_model")
    if vegetation in {"bare", "bare_soil", "soil", "barren"}:
        return "bare"
    if vegetation in {"hom", "homogeneous"}:
        return "homogeneous"
    if vegetation in {"row", "raw"}:
        return "row"
    if vegetation in {"crown", "forest"}:
        return "crown"
    if settings["distribution"] == 2 or settings["surface_type"] in (12, 14):
        return "row"
    if settings["surface_type"] in (1, 2, 3, 4, 5):
        return "crown"
    return "homogeneous"


def _config_path(config, base, key, default):
    value = config.get(key, default)
    path = Path(value)
    return path if path.is_absolute() else base / path


def _spectral_values(settings, wavelength, model=""):
    optics = band_optics(settings, wavelength)
    # Urban streets replace the soil background in the shortwave.
    background = "street_reflectance" if model == "urban" and float(wavelength) < 2.5 else "soil_reflectance"
    return float(optics["leaf_reflectance"]), float(optics["leaf_transmittance"]), float(optics[background])


_tir_band_settings = tir_band_settings


def _band_observation(settings, meteo, state, name, wavelength, vza, vaa, tir_reflectance=False):
    # A solar reflectance is undefined without incident daylight. In particular,
    # never pass below-horizon solar angles into the crown hotspot geometry.
    # Thermal emission and the energy-balance time stepping remain active.
    if float(wavelength) < 2.5 and (
        float(meteo["shortwave"]) <= 0.0
        or float(state.get("sza_global", state["sza"])) >= 90.0
    ):
        return RTStructure.from_mapping({
            "wavelength_um": wavelength,
            "reflectance": float("nan"),
            "transmittance": float("nan"),
            "reflected_shortwave": 0.0 if float(meteo["shortwave"]) <= 0.0 else float("nan"),
            "radiance": float("nan"),
            "brightness_temperature_K": float("nan"),
        })
    if settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 2)) == 1:
        local_vza, local_vaa, visible = slope_direction(
            vza,
            vaa,
            settings.get("terrain_slope", 0.0),
            settings.get("terrain_aspect", 0.0),
        )
        # A direction behind the slope (or exactly grazing it) sees no surface.
        # TiRT leaves these rows empty; keep the same convention here.
        if not visible or local_vza >= 90.0:
            return RTStructure.from_mapping({
                "wavelength_um": wavelength,
                "reflectance": float("nan"),
                "transmittance": float("nan"),
                "reflected_shortwave": float("nan"),
                "radiance": float("nan"),
                "brightness_temperature_K": float("nan"),
            })
    else:
        local_vza, local_vaa = float(vza), float(vaa)
    sunlit = bool(state.get("sunlit_surface", True))
    leaf_rho, leaf_tau, soil_rho = _spectral_values(settings, wavelength, name)
    radiance = brightness_temperature = reflected = float("nan")
    if float(wavelength) >= 2.5:
        band_settings = _tir_band_settings(settings, wavelength)
        # Thermal observations use the TiRT directional model at the energy-balance temperatures;
        # the thermal reflectance is 1 - directional emissivity (Kirchhoff).
        radiance, brightness_temperature = tirt_brightness(band_settings, state, wavelength, vza, vaa)
        reflectance = (1.0 - tirt_emissivity(band_settings, state, wavelength, vza, vaa)
                       if tir_reflectance else float("nan"))
    else:
        # The state holds the meteorology-derived diffuse fraction; settings
        # may still contain -1 (automatic partition), which is not a fraction.
        from base.radiation import _shortwave_diffuse_fraction
        diffuse = state.get("shortwave_diffuse_fraction")
        if diffuse is None:
            diffuse = _shortwave_diffuse_fraction(settings, meteo["shortwave"], state["sza"])
        if name == "urban":
            optics = band_optics(settings, wavelength)
            urban_settings = dict(
                settings, shortwave_diffuse_fraction=diffuse,
                urban_reflectance_roof=optics["roof_reflectance"], urban_reflectance_wall=optics["wall_reflectance"],
                urban_reflectance_street=optics["street_reflectance"],
                leaf_reflectance=leaf_rho, leaf_transmittance=leaf_tau)
            for key in ("urban_reflectance_roof_spectrum", "urban_reflectance_wall_spectrum",
                        "urban_reflectance_street_spectrum", "leaf_reflectance_spectrum", "leaf_transmittance_spectrum"):
                urban_settings.pop(key, None)
            reflectance = urban_directional_reflectance(
                urban_settings, wavelength, local_vza, local_vaa, state["sza"], state["saa"])
        elif settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 2)) == 2:
            reflectance = directional_reflectance(
                settings, state["sza"], state["saa"], local_vza, local_vaa, leaf_rho, leaf_tau, soil_rho,
                sunlit=sunlit, use_hapke=False, diffuse_fraction=diffuse)
        else:
            soil_direct = None
            if int(settings.get("hapke_enable", settings.get("hapke_optical_enable", 0))) == 1:
                from base.soilopt import hapke_reflectance
                soil_direct = hapke_reflectance(local_vza, settings.get("hapke_K", 1.0), soil_rho)
            if settings.get("surface_model") == "terrain":
                reflectance = tirteb_radiation.slope_reflectance(
                    settings, state.get("sza_global", state["sza"]), state.get("saa_global", state["saa"]),
                    state["sza"], state["saa"], local_vza, local_vaa, leaf_rho, leaf_tau, soil_rho,
                    diffuse, sunlit, soil_direct)
            else:
                reflectance = tirteb_radiation.reflectance(
                    settings, state["sza"], state["saa"], local_vza, local_vaa, leaf_rho, leaf_tau, soil_rho,
                    diffuse, sunlit, soil_direct)
        reflected = max(float(meteo["shortwave"]), 0.0) * reflectance
    fractions = tirteb_radiation.view_fractions(settings, local_vza, local_vaa, state["sza"], state["saa"], sunlit)
    leaf_fraction = float(fractions[2] + fractions[3])
    if name == "urban":
        urban_fractions = urban_component_fractions(settings, local_vza, local_vaa, state["sza"], state["saa"])
        street_fraction = float(urban_fractions[2] + urban_fractions[3])
        transmittance = street_fraction * urban_canopy_fraction(settings, local_vza, local_vaa) * leaf_fraction * leaf_tau
    else:
        transmittance = leaf_fraction * leaf_tau
    return RTStructure.from_mapping({
        "wavelength_um": wavelength,
        "reflectance": reflectance,
        "transmittance": transmittance,
        "reflected_shortwave": reflected,
        "radiance": radiance,
        "brightness_temperature_K": brightness_temperature,
    })


def _band_name(wavelength):
    if wavelength < 0.7:
        return "optical"
    if wavelength < 2.5:
        return "infrared"
    return "thermal_infrared"


def _observation_band_name(wavelength):
    if abs(float(wavelength) - 0.55) < 1e-6:
        return "green"
    if abs(float(wavelength) - 0.86) < 1e-6:
        return "nir"
    return _band_name(wavelength)


def _wavelength_tag(wavelength):
    return f"{float(wavelength):g}".replace(".", "p")


def _urban_surface_totals(settings, state):
    from base.scene_fluxes import scene_totals
    return scene_totals(settings, state, urban=True)


def _flux_weights(settings, state):
    """Scene weights of the soil (per soil area) and big-leaf (per ground) sunlit and shaded fluxes."""
    lai = float(settings.get("lai", 0.0))
    leaf = (state["lai_sunlit"] / lai, state["lai_shaded"] / lai) if lai > 0.0 else (0.0, 0.0)
    return state["soil_sunlit_fraction"], state["soil_shaded_fraction"], *leaf


def _diagnostics(state, model):
    nan = float("nan")
    row = {key: int(bool(state.get(key, False))) for key in BOOLEAN_DIAGNOSTICS}
    row.update({key: state.get(key, nan) for key in NUMERIC_DIAGNOSTICS})
    row.update({key: nan if model == "urban" else state.get(key, nan) for key in AERODYNAMIC_DIAGNOSTICS})
    row["temperature_solver"] = state.get("temperature_solver", "relaxation")
    return row


def _sunlit_shares(settings, state, model):
    """Sunlit share of each surface's energy-balance area weight (shaded share = 1 - value); NaN if absent."""
    def share(sunlit, shaded):
        total = float(sunlit) + float(shaded)
        return float(sunlit) / total if total > 0.0 else float("nan")

    nan = float("nan")
    if model == "urban":
        from base.scene_fluxes import component_weights
        w = component_weights(settings, state, urban=True)
        street = share(w["streetsunlit"], w["streetshaded"])
        shares = {"frac_soil_sunlit": street, "frac_leaf_sunlit": share(w["leafsunlit"], w["leafshaded"]),
                  "frac_roof_sunlit": share(w["roofsunlit"], w["roofshaded"]),
                  "frac_wall_sunlit": share(w["wallsunlit"], w["wallshaded"]), "frac_street_sunlit": street}
    else:
        soil_sun, soil_shade, leaf_sun, leaf_shade = _flux_weights(settings, state)
        shares = {"frac_soil_sunlit": share(soil_sun, soil_shade), "frac_leaf_sunlit": share(leaf_sun, leaf_shade),
                  "frac_roof_sunlit": nan, "frac_wall_sunlit": nan, "frac_street_sunlit": nan}
    if float(settings.get("lai", 0.0)) <= 0.0:
        shares["frac_leaf_sunlit"] = nan
    return shares


def _observation_quantity(config, wavelength):
    if float(wavelength) < 2.5:
        return "reflectance"
    quantity = str(
        config.get("observation_quantity", config.get("longwave_quantity", "temperature"))
    ).strip().lower()
    if quantity in {"reflectance", "rho", "refl"}:
        return "reflectance"
    if quantity in {"temperature", "temp", "brightness_temperature", "bt"}:
        return "temperature"
    raise ValueError("observation_quantity must be temperature or reflectance")


def _observation_columns(wavelength, values, quantity):
    prefix = _wavelength_tag(wavelength)
    if quantity == "temperature":
        return {prefix: values.get("brightness_temperature_K", float("nan"))}
    return {prefix: values.get("reflectance", float("nan"))}


# Half a minute in days: meteorological DOY values carry only ~4 decimals (about 9 s).
OBSERVATION_TIME_TOLERANCE_DAYS = 0.5 / 1440.0


def _observation_hour_fractions(config):
    """Parse observation_hours as local clock hours 0-24 (decimals allowed) and return day fractions."""
    raw = config.get("observation_hours", "")
    fractions = []
    for item in str(raw).replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        try:
            hour = float(item)
        except ValueError:
            hour = float("nan")
        if not 0.0 <= hour <= 24.0:
            raise ValueError(f"observation_hours must be clock hours between 0 and 24 (e.g. 12 or 13.5); got {item!r}")
        fractions.append((hour / 24.0) % 1.0)
    return tuple(fractions)


def _matches_observation_hour(time_value, fractions, tolerance=OBSERVATION_TIME_TOLERANCE_DAYS):
    # Meteorological DOY values carry only ~4 decimals (about 9 s), so match to the half minute.
    fraction = float(time_value) % 1.0
    return any(min(abs(fraction - target), 1.0 - abs(fraction - target)) <= tolerance for target in fractions)


def _matches_observation_time(time_value, times, tolerance=OBSERVATION_TIME_TOLERANCE_DAYS):
    """Absolute DOY times match within the same half-minute tolerance as observation_hours."""
    value = float(time_value)
    return any(abs(value - float(target)) <= tolerance for target in times)


def _csv_value(value):
    if isinstance(value, np.ndarray):
        return json.dumps(value.tolist(), separators=(",", ":"))
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, separators=(",", ":"), default=str)
    return value


def _intermediate_row(state, node, time_value, model, totals):
    row = {"node": node, "time": time_value, "model": model}
    row.update({key: _csv_value(value) for key, value in state.as_dict().items()})
    row.update({f"{key}_total": value for key, value in totals.items()})
    return row


SNAPSHOT_MARKER = "# Effective TIRTEB input"


def _is_input_snapshot(path):
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    return text.startswith(SNAPSHOT_MARKER) or any(line.strip().startswith("[resolved_inputs]") for line in text.splitlines())


def _snapshot_path(input_path, output_dir):
    """input.csv beside the results, unless that would overwrite the source or a user file."""
    snapshot = output_dir / "input.csv"
    if snapshot.resolve() == input_path.resolve() or (snapshot.is_file() and not _is_input_snapshot(snapshot)):
        return output_dir / "input_used.csv"
    return snapshot


def _write_input_snapshot(input_path, snapshot, config, overrides, input_data):
    """Keep the source, effective overrides and local copies of referenced files."""
    original = input_path.read_bytes()
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    source_snapshot = snapshot.with_name("source_input.csv")
    default_snapshot = snapshot.with_name("input.csv")
    if source_snapshot.resolve() != input_path.resolve() and (
        input_path.resolve() != default_snapshot.resolve() or not source_snapshot.exists()
    ):
        source_snapshot.write_bytes(original)
    references = {key: config[key] for key in INPUT_FILE_KEYS if config.get(key) not in (None, "")}
    references.setdefault("meteo", str(input_data.meteo_path))
    for key, alias, default in (
        ("esun_file", "solar_direct_file", "data/defined/Esun_.dat"),
        ("esky_file", "solar_sky_file", "data/defined/Esky_.dat"),
        ("optical_coefficients_file", "optipar_file", "data/defined/optipar_new.csv"),
    ):
        if key not in references and alias not in references:
            references[key] = default
    if input_data.observation_direction_mode in {"file", "txt", "read"}:
        references.setdefault("directions_file", str(input_data.directions_file))
    if input_data.observation_spectral_mode in {"file", "txt", "read"}:
        references.setdefault("spectrals_file", str(input_data.spectrals_file))

    frozen = {}
    resolved_references = {}
    for key, value in sorted(references.items()):
        source = Path(str(value))
        if not source.is_absolute():
            source = input_path.parent / source
        source = source.resolve()
        if source.is_file():
            if source not in frozen:
                filename = source.name if source.name.startswith(key + "__") else f"{key}__{source.name}"
                target = snapshot.parent / "inputs" / filename
                target.parent.mkdir(parents=True, exist_ok=True)
                if source != target.resolve():
                    shutil.copy2(source, target)
                frozen[source] = target.resolve()
            resolved_references[key] = frozen[source].relative_to(snapshot.parent.resolve()).as_posix()
        else:
            # Inactive optional sources can remain absent, but never relative.
            resolved_references[key] = str(source)

    text = original.decode("utf-8-sig")
    if text.startswith(SNAPSHOT_MARKER):
        text = text.split("\n", 1)[1] if "\n" in text else ""
    with snapshot.open("w", newline="", encoding="utf-8") as file:
        file.write(f"{SNAPSHOT_MARKER} written by tirteb_run.py; overwritten on rerun,,,\n")
        file.write(text)
        if text and not text.endswith("\n"):
            file.write("\n")
        writer = csv.writer(file)
        if overrides:
            file.write("\n[effective_overrides]\n")
            for key, value in overrides.items():
                if value is not None:
                    writer.writerow(["effective_overrides", key, value, "applied for this case"])
        file.write("\n[resolved_inputs]\n")
        writer.writerow(["resolved_inputs", "case_result_dir", ".", "write replay results beside this snapshot"])
        for key, value in resolved_references.items():
            writer.writerow(["resolved_inputs", key, value, "resolved input reference"])


def _run_csv(input_path=DEFAULT_INPUT, output_path=None, **overrides):
    input_path = Path(input_path).resolve()
    config = read_input(input_path)
    base = input_path.parent
    config["_input_base"] = str(base)
    config.update({key: str(value) for key, value in overrides.items() if value is not None})
    input_data = InputStructure.from_config(config, base)

    meteo_path = input_data.meteo_path
    if output_path is None:
        output_path = input_data.case_dir / "output.csv"
    else:
        output_path = Path(output_path)
        if not output_path.is_absolute():
            output_path = base / output_path

    year, doy = input_data.year, input_data.doy
    latitude, longitude = input_data.latitude, input_data.longitude
    dem, max_iter = input_data.dem, input_data.max_iter
    observation_vza = input_data.observation_vza or DEFAULT_OBSERVATION_VZA
    observation_vaa = input_data.observation_vaa or DEFAULT_OBSERVATION_VAA
    observation_times = input_data.observation_times
    if input_data.observation_direction_mode in {"file", "txt", "read"}:
        observation_directions = read_directions(input_data.directions_file)
    else:
        observation_directions = tuple(
            (view_zenith, view_azimuth)
            for view_zenith in observation_vza
            for view_azimuth in observation_vaa
        )
    observation_time_mode = input_data.observation_time_mode
    if observation_time_mode in {"all", "every", "full"}:
        selected_times = None
        selected_fractions = None
    elif observation_time_mode in {"time", "selected", "specific", "some", "exact"}:
        selected_times = tuple(observation_times)
        if not selected_times:
            raise ValueError("observation_time_mode=time requires observation_times")
        selected_fractions = None
    elif observation_time_mode in {"hour", "hours", "fraction", "fractional"}:
        selected_times = None
        selected_fractions = _observation_hour_fractions(config)
        if not selected_fractions:
            raise ValueError("observation_time_mode=hour requires observation_hours (clock hours 0-24, e.g. 12;13.5)")
    else:
        raise ValueError("observation_time_mode must be all, time, or hour")

    all_records = read_meteo(meteo_path)
    records = [
        record for record in all_records
        if input_data.start_time <= float(record["time"]) <= input_data.end_time
    ]
    if not records:
        raise ValueError(
            f"start_time={input_data.start_time:g} and end_time={input_data.end_time:g} select no meteorological data"
        )
    observation_selected = [
        (selected_times is None and selected_fractions is None)
        or (selected_times is not None and _matches_observation_time(record["time"], selected_times))
        or (selected_fractions is not None and _matches_observation_hour(record["time"], selected_fractions))
        for record in records
    ]
    if not any(observation_selected):
        selector = config.get("observation_hours", "") if selected_fractions is not None else config.get("observation_times", "")
        available_times = [float(record["time"]) for record in records]
        if selected_fractions is not None:
            nearest = sorted({
                min(available_times, key=lambda value: min(abs(value % 1.0 - target), 1.0 - abs(value % 1.0 - target))) % 1.0
                for target in selected_fractions
            })
            raise ValueError(
                f"observation_hours {selector!r} match no meteorological record; nearest record clock hours: "
                + ";".join(f"{value * 24.0:.4g}" for value in nearest)
            )
        suggestions = sorted({
            min(available_times, key=lambda value: abs(value - target))
            for target in selected_times
        })[:5] if selected_times is not None else available_times[:5]
        available = ";".join(f"{value:.12g}" for value in suggestions)
        raise ValueError(
            f"Observation time {selector!r} matches no meteorological record in the simulation interval "
            f"({float(records[0]['time']):.12g}–{float(records[-1]['time']):.12g}); "
            f"available times: {available}. Specify an existing meteorological time or use observation_time_mode=all"
        )
    settings = read_defined(None, config)
    observation_wavelengths = tuple(settings["observation_wavelengths"])
    if model_name(settings) == "urban":
        scene_height = max(
            float(settings.get("urban_building_height", 0.0)),
            float(settings.get("aerodynamic_height", settings.get("height", 0.0))),
        )
        # Keep the user's reference height; only lift it above the roughness elements
        # (aerodynamic_resistance applies the same lower bound).
        settings["z"] = max(float(settings.get("z", 15.0)), scene_height + 0.1)
    # Initialize soil and urban surfaces from the full-day air-temperature
    # mean of the first simulated day. Legacy input temperatures remain
    # readable for compatibility but are not used for initialization.
    start_day = int(float(records[0]["time"]))
    daily_air = [
        float(record["air_temperature_c"])
        for record in all_records
        if int(float(record["time"])) == start_day
    ]
    if not daily_air:
        daily_air = [float(record["air_temperature_c"]) for record in records]
    initial_surface_temperature = sum(daily_air) / len(daily_air)
    settings["soil_temperature_c"] = initial_surface_temperature
    settings["urban_temperature_c"] = initial_surface_temperature
    name = model_name(settings)
    previous = None
    output = []
    intermediate = []
    final_states = []
    vertical_vza, vertical_vaa = 0.0, 0.0

    for node, meteo in enumerate(records):
        if settings.get("surface_model") == "terrain":
            if int(settings.get("terrain_type", 2)) == 1:
                terrain_slope = float(settings.get("terrain_slope", 0.0))
                terrain_aspect = float(settings.get("terrain_aspect", 0.0))
            else:
                radius = max(float(settings.get("terrain_radius", 100.0)), 1e-6)
                height = float(settings.get("terrain_height", 20.0))
                terrain_slope = float(np.degrees(np.arctan2(height, radius)))
                terrain_aspect = float(settings.get("terrain_aspect", 0.0))
        else:
            terrain_slope, terrain_aspect = 0.0, 0.0
        state = run_node(
            settings, meteo, node, year, doy, latitude, longitude,
            dem, terrain_slope, vertical_vza, vertical_vaa, previous, max_iter, terrain_aspect,
            input_data.utc_offset_hours,
        )
        vertical_rt = _band_observation(settings, meteo, state, name, 10.5, vertical_vza, vertical_vaa)
        green_rt = _band_observation(settings, meteo, state, name, 0.55, vertical_vza, vertical_vaa)
        nir_rt = _band_observation(settings, meteo, state, name, 0.86, vertical_vza, vertical_vaa)
        if name == "urban":
            urban_totals = _urban_surface_totals(settings, state)
            rn_total, le_total, h_total, g_soil = (urban_totals[key] for key in ("Rn", "LE", "H", "G"))
        else:
            # Soil fluxes are per unit soil area, leaf fluxes per unit ground of the big leaf.
            weights = _flux_weights(settings, state)
            totals_by_prefix = {
                prefix: sum(weight * state.get(prefix + key, 0.0) for weight, key in zip(
                    weights, ("soilsunlit", "soilshaded", "leafsunlit", "leafshaded")))
                for prefix in ("N_", "LE_", "H_")
            }
            rn_total, le_total, h_total = totals_by_prefix["N_"], totals_by_prefix["LE_"], totals_by_prefix["H_"]
            g_soil = weights[0] * state.get("G_sunlit", 0.0) + weights[1] * state.get("G_shaded", 0.0)
        totals = {"Rn": rn_total, "LE": le_total, "H": h_total, "G": g_soil}
        # No leaves exist without LAI (bare soil, urban bare): leave their temperatures empty.
        has_leaves = float(settings.get("lai", 0.0)) > 0.0
        leaf_temperature_sunlit = state["Tls"] - 273.15 if has_leaves else float("nan")
        leaf_temperature_shaded = state["Tlh"] - 273.15 if has_leaves else float("nan")
        component_temperature_keys = (
            "Troof_sunlit", "Troof_shaded", "Twall_sunlit", "Twall_shaded",
            "Tstreet_sunlit", "Tstreet_shaded",
        )
        urban_temperatures = (
            {key: state[key] - 273.15 for key in component_temperature_keys}
            if name == "urban" else
            {key: float("nan") for key in component_temperature_keys}
        )
        output.append({
            "node": node, "time": meteo["time"], "veg": name,
            "surface": settings.get("surface_model", "legacy"),
            "Ta_C": meteo["air_temperature_c"],
            "Tss_C": state["Tss"] - 273.15, "Tsh_C": state["Tsh"] - 273.15,
            "Tls_C": leaf_temperature_sunlit, "Tlh_C": leaf_temperature_shaded,
            "Troof_sunlit_C": urban_temperatures["Troof_sunlit"],
            "Troof_shaded_C": urban_temperatures["Troof_shaded"],
            "Twall_sunlit_C": urban_temperatures["Twall_sunlit"],
            "Twall_shaded_C": urban_temperatures["Twall_shaded"],
            "Tstreet_sunlit_C": urban_temperatures["Tstreet_sunlit"],
            "Tstreet_shaded_C": urban_temperatures["Tstreet_shaded"],
            "Tleaf_sunlit_C": leaf_temperature_sunlit if name == "urban" else float("nan"),
            "Tleaf_shaded_C": leaf_temperature_shaded if name == "urban" else float("nan"),
            **_sunlit_shares(settings, state, name),
            "Rn": rn_total, "LE": le_total, "H": h_total, "G": g_soil,
            "Tv_C": vertical_rt["brightness_temperature_K"] - 273.15,
            "green": green_rt["reflectance"],
            "nir": nir_rt["reflectance"],
            "closed": state["closed"],
            "L": state["L"], "ustar": state["ustar"],
            "stability_closed": state["stability_closed"],
            "stability_iterations": state["stability_iterations"],
            "stability_residual_m_inv": state["stability_residual_m_inv"],
            "iterations": state["iterations"],
            "energy_balance_residual_W_m2": state["energy_balance_residual_W_m2"],
            **_diagnostics(state, name),
        })
        intermediate.append(_intermediate_row(state, node, meteo["time"], name, totals))
        final_states.append(state)
        previous = state

    unclosed = [row for row in output if not row["closed"]]
    if unclosed:
        maximum_residual = max(row["energy_balance_residual_W_m2"] for row in unclosed)
        warnings.warn(
            f"Energy balance did not close at {len(unclosed)}/{len(output)} time steps; "
            f"maximum component residual {maximum_residual:g} W/m²; see the closed and energy_balance_residual_W_m2 columns in output.csv",
            RuntimeWarning,
            stacklevel=2,
        )
    message = clearness_index_warning(final_states)
    if message:
        warnings.warn(message, RuntimeWarning, stacklevel=2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_snapshot = _snapshot_path(input_path, output_path.parent)
    _write_input_snapshot(input_path, input_snapshot, config, overrides, input_data)
    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(output)

    intermediate_path = output_path.with_name("intermediate.csv")
    if input_data.output_intermediate:
        intermediate_fields = ["node", "time", "model"]
        for row in intermediate:
            for key in row:
                if key not in intermediate_fields:
                    intermediate_fields.append(key)
        with intermediate_path.open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=intermediate_fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(intermediate)
    elif intermediate_path.exists():
        intermediate_path.unlink()

    observations_path = output_path.with_name("observations.csv")
    observations = []
    for node, meteo in enumerate(records):
        selected_time = observation_selected[node]
        state = final_states[node]
        for view_zenith, view_azimuth in observation_directions:
            if not selected_time:
                continue
            results = {}
            for wavelength in observation_wavelengths:
                quantity = _observation_quantity(config, wavelength)
                values = _band_observation(settings, meteo, state, name, wavelength, view_zenith, view_azimuth,
                                           tir_reflectance=quantity == "reflectance")
                results.update(_observation_columns(wavelength, values, quantity))
            observations.append(ObservationStructure(
                node=node, time=meteo["time"], model=name,
                surface_model=settings.get("surface_model", "legacy"),
                vza=view_zenith, vaa=view_azimuth,
                sza=state["sza"], saa=state["saa"], results=results,
            ).as_row())
    with observations_path.open("w", newline="", encoding="utf-8") as file:
        observation_fields = list(OBSERVATION_BASE_FIELDS)
        for wavelength in observation_wavelengths:
            quantity = _observation_quantity(config, wavelength)
            observation_fields.extend(_observation_columns(wavelength, {}, quantity).keys())
        writer = csv.DictWriter(file, fieldnames=observation_fields)
        writer.writeheader()
        writer.writerows(observations)
    return output_path


def run(input_path=DEFAULT_INPUT, output_path=None, *, output_dir=None,
        plots=True, show=False, **overrides):
    """JSON jobs export two time-series figures/CSVs; legacy CSV calls are preserved."""
    source = Path(input_path).resolve()
    if source.suffix.lower() == ".json":
        if output_path is not None:
            raise ValueError("JSON jobs use output_dir/--output-dir, not --output")
        from base.energy_runner import run_job
        return run_job(source, _run_csv, INPUT_FILE_KEYS, output_dir=output_dir,
                       plots=plots, show=show, overrides=overrides)
    if output_dir is not None or show:
        raise ValueError("CSV inputs use output_path/--output; JSON jobs support output_dir and plots")
    return _run_csv(source, output_path, **overrides)


SOIL_METHODS = ("direct", "surface", "depth")
STOMATAL_METHODS = ("farquhar", "ballberry")
SOIL_TEMPERATURE_METHODS = ("direct", "robust")
COMPARISON_FIELDS = ("soil_heat_method", "stomatal_method", "soil_temperature_method")


def run_comparison(input_path=DEFAULT_INPUT, output_root=None, start_time=None, end_time=None, max_iter=None):
    """Compare soil heat, stomatal and soil temperature methods for this scene."""
    input_path = Path(input_path).resolve()
    if input_path.suffix.lower() == ".json":
        raise ValueError("--compare uses CSV: first run the JSON, then pass its generated input.csv")
    config = read_input(input_path)
    input_data = InputStructure.from_config(config, input_path.parent)
    output_root = Path(output_root).resolve() if output_root else input_data.case_dir.resolve()
    summary = []
    for soil_heat_method in SOIL_METHODS:
        for stomatal_method in STOMATAL_METHODS:
            for soil_temperature_method in SOIL_TEMPERATURE_METHODS:
                case_dir = output_root / f"{soil_heat_method}_{stomatal_method}_{soil_temperature_method}"
                overrides = {
                    "case_result_dir": str(case_dir),
                    "soil_heat_method": soil_heat_method,
                    "stomatal_method": stomatal_method,
                    "soil_temperature_method": soil_temperature_method,
                    "observation_time_mode": "all",
                }
                if start_time is not None:
                    overrides["start_time"] = start_time
                if end_time is not None:
                    overrides["end_time"] = end_time
                if max_iter is not None:
                    overrides["max_iter"] = max_iter
                result = run(input_path, **overrides)
                methods = {"soil_heat_method": soil_heat_method, "stomatal_method": stomatal_method,
                           "soil_temperature_method": soil_temperature_method}
                with Path(result).open(newline="", encoding="utf-8") as file:
                    summary.extend(dict(methods, **row) for row in csv.DictReader(file))
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=COMPARISON_FIELDS + tuple(OUTPUT_FIELDS))
        writer.writeheader()
        writer.writerows(summary)
    return summary_path


def main(argv=None, *, default_input=None, default_output_dir=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--max-iter", type=int)
    parser.add_argument("--output-dir", type=Path, help="JSON output directory")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args(argv)
    # Preserve the old CSV default for --compare; normal runs use input_tirteb.csv
    # unless --input points to a JSON scene.
    input_path = args.input or default_input or DEFAULT_INPUT
    if args.show and args.no_plots:
        parser.error("--show cannot be combined with --no-plots")
    if args.compare:
        if args.output_dir is not None or args.show:
            parser.error("--compare uses --output-root and does not plot")
        print(run_comparison(input_path, args.output_root, args.start_time, args.end_time, args.max_iter))
    else:
        if args.output_dir is None and Path(input_path).suffix.lower() == ".json":
            args.output_dir = default_output_dir
        print(run(input_path, args.output, output_dir=args.output_dir,
                  plots=not args.no_plots, show=args.show,
                  start_time=args.start_time, end_time=args.end_time, max_iter=args.max_iter))


if __name__ == "__main__":
    # 默认读取 input_tirteb.csv；JSON 场景用 --input 指定。
    main()
