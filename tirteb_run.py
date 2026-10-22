import argparse
import csv
import json
import shutil
import warnings
from pathlib import Path

import numpy as np

from base.brightness import urban_canopy_fraction, urban_component_fractions, urban_directional_reflectance
from base.tirteb_tirt import tir_band_settings, tir_spectral_value, tirt_brightness
from base.modelio import read_defined, read_directions, read_input, read_meteo, read_spectrals
from base.model import run_node
from base.solar import slope_direction
from base.structs import InputStructure, ObservationStructure, RTStructure
from base.radiation import component_proportions, directional_reflectance, legacy_directional_reflectance


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "input_tirteb.csv"
FEATURE_BANDS = ((0.55, "optical"), (0.86, "infrared"), (10.5, "thermal_infrared"))
DEFAULT_OBSERVATION_VZA = (0.0, 30.0, 60.0)
DEFAULT_OBSERVATION_VAA = (0.0, 90.0, 180.0)
OUTPUT_FIELDS = (
    "node", "time", "veg", "surface", "Ta_C",
    "Tss_C", "Tsh_C", "Tls_C", "Tlh_C", "LE", "H", "G",
    "Troof_sunlit_C", "Troof_shaded_C", "Twall_sunlit_C", "Twall_shaded_C",
    "Tstreet_sunlit_C", "Tstreet_shaded_C", "Tleaf_sunlit_C", "Tleaf_shaded_C",
    "Tv_C", "red", "nir",
    "closed", "iterations", "energy_balance_residual_W_m2",
)
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


def _config_numbers(config, key, default):
    value = config.get(key)
    if value in (None, ""):
        return tuple(default)
    text = str(value).replace(";", ",")
    return tuple(float(item.strip()) for item in text.split(",") if item.strip())


def _spectral_values(settings, wavelength, model=""):
    wavelength_nm = float(wavelength) * 1000.0
    index = int(round(wavelength_nm - 400.0))
    leaf_rho = settings.get("leaf_reflectance", 0.0)
    leaf_tau = settings.get("leaf_transmittance", 0.0)
    soil_rho = settings.get("soil_reflectance", 0.0)
    leaf_spectrum = settings.get("leaf_reflectance_spectrum", [])
    tau_spectrum = settings.get("leaf_transmittance_spectrum", [])
    soil_spectrum = settings.get("soil_reflectance_spectrum", [])
    if 0 <= index < 2001:
        if len(leaf_spectrum) == 2001:
            leaf_rho = leaf_spectrum[index]
        if len(tau_spectrum) == 2001:
            leaf_tau = tau_spectrum[index]
        if len(soil_spectrum) == 2001:
            soil_rho = soil_spectrum[index]
        if model == "urban":
            urban_street = settings.get("urban_reflectance_street_spectrum", [])
            if len(urban_street) == len(settings.get("wavelengths", [])):
                urban_index = int(round(wavelength_nm - 400.0))
                if 0 <= urban_index < len(urban_street):
                    soil_rho = urban_street[urban_index]
    if float(wavelength) >= 2.5:
        leaf_rho = _tir_spectral_value(settings, "leaf_reflectance_tir_spectrum", "leaf_reflectance", wavelength)
        leaf_tau = _tir_spectral_value(settings, "leaf_transmittance_tir_spectrum", "leaf_transmittance", wavelength)
        soil_rho = _tir_spectral_value(settings, "soil_reflectance_tir_spectrum", "soil_reflectance", wavelength)
    return float(leaf_rho), float(leaf_tau), float(soil_rho)


_tir_spectral_value = tir_spectral_value
_tir_band_settings = tir_band_settings


def _band_observation(settings, meteo, state, name, wavelength, vza, vaa):
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
    if settings.get("surface_model") == "terrain" and int(settings.get("terrain_type", 1)) == 1:
        local_vza, local_vaa, _ = slope_direction(
            vza,
            vaa,
            state.get("terrain_slope", 0.0),
            state.get("terrain_aspect", 0.0),
        )
    else:
        local_vza, local_vaa = float(vza), float(vaa)
    sunlit = state.get("sunlit_surface", True)
    leaf_rho, leaf_tau, soil_rho = _spectral_values(settings, wavelength, name)
    if name == "urban" and float(wavelength) < 2.5:
        # The state holds the meteorology-derived diffuse fraction; settings
        # may still contain -1 (automatic partition), which is not a fraction.
        from base.radiation import _shortwave_diffuse_fraction
        diffuse = state.get("shortwave_diffuse_fraction")
        if diffuse is None:
            diffuse = _shortwave_diffuse_fraction(settings, meteo["shortwave"], state["sza"])
        reflectance = urban_directional_reflectance(
            dict(settings, shortwave_diffuse_fraction=diffuse),
            wavelength,
            local_vza,
            local_vaa,
            state["sza"],
            state["saa"],
        )
    else:
        observation_reflectance = directional_reflectance if float(wavelength) < 2.5 else legacy_directional_reflectance
        illumination = ({'diffuse_fraction': state.get('shortwave_diffuse_fraction',0.0)}
                        if float(wavelength) < 2.5 else {})
        reflectance = observation_reflectance(
            settings,
            state["sza"],
            state["saa"],
            local_vza,
            local_vaa,
            leaf_rho,
            leaf_tau,
            soil_rho,
            sunlit=sunlit,
            use_hapke=(
                float(wavelength) < 2.5
                and int(settings.get("hapke_enable", settings.get("hapke_optical_enable", 0))) == 1
            ),
            **illumination,
        )
    proportions = component_proportions(
        settings, state["sza"], state["saa"], local_vza, local_vaa, sunlit=sunlit
    )
    leaf_fraction = float(proportions[2] + proportions[3])
    if name == "urban":
        urban_fractions = urban_component_fractions(settings, local_vza, local_vaa, state["sza"], state["saa"])
        street_fraction = float(urban_fractions[2] + urban_fractions[3])
        transmittance = street_fraction * urban_canopy_fraction(settings, local_vza, local_vaa) * leaf_fraction * leaf_tau
    else:
        transmittance = leaf_fraction * leaf_tau
    if wavelength >= 2.5:
        band_settings = _tir_band_settings(settings, wavelength)
        # Thermal observations use the TiRT directional model at the energy-balance temperatures.
        radiance, brightness_temperature = tirt_brightness(band_settings, state, wavelength, vza, vaa)
        reflected = float("nan")
    else:
        radiance = float("nan")
        brightness_temperature = float("nan")
        reflected = max(float(meteo["shortwave"]), 0.0) * reflectance
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
        return "red"
    if abs(float(wavelength) - 0.86) < 1e-6:
        return "nir"
    return _band_name(wavelength)


def _wavelength_tag(wavelength):
    return f"{float(wavelength):g}".replace(".", "p")


def _urban_surface_totals(settings, state):
    density = max(float(settings.get("urban_density", 0.0015)), 0.0)
    length = max(float(settings.get("urban_building_length", 10.0)), 0.0)
    width = max(float(settings.get("urban_building_width", 10.0)), 0.0)
    height = max(float(settings.get("urban_building_height", 20.0)), 0.0)
    roof = density * length * width
    wall = density * 2.0 * height * (length + width)
    street = max(1.0 - roof, 0.0)
    total = max(roof + wall + street, 1.0e-12)
    weights = {"roof": roof / total, "wall": wall / total, "street": street / total}
    canopy = urban_canopy_fraction(settings, 0.0)
    canopy_weight = weights["street"] * canopy
    weights["street"] *= 1.0 - canopy
    totals = {"LE": 0.0, "H": 0.0, "G": 0.0}
    for component, weight in weights.items():
        for suffix in ("sunlit", "shaded"):
            totals["LE"] += weight * 0.5 * state.get("LE_" + component + suffix, 0.0)
            totals["H"] += weight * 0.5 * state.get("H_" + component + suffix, 0.0)
            totals["G"] += weight * 0.5 * state.get("G_" + component + suffix, 0.0)
    for suffix in ("sunlit", "shaded"):
        totals["LE"] += canopy_weight * 0.5 * state.get("LE_leaf" + suffix, 0.0)
        totals["H"] += canopy_weight * 0.5 * state.get("H_leaf" + suffix, 0.0)
        totals["G"] += canopy_weight * 0.5 * state.get("G_leaf" + suffix, 0.0)
    return totals


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


def _observation_hour_fractions(config):
    """Parse hour selectors as fractions of a day; *.5 means every time ending in .5."""
    raw = config.get("observation_hours", config.get("observation_times", ""))
    fractions = []
    for item in str(raw).replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        if item.startswith("*"):
            item = item[1:]
        try:
            value = float(item)
        except ValueError:
            continue
        fraction = value if abs(value) <= 1.0 else value / 24.0
        fractions.append(fraction % 1.0)
    return tuple(fractions)


def _matches_observation_hour(time_value, fractions, tolerance=1.0e-7):
    fraction = float(time_value) % 1.0
    return any(abs(fraction - target) <= tolerance for target in fractions)


def _csv_value(value):
    if isinstance(value, np.ndarray):
        return json.dumps(value.tolist(), separators=(",", ":"))
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, separators=(",", ":"), default=str)
    return value


def _intermediate_row(state, node, time_value, model, proportions, totals):
    row = {"node": node, "time": time_value, "model": model}
    row.update({key: _csv_value(value) for key, value in state.as_dict().items()})
    for index, value in enumerate(proportions):
        row[f"component_fraction_{index}"] = value
    row.update({f"{key}_total": value for key, value in totals.items()})
    return row


def _write_input_snapshot(input_path, snapshot, config, overrides, input_data):
    """Keep the source, effective overrides and local copies of referenced files."""
    original = input_path.read_bytes()
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    source_snapshot = snapshot.with_name("source_input.csv")
    if source_snapshot.resolve() != input_path.resolve() and (
        input_path.resolve() != snapshot.resolve() or not source_snapshot.exists()
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
    with snapshot.open("w", newline="", encoding="utf-8") as file:
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


def run(input_path=DEFAULT_INPUT, output_path=None, **overrides):
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
    observation_wavelengths = _config_numbers(
        config, "observation_wavelengths", [band[0] for band in FEATURE_BANDS]
    )
    observation_wavelengths = input_data.observation_wavelengths or observation_wavelengths
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
    if input_data.observation_spectral_mode in {"file", "txt", "read"}:
        observation_wavelengths = read_spectrals(input_data.spectrals_file)
    observation_time_mode = input_data.observation_time_mode
    if observation_time_mode in {"all", "every", "full"}:
        selected_times = None
        selected_fractions = None
    elif observation_time_mode in {"time", "selected", "specific", "some", "exact"}:
        selected_times = {round(value, 8) for value in observation_times}
        if not selected_times:
            raise ValueError("observation_time_mode=time requires observation_times")
        selected_fractions = None
    elif observation_time_mode in {"hour", "hours", "fraction", "fractional"}:
        selected_times = None
        selected_fractions = _observation_hour_fractions(config)
        if not selected_fractions:
            raise ValueError("observation_time_mode=hour requires observation_hours or *.5 in observation_times")
    else:
        raise ValueError("observation_time_mode must be all, time, or hour")

    all_records = read_meteo(meteo_path)
    records = [
        record for record in all_records
        if input_data.start_time <= float(record["time"]) <= input_data.end_time
    ]
    if not records:
        raise ValueError(
            f"start_time={input_data.start_time:g} 和 end_time={input_data.end_time:g} 未选中气象数据"
        )
    observation_selected = [
        (selected_times is None and selected_fractions is None)
        or (selected_times is not None and round(float(record["time"]), 8) in selected_times)
        or (selected_fractions is not None and _matches_observation_hour(record["time"], selected_fractions))
        for record in records
    ]
    if not any(observation_selected):
        selector = config.get("observation_hours", config.get("observation_times", "")) if selected_fractions is not None else config.get("observation_times", "")
        available_times = [float(record["time"]) for record in records]
        suggestions = sorted({
            min(available_times, key=lambda value: abs(value - target))
            for target in selected_times
        })[:5] if selected_times is not None else available_times[:5]
        available = ";".join(f"{value:.12g}" for value in suggestions)
        raise ValueError(
            f"观测时间 {selector!r} 未命中模拟区间内的气象记录 "
            f"({float(records[0]['time']):.12g}–{float(records[-1]['time']):.12g})；"
            f"可选实际时间：{available}。请指定已有气象时间或使用 observation_time_mode=all"
        )
    settings = read_defined(None, config)
    if model_name(settings) == "urban":
        scene_height = max(
            float(settings.get("urban_building_height", 0.0)),
            float(settings.get("height", 0.0)),
        )
        settings["z"] = max(15.0, scene_height)
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
            if int(settings.get("terrain_type", 1)) == 1:
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
        red_rt = _band_observation(settings, meteo, state, name, 0.55, vertical_vza, vertical_vaa)
        nir_rt = _band_observation(settings, meteo, state, name, 0.86, vertical_vza, vertical_vaa)
        proportions = component_proportions(
            settings,
            state["sza"],
            state["saa"],
            state["vza"],
            state["vaa"],
            sunlit=state.get("sunlit_surface", True),
        )
        if name == "urban":
            urban_totals = _urban_surface_totals(settings, state)
            le_total, h_total, g_soil = urban_totals["LE"], urban_totals["H"], urban_totals["G"]
        else:
            le_total = sum(weight * state.get(key, 0.0) for weight, key in zip(
                proportions, ("LE_soilsunlit", "LE_soilshaded", "LE_leafsunlit", "LE_leafshaded"),
            ))
            h_total = sum(weight * state.get(key, 0.0) for weight, key in zip(
                proportions, ("H_soilsunlit", "H_soilshaded", "H_leafsunlit", "H_leafshaded"),
            ))
            g_soil = proportions[0] * state.get("G_sunlit", 0.0) + proportions[1] * state.get("G_shaded", 0.0)
        totals = {"LE": le_total, "H": h_total, "G": g_soil}
        leaf_temperature_sunlit = (
            float("nan") if name == "bare" else state["Tls"] - 273.15
        )
        leaf_temperature_shaded = (
            float("nan") if name == "bare" else state["Tlh"] - 273.15
        )
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
            "LE": le_total, "H": h_total, "G": g_soil,
            "Tv_C": vertical_rt["brightness_temperature_K"] - 273.15,
            "red": red_rt["reflectance"],
            "nir": nir_rt["reflectance"],
            "closed": state["closed"],
            "iterations": state["iterations"],
            "energy_balance_residual_W_m2": state["energy_balance_residual_W_m2"],
        })
        intermediate.append(_intermediate_row(state, node, meteo["time"], name, proportions, totals))
        final_states.append(state)
        previous = state

    unclosed = [row for row in output if not row["closed"]]
    if unclosed:
        maximum_residual = max(row["energy_balance_residual_W_m2"] for row in unclosed)
        warnings.warn(
            f"{len(unclosed)}/{len(output)} 个时间点的能量平衡未闭合；"
            f"最大组分残差 {maximum_residual:g} W/m²，详见 output.csv 的 closed 和 energy_balance_residual_W_m2",
            RuntimeWarning,
            stacklevel=2,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_snapshot = output_path.parent / "input.csv"
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
                values = _band_observation(settings, meteo, state, name, wavelength, view_zenith, view_azimuth)
                quantity = _observation_quantity(config, wavelength)
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


SOIL_METHODS = ("direct", "surface", "depth")
STOMATAL_METHODS = ("farquhar", "ballberry")
SOIL_TEMPERATURE_METHODS = ("direct", "robust")


def run_comparison(input_path=DEFAULT_INPUT, output_root=None, start_time=None, end_time=None, max_iter=None):
    """Compare soil heat, stomatal and soil temperature methods for this scene."""
    input_path = Path(input_path).resolve()
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
                with Path(result).open(newline="", encoding="utf-8") as file:
                    summary.extend(csv.DictReader(file))
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(summary)
    return summary_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--max-iter", type=int)
    args = parser.parse_args()
    if args.compare:
        print(run_comparison(args.input, args.output_root, args.start_time, args.end_time, args.max_iter))
    else:
        run(args.input, args.output)


if __name__ == "__main__":
    main()
