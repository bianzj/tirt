import math
from datetime import date, timedelta

from .aerodynamic import aerodynamic_resistance
from .biochemical import suresist, suresist_urban
from .eb import rebalance, rebalance_urban
from .evapo import evapotranspiration, urban_evapotranspiration, update_thermal_history
from .radiation import (
    merge_radiation,
    prepare_longwave,
    prepare_shortwave,
    radiative_transfer_diagnostics,
)
from .solar import slope_direction, solar_angles
from .structs import FluxStructure


def run_node(settings, meteo, node, year, doy, latitude, longitude, dem=0.0, slope=0.0, vza=0.0, vaa=0.0, initial_state=None, max_iter=50, aspect=0.0, utc_offset_hours=0.0):
    # Meteorological time is local time; solar geometry is evaluated in UTC.
    local_time = float(meteo["time"])
    utc_time = local_time - float(utc_offset_hours) / 24.0
    local_date = date(int(year), 1, 1) + timedelta(days=int(doy) - 1 + math.floor(local_time) - int(doy))
    utc_date = local_date + timedelta(days=math.floor(utc_time) - math.floor(local_time))
    utc_doy = (utc_date - date(utc_date.year, 1, 1)).days + 1
    global_sza, global_saa = solar_angles(
        utc_date.year, utc_doy, utc_time, latitude, longitude, timezone=0.0
    )
    global_vza, global_vaa = float(vza), float(vaa)
    sunlit_surface = True
    if settings.get("surface_model") == "terrain":
        if int(settings.get("terrain_type", 1)) == 1:
            slope = float(settings.get("terrain_slope", settings.get("slope", slope)))
            aspect = float(settings.get("terrain_aspect", settings.get("aspect", aspect)))
            sza, saa, sunlit_surface = slope_direction(global_sza, global_saa, slope, aspect)
            vza, vaa, _ = slope_direction(global_vza, global_vaa, slope, aspect)
        else:
            slope = 0.0
            aspect = 0.0
            sza, saa = global_sza, global_saa
            vza, vaa = global_vza, global_vaa
            sunlit_surface = global_sza < 90.0
    else:
        slope, aspect = 0.0, 0.0
        sza, saa = global_sza, global_saa
    is_urban = (
        str(settings.get("surface_model", "")).strip().lower() == "urban"
        or int(settings.get("surface_type", -1)) == 13
    )
    soil_initial = float(settings.get("soil_temperature_c", 0.02)) + 273.15
    urban_initial = float(settings.get("urban_temperature_c", settings.get("soil_temperature_c", 0.02))) + 273.15
    wall_layers = max(int(settings.get("urban_wall_layers", 3)), 1)
    wall_boundary = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
    wall_profile = [wall_boundary] * (wall_layers + 1)
    surface_initial = urban_initial if is_urban else soil_initial
    state = {
        "Tss": soil_initial, "Tsh": soil_initial, "Tls": surface_initial, "Tlh": surface_initial,
        "Troof_sunlit": surface_initial, "Troof_shaded": surface_initial,
        "Twall_sunlit": surface_initial, "Twall_shaded": surface_initial,
        "L": 1e9, "dG_term": math.sqrt(settings["soil_cs"] * settings["soil_rhos"] * settings["soil_lambda"]) * 2 / (math.sqrt(math.pi) * math.sqrt(3600)),
        "dt_seconds": 3600.0,
        "surface_history_sunlit": [soil_initial], "surface_history_shaded": [soil_initial],
        "surface_history_roof_sunlit": [urban_initial], "surface_history_roof_shaded": [urban_initial],
        "surface_history_wall_sunlit": [urban_initial], "surface_history_wall_shaded": [urban_initial],
        "surface_history_street_sunlit": [soil_initial], "surface_history_street_shaded": [soil_initial],
        "depth_profile_sunlit": [soil_initial] * 4, "depth_profile_shaded": [soil_initial] * 4,
        "depth_profile_roof_sunlit": [urban_initial] * 4, "depth_profile_roof_shaded": [urban_initial] * 4,
        "depth_profile_wall_sunlit": wall_profile[:], "depth_profile_wall_shaded": wall_profile[:],
        "depth_profile_street_sunlit": [soil_initial] * 4, "depth_profile_street_shaded": [soil_initial] * 4,
        "cs": settings.get("ca", 380.0), "ci": settings.get("ca", 380.0),
        "es": settings.get("ea", meteo["vapor_pressure"]),
    }
    if initial_state:
        for key in ("Tss", "Tsh", "Tls", "Tlh", "Troof_sunlit", "Troof_shaded", "Twall_sunlit", "Twall_shaded", "L", "dt_seconds",
                    "surface_history_sunlit", "surface_history_shaded",
                    "surface_history_roof_sunlit", "surface_history_roof_shaded",
                    "surface_history_wall_sunlit", "surface_history_wall_shaded",
                    "surface_history_street_sunlit", "surface_history_street_shaded",
                    "depth_profile_sunlit", "depth_profile_shaded",
                    "depth_profile_roof_sunlit", "depth_profile_roof_shaded",
                    "depth_profile_wall_sunlit", "depth_profile_wall_shaded",
                    "depth_profile_street_sunlit", "depth_profile_street_shaded"):
            if key in initial_state:
                state[key] = list(initial_state[key]) if isinstance(initial_state[key], (list, tuple)) else initial_state[key]
        if "time" in initial_state:
            state["dt_seconds"] = max(abs(float(meteo["time"]) - float(initial_state["time"])) * 86400.0, 1.0)
    # Keep the wall profile shape and its indoor boundary fixed when carrying state forward.
    for profile_key in ("depth_profile_wall_sunlit", "depth_profile_wall_shaded"):
        profile = list(state.get(profile_key, ()))
        state[profile_key] = (profile[:wall_layers + 1] + wall_profile)[:wall_layers + 1]
        state[profile_key][-1] = wall_boundary
    # Wall illumination must be available during the first radiation iteration.
    state.update(sza=sza, saa=saa, vza=vza, vaa=vaa)
    shortwave = prepare_shortwave(
        settings, meteo, sza, saa, vza, vaa, sunlit=sunlit_surface
    )
    iterations = 0
    closed = False
    netrad = {}

    def evaluate_flux(rebalance_state=True):
        nonlocal netrad
        longwave = prepare_longwave(settings, meteo, state, sza, saa, vza, vaa)
        netrad = merge_radiation(shortwave, longwave)
        aerodynamic_resistance(settings, meteo, state, is_urban)
        if is_urban:
            # Urban aerodynamic_resistance computes separate canyon and leaf
            # transfer paths; only the asphalt street shares its resistance
            # with the legacy soil state used by the radiation/flux helpers.
            state["raa_soil"] = state["raa_street"]
            suresist_urban(settings, meteo, state, netrad)
            _ensure_urban_state(state)
            urban_evapotranspiration(settings, meteo, state, netrad)
            return rebalance_urban(settings, meteo, state) if rebalance_state else False
        suresist(settings, meteo, state, netrad)
        evapotranspiration(settings, meteo, state, netrad)
        return rebalance(settings, meteo, state) if rebalance_state else False

    # The input surface option is authoritative; surface_type=13 is kept for
    # Keep the input surface option authoritative for all scenarios.
    for kiter in range(max_iter):
        closed = evaluate_flux()
        if is_urban:
            _ensure_urban_state(state)
        iterations = kiter
        if closed:
            break

    if is_urban:
        _ensure_urban_state(state)
    update_thermal_history(settings, state)
    state.update(
        time=meteo["time"], sza=sza, saa=saa, vza=vza, vaa=vaa,
        sza_global=global_sza, saa_global=global_saa,
        vza_global=global_vza, vaa_global=global_vaa,
        terrain_slope=slope, terrain_aspect=aspect,
        sunlit_surface=sunlit_surface,
        iterations=iterations, closed=closed,
    )
    state.update(radiative_transfer_diagnostics(settings, state, vza, vaa))
    state.update(
        shortwave_direct=netrad.get("shortwave_direct", 0.0),
        shortwave_diffuse=netrad.get("shortwave_diffuse", 0.0),
        shortwave_global=netrad.get("shortwave_global", 0.0),
        shortwave_direct_incident=netrad.get("shortwave_direct_incident", 0.0),
        shortwave_diffuse_incident=netrad.get("shortwave_diffuse_incident", 0.0),
        shortwave_diffuse_fraction=netrad.get("shortwave_diffuse_fraction", 0.0),
        longwave_down=netrad.get("longwave_down", 0.0),
        tir_scatter_soil=netrad.get("tir_scatter_soil", 0.0),
        tir_scatter_leaf=netrad.get("tir_scatter_leaf", 0.0),
    )
    if not is_urban:
        for component in ("soil", "leaf"):
            for suffix, sunlit in (("sunlit", True), ("shaded", False)):
                direct = float(netrad.get(f"directVrad_{component}", 0.0)) if sunlit else 0.0
                diffuse = float(netrad.get(f"diffuseVrad_{component}", 0.0))
                longwave = float(netrad.get(f"diffuseTrad_{component}", 0.0))
                state[f"SW_direct_{component}_{suffix}"] = direct
                state[f"SW_diffuse_{component}_{suffix}"] = diffuse
                state[f"LW_input_{component}_{suffix}"] = longwave
                state[f"radiation_input_{component}_{suffix}"] = direct + diffuse + longwave
    return FluxStructure.from_mapping(state)


def _ensure_urban_state(state):
    state.setdefault("Troof_sunlit", state["Tls"])
    state.setdefault("Troof_shaded", state["Tlh"])
    state.setdefault("Twall_sunlit", state["Tls"])
    state.setdefault("Twall_shaded", state["Tlh"])
    state["Trs"] = state["Troof_sunlit"]
    state["Trh"] = state["Troof_shaded"]
    state["Tws"] = state["Twall_sunlit"]
    state["Twh"] = state["Twall_shaded"]
    state["Tstreet_sunlit"] = state["Tss"]
    state["Tstreet_shaded"] = state["Tsh"]
