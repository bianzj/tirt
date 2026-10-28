import math

from .eb import CP, MAIR, MH20, SIGMASB, air_density, latent_heat, saturation_vapor_pressure
from .heat_conduction import (
    HISTORY_SECONDS, column_surface_flux, conduction_profile, history_intervals, surface_memory_flux,
    thermal_properties, time_step,
)


def _method(value, default="direct"):
    text = str(value if value is not None else default).strip().lower().replace("-", "_")
    if text in {"0", "direct", "net", "netrad", "radiation", "ratio"}:
        return "direct"
    if text in {"1", "surface", "surface_diff", "surface_differential", "inertia"}:
        return "surface"
    if text in {"2", "depth", "depth_diff", "depth_differential", "profile"}:
        return "depth"
    raise ValueError("heat flux method must be direct, surface or depth")


def _component_flux(state, name, temp, ra, rs, incoming, factor, qa, e_to_q, ta, rho, emissivity=1.0):
    ra = max(float(ra), 1.0)
    rs = max(float(rs), 0.0)
    ei = saturation_vapor_pressure(temp - 273.15)
    qi = ei * e_to_q
    lam = latent_heat(temp)
    state["LE_" + name] = rho / (ra + rs) * lam * (qi - qa)
    state["H_" + name] = rho * CP / ra * (temp - ta)
    outgoing_longwave = max(float(emissivity), 0.0) * SIGMASB * temp ** 4 * factor
    state["LW_out_" + name] = outgoing_longwave
    state["N_" + name] = incoming - outgoing_longwave


def _thermal_inertia(settings, urban=False):
    conductivity, capacity = thermal_properties(settings, urban=urban)
    return math.sqrt(conductivity * capacity)


def _surface_differential(settings, state, temp, history_key, urban=False):
    history = list(state.get(history_key, ()))
    if not history:
        return 0.0
    dt = time_step(settings, state)
    intervals = history_intervals(state, len(history) - 1, dt)
    return surface_memory_flux(temp, history, intervals, dt,
                               _thermal_inertia(settings, urban=urban))


def _depth_differential(settings, state, temp, profile_key, urban=False):
    # Uses the same end-of-step profile that update_thermal_history commits.
    return column_surface_flux(settings, state, temp, profile_key, urban=urban)


def _envelope_differential(settings, state, temp, profile_key):
    """Finite-difference conduction through the roof/wall to the fixed indoor
    temperature (urban_internal_temperature_c = 24 C).

    The surface heat flux is the first-layer gradient
        G = urban_lambda * (T_surface - T_node1) / dz
    where dz = urban_wall_thickness / urban_wall_layers.  At steady state the
    profile is linear and any layer count gives the same value
        G = urban_lambda * (T_surface - T_inside) / urban_wall_thickness.
    """
    thickness = max(float(settings.get("urban_wall_thickness", 0.5)), 1.0e-3)
    internal = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
    component = "sunlit" if "sunlit" in profile_key else "shaded"
    nodes, dz, conductivity = conduction_profile(
        settings, state, temp, profile_key, urban=True, envelope=True,
    )
    layers = len(nodes) - 1
    layer_fluxes = [
        conductivity * (nodes[index] - nodes[index + 1]) / dz
        for index in range(layers)
    ]
    for index, flux in enumerate(layer_fluxes, 1):
        state[f"G_wall_layer{index}_{component}"] = flux
    state[f"G_wall_surface_{component}"] = layer_fluxes[0]
    state[f"G_wall_inner_{component}"] = layer_fluxes[-1]
    state[f"G_wall_equivalent_{component}"] = conductivity * (float(temp) - internal) / thickness
    state[f"G_wall_balance_{component}"] = layer_fluxes[0]
    return layer_fluxes[0]


def _heat_flux(settings, state, temp, net, method, history_key, profile_key, fraction, urban=False):
    method = _method(method)
    if method == "direct":
        return float(net) * float(fraction)
    if method == "surface":
        return _surface_differential(settings, state, temp, history_key, urban=urban)
    if urban and ("roof" in profile_key or "wall" in profile_key):
        return _envelope_differential(settings, state, temp, profile_key)
    return _depth_differential(settings, state, temp, profile_key, urban=urban)


def _incoming(netrad, component, sunlit):
    return float(netrad.get(f"{component}_{'sunlit' if sunlit else 'shaded'}", 0.0))


def _urban_absorption_spectrum(settings, component, default):
    reflectance = settings.get("urban_reflectance_" + component + "_spectrum", ())
    if len(reflectance) >= 2001:
        return [
            1.0 - min(max(float(value), 0.0), 0.99)
            for value in reflectance[:2001]
        ]
    return [max(min(float(default), 1.0), 0.0)] * 2001


def _urban_longwave_inputs(settings, state, sky):
    """Reciprocal city radiosity, including the complete canopy/soil boundary."""
    from .urban_energy import longwave_inputs
    return longwave_inputs(settings, state, sky)


def _urban_component_incoming(settings, state, netrad):
    """Absorbed city inputs from a common shortwave/longwave exchange geometry."""
    from .urban_energy import shortwave_inputs
    shortwave = netrad.get("_urban_shortwave")
    if shortwave is None:
        shortwave = shortwave_inputs(settings, state, netrad)
    state.update(shortwave)
    longwave = _urban_longwave_inputs(settings, state, netrad.get("longwave_down", 0.0))
    for key, value in longwave.items():
        if key.startswith("urban_") or key.endswith("_fraction"):
            state[key] = value
    result = {}
    for component in ("roof", "wall", "street", "leaf"):
        inputs = []
        for suffix in ("sunlit", "shaded"):
            lw = float(longwave[f"{component}_{suffix}"])
            state[f"LW_input_{component}_{suffix}"] = lw
            state[f"LW_down_{component}_{suffix}"] = float(longwave[f"LW_down_{component}_{suffix}"])
            total = (float(shortwave.get(f"SW_direct_{component}_{suffix}", 0.0))
                     + float(shortwave.get(f"SW_diffuse_{component}_{suffix}", 0.0)) + lw)
            state[f"radiation_input_{component}_{suffix}"] = total
            state[f"raa_{component}_{suffix}"] = state.get(f"raa_{component}", 0.0)
            state[f"rss_{component}_{suffix}"] = state.get("rssunlit" if suffix == "sunlit" else "rsshaded", 0.0) if component == "leaf" else state.get(f"rss{component}", 0.0)
            inputs.append(total)
        result[component] = tuple(inputs)
    return result


def evapotranspiration(settings, meteo, state, netrad):
    p = max(float(meteo["pressure"]), 1.0e-6)
    e_to_q = MH20 / MAIR / p
    qa = meteo["vapor_pressure"] * e_to_q
    ta = meteo["air_temperature_c"] + 273.15
    rho = air_density(meteo)
    for key in ("LE_leafsunlit", "LE_leafshaded", "H_leafsunlit", "H_leafshaded", "N_leafsunlit", "N_leafshaded"):
        state.setdefault(key, 0.0)
    if settings["lai"] > 0.0:
        # One big leaf per unit ground (bulk resistances): LAI times the per-leaf-area radiation, two faces.
        lai = float(settings["lai"])
        emissivity = float(settings["leaf_emissivity"])
        _component_flux(state, "leafsunlit", state["Tls"], state["raa_leaf"], state["rssunlit"], lai * _incoming(netrad, "leaf", True), 2.0 * lai, qa, e_to_q, ta, rho, emissivity)
        _component_flux(state, "leafshaded", state["Tlh"], state["raa_leaf"], state["rsshaded"], lai * _incoming(netrad, "leaf", False), 2.0 * lai, qa, e_to_q, ta, rho, emissivity)
    soil_emissivity = float(settings["soil_emissivity"])
    _component_flux(state, "soilsunlit", state["Tss"], state["raa_soil"], state["rss"], _incoming(netrad, "soil", True), 1.0, qa, e_to_q, ta, rho, soil_emissivity)
    _component_flux(state, "soilshaded", state["Tsh"], state["raa_soil"], state["rss"], _incoming(netrad, "soil", False), 1.0, qa, e_to_q, ta, rho, soil_emissivity)
    state["G_sunlit"] = _heat_flux(
        settings, state, state["Tss"], state["N_soilsunlit"], settings.get("soil_heat_method", "direct"),
        "surface_history_sunlit", "depth_profile_sunlit", settings.get("soil_heat_fraction", 0.25),
    )
    state["G_shaded"] = _heat_flux(
        settings, state, state["Tsh"], state["N_soilshaded"], settings.get("soil_heat_method", "direct"),
        "surface_history_shaded", "depth_profile_shaded", settings.get("soil_heat_fraction", 0.25),
    )
    return state


def urban_evapotranspiration(settings, meteo, state, netrad):
    p = max(float(meteo["pressure"]), 1.0e-6)
    e_to_q = MH20 / MAIR / p
    qa = meteo["vapor_pressure"] * e_to_q
    ta = meteo["air_temperature_c"] + 273.15
    rho = air_density(meteo)
    incoming = _urban_component_incoming(settings, state, netrad)
    emissivity = {
        "roof": float(settings.get("urban_emissivity_roof", 0.95)),
        "wall": float(settings.get("urban_emissivity_wall", 0.92)),
        "street": float(settings.get("urban_emissivity_street", 0.955)),
    }
    components = (
        ("roof", state["Troof_sunlit"], state["Troof_shaded"], state.get("raa_roof", 1.0), state.get("rssroof", 3662.8562)),
        ("wall", state["Twall_sunlit"], state["Twall_shaded"], state.get("raa_wall", 1.0), state.get("rsswall", 3662.8562)),
        ("street", state["Tss"], state["Tsh"], state.get("raa_street", 1.0), state.get("rssstreet", 3662.8562)),
    )
    for component, temp_sun, temp_shade, ra, rs in components:
        for sunlit, temp in ((True, temp_sun), (False, temp_shade)):
            name = component + ("sunlit" if sunlit else "shaded")
            _component_flux(
                state, name, temp, ra, rs, incoming[component][0 if sunlit else 1],
                1.0, qa, e_to_q, ta, rho, emissivity=emissivity[component],
            )
            state["LE_" + name] = 0.0
            heat_method = (
                settings.get("urban_heat_method", "depth")
                if component in {"roof", "wall"}
                else settings.get("soil_heat_method", "surface")
            )
            state["G_" + name] = _heat_flux(
                settings, state, temp, state["N_" + name], heat_method,
                "surface_history_" + component + ("_sunlit" if sunlit else "_shaded"),
                "depth_profile_" + component + ("_sunlit" if sunlit else "_shaded"),
                (
                    settings.get("urban_heat_fraction", 0.4)
                    if component in {"roof", "wall"}
                    else settings.get("soil_heat_fraction", 0.25)
                ),
                urban=component in {"roof", "wall"},
            )
    lai = float(settings.get("lai", 0.0))
    if lai > 0.0:
        for sunlit, temp, rs in ((True, state["Tls"], state.get("rssunlit", 625000.0)), (False, state["Tlh"], state.get("rsshaded", 625000.0))):
            suffix = "sunlit" if sunlit else "shaded"
            name = "leafsunlit" if sunlit else "leafshaded"
            # Big leaf per unit street ground, as for the plane canopies.
            _component_flux(
                state, name, temp, state.get("raa_leaf", state.get("raa_street", 1.0)), rs,
                lai * state[f"radiation_input_leaf_{suffix}"], 2.0 * lai, qa, e_to_q, ta, rho,
                emissivity=float(settings.get("leaf_emissivity", 0.975)),
            )
            state["G_" + name] = 0.0
    else:
        for suffix in ("sunlit", "shaded"):
            for key in ("LE_", "H_", "N_", "G_"):
                state[key + "leaf" + suffix] = 0.0
    return state


def update_thermal_history(settings, state):
    dt = time_step(settings, state)
    mappings = (
        ("Tss", "sunlit"), ("Tsh", "shaded"),
        ("Troof_sunlit", "roof_sunlit"), ("Troof_shaded", "roof_shaded"),
        ("Twall_sunlit", "wall_sunlit"), ("Twall_shaded", "wall_shaded"),
        ("Tss", "street_sunlit"), ("Tsh", "street_shaded"),
    )
    longest = max(len(state.get("surface_history_" + suffix, ()))
                  for _, suffix in mappings)
    intervals = [dt] + history_intervals(state, max(longest - 1, 0), dt)
    elapsed = 0.0
    retained = 0
    for interval in intervals:
        elapsed += interval
        retained += 1
        if elapsed >= HISTORY_SECONDS:
            break
    # Keep the oldest bracketing sample, so an interval intersecting the 48 h
    # cutoff can be integrated with its original slope and duration.
    for temp_key, suffix in mappings:
        history_key = "surface_history_" + suffix
        state[history_key] = ([float(state[temp_key])]
                              + list(state.get(history_key, ())))[:retained + 1]
        profile_key = "depth_profile_" + suffix
        # Street G uses the soil thermal material; advance that exact material.
        urban = suffix.startswith(("roof_", "wall_"))
        envelope = suffix.startswith(("roof_", "wall_"))
        state[profile_key] = conduction_profile(
            settings, state, float(state[temp_key]), profile_key,
            urban=urban, envelope=envelope,
        )[0]
    state["thermal_history_dt_seconds"] = intervals[:retained]


__all__ = ["evapotranspiration", "urban_evapotranspiration", "update_thermal_history"]
