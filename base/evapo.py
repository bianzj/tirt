import math

from .eb import CP, MAIR, MH20, RHOA, SIGMASB, latent_heat, saturation_vapor_pressure


DT_SECONDS = 3600.0
HISTORY_LENGTH = 12


def _method(value, default="direct"):
    text = str(value if value is not None else default).strip().lower().replace("-", "_")
    if text in {"0", "direct", "net", "netrad", "radiation", "ratio"}:
        return "direct"
    if text in {"1", "surface", "surface_diff", "surface_differential", "inertia"}:
        return "surface"
    if text in {"2", "depth", "depth_diff", "depth_differential", "profile"}:
        return "depth"
    raise ValueError("heat flux method must be direct, surface or depth")


def _component_flux(state, name, temp, ra, rs, incoming, factor, qa, e_to_q, ta, emissivity=1.0):
    ra = max(float(ra), 1.0)
    rs = max(float(rs), 0.0)
    ei = saturation_vapor_pressure(temp - 273.15)
    qi = ei * e_to_q
    lam = latent_heat(temp)
    state["LE_" + name] = RHOA / (ra + rs) * lam * (qi - qa)
    state["H_" + name] = RHOA * CP / ra * (temp - ta)
    outgoing_longwave = max(float(emissivity), 0.0) * SIGMASB * temp ** 4 * factor
    state["LW_out_" + name] = outgoing_longwave
    state["N_" + name] = incoming - outgoing_longwave


def _thermal_inertia(settings, urban=False):
    prefix = "urban_" if urban else "soil_"
    return math.sqrt(max(
        float(settings.get(prefix + "cs", 1180.0 if not urban else 1000000.0))
        * float(settings.get(prefix + "rhos", 1100.0 if not urban else 1800.0))
        * float(settings.get(prefix + "lambda", 1.55 if not urban else 1.5)),
        0.0,
    ))


def _surface_differential(settings, state, temp, history_key, urban=False):
    history = list(state.get(history_key, ()))
    if not history:
        return 0.0
    dt = max(float(state.get("dt_seconds", DT_SECONDS)), 1.0)
    current_c = temp - 273.15
    ssum = 0.0
    for index, previous in enumerate(history[:HISTORY_LENGTH]):
        previous_c = float(previous) - 273.15
        elapsed = dt * (index + 1.0)
        weight = math.sqrt(elapsed) - math.sqrt(max(elapsed - dt, 0.0))
        if index == 0:
            difference = current_c - previous_c
        else:
            older = float(history[index - 1]) - 273.15
            difference = older - previous_c
        ssum += difference / dt * weight
    return 2.0 * _thermal_inertia(settings, urban=urban) / math.sqrt(math.pi) * ssum


def _depth_differential(settings, state, temp, profile_key, urban=False):
    profile = list(state.get(profile_key, ()))
    if len(profile) < 2:
        return 0.0
    prefix = "urban_" if urban else "soil_"
    dz = max(float(settings.get(prefix + "depth_step", 0.05)), 1.0e-3)
    conductivity = max(float(settings.get(prefix + "lambda", 1.55 if not urban else 1.5)), 1.0e-8)
    # Node 0 is the surface boundary; node 1 is the first subsurface node.
    return conductivity * (temp - float(profile[1])) / dz


def _envelope_differential(settings, state, temp, profile_key):
    """Finite-difference conduction through the roof/wall to the fixed indoor
    temperature (urban_internal_temperature_c = 24 C).

    The surface heat flux is the first-layer gradient
        G = urban_lambda * (T_surface - T_node1) / dz
    where dz = urban_wall_thickness / urban_wall_layers.  At steady state the
    profile is linear and any layer count gives the same value
        G = urban_lambda * (T_surface - T_inside) / urban_wall_thickness.
    """
    layers = max(int(settings.get("urban_wall_layers", 3)), 1)
    thickness = max(float(settings.get("urban_wall_thickness", 0.5)), 1.0e-3)
    dz = thickness / layers
    conductivity = max(float(settings.get("urban_lambda", 1.5)), 1.0e-8)
    internal = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
    component = "sunlit" if "sunlit" in profile_key else "shaded"
    profile = list(state.get(profile_key, ()))
    nodes = [float(temp)]
    nodes.extend(float(value) for value in profile[1:layers])
    nodes.append(internal)
    nodes = (nodes[:layers + 1] + [internal] * (layers + 1))[:layers + 1]
    while len(nodes) < layers + 1:
        nodes.append(internal)
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
    direct = netrad.get("directVrad_" + component, 0.0) + netrad.get("directTrad_" + component, 0.0)
    diffuse = netrad.get("diffuseVrad_" + component, 0.0) + netrad.get("diffuseTrad_" + component, 0.0)
    return direct + diffuse if sunlit else diffuse


def _urban_absorption_spectrum(settings, component, default):
    reflectance = settings.get("urban_reflectance_" + component + "_spectrum", ())
    if len(reflectance) >= 2001:
        return [
            1.0 - min(max(float(value), 0.0), 0.99)
            for value in reflectance[:2001]
        ]
    return [max(min(float(default), 1.0), 0.0)] * 2001


def _urban_longwave_inputs(settings, state, sky):
    """Return absorbed urban longwave input using sky and canyon surface emission."""
    sky = max(float(sky), 0.0)
    eps = {
        "roof": min(max(float(settings.get("urban_emissivity_roof", 0.95)), 0.0), 1.0),
        "wall": min(max(float(settings.get("urban_emissivity_wall", 0.92)), 0.0), 1.0),
        "street": min(max(float(settings.get("urban_emissivity_street", 0.955)), 0.0), 1.0),
        "leaf": min(max(float(settings.get("leaf_emissivity", 0.975)), 0.0), 1.0),
    }
    length = max(float(settings.get("urban_building_length", 10.0)), 1.0e-3)
    width = max(float(settings.get("urban_building_width", 10.0)), 1.0e-3)
    height = max(float(settings.get("urban_building_height", 20.0)), 0.0)
    plan_scale = math.sqrt(length * width)
    aspect = height / plan_scale
    street_svf = min(max(1.0 / (1.0 + 2.0 * aspect), 0.05), 0.95)
    # A vertical wall sees at most half of its hemisphere as sky. The
    # street-derived term describes how open that upper half is.
    wall_svf = min(max(0.25 * (1.0 + street_svf), 0.1), 0.45)
    leaf_svf = min(max(0.5 * street_svf + 0.15, 0.1), 0.8)

    temperatures = {
        "roof": (float(state.get("Troof_sunlit", 293.15)), float(state.get("Troof_shaded", 293.15))),
        "wall": (float(state.get("Twall_sunlit", 293.15)), float(state.get("Twall_shaded", 293.15))),
        "street": (float(state.get("Tss", 293.15)), float(state.get("Tsh", 293.15))),
    }
    emission = {
        component: tuple(eps[component] * SIGMASB * temperature ** 4 for temperature in values)
        for component, values in temperatures.items()
    }
    mean_emission = {
        component: sum(values) / len(values)
        for component, values in emission.items()
    }
    result = {}
    for suffix in ("sunlit", "shaded"):
        roof_incident = sky
        wall_surroundings = (
            0.5 * mean_emission["wall"]
            + 0.3 * mean_emission["street"]
            + 0.2 * mean_emission["roof"]
        )
        wall_incident = wall_svf * sky + (1.0 - wall_svf) * wall_surroundings
        street_incident = street_svf * sky + (1.0 - street_svf) * mean_emission["wall"]
        leaf_surroundings = 0.6 * mean_emission["street"] + 0.4 * mean_emission["wall"]
        leaf_incident = leaf_svf * sky + (1.0 - leaf_svf) * leaf_surroundings
        result[f"roof_{suffix}"] = eps["roof"] * roof_incident
        result[f"wall_{suffix}"] = eps["wall"] * wall_incident
        result[f"street_{suffix}"] = eps["street"] * street_incident
        result[f"leaf_{suffix}"] = 2.0 * eps["leaf"] * leaf_incident
        for component, incident in (
            ("roof", roof_incident), ("wall", wall_incident),
            ("street", street_incident), ("leaf", leaf_incident),
        ):
            result[f"LW_down_{component}_{suffix}"] = incident
    result.update(
        urban_street_svf=street_svf,
        urban_wall_svf=wall_svf,
        urban_leaf_svf=leaf_svf,
    )
    return result


def _urban_component_incoming(settings, state, netrad):
    """Calculate absorbed shortwave and longwave input for each urban component."""
    diffuse_soil = float(netrad.get("diffuseVrad_soil", 0.0))
    direct_incident = float(netrad.get("directVrad_incident", 0.0))
    direct_spectrum = list(netrad.get("directVrad_incident_spectrum", ()))
    diffuse_spectrum = list(netrad.get("diffuseVrad_incident_spectrum", ()))
    if len(diffuse_spectrum) != 2001:
        diffuse_spectrum = list(netrad.get("diffuseVrad_soil_incident_spectrum", ()))
    solar_zenith = max(float(state.get("sza", 90.0)), 0.0)
    saa = float(state.get("saa", 0.0))
    # Average the positive solar projection over illuminated facades only.
    # Shaded facades are represented by the separate shaded-wall component.
    wall_length = max(float(settings.get("urban_building_length", 10.0)), 1.0e-3)
    wall_width = max(float(settings.get("urban_building_width", 10.0)), 1.0e-3)
    wall_orientations = (0.0, 90.0, 180.0, 270.0)
    wall_areas = (wall_length, wall_width, wall_length, wall_width)
    wall_area_total = sum(wall_areas)
    illuminated = [
        (area, max(math.cos(math.radians(saa - orientation)), 0.0))
        for orientation, area in zip(wall_orientations, wall_areas)
    ]
    illuminated_area = sum(area for area, projection in illuminated if projection > 1.0e-12)
    wall_orientation = (
        sum(area * projection for area, projection in illuminated) / illuminated_area
        if illuminated_area > 0.0 else 0.0
    )
    wall_sunlit_fraction = illuminated_area / wall_area_total
    wall_direct_factor = 0.0
    if solar_zenith < 90.0:
        sza = min(solar_zenith, 89.0)
        wall_direct_factor = min(3.0, math.tan(math.radians(sza)) * wall_orientation)
    wall_diffuse_factor = min(
        max(float(settings.get("urban_wall_diffuse_factor", 0.5)), 0.0), 1.0
    )
    absorption = {
        "roof": max(min(float(settings.get("urban_absorption_roof", 0.95)), 1.0), 0.0),
        "wall": max(min(float(settings.get("urban_absorption_wall", 0.92)), 1.0), 0.0),
        "street": max(min(float(settings.get("urban_absorption_street", 0.955)), 1.0), 0.0),
    }
    direct_factor = {"roof": 1.0, "wall": wall_direct_factor, "street": 1.0}
    diffuse_factor = {"roof": 1.0, "wall": wall_diffuse_factor, "street": 1.0}
    street_factor = max(float(settings.get("urban_street_radiation_factor", 1.0)), 0.0)
    longwave = _urban_longwave_inputs(settings, state, netrad.get("longwave_down", 0.0))
    state["wall_orientation"] = wall_orientation
    state["wall_direct_factor"] = wall_direct_factor
    state["wall_sunlit_fraction"] = wall_sunlit_fraction
    state["wall_diffuse_factor"] = wall_diffuse_factor
    for key in ("urban_street_svf", "urban_wall_svf", "urban_leaf_svf"):
        state[key] = longwave[key]

    for suffix in ("sunlit", "shaded"):
        state[f"LW_input_leaf_{suffix}"] = longwave[f"leaf_{suffix}"]
        state[f"LW_down_leaf_{suffix}"] = longwave[f"LW_down_leaf_{suffix}"]

    result = {}
    for component in ("roof", "wall", "street"):
        absorption_spectrum = _urban_absorption_spectrum(settings, component, absorption[component])
        if len(direct_spectrum) == 2001 and len(diffuse_spectrum) == 2001:
            direct = direct_factor[component] * sum(
                source * band_absorption
                for source, band_absorption in zip(direct_spectrum, absorption_spectrum)
            )
            diffuse_sw = diffuse_factor[component] * sum(
                source * band_absorption
                for source, band_absorption in zip(diffuse_spectrum, absorption_spectrum)
            )
        else:
            direct = direct_incident * direct_factor[component] * absorption[component]
            diffuse_sw = diffuse_soil * diffuse_factor[component] * absorption[component]
        overall_factor = street_factor if component == "street" else 1.0
        direct_input = direct * overall_factor
        diffuse_input = diffuse_sw * overall_factor
        for suffix in ("sunlit", "shaded"):
            longwave_input = longwave[f"{component}_{suffix}"]
            state[f"SW_direct_{component}_{suffix}"] = direct_input if suffix == "sunlit" else 0.0
            state[f"SW_diffuse_{component}_{suffix}"] = diffuse_input
            state[f"LW_input_{component}_{suffix}"] = longwave_input
            state[f"LW_down_{component}_{suffix}"] = longwave[f"LW_down_{component}_{suffix}"]
            state[f"radiation_input_{component}_{suffix}"] = (
                state[f"SW_direct_{component}_{suffix}"] + diffuse_input + longwave_input
            )
            state[f"raa_{component}_{suffix}"] = state.get(f"raa_{component}", 0.0)
            state[f"rss_{component}_{suffix}"] = state.get(f"rss{component}", 0.0)
        result[component] = (
            direct_input + diffuse_input + longwave["%s_sunlit" % component],
            diffuse_input + longwave["%s_shaded" % component],
        )
    return result

def evapotranspiration(settings, meteo, state, netrad):
    p = max(float(meteo["pressure"]), 1.0e-6)
    e_to_q = MH20 / MAIR / p
    qa = meteo["vapor_pressure"] * e_to_q
    ta = meteo["air_temperature_c"] + 273.15
    for key in ("LE_leafsunlit", "LE_leafshaded", "H_leafsunlit", "H_leafshaded", "N_leafsunlit", "N_leafshaded"):
        state.setdefault(key, 0.0)
    if settings["lai"] > 0.0:
        _component_flux(state, "leafsunlit", state["Tls"], state["raa_leaf"], state["rssunlit"], _incoming(netrad, "leaf", True), 2.0, qa, e_to_q, ta)
        _component_flux(state, "leafshaded", state["Tlh"], state["raa_leaf"], state["rsshaded"], _incoming(netrad, "leaf", False), 2.0, qa, e_to_q, ta)
    _component_flux(state, "soilsunlit", state["Tss"], state["raa_soil"], state["rss"], _incoming(netrad, "soil", True), 1.0, qa, e_to_q, ta)
    _component_flux(state, "soilshaded", state["Tsh"], state["raa_soil"], state["rss"], _incoming(netrad, "soil", False), 1.0, qa, e_to_q, ta)
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
                1.0, qa, e_to_q, ta, emissivity=emissivity[component],
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
    for suffix, sunlit in (("sunlit", True), ("shaded", False)):
        state[f"SW_direct_leaf_{suffix}"] = float(netrad.get("directVrad_leaf", 0.0)) if sunlit else 0.0
        state[f"SW_diffuse_leaf_{suffix}"] = float(netrad.get("diffuseVrad_leaf", 0.0))
        state[f"LW_input_leaf_{suffix}"] = float(state.get(f"LW_input_leaf_{suffix}", 0.0))
        state[f"LW_down_leaf_{suffix}"] = float(state.get(f"LW_down_leaf_{suffix}", 0.0))
        state[f"radiation_input_leaf_{suffix}"] = (
            state[f"SW_direct_leaf_{suffix}"]
            + state[f"SW_diffuse_leaf_{suffix}"]
            + state[f"LW_input_leaf_{suffix}"]
        )
        state[f"raa_leaf_{suffix}"] = state.get("raa_leaf", 0.0)
        state[f"rss_leaf_{suffix}"] = state.get("rssunlit" if sunlit else "rsshaded", 0.0)
    if float(settings.get("lai", 0.0)) > 0.0:
        for sunlit, temp, rs in ((True, state["Tls"], state.get("rssunlit", 625000.0)), (False, state["Tlh"], state.get("rsshaded", 625000.0))):
            suffix = "sunlit" if sunlit else "shaded"
            name = "leafsunlit" if sunlit else "leafshaded"
            leaf_incoming = (
                state[f"SW_direct_leaf_{suffix}"]
                + state[f"SW_diffuse_leaf_{suffix}"]
                + state[f"LW_input_leaf_{suffix}"]
            )
            _component_flux(
                state, name, temp, state.get("raa_leaf", state.get("raa_street", 1.0)), rs,
                leaf_incoming, 2.0, qa, e_to_q, ta,
                emissivity=float(settings.get("leaf_emissivity", 0.975)),
            )
            state["G_" + name] = 0.0
    else:
        for suffix in ("sunlit", "shaded"):
            for key in ("LE_", "H_", "N_", "G_"):
                state[key + "leaf" + suffix] = 0.0
    return state


def update_thermal_history(settings, state):
    for temp_key, history_key in (("Tss", "surface_history_sunlit"), ("Tsh", "surface_history_shaded")):
        state[history_key] = [float(state[temp_key])] + list(state.get(history_key, ()))[:HISTORY_LENGTH - 1]
    for component, temp_sun, temp_shade in (("roof", "Troof_sunlit", "Troof_shaded"), ("wall", "Twall_sunlit", "Twall_shaded"), ("street", "Tss", "Tsh")):
        for suffix, key in (("sunlit", temp_sun), ("shaded", temp_shade)):
            history_key = "surface_history_" + component + "_" + suffix
            state[history_key] = [float(state[key])] + list(state.get(history_key, ()))[:HISTORY_LENGTH - 1]

    profiles = (
        ("Tss", "depth_profile_sunlit"), ("Tsh", "depth_profile_shaded"),
        ("Troof_sunlit", "depth_profile_roof_sunlit"), ("Troof_shaded", "depth_profile_roof_shaded"),
        ("Twall_sunlit", "depth_profile_wall_sunlit"), ("Twall_shaded", "depth_profile_wall_shaded"),
        ("Tss", "depth_profile_street_sunlit"), ("Tsh", "depth_profile_street_shaded"),
    )
    dt = max(float(state.get("dt_seconds", DT_SECONDS)), 1.0)
    for temp_key, profile_key in profiles:
        surface = float(state[temp_key])
        if "roof" in profile_key or "wall" in profile_key:
            layers = max(int(settings.get("urban_wall_layers", 3)), 1)
            thickness = max(float(settings.get("urban_wall_thickness", 0.5)), 1.0e-3)
            dz = thickness / layers
            internal = float(settings.get("urban_internal_temperature_c", 24.0)) + 273.15
            profile = list(state.get(profile_key, (internal,) * (layers + 1)))
            profile = (profile[:layers + 1] + [internal] * (layers + 1))[:layers + 1]
            profile[-1] = internal
            old = profile[:]
            profile[0] = surface
            diffusivity = max(float(settings.get("urban_lambda", 1.5)), 0.0) / max(float(settings.get("urban_cs", 1000000.0)), 1.0e-6)
            factor = min(0.45, diffusivity * dt / (dz * dz))
            for index in range(1, layers):
                profile[index] = old[index] + factor * (old[index - 1] - 2.0 * old[index] + old[index + 1])
            profile[-1] = internal
            state[profile_key] = profile
            continue
        urban = "roof" in profile_key or "wall" in profile_key
        prefix = "urban_" if urban else "soil_"
        dz = max(float(settings.get(prefix + "depth_step", 0.05)), 1.0e-3)
        diffusivity = max(float(settings.get(prefix + "lambda", 1.5 if urban else 1.55)), 0.0) / max(float(settings.get(prefix + "cs", 1000000.0 if urban else 1180.0)), 1.0e-6)
        profile = list(state.get(profile_key, (surface,) * 4))
        if len(profile) < 4:
            profile.extend([profile[-1] if profile else surface] * (4 - len(profile)))
        old = profile[:]
        profile[0] = surface
        factor = min(0.45, diffusivity * dt / (dz * dz))
        for index in range(1, len(profile) - 1):
            profile[index] = old[index] + factor * (old[index - 1] - 2.0 * old[index] + old[index + 1])
        profile[-1] = old[-1]
        state[profile_key] = profile


__all__ = ["evapotranspiration", "urban_evapotranspiration", "update_thermal_history"]
