import math


CP = 1004.0
RD = 287.05
MAIR = 28.96
MH20 = 18.0
SIGMASB = 5.67e-8
RAD_THRESHOLD = 5.0
# Keep the iterative energy-balance solve in a physically usable range.  These
# are numerical guards, not expected operating limits for normal cases.
TMIN_THRESHOLD = 240.0
TMAX_THRESHOLD = 380.0
WC = 0.8
MAX_TEMPERATURE_STEP = 5.0


def air_density(meteo):
    """Moist-air density [kg m-3] = p / (Rd*Tv); pressure and vapour pressure are in hPa."""
    pressure = float(meteo["pressure"]) * 100.0
    vapor = float(meteo["vapor_pressure"]) * 100.0
    virtual = (float(meteo["air_temperature_c"]) + 273.15) / (1.0 - 0.378 * vapor / pressure)
    return pressure / (RD * virtual)


def saturation_vapor_pressure(temp_c):
    return 6.107 * 10.0 ** (7.5 * temp_c / (237.3 + temp_c))


def saturation_slope(temp_c):
    es = saturation_vapor_pressure(temp_c)
    return es * 2.3026 * 7.5 * 237.3 / (237.3 + temp_c) ** 2


def latent_heat(temp_k):
    return (2.501 - 0.002361 * (temp_k - 273.15)) * 1.0e6


def _temperature_method(value):
    text = str(value if value is not None else "direct").strip().lower().replace("-", "_")
    if text in {"direct", "0", "linear"}:
        return "direct"
    if text in {"robust", "1", "stable", "bounded"}:
        return "robust"
    raise ValueError("temperature_method must be direct or robust")


def _update_temperature(state, key, ra, rs, gterm, emissivity, pressure, closed, name=None, gkey=None, method="direct", lai=0.0, soil=False, faces=1.0, rho=None):
    name = name or {"Tls": "leafsunlit", "Tlh": "leafshaded", "Tss": "soilsunlit", "Tsh": "soilshaded"}[key]
    old = float(state[key])
    ra = max(float(ra), 1.0)
    rs = max(float(rs), 0.0)
    pressure = max(float(pressure), 1.0e-6)
    rho = float(rho)
    old_for_derivative = min(max(old, 180.0), 380.0)
    residual = float(state.get("N_" + name, 0.0)) - float(state.get("LE_" + name, 0.0)) - float(state.get("H_" + name, 0.0))
    if gkey:
        residual -= float(state.get(gkey, 0.0))
    if method == "robust":
        # ebal.m: T = T + Wc*EBres / (heat-conduction + 4*sigma*T^3).
        # The MATLAB soil expression contains two air heat-conduction terms.
        effective_resistance = max(ra + rs, 1.0)
        conduction = rho * CP / effective_resistance
        if soil:
            conduction *= 2.0
        else:
            conduction /= max(float(lai) + 1.0, 1.0)
        denominator = conduction + 4.0 * SIGMASB * old_for_derivative ** 3
    else:
        denominator = (
            rho * CP / ra
            + rho * latent_heat(old_for_derivative) * (MH20 / MAIR / pressure)
            * saturation_slope(old_for_derivative - 273.15) / (ra + rs)
            + 4.0 * faces * emissivity * SIGMASB * old_for_derivative ** 3
            + max(float(gterm), 0.0)
        )
    if not math.isfinite(residual) or not math.isfinite(denominator) or denominator <= 0.0:
        return False
    raw_delta = WC * residual / denominator
    if not math.isfinite(raw_delta):
        return False
    delta = max(-MAX_TEMPERATURE_STEP, min(MAX_TEMPERATURE_STEP, raw_delta))
    new = min(max(old + delta, TMIN_THRESHOLD), TMAX_THRESHOLD)
    if not math.isfinite(new):
        return False
    state[key] = new
    # A clipped Newton step is deliberately reported as open so the caller
    # gets another iteration instead of passing a partially solved state on.
    clipped = abs(delta - raw_delta) > 1.0e-12 or abs(new - (old + delta)) > 1.0e-12
    return closed and not clipped and abs(residual) <= RAD_THRESHOLD


def rebalance(settings, meteo, state):
    ta = meteo["air_temperature_c"] + 273.15
    leaf_method = _temperature_method(settings.get("temperature_method", "direct"))
    # A soil override is independent; other components retain temperature_method.
    soil_method = _temperature_method(settings.get("soil_temperature_method", leaf_method))
    rho = air_density(meteo)
    closed = True
    if settings["lai"] > 0.0:
        closed = _update_temperature(state, "Tls", state["raa_leaf"], state["rssunlit"], 0.0, settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, lai=settings["lai"], faces=2.0 * settings["lai"], rho=rho)
        closed = _update_temperature(state, "Tlh", state["raa_leaf"], state["rsshaded"], 0.0, settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, lai=settings["lai"], faces=2.0 * settings["lai"], rho=rho)
    closed = _update_temperature(state, "Tss", state["raa_soil"], state["rss"], state["dG_term"], settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_sunlit", method=soil_method, lai=settings["lai"], soil=True, rho=rho)
    closed = _update_temperature(state, "Tsh", state["raa_soil"], state["rss"], state["dG_term"], settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_shaded", method=soil_method, lai=settings["lai"], soil=True, rho=rho)
    if soil_method == "robust":
        for key in ("Tss", "Tsh"):
            if abs(state[key] - 273.15) > 100.0:
                state[key] = ta
                closed = False
    if state["Tss"] > TMAX_THRESHOLD or state.get("Tls", ta) > TMAX_THRESHOLD:
        state["Tss"], state["Tls"] = ta + 2.0, ta + 1.0
        closed = False
    if state.get("Tlh", ta) < TMIN_THRESHOLD or state["Tsh"] < TMIN_THRESHOLD:
        state["Tlh"], state["Tsh"] = ta + 2.0, ta + 1.0
        closed = False
    for key, offset in (("Tls", 2.0), ("Tlh", 1.0), ("Tss", 2.0), ("Tsh", 1.0)):
        if not math.isfinite(state[key]):
            state[key] = ta + offset
            closed = False
    return closed


def _urban_temperature(state, temp_key, components, resistances, emissivity, settings, pressure, closed, method, rho):
    residual = 0.0
    denominator = 0.0
    for component, suffix in components:
        name = component + suffix
        ra, rs = resistances[component]
        ra = max(float(ra), 1.0)
        rs = max(float(rs), 0.0)
        residual += float(state.get("N_" + name, 0.0)) - float(state.get("LE_" + name, 0.0)) - float(state.get("H_" + name, 0.0)) - float(state.get("G_" + name, 0.0))
        old = min(max(float(state[temp_key]), 180.0), 380.0)
        if method == "robust":
            denominator += rho * CP / (ra + rs) + 4.0 * emissivity[component] * SIGMASB * old ** 3
        else:
            denominator += rho * CP / ra + rho * latent_heat(old) * (MH20 / MAIR / max(float(pressure), 1.0e-6)) * saturation_slope(old - 273.15) / (ra + rs) + 4.0 * emissivity[component] * SIGMASB * old ** 3
    if not math.isfinite(residual) or not math.isfinite(denominator) or denominator <= 0.0:
        return False
    raw_delta = WC * residual / denominator
    if not math.isfinite(raw_delta):
        return False
    delta = max(-MAX_TEMPERATURE_STEP, min(MAX_TEMPERATURE_STEP, raw_delta))
    old_temperature = float(state[temp_key])
    new_temperature = min(max(old_temperature + delta, TMIN_THRESHOLD), TMAX_THRESHOLD)
    state[temp_key] = new_temperature
    clipped = abs(delta - raw_delta) > 1.0e-12 or abs(new_temperature - (old_temperature + delta)) > 1.0e-12
    return closed and not clipped and abs(residual) <= RAD_THRESHOLD


def rebalance_urban(settings, meteo, state):
    urban_method = _temperature_method(settings.get("temperature_method", "direct"))
    rho = air_density(meteo)
    resistance = {
        "roof": (state.get("raa_roof", 1.0), state.get("rssroof", 3662.8562)),
        "wall": (state.get("raa_wall", 1.0), state.get("rsswall", 3662.8562)),
        "street": (state.get("raa_street", 1.0), state.get("rssstreet", 3662.8562)),
    }
    emissivity = {
        "roof": float(settings.get("urban_emissivity_roof", 0.950)),
        "wall": float(settings.get("urban_emissivity_wall", 0.920)),
        "street": float(settings.get("urban_emissivity_street", 0.955)),
    }
    closed = True
    for temp_key, component, suffix in (
        ("Troof_sunlit", "roof", "sunlit"),
        ("Troof_shaded", "roof", "shaded"),
        ("Twall_sunlit", "wall", "sunlit"),
        ("Twall_shaded", "wall", "shaded"),
        ("Tss", "street", "sunlit"),
        ("Tsh", "street", "shaded"),
    ):
        closed = _urban_temperature(
            state, temp_key, ((component, suffix),), resistance, emissivity,
            settings, meteo["pressure"], closed, urban_method, rho,
        )
    if float(settings.get("lai", 0.0)) > 0.0:
        closed = _update_temperature(
            state, "Tls", state.get("raa_leaf", state.get("raa_street", 1.0)),
            state.get("rssunlit", 625000.0), 0.0, settings["leaf_emissivity"],
            meteo["pressure"], closed, name="leafsunlit", method=urban_method,
            lai=settings["lai"], faces=2.0 * float(settings["lai"]), rho=rho,
        )
        closed = _update_temperature(
            state, "Tlh", state.get("raa_leaf", state.get("raa_street", 1.0)),
            state.get("rsshaded", 625000.0), 0.0, settings["leaf_emissivity"],
            meteo["pressure"], closed, name="leafshaded", method=urban_method,
            lai=settings["lai"], faces=2.0 * float(settings["lai"]), rho=rho,
        )
    else:
        state["Tls"] = 0.5 * (state["Troof_sunlit"] + state["Twall_sunlit"])
        state["Tlh"] = 0.5 * (state["Troof_shaded"] + state["Twall_shaded"])
    ta = meteo["air_temperature_c"] + 273.15
    for key, offset in (
        ("Troof_sunlit", 2.0), ("Troof_shaded", 2.0),
        ("Twall_sunlit", 2.0), ("Twall_shaded", 2.0),
        ("Tss", 1.0), ("Tsh", 1.0), ("Tls", 2.0), ("Tlh", 2.0),
    ):
        if not math.isfinite(state[key]) or state[key] < TMIN_THRESHOLD or state[key] > TMAX_THRESHOLD:
            state[key] = ta + offset
            closed = False
    return closed


__all__ = ["air_density", "saturation_vapor_pressure", "saturation_slope", "latent_heat", "rebalance", "rebalance_urban"]
