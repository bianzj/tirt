import math


RHOA = 1.2047
CP = 1004.0
MAIR = 28.96
MH20 = 18.0
SIGMASB = 5.67e-8
RAD_THRESHOLD = 5.0
# Keep the iterative energy-balance solve in a physically usable range.  These
# are numerical guards, not expected operating limits for normal cases.
TMIN_THRESHOLD = 240.0
TMAX_THRESHOLD = 380.0
WC = 0.9
MAX_TEMPERATURE_STEP = 1.5


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


def _jacobian(state, name, temp, ra, rs, emissivity, pressure, factor, latent=True):
    # Full -dR/dT of R = N - LE - H - G (Newton-Raphson as in SCOPE ebal, Van der Tol et al. 2009); factor = emitting sides.
    jacobian = (
        factor * 4.0 * emissivity * SIGMASB * temp ** 3
        + RHOA * CP / ra
        + float(state.get("dGdT_" + name, 0.0))
    )
    if latent:
        jacobian += RHOA * latent_heat(temp) * (MH20 / MAIR / pressure) * saturation_slope(temp - 273.15) / (ra + rs)
    return jacobian


def _newton_step(state, key, name, residual, jacobian, closed, method):
    if not math.isfinite(residual) or not math.isfinite(jacobian) or jacobian <= 0.0:
        return False
    relax = WC
    if method == "robust":
        # Halve this component's relaxation whenever its residual changes sign (damps oscillation).
        relax = float(state.get("_eb_relax_" + name, WC))
        previous = state.get("_eb_res_" + name)
        if previous is not None and previous * residual < 0.0 and abs(residual) > RAD_THRESHOLD:
            relax = max(0.5 * relax, 0.1)
        state["_eb_relax_" + name] = relax
        state["_eb_res_" + name] = residual
    raw_delta = relax * residual / jacobian
    if not math.isfinite(raw_delta):
        return False
    delta = max(-MAX_TEMPERATURE_STEP, min(MAX_TEMPERATURE_STEP, raw_delta))
    old = float(state[key])
    new = min(max(old + delta, TMIN_THRESHOLD), TMAX_THRESHOLD)
    state[key] = new
    # A clipped Newton step is deliberately reported as open so the caller
    # gets another iteration instead of passing a partially solved state on.
    clipped = abs(delta - raw_delta) > 1.0e-12 or abs(new - (old + delta)) > 1.0e-12
    return closed and not clipped and abs(residual) <= RAD_THRESHOLD


def _update_temperature(state, key, ra, rs, emissivity, pressure, closed, name=None, gkey=None, method="direct", factor=1.0):
    name = name or {"Tls": "leafsunlit", "Tlh": "leafshaded", "Tss": "soilsunlit", "Tsh": "soilshaded"}[key]
    ra = max(float(ra), 1.0)
    rs = max(float(rs), 0.0)
    pressure = max(float(pressure), 1.0e-6)
    temp = min(max(float(state[key]), 180.0), 380.0)
    residual = float(state.get("N_" + name, 0.0)) - float(state.get("LE_" + name, 0.0)) - float(state.get("H_" + name, 0.0))
    if gkey:
        residual -= float(state.get(gkey, 0.0))
    jacobian = _jacobian(state, name, temp, ra, rs, emissivity, pressure, factor)
    return _newton_step(state, key, name, residual, jacobian, closed, method)


def rebalance(settings, meteo, state):
    ta = meteo["air_temperature_c"] + 273.15
    leaf_method = _temperature_method(settings.get("temperature_method", "direct"))
    soil_method = _temperature_method(settings.get("soil_temperature_method", leaf_method))
    closed = True
    if settings["lai"] > 0.0:
        closed = _update_temperature(state, "Tls", state["raa_leaf"], state["rssunlit"], settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, factor=2.0)
        closed = _update_temperature(state, "Tlh", state["raa_leaf"], state["rsshaded"], settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, factor=2.0)
    closed = _update_temperature(state, "Tss", state["raa_soil"], state["rss"], settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_sunlit", method=soil_method)
    closed = _update_temperature(state, "Tsh", state["raa_soil"], state["rss"], settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_shaded", method=soil_method)
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


def _urban_temperature(state, temp_key, name, ra, emissivity, closed, method):
    # Urban LE is fixed at zero (evapo.urban_evapotranspiration), so the Jacobian has no latent term.
    ra = max(float(ra), 1.0)
    temp = min(max(float(state[temp_key]), 180.0), 380.0)
    residual = float(state.get("N_" + name, 0.0)) - float(state.get("LE_" + name, 0.0)) - float(state.get("H_" + name, 0.0)) - float(state.get("G_" + name, 0.0))
    jacobian = _jacobian(state, name, temp, ra, 0.0, emissivity, 1.0, 1.0, latent=False)
    return _newton_step(state, temp_key, name, residual, jacobian, closed, method)


def rebalance_urban(settings, meteo, state):
    urban_method = _temperature_method(settings.get("temperature_method", "direct"))
    resistance = {
        "roof": state.get("raa_roof", 1.0),
        "wall": state.get("raa_wall", 1.0),
        "street": state.get("raa_street", 1.0),
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
            state, temp_key, component + suffix, resistance[component],
            emissivity[component], closed, urban_method,
        )
    if float(settings.get("lai", 0.0)) > 0.0:
        closed = _update_temperature(
            state, "Tls", state.get("raa_leaf", state.get("raa_street", 1.0)),
            state.get("rssunlit", 625000.0), settings["leaf_emissivity"],
            meteo["pressure"], closed, name="leafsunlit", method=urban_method,
            factor=2.0,
        )
        closed = _update_temperature(
            state, "Tlh", state.get("raa_leaf", state.get("raa_street", 1.0)),
            state.get("rsshaded", 625000.0), settings["leaf_emissivity"],
            meteo["pressure"], closed, name="leafshaded", method=urban_method,
            factor=2.0,
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


__all__ = ["saturation_vapor_pressure", "saturation_slope", "latent_heat", "rebalance", "rebalance_urban"]
