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


def heat_derivative(settings, state, method_key, default_method, urban=False, envelope=False):
    """Upper estimate of dG/dT [W m-2 K-1] for the configured ground/envelope heat method.

    Used only in relaxation denominators; overestimating it damps but never
    destabilises the update. The direct method (G = fraction*N) contributes a
    small negative slope and is taken as 0.
    """
    from .evapo import _method
    from .heat_conduction import column_layers, thermal_properties, time_step

    method = _method(settings.get(method_key, default_method))
    if method == "direct":
        return 0.0
    conductivity, capacity = thermal_properties(settings, urban=urban)
    dt = time_step(settings, state)
    if method == "surface":
        return 2.0 * math.sqrt(conductivity * capacity) / math.sqrt(math.pi * dt)
    if envelope:
        layers = max(int(settings.get("urban_wall_layers", 3)), 1)
        dz = max(float(settings.get("urban_wall_thickness", 0.5)), 1.0e-3) / layers
        return conductivity / dz
    dz = column_layers()[0]
    return conductivity / dz + capacity * 0.5 * dz / dt


def _jacobian(rho, ra, rs, pressure, temp, emissivity, faces, gterm, latent=True):
    """d(N - LE - H - G)/dT magnitude for one component (per unit ground area)."""
    value = rho * CP / ra + 4.0 * faces * emissivity * SIGMASB * temp ** 3 + max(float(gterm), 0.0)
    if latent:
        value += (rho * latent_heat(temp) * (MH20 / MAIR / pressure)
                  * saturation_slope(temp - 273.15) / (ra + rs))
    return value


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
    # Relaxed Newton step on the complete per-ground-area Jacobian for both
    # methods. The former "robust" ebal.m-style denominator (rho*cp/(ra+rs)
    # without the 2*LAI emission, dLE/dT and dG/dT terms) under-estimated the
    # slope several fold, so Wc*J/denominator exceeded 2 and the update
    # oscillated. Both methods keep the 5 K step limit and temperature bounds.
    denominator = _jacobian(rho, ra, rs, pressure, old_for_derivative, emissivity,
                            faces, gterm, latent=True)
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
    # dG/dT of the configured soil heat method (surface memory, column or direct).
    soil_gterm = heat_derivative(settings, state, "soil_heat_method", "direct")
    if settings["lai"] > 0.0:
        closed = _update_temperature(state, "Tls", state["raa_leaf"], state["rssunlit"], 0.0, settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, lai=settings["lai"], faces=2.0 * settings["lai"], rho=rho)
        closed = _update_temperature(state, "Tlh", state["raa_leaf"], state["rsshaded"], 0.0, settings["leaf_emissivity"], meteo["pressure"], closed, method=leaf_method, lai=settings["lai"], faces=2.0 * settings["lai"], rho=rho)
    closed = _update_temperature(state, "Tss", state["raa_soil"], state["rss"], soil_gterm, settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_sunlit", method=soil_method, lai=settings["lai"], soil=True, rho=rho)
    closed = _update_temperature(state, "Tsh", state["raa_soil"], state["rss"], soil_gterm, settings["soil_emissivity"], meteo["pressure"], closed, gkey="G_shaded", method=soil_method, lai=settings["lai"], soil=True, rho=rho)
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
        # H uses ra alone (evapo._component_flux); urban LE is forced to 0, so it has
        # no slope; G is the roof/wall envelope or the street (soil-material) flux.
        if component in {"roof", "wall"}:
            gterm = heat_derivative(settings, state, "urban_heat_method", "depth", urban=True, envelope=True)
        else:
            gterm = heat_derivative(settings, state, "soil_heat_method", "surface")
        denominator += _jacobian(rho, ra, rs, max(float(pressure), 1.0e-6), old,
                                 emissivity[component], 1.0, gterm, latent=False)
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
    # urban_temperature_method is parsed by read_defined (falls back to temperature_method).
    urban_method = _temperature_method(
        settings.get("urban_temperature_method") or settings.get("temperature_method", "direct"))
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
