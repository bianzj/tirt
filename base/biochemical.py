import math

from .eb import saturation_vapor_pressure


RHOA = 1.2047
MAIR = 28.96
RSS = 50.0


def _root(a, b, c, smaller=True):
    if a == 0.0:
        return -c / b if b else float("nan")
    discriminant = max(b * b - 4.0 * a * c, 0.0)
    roots = ((-b - math.sqrt(discriminant)) / (2.0 * a), (-b + math.sqrt(discriminant)) / (2.0 * a))
    return min(roots) if smaller else max(roots)


def farquhar_resistance(settings, meteo, state, netrad):
    """Biochemical::suresist -> biochemical_Farquhar, in C++ order and units."""
    state["rss"] = max(float(settings.get("soil_rss", RSS)), 0.0)
    if settings["surface_type"] == 16 or settings["lai"] <= 0.0:
        return state

    leaf_type = float(settings["leaf_type"])
    vcmax25 = float(settings["vcmax"])
    ballberry = float(settings["ballberry_m"])
    rdparam = float(settings["rdparam"])
    tcor = float(settings.get("tcor", 1.0))
    stress = float(settings.get("stressfactor", 1.0))
    tparam = list(settings.get("tparam", (0.2, 0.3, 281.0, 313.0, 328.0)))
    while len(tparam) < 5:
        tparam.append(0.0)
    tref = 25.0 + 273.15
    p = float(meteo["pressure"])
    cs0 = float(state.get("cs", settings.get("ca", 380.0)))
    eb = float(state.get("es", settings.get("ea", meteo["vapor_pressure"])))
    oxygen = float(settings.get("oa", 209.0))
    q = float(netrad.get("directPrad_leaf", 0.0) + netrad.get("diffusePrad_leaf", 0.0))
    effcon = 1.0 / 6.0 if leaf_type >= 3.5 else 1.0 / 5.0

    # The source has this historical assignment for both components.
    for issunlit in (0, 1):
        tc = state["Tlh"]
        cs = cs0 * (1.0e-6 * (p * 1.0e-3))
        o = (oxygen * 1.0e-3) * (p * 1.0e-3)
        if leaf_type > 3.5:
            o = 0.0
        kc25 = 350.0e-6
        ko25 = 450.0e-3
        kpep25 = (vcmax25 / 56.0) * 1.0e6
        atheta = 0.8
        kf = 0.05
        kd = max(0.8738, 0.0301 * (tc - 273.15) + 0.0773)
        kp = 4.0
        qt = 0.1 * (tc - tref) * tcor
        th = 1.0 + tcor * math.exp(tparam[1] * (tc - tparam[3]))
        tl = 1.0 + tcor * math.exp(tparam[0] * (tparam[2] - tc))
        kc = kc25 * math.exp(math.log(2.1) * qt)
        ko = ko25 * math.exp(math.log(1.2) * qt)
        kpep = kpep25 * math.exp(math.log(1.8) * qt)
        if leaf_type < 3.5:
            vcmax = vcmax25 * math.exp(math.log(2.1) * qt) / th * stress
            mm = kc * (1.0 + o / ko)
            vs_c3 = vcmax25 / 2.0 * math.exp(math.log(1.8) * qt)
            min_ci = 0.3
        else:
            vcmax = vcmax25 * math.exp(math.log(2.1) * qt) / (th * tl) * stress
            mm = 0.0
            vs_c3 = 0.0
            min_ci = 0.1
        specificity = 2600.0 * math.exp(math.log(0.75) * qt)
        rd = rdparam * vcmax25 * math.exp(math.log(1.8) * qt) / (1.0 + math.exp(1.3 * (tc - tparam[4])))
        po0 = kp / (kf + kd + kp)
        je = 0.5 * po0 * q
        gamma = 0.5 * o / specificity
        rh = min(1.0, eb / saturation_vapor_pressure(tc - 273.15))
        ci = max(min_ci * cs, cs * (1.0 - 1.6 / (ballberry * rh)))
        if leaf_type < 3.5:
            vs = vs_c3
            vc = vcmax * (ci - gamma) / (mm + ci)
            co2_per_electron = (ci - gamma) / (ci + 2.0 * gamma) * effcon
            ve = je * co2_per_electron
        else:
            vc = vcmax
            vs = kpep * ci
            co2_per_electron = effcon
            ve = je * co2_per_electron
        v = _root(atheta, -(vc + ve), vc * ve, True)
        ag = _root(0.98, -(v + vs), v * vs, True)
        assimilation = ag - rd
        gs = 1.6 * assimilation * (1.0e-6 * (p * 1.0e-3)) / (cs - ci)
        rcw = (RHOA / (MAIR * 1.0e-3)) / gs if gs else float("nan")
        if assimilation < 0.0 or not math.isfinite(rcw):
            rcw = 0.625 * 1.0e6
        if issunlit == 0:
            state["rsshaded"] = rcw
        else:
            state["rssunlit"] = rcw
        if math.isfinite(ci):
            state["ci"] = ci / (1.0e-6 * (p * 1.0e-3))
    return state


def suresist(settings, meteo, state, netrad):
    method = str(settings.get("stomatal_method", "farquhar")).strip().lower().replace("-", "_")
    if method in {"farquhar", "mechanistic", "photosynthesis", "0"}:
        return farquhar_resistance(settings, meteo, state, netrad)
    if method in {"ballberry", "ball_berry", "empirical", "light_response", "1"}:
        return ballberry_resistance(settings, meteo, state, netrad)
    raise ValueError("stomatal_method must be farquhar or ballberry")


def ballberry_resistance(settings, meteo, state, netrad):
    """Empirical Ball-Berry/light-response stomatal resistance in s m-1."""
    state["rss"] = max(float(settings.get("soil_rss", RSS)), 0.0)
    if settings["lai"] <= 0.0:
        state["rssunlit"] = state["rsshaded"] = 625000.0
        return state
    pressure = max(float(meteo["pressure"]), 1.0e-6)
    air_temp = float(meteo["air_temperature_c"]) + 273.15
    ambient_es = saturation_vapor_pressure(air_temp - 273.15)
    rh = min(max(float(meteo["vapor_pressure"]) / max(ambient_es, 1.0e-6), 0.05), 1.2)
    ca = max(float(state.get("cs", settings.get("ca", 380.0))), 1.0)
    par = max(float(netrad.get("directPrad_leaf", 0.0)) + float(netrad.get("diffusePrad_leaf", 0.0)), 0.0)
    light = par / (par + 100.0)
    temperature = math.exp(-((float(state["Tlh"]) - 298.15) / 18.0) ** 2)
    vcmax = max(float(settings.get("vcmax", 80.0)), 0.0)
    stress = min(max(float(settings.get("stressfactor", 1.0)), 0.0), 1.0)
    assimilation = vcmax * light * temperature * stress
    slope = max(float(settings.get("ballberry_m", 9.0)), 0.0)
    g0 = max(float(settings.get("stomatal_g0", 0.01)), 1.0e-5)
    gs = g0 + (1.0 + slope * rh) * assimilation / ca
    gs = min(max(gs, 1.0e-5), 2.0)
    resistance = (RHOA / (MAIR * 1.0e-3)) / gs
    state["rsshaded"] = resistance
    state["rssunlit"] = resistance
    state["ci"] = ca
    return state


def suresist_urban(settings, meteo, state, netrad):
    # Aerodynamic resistance is supplied by the urban canopy branch. Only the
    # leaf surface resistance changes with vegetation physiology.
    resistance = max(float(settings.get("urban_rss", 3662.8562)), 0.0)
    state["rssroof"] = state["rsswall"] = state["rssstreet"] = resistance
    state["rss"] = resistance
    if float(settings.get("lai", 0.0)) <= 0.0:
        state["rssunlit"] = state["rsshaded"] = resistance
        return state
    method = str(settings.get("stomatal_method", "farquhar")).strip().lower()
    if method in {"ballberry", "ball_berry", "empirical", "light_response", "1"}:
        ballberry_resistance(settings, meteo, state, netrad)
    else:
        farquhar_resistance(settings, meteo, state, netrad)
    state["rss"] = resistance
    return state


__all__ = ["farquhar_resistance", "ballberry_resistance", "suresist", "suresist_urban"]
