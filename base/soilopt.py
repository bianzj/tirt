import math
from functools import lru_cache

import numpy as np

from .leafopt import calctav


def _soil_values(soil_optical):
    if isinstance(soil_optical, dict):
        return soil_optical["BSMBrightness"], soil_optical["BSMlat"], soil_optical["BSMlon"]
    return soil_optical[0], soil_optical[1], soil_optical[2]


def bsm(coeff, soil_optical, soil_moisture):
    """Calculate wet-soil reflectance as SoilOpt::bsm."""
    brightness, lat, lon = _soil_values(soil_optical)
    smc = soil_moisture[5] if isinstance(soil_moisture, (list, tuple)) else soil_moisture
    smcp = 0.25
    film = 0.015
    rd = 3.1415926 / 180.0
    f1 = brightness * math.sin(lat * rd)
    f2 = brightness * math.cos(lat * rd) * math.sin(lon * rd)
    f3 = brightness * math.cos(lat * rd) * math.cos(lon * rd)

    reflectance = []
    for row in coeff[:2001]:
        gsv1, gsv2, gsv3 = row[14], row[15], row[16]
        kw, nw = row[5], row[17]
        rdry = gsv1 * f1 + gsv2 * f2 + gsv3 * f3
        tw = math.exp(-kw * film)
        tav90_nw = calctav(90.0, nw)
        rbac = 1.0 - (1.0 - rdry) * (rdry * calctav(90.0, 2.0) / tav90_nw + 1.0 - rdry)
        p = 1.0 - tav90_nw / (nw * nw)
        # Fresnel reflectance of the water film surface (SCOPE BSM: Rw = 1 - tav(40, nw)).
        rw = 1.0 - calctav(40.0, nw)
        radd = (1.0 - rw) * (1.0 - p) * rbac / (1.0 - p * rbac)
        mu = (smc - 0.05) / smcp
        fdry = math.exp(-mu)
        fmul = (math.exp(tw * mu) - 1.0) * fdry
        reflectance.append(rdry * fdry + rw * (1.0 - fdry) + radd * fmul)
    return reflectance


def soil_emissivity(reflectance):
    return 1.0 - reflectance


def hapke_emissivity(view_zenith, K=1.0, albedo=0.05):
    """TIRT's simplified Hapke/Thapke directional soil emissivity; ``albedo`` is the single-scattering albedo w."""
    angle = np.asarray(view_zenith, dtype=float)
    K = max(float(K), 1.0e-8)
    albedo = np.clip(np.asarray(albedo, dtype=float), 0.0, 0.999999)
    base_emissivity = np.sqrt(1.0 - albedo)
    term = 2.0 * np.cos(np.deg2rad(angle)) / K
    denominator = 1.0 + term * base_emissivity
    result = base_emissivity * (1.0 + term) / np.maximum(denominator, 1.0e-12)
    return np.clip(result, 0.0, 1.0)


def _hapke_hemispherical_emissivity(gamma, K):
    """Closed-form 2*integral(mu * eps(mu)) over [0, 1] for eps = gamma (1 + a mu) / (1 + a gamma mu), a = 2/K."""
    b = 2.0 * gamma / K
    integral = 0.5 - b / 3.0 if b < 1.0e-6 else 1.0 / b - math.log1p(b) / (b * b)
    return 1.0 + 2.0 * (gamma - 1.0) * integral


@lru_cache(maxsize=4096)
def single_scattering_albedo(hemispherical_emissivity, K=1.0):
    """Hapke single-scattering albedo w whose directional emissivity integrates to the given hemispherical emissivity."""
    target = min(max(float(hemispherical_emissivity), 0.0), 1.0)
    K = max(float(K), 1.0e-8)
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if _hapke_hemispherical_emissivity(middle, K) < target:
            low = middle
        else:
            high = middle
    gamma = 0.5 * (low + high)
    return 1.0 - gamma * gamma


def hapke_reflectance(view_zenith, K=1.0, albedo=0.05):
    """Simple Hapke/Thapke directional reflectance from soil albedo."""
    return float(1.0 - float(hapke_emissivity(view_zenith, K, albedo)))


def soil_tir_emissivity(
    soil_moisture,
    enabled=0,
    dry_moisture=0.05,
    wet_moisture=0.45,
    dry_emissivity=0.955,
    wet_emissivity=0.99,
):
    """Linear soil-water correction for thermal-infrared emissivity."""
    dry_moisture = float(dry_moisture)
    wet_moisture = max(float(wet_moisture), dry_moisture + 1.0e-8)
    if int(enabled) != 1:
        return float(np.clip(dry_emissivity, 0.0, 1.0))
    fraction = (float(soil_moisture) - dry_moisture) / (wet_moisture - dry_moisture)
    fraction = min(max(fraction, 0.0), 1.0)
    return float(np.clip(
        float(dry_emissivity) + fraction * (float(wet_emissivity) - float(dry_emissivity)),
        0.0,
        1.0,
    ))


def hapke_thermal_emissivity(
    view_zenith,
    K=1.0,
    albedo=0.05,
    soil_moisture=0.0,
    water_enabled=0,
    dry_moisture=0.05,
    wet_moisture=0.45,
    dry_emissivity=0.955,
    wet_emissivity=0.99,
):
    """Directional Hapke soil emissivity; ``albedo`` is the hemispherical TIR reflectance (1 - emissivity)."""
    hemispherical_emissivity = 1.0 - float(albedo)
    if int(water_enabled) == 1:
        hemispherical_emissivity = soil_tir_emissivity(
            soil_moisture,
            enabled=1,
            dry_moisture=dry_moisture,
            wet_moisture=wet_moisture,
            dry_emissivity=dry_emissivity,
            wet_emissivity=wet_emissivity,
        )
    return hapke_emissivity(view_zenith, K, single_scattering_albedo(hemispherical_emissivity, float(K)))


__all__ = [
    "bsm", "soil_emissivity", "hapke_emissivity", "hapke_reflectance",
    "soil_tir_emissivity", "hapke_thermal_emissivity",
]
