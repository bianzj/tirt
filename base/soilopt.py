import math
from functools import lru_cache

import numpy as np

from .leafopt import calctav

#: SCOPE's BSM soil-moisture capacity (reference value 0.25) in m3 m-3.
BSM_WATER_CAPACITY = 0.25
#: BSM residual (dry) and maximum (saturated) soil water, in m3 m-3.
BSM_WATER_RESIDUAL = 0.05
BSM_WATER_MAX = 0.55


def volumetric_water(value):
    """Normalize soil water to the model's m3 m-3 unit.

    TIRTEB uses volume fractions everywhere, with 0.25 as the BSM reference
    capacity. Legacy percent inputs are still accepted: a value above 1 is
    divided by 100.
    """
    fraction = float(value)
    return fraction / 100.0 if fraction > 1.0 else fraction


def bsm_water_mu(soil_moisture):
    """Poisson mu of the BSM water-film model, computed in m3 m-3."""
    water = min(max(volumetric_water(soil_moisture), BSM_WATER_RESIDUAL), BSM_WATER_MAX)
    return (water - BSM_WATER_RESIDUAL) / BSM_WATER_CAPACITY


def _soil_values(soil_optical):
    if isinstance(soil_optical, dict):
        return soil_optical["BSMBrightness"], soil_optical["BSMlat"], soil_optical["BSMlon"]
    return soil_optical[0], soil_optical[1], soil_optical[2]


def bsm(coeff, soil_optical, soil_moisture):
    """Calculate wet-soil reflectance as SCOPE's BSM; moisture is m3 m-3."""
    brightness, lat, lon = _soil_values(soil_optical)
    theta = soil_moisture[5] if isinstance(soil_moisture, (list, tuple)) else soil_moisture
    # SCOPE's soilwat uses mu = (SMp - 5) / SMC with SMp in volume percent.
    # TIRTEB keeps both theta and the 0.25 capacity in m3 m-3, so the
    # equivalent expression is mu = (theta - 0.05) / 0.25.
    mu = bsm_water_mu(theta)
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
        # SCOPE: rbac = 1 - (1-rdry) * (rdry * tav(90, 2/nw) / tav(90, 2) + 1-rdry)
        rbac = 1.0 - (1.0 - rdry) * (
            rdry * calctav(90.0, 2.0 / nw) / calctav(90.0, 2.0) + 1.0 - rdry)
        p = 1.0 - tav90_nw / (nw * nw)
        # SCOPE: Rw = 1 - tav(40, nw); calctav returns transmittance, not reflectance.
        rw = 1.0 - calctav(40.0, nw)
        radd = (1.0 - rw) * (1.0 - p) * rbac / (1.0 - p * rbac)
        fdry = math.exp(-mu)
        fmul = (math.exp(tw * mu) - 1.0) * fdry
        reflectance.append(rdry * fdry + rw * (1.0 - fdry) + radd * fmul)
    return reflectance


def soil_emissivity(reflectance):
    return 1.0 - reflectance


def hapke_emissivity(view_zenith, K=1.0, albedo=0.05):
    """TIRT's simplified Hapke/Thapke directional soil emissivity."""
    angle = np.asarray(view_zenith, dtype=float)
    K = max(float(K), 1.0e-8)
    albedo = np.clip(np.asarray(albedo, dtype=float), 0.0, 0.999999)
    base_emissivity = np.sqrt(1.0 - albedo)
    term = 2.0 * np.cos(np.deg2rad(angle)) / K
    denominator = 1.0 + term * base_emissivity
    result = base_emissivity * (1.0 + term) / np.maximum(denominator, 1.0e-12)
    return np.clip(result, 0.0, 1.0)


def hapke_reflectance(view_zenith, K=1.0, albedo=0.05):
    """Simple Hapke/Thapke directional reflectance from soil albedo."""
    return float(1.0 - float(hapke_emissivity(view_zenith, K, albedo)))


_HAPKE_MU, _HAPKE_WEIGHTS = np.polynomial.legendre.leggauss(16)
_HAPKE_MU, _HAPKE_WEIGHTS = 0.5 * (_HAPKE_MU + 1.0), 0.5 * _HAPKE_WEIGHTS


def hapke_hemispherical_emissivity(K=1.0, albedo=0.05):
    """Cosine-weighted hemispherical mean of hapke_emissivity."""
    zenith = np.degrees(np.arccos(_HAPKE_MU))
    return float(2.0 * np.sum(_HAPKE_WEIGHTS * _HAPKE_MU * hapke_emissivity(zenith, K, albedo)))


@lru_cache(maxsize=256)
def hapke_single_scattering_albedo(hemispherical_emissivity, K=1.0):
    """Single-scattering albedo whose Hapke hemispherical emissivity equals the input."""
    target = float(hemispherical_emissivity)
    if not 0.0 < target <= 1.0:
        raise ValueError("Hapke soil emissivity must be in (0, 1]")
    if target >= hapke_hemispherical_emissivity(K, 0.0):
        return 0.0
    from scipy.optimize import brentq
    return float(brentq(lambda w: hapke_hemispherical_emissivity(K, w) - target, 0.0, 0.999999, xtol=1.0e-12))


def hapke_directional_emissivity(view_zenith, K=1.0, hemispherical_emissivity=0.955):
    """Directional Hapke soil emissivity that preserves the given hemispherical emissivity."""
    albedo = hapke_single_scattering_albedo(float(hemispherical_emissivity), float(K))
    return hapke_emissivity(view_zenith, K, albedo)


def soil_tir_emissivity(
    soil_moisture,
    enabled=0,
    dry_moisture=0.05,
    wet_moisture=0.45,
    dry_emissivity=0.955,
    wet_emissivity=0.99,
):
    """Linear soil-water correction for TIR emissivity; water is m3 m-3."""
    dry_moisture = volumetric_water(dry_moisture)
    wet_moisture = max(volumetric_water(wet_moisture), dry_moisture + 1.0e-8)
    if int(enabled) != 1:
        return float(np.clip(dry_emissivity, 0.0, 1.0))
    fraction = (volumetric_water(soil_moisture) - dry_moisture) / (wet_moisture - dry_moisture)
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
    """Directional Hapke soil emissivity with an optional water correction."""
    if int(water_enabled) == 1:
        water_emissivity = soil_tir_emissivity(
            soil_moisture,
            enabled=1,
            dry_moisture=dry_moisture,
            wet_moisture=wet_moisture,
            dry_emissivity=dry_emissivity,
            wet_emissivity=wet_emissivity,
        )
        # Kirchhoff's relation supplies the current TIR albedo independently
        # of the optical BSM spectrum.
        albedo = 1.0 - water_emissivity
    return hapke_emissivity(view_zenith, K, albedo)


__all__ = [
    "BSM_WATER_CAPACITY", "BSM_WATER_RESIDUAL", "BSM_WATER_MAX",
    "volumetric_water", "bsm_water_mu", "bsm",
    "soil_emissivity", "hapke_emissivity", "hapke_reflectance",
    "hapke_hemispherical_emissivity", "hapke_single_scattering_albedo", "hapke_directional_emissivity",
    "soil_tir_emissivity", "hapke_thermal_emissivity",
]
