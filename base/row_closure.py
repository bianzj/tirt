"""Row closures on one footing: legacy, conditioned and the analytical reference.

legacy      : CI transmittance, CI layer hotspot, equivalent-homogeneous multiple scattering.
conditioned : legacy + hotspot correlation compensation; multiple scattering as legacy.
reference   : parametric geometry (exact cross-section integrals) + cross-section radiosity.

The legacy row gap P(theta, phi) gives the clumping index Omega = -mu ln P / (G LAI).
The layer hotspot puts Omega into the homogeneous hotspot equations, i.e. it treats the
sun and view paths as independent equivalent layers. The compensation restores their
row-scale correlation: each direction becomes a two-phase surrogate (clear fraction c,
covered depth tau_c) with the same gap and mean depth G LAI/mu; the sun and view covered
strips overlap by the exact periodic strip overlap at their across-row displacement, and
inside the shared row the layer (leaf-scale) overlap sqrt(tau_s tau_v) w applies. Visible
leaves lie inside a row, so their sun path is always covered. The compensation vanishes
when Omega = 1. Equivalent-homogeneous scattering uses the rows' view gap, hemispherical
gap and the sunlit-leaf fraction (1 - P_s) mu_s / (G LAI) implied by the CI transmittance.
No fitted coefficient enters.
"""
import functools
import math

import numpy as np

from .gap import gap_probability_row_analytical

G = 0.5
CLOSURES = ("legacy", "conditioned", "reference")
_DEPTH, _DEPTH_W = np.polynomial.legendre.leggauss(8)
_DEPTH, _DEPTH_W = 0.5 * (_DEPTH + 1.0), 0.5 * _DEPTH_W
_STEP = 1.0e-4
_MU, _MU_W = np.polynomial.legendre.leggauss(12)
_MU, _MU_W = 0.5 * (_MU + 1.0), 0.5 * _MU_W
_AZIMUTHS = (np.arange(12) + 0.5) * 7.5


def _check(kind):
    if kind not in CLOSURES:
        raise ValueError(f"row_closure must be one of {CLOSURES}")


def _zenith(angle):
    return min(abs(float(angle)), 89.999)


def _gap(lai, width, blank, height, zenith, relative):
    return float(gap_probability_row_analytical(lai, width, blank, height, np.array([zenith]),
                                                np.array([abs(relative) % 180.0]))[0])


@functools.lru_cache(maxsize=64)
def hemispherical_gap(lai, width, blank, height):
    """Cosine-weighted hemispherical gap of the rows."""
    zenith = np.degrees(np.arccos(_MU))
    gaps = np.array([[_gap(lai, width, blank, height, z, a) for a in _AZIMUTHS] for z in zenith]).mean(axis=1)
    return float(2.0 * np.sum(_MU_W * _MU * gaps))


def _phases(gap, depth):
    """Clear fraction c and covered depth t with (1 - c) t = depth and c + (1 - c) exp(-t) = gap."""
    if depth < 1.0e-12:
        return 1.0, 0.0
    gap = min(max(gap, math.exp(-depth)), 1.0 - 1.0e-12)
    if gap - math.exp(-depth) < 1.0e-12:
        return 0.0, depth
    low, high = depth, 2.0 * depth
    while depth / high * -math.expm1(-high) > 1.0 - gap:
        high *= 2.0
    for _ in range(80):
        middle = 0.5 * (low + high)
        if depth / middle * -math.expm1(-middle) > 1.0 - gap:
            low = middle
        else:
            high = middle
    covered = 0.5 * (low + high)
    return 1.0 - depth / covered, covered


def _strip_overlap(a, b, offset, period):
    """Overlap length of two periodic strips of lengths a, b whose centres are offset."""
    d = abs(offset) % period
    d = min(d, period - d)
    return min(a, b, max(0.0, 0.5 * (a + b) - d) + max(0.0, 0.5 * (a + b) - (period - d)))


def _layer(lai, width, blank, height, sun, view, correlation, fraction):
    """Compensated joint gap, view gap and covered sun path of the top `fraction` of the rows."""
    (zs, rs), (zv, rv) = sun, view
    layer_lai, layer_height = lai * fraction, height * fraction
    gs, gv = _gap(layer_lai, width, blank, layer_height, zs, rs), _gap(layer_lai, width, blank, layer_height, zv, rv)
    cs, ts = _phases(gs, G * layer_lai / math.cos(math.radians(zs)))
    cv, tv = _phases(gv, G * layer_lai / math.cos(math.radians(zv)))
    period = width + blank
    shift_s = layer_height * math.tan(math.radians(zs)) * math.sin(math.radians(rs))
    shift_v = layer_height * math.tan(math.radians(zv)) * math.sin(math.radians(rv))
    a, b = 1.0 - cs, 1.0 - cv
    both = _strip_overlap(a * period, b * period, 0.5 * (shift_s - shift_v), period) / period
    leaf = 1.0 if correlation == 0.0 else -math.expm1(-correlation * fraction) / (correlation * fraction)
    shared = min(ts, tv, math.sqrt(ts * tv) * leaf)
    joint = (1.0 - a - b + both) + (a - both) * math.exp(-ts) + (b - both) * math.exp(-tv) + \
        both * math.exp(-(ts + tv - shared))
    return min(joint, gs, gv), gv, math.exp(-(ts - shared))


def leaf_correlation(hspot, vza, vaa, sza, saa):
    """Full-depth leaf-scale hotspot correlation w = (1 - exp(-a)) / a of the layer hotspot."""
    zs, zv = _zenith(sza), _zenith(vza)
    tan_s, tan_v = math.tan(math.radians(zs)), math.tan(math.radians(zv))
    separation = math.sqrt(max(0.0, (tan_v - tan_s) ** 2 +
                               4.0 * tan_v * tan_s * math.sin(math.radians(float(vaa) - float(saa)) / 2.0) ** 2))
    if separation == 0.0:
        return 1.0
    if hspot <= 0.0:
        return 0.0
    a = separation * 2.0 / (G / math.cos(math.radians(zs)) + G / math.cos(math.radians(zv))) / hspot
    return -math.expm1(-a) / a


def _conditioned_setup(hspot, vza, vaa, sza, saa, row_azimuth):
    zs, zv = _zenith(sza), _zenith(vza)
    sun, view = (zs, float(saa) - row_azimuth), (zv, float(vaa) - row_azimuth)
    tan_s, tan_v = math.tan(math.radians(zs)), math.tan(math.radians(zv))
    separation = math.sqrt(max(0.0, (tan_v - tan_s) ** 2 +
                               4.0 * tan_v * tan_s * math.sin(math.radians(float(vaa) - float(saa)) / 2.0) ** 2))
    correction = 2.0 / (G / math.cos(math.radians(zs)) + G / math.cos(math.radians(zv)))
    correlation = (0.0 if separation == 0.0 else math.inf) if hspot <= 0.0 else separation * correction / hspot
    return sun, view, correlation


def conditioned_soil(lai, hspot, width, blank, height, vza, vaa, sza, saa, row_azimuth):
    """Sunlit and shaded soil of the compensated (dual-scale) hotspot."""
    if lai <= 0.0:
        return np.array([1.0, 0.0])
    soil_sun, soil, _ = _layer(lai, width, blank, height, *_conditioned_setup(hspot, vza, vaa, sza, saa, row_azimuth), 1.0)
    soil = min(max(soil, 0.0), 1.0)
    soil_sun = min(max(soil_sun, 0.0), soil)
    return np.array([soil_sun, soil - soil_sun])


def _conditioned_fractions(lai, hspot, width, blank, height, vza, vaa, sza, saa, row_azimuth):
    sun, view, correlation = _conditioned_setup(hspot, vza, vaa, sza, saa, row_azimuth)
    soil_sun, soil, _ = _layer(lai, width, blank, height, sun, view, correlation, 1.0)
    leaf_sun = 0.0
    for depth, weight in zip(_DEPTH, _DEPTH_W):
        lower, upper = max(depth - _STEP, 1.0e-6), min(depth + _STEP, 1.0)
        visible = (_layer(lai, width, blank, height, sun, view, correlation, lower)[1] -
                   _layer(lai, width, blank, height, sun, view, correlation, upper)[1]) / (upper - lower)
        leaf_sun += weight * visible * _layer(lai, width, blank, height, sun, view, correlation, depth)[2]
    return np.array([max(soil_sun, 0.0), soil, leaf_sun])


def _legacy_fractions(lai, hspot, width, blank, height, vza, vaa, sza, saa, row_azimuth):
    from .hotspot import hotspot_analytical, hotspot_layer
    zs, zv = _zenith(sza), _zenith(vza)
    gv = max(_gap(lai, width, blank, height, zv, float(vaa) - row_azimuth), 1.0e-12)
    gs = max(_gap(lai, width, blank, height, zs, float(saa) - row_azimuth), 1.0e-12)
    omega_v = -math.log(gv) * math.cos(math.radians(zv)) / (G * lai)
    omega_s = -math.log(gs) * math.cos(math.radians(zs)) / (G * lai)
    relative = np.array([abs((float(vaa) - float(saa) + 180.0) % 360.0 - 180.0)])
    angles = (np.array([zv]), np.array([zs]), relative)
    soil_sun = float(np.ravel(hotspot_layer(lai, max(hspot, 1.0e-12), *angles, omega_s, omega_v))[0])
    leaf_sun = float(np.ravel(hotspot_analytical(lai, max(hspot, 1.0e-12), *angles, omega_s, omega_v))[0])
    return np.array([soil_sun, gv, leaf_sun])


def fractions(kind, lai, hspot, width, blank, height, vza, vaa, sza, saa, row_azimuth):
    """Sunlit soil, shaded soil, sunlit leaf, shaded leaf for one direction (degrees)."""
    _check(kind)
    if lai <= 0.0:
        return np.array([1.0, 0.0, 0.0, 0.0])
    if kind == "reference":
        from .row_parametric import row_fractions
        return row_fractions(lai, max(hspot, 1.0e-12), width, blank, height, vza, vaa, sza, saa, row_azimuth)
    build = _legacy_fractions if kind == "legacy" else _conditioned_fractions
    soil_sun, soil, leaf_sun = build(lai, hspot, width, blank, height, vza, vaa, sza, saa, row_azimuth)
    soil = min(max(soil, 0.0), 1.0)
    soil_sun = min(max(soil_sun, 0.0), soil)
    leaf_sun = min(max(leaf_sun, 0.0), 1.0 - soil)
    return np.array([soil_sun, soil - soil_sun, leaf_sun, 1.0 - soil - leaf_sun])


def sun_gap(lai, width, blank, height, sza, saa, row_azimuth):
    return _gap(lai, width, blank, height, _zenith(sza), float(saa) - row_azimuth)


def thermal_weights(hotspot, scatter, lai, hspot, width, blank, height, emis_leaf, emis_soil,
                    vza, vaa, sza, saa, row_azimuth):
    """Effective component emissivities [soil sun, soil shade, leaf sun, leaf shade] (black sky)."""
    _check(hotspot)
    _check(scatter)
    shape = (float(lai), float(width), float(blank), float(height))
    direct = fractions(hotspot, *shape[:1], hspot, *shape[1:], vza, vaa, sza, saa, row_azimuth) * \
        np.array([emis_soil, emis_soil, emis_leaf, emis_leaf])
    if lai <= 0.0:
        return direct
    if scatter == "reference":
        from .row_radiosity import scatter_emissivity
        seen = direct[0] / emis_soil + direct[1] / emis_soil
        return direct + scatter_emissivity(*shape, emis_leaf, emis_soil, vza, vaa, sza, saa, row_azimuth, soil_view=seen)
    from .emissivity import _equivalent_hemispherical_gap, emissivity_scatter_analytical
    zv, zs = _zenith(vza), _zenith(sza)
    view_gap = _gap(*shape, zv, float(vaa) - row_azimuth)
    hemi_gap = float(_equivalent_hemispherical_gap(hemispherical_gap(*shape)))
    sunlit = (1.0 - sun_gap(*shape, zs, saa, row_azimuth)) * math.cos(math.radians(zs)) / (G * lai)
    extra = np.asarray(emissivity_scatter_analytical(lai, emis_leaf, emis_soil, np.array([zv]), np.array([zs]),
                                                     view_gap, hemi_gap, sunlit)[0], float).reshape(4)
    budget = max(0.0, 1.0 - float(direct.sum()))
    if extra.sum() > budget:
        extra *= budget / extra.sum()
    return direct + extra


def optical_reflectance(hotspot, scatter, settings, sza, saa, vza, vaa, rho, tau, soil_rho,
                        diffuse_fraction, sunlit=True, soil_direct=None):
    """BRF/HDRF of a plane row canopy with the chosen hotspot and scattering closures."""
    from .radiation import relative_azimuth
    from .vnir import canopy_reflectance, spherical_coefficients
    _check(hotspot)
    _check(scatter)
    lai, hspot = float(settings.get("lai", 0.0)), float(settings.get("hspot", 0.1))
    width, blank, height, azimuth = (float(settings[key]) for key in
                                     ("row_width", "row_blank", "row_height", "row_azimuth"))
    soil_direct = float(soil_rho) if soil_direct is None else float(soil_direct)
    raa = float(relative_azimuth(saa, vaa))
    parts = fractions(hotspot, lai, hspot, width, blank, height, vza, vaa, sza, saa, azimuth)
    if not sunlit:
        parts = np.array([0.0, parts[0] + parts[1], 0.0, parts[2] + parts[3]])
    if scatter == "reference":
        from .row_radiosity import shortwave_multiple
        direct = (1.0 - diffuse_fraction) if sunlit else 0.0
        _, ko, sob, sof = spherical_coefficients(float(sza), float(vza), raa)
        single = direct * ((sob * rho + sof * tau) / ko * parts[2] + soil_direct * parts[0])
        return float(single + shortwave_multiple(lai, width, blank, height, rho, tau, soil_rho, vza, vaa, sza, saa,
                                                 azimuth, direct, diffuse_fraction, soil_view=parts[0] + parts[1]))
    return float(canopy_reflectance(lai, hspot, float(rho), float(tau), float(soil_rho), float(sza), float(vza), raa,
                                    parts, sun_gap(lai, width, blank, height, sza, saa, azimuth), False,
                                    float(diffuse_fraction), sunlit, soil_direct))
