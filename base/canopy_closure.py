"""Canopy closures shared by thermal and optical bands: four components + node radiosity.

Every level keeps the same form - first order from the four components (phase coefficients / emissivities),
multiple scattering from geometric-optics node radiosity with bi-Lambertian leaves (hom3 / row3 / crown2),
and the same terrain / urban macro structures. Only the four components differ:
  standard : CI transmittance + dual-scale hotspot. Sunlit soil: hom = layer hotspot, rows = CI two-phase strip
             overlap, crowns = Poisson-crown lens overlap; sunlit leaves: covered-phase two-layer approximation.
  upgrade  : exact geometric optics. hom = Beer-Lambert + K-integral hotspot (4SAIL/Kuusk joint gap), rows =
             cross-section integral (row parametric), crowns = FRT-type Poisson spheroids (crown parametric).
  legacy   : the earlier closures everywhere (not handled here).
"""
import functools
import math

import numpy as np

from . import node_radiosity, row_closure
from .crown_parametric import crown_fractions
from .hotspot import two_layer_sunlit_leaf_covered
from .vnir import joint_transmission, spherical_coefficients

G = 0.5
NODES = {"hom": "hom3", "row": "row3", "crown": "crown2"}
HOMOGENEOUS = (0.5, 0.0, 1.0)
LEVELS = ("standard", "upgrade", "legacy")


def closure(settings):
    value = str(settings.get("canopy_closure") or "standard").strip().lower()
    if value not in LEVELS:
        raise ValueError("canopy_closure must be standard, upgrade or legacy")
    return value


def _cosine(angle):
    return math.cos(math.radians(min(abs(float(angle)), 89.9)))


def _relative(vaa, saa):
    return abs((float(vaa) - float(saa) + 180.0) % 360.0 - 180.0)


def fractions(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth=0.0, base=None, level="standard"):
    """Sunlit soil, shaded soil, sunlit leaf, shaded leaf for one direction (degrees).

    shape: rows (width, blank, height); crowns (density, vertical semiaxis, horizontal radius); hom ignored.
    """
    if lai <= 0.0 or (kind == "crown" and shape[0] <= 0.0):
        return np.array([1.0, 0.0, 0.0, 0.0])
    zv, zs = min(abs(float(vza)), 89.9), min(abs(float(sza)), 89.9)
    if kind == "crown":
        density, vertical, horizontal = shape
        args = (lai, density, hspot, vertical, horizontal, vertical if base is None else base)
        parts = crown_fractions(*args, zv, zs, _relative(vaa, saa))[:, 0]
        if level == "upgrade":
            return parts
        soil = parts[:2]
        sun = crown_fractions(*args, zs, zs, 0.0)[:, 0]
        sun_gap = sun[0] + sun[1]
    elif level == "upgrade":
        if kind == "row":
            return np.asarray(row_closure.fractions("reference", lai, hspot, *shape, vza, vaa, sza, saa, row_azimuth))
        ks, kv = G / math.cos(math.radians(zs)), G / math.cos(math.radians(zv))
        joint, path = joint_transmission(lai, ks, kv, hspot, zs, zv, _relative(vaa, saa))
        gap = math.exp(-kv * lai)
        return np.array([joint, gap - joint, kv * path, 1.0 - gap - kv * path])
    else:
        shape = HOMOGENEOUS if kind == "hom" else tuple(shape)
        soil = row_closure.conditioned_soil(lai, hspot, *shape, vza, vaa, sza, saa, row_azimuth)
        sun_gap = row_closure.sun_gap(lai, *shape, sza, saa, row_azimuth)
    view_gap = soil[0] + soil[1]
    leaf = two_layer_sunlit_leaf_covered(lai, view_gap, sun_gap, _cosine(vza), _cosine(sza),
                                         row_closure.leaf_correlation(hspot, vza, vaa, sza, saa))
    return np.array([soil[0], soil[1], leaf, 1.0 - view_gap - leaf])


def _shape(kind, shape):
    return HOMOGENEOUS if kind == "hom" else tuple(float(x) for x in shape)


def thermal_weights(kind, lai, hspot, shape, emis_leaf, emis_soil, vza, vaa, sza, saa, row_azimuth=0.0,
                    base=None, emis_soil_direct=None, level="standard"):
    """Effective component emissivities [soil sun, soil shade, leaf sun, leaf shade] (black sky)."""
    es_direct = emis_soil if emis_soil_direct is None else emis_soil_direct
    parts = fractions(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth, base, level)
    direct = parts * np.array([es_direct, es_direct, emis_leaf, emis_leaf])
    if lai <= 0.0:
        return direct
    return direct + node_radiosity.thermal_scatter(NODES[kind], lai, _shape(kind, shape), emis_leaf, emis_soil,
                                                   vza, vaa, sza, saa, row_azimuth)


def reflectance(kind, lai, hspot, shape, rho, tau, soil, vza, vaa, sza, saa, row_azimuth=0.0, base=None,
                diffuse=0.0, sunlit=True, soil_direct=None, level="standard"):
    """BRF / HDRF: first-order four components + node-radiosity higher orders and sky."""
    soil_first = soil if soil_direct is None else soil_direct
    direct = (1.0 - diffuse) if sunlit else 0.0
    if lai <= 0.0 or (kind == "crown" and shape[0] <= 0.0):
        return float(direct * soil_first + diffuse * soil)
    parts = fractions(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth, base, level)
    if not sunlit:
        parts = np.array([0.0, parts[0] + parts[1], 0.0, parts[2] + parts[3]])
    _, ko, sob, sof = spherical_coefficients(min(float(sza), 89.9), min(float(vza), 89.9), _relative(vaa, saa))
    single = direct * ((sob * rho + sof * tau) / ko * parts[2] + soil_first * parts[0])
    return float(single + node_radiosity.shortwave_multiple(NODES[kind], lai, _shape(kind, shape), rho, tau, soil,
                                                            vza, vaa, sza, saa, row_azimuth, direct, diffuse))


# ---- canopies inside terrain / urban models -----------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def fraction_hook(level):
    """canopy_fraction_model(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth, base) of one level."""
    def hook(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth=0.0, base=None):
        return fractions(kind, lai, hspot, shape, vza, vaa, sza, saa, row_azimuth, base, level)
    hook.level = level
    return hook


def canopy_scatter_hook(kind, lai, emis_leaf, emis_soil, vza, sza, density, vertical, horizontal,
                        view_azimuth=None, solar_azimuth=None, row_shape=None):
    """Node-radiosity canopy scattering (4, N) with the canopy_scatter signature (same for every level).

    Rows need the view / solar azimuths and row_shape = (width, blank, height, row azimuth); hom and crown2
    schedules are azimuth-free.
    """
    el, es = float(np.ravel(emis_leaf)[0]), float(np.ravel(emis_soil)[0])
    vza, sza = (a.ravel() for a in np.broadcast_arrays(np.atleast_1d(vza).astype(float), np.atleast_1d(sza).astype(float)))
    kind = "hom" if kind in ("hom", "homogeneous") else kind
    if lai <= 0.0 or (kind == "crown" and density <= 0.0):
        return np.zeros((4, vza.size))
    vaa = np.zeros_like(vza) if view_azimuth is None else np.broadcast_to(np.atleast_1d(view_azimuth).astype(float), vza.shape)
    saa = np.zeros_like(vza) if solar_azimuth is None else np.broadcast_to(np.atleast_1d(solar_azimuth).astype(float), vza.shape)
    if kind == "crown":
        shape, row_azimuth = (float(density), float(vertical), float(horizontal)), 0.0
    elif kind == "row":
        shape, row_azimuth = tuple(float(x) for x in row_shape[:3]), float(row_shape[3])
    else:
        shape, row_azimuth = HOMOGENEOUS, 0.0
    return np.asarray([node_radiosity.thermal_scatter(NODES[kind], float(lai), shape, el, es, v, a, s, sa, row_azimuth)
                       for v, a, s, sa in zip(vza, vaa, sza, saa)]).T


def install(model, level):
    """Attach the canopy closure of `level` to a TerrainForest / Urban_Veg model."""
    model.canopy_fraction_model = fraction_hook(level)
    model.canopy_scatter_model = canopy_scatter_hook


def structure(settings):
    """(kind, lai, hspot, shape, row azimuth, crown base) of an optical / thermal settings dict."""
    kind = str(settings.get("vegetation_model", "hom")).strip().lower()
    kind = {"homogeneous": "hom", "raw": "row", "forest": "crown"}.get(kind, kind)
    lai, hspot = float(settings.get("lai", 0.0)), float(settings.get("hspot", 0.1))
    if kind == "crown":
        lai, hspot, shape, base = _crown(settings)
        return kind, lai, hspot, shape, 0.0, base
    if kind == "row":
        shape = tuple(float(settings.get(k, d)) for k, d in (("row_width", 0.5), ("row_blank", 0.5), ("row_height", 1.0)))
        return kind, lai, hspot, shape, float(settings.get("row_azimuth", 0.0)), None
    return "hom", lai, hspot, HOMOGENEOUS, 0.0, None


def beam_albedo(kind, lai, shape, rho, tau, soil, sza, saa=0.0, row_azimuth=0.0):
    if lai <= 1.0e-8 or (kind == "crown" and shape[0] <= 0.0):
        return float(soil)
    return node_radiosity.beam_response(NODES[kind], lai, _shape(kind, shape), rho, tau, soil, 0.0, 0.0,
                                        sza, saa, row_azimuth)[1]


def sky_albedo(kind, lai, shape, rho, tau, soil, row_azimuth=0.0):
    if lai <= 1.0e-8 or (kind == "crown" and shape[0] <= 0.0):
        return float(soil)
    return node_radiosity.sky_response(NODES[kind], lai, _shape(kind, shape), rho, tau, soil, 0.0, 0.0, row_azimuth)[1]


def _crown(settings):
    lai, density = float(settings["lai"]), float(settings["stand_density"])
    vertical, horizontal = float(settings["height"]), float(settings["width"])
    base = next((float(settings[k]) for k in ("crown_bottom", "crown_base_height", "crown_height")
                 if settings.get(k) not in (None, "")), vertical)
    return lai, float(settings.get("hspot", 0.1)), (density, vertical, horizontal), base


def crown_reflectance(settings, sza, vza, raa, rho, tau, soil, diffuse=0.0, sunlit=True):
    """Drop-in for optical_crown.crown_reflectance (relative azimuth only)."""
    lai, hspot, shape, base = _crown(settings)
    return reflectance("crown", lai, hspot, shape, rho, tau, soil, vza, raa, sza, 0.0, 0.0, base, diffuse, sunlit,
                       level=closure(settings))
