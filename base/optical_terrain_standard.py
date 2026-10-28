"""Canopy optics of compound terrain for the standard and upgrade closures, in the same structure as the thermal
TerrainForest path.

First order   : the thermal model's visible four components (plane + cone facets, local geometry, hooks of
                canopy_closure) weighted by the spherical-leaf phase coefficients.
Canopy orders : node radiosity of every visible plane / facet canopy in its local frame (beam: orders >= 2,
                sky: all orders), as the thermal canopy scatter term.
Terrain       : single-bounce exchange with the thermal view factors (Terrain.scatter_terms); reflectors use the
                canopy diffuse albedo, sources the sunlit / sky-lit canopy radiosity.
Covers: homogeneous, row and crown canopies.
"""
import math
from functools import lru_cache

import numpy as np

from . import canopy_closure, node_radiosity
from .solar import slope_direction
from .terrain_exchange import TerrainExchange
from .terrain_forest import TerrainForest
from .terrain_rt import _shape_values


@lru_cache(maxsize=8192)
def _beam(kind, lai, shape, rho, tau, soil, vza, vaa, sza, saa, row_azimuth):
    return node_radiosity.beam_response(canopy_closure.NODES[kind], lai, shape, rho, tau, soil, vza, vaa, sza, saa,
                                        row_azimuth)


@lru_cache(maxsize=8192)
def _sky(kind, lai, shape, rho, tau, soil, vza, vaa, row_azimuth):
    return node_radiosity.sky_response(canopy_closure.NODES[kind], lai, shape, rho, tau, soil, vza, vaa, row_azimuth)


def _facet_canopy(kind, lai, shape, cosine):
    """Canopy on a facet tilted by beta: area densities per facet area scale with cos(beta)."""
    if kind == "crown":
        return lai * cosine, (shape[0] * cosine, shape[1], shape[2])
    return lai * cosine, shape


@lru_cache(maxsize=1024)
def _facet_sources(kind, lai, shape, rho, tau, soil, slope, sza, saa, row_azimuth, azimuths):
    """Mean beam radiosity of a sunlit facet canopy (per unit horizontal irradiance) and its diffuse albedo."""
    values = []
    for aspect in np.arange(azimuths) * 360.0 / azimuths:
        local_sun, local_saa, lit = slope_direction(sza, saa, slope, aspect)
        if lit:
            incident = math.cos(math.radians(local_sun)) / math.cos(math.radians(sza))
            values.append(incident * _beam(kind, lai, shape, rho, tau, soil, 0.0, 0.0, min(local_sun, 89.0),
                                           float(local_saa), row_azimuth)[1])
    return (float(np.mean(values)) if values else 0.0, _sky(kind, lai, shape, rho, tau, soil, 0.0, 0.0, row_azimuth)[1])


def reflectance(settings, sza, saa, vza, vaa, rho, tau, soil, diffuse_fraction=0.0, sunlit=True):
    kind, lai, hspot, shape, row_azimuth, base = canopy_closure.structure(settings)
    shape = tuple(float(x) for x in shape)
    shapes = np.asarray(_shape_values(settings), float)
    if not np.any(shapes[:, 2] > 0.):
        return canopy_closure.reflectance(kind, lai, hspot, shape, rho, tau, soil, vza, vaa, sza, saa,
                                          row_azimuth, base, diffuse_fraction, sunlit,
                                          level=canopy_closure.closure(settings))
    model = TerrainForest()
    canopy_closure.install(model, canopy_closure.closure(settings))
    if kind == "crown":
        model.forestshape = np.array([lai, shape[0], shape[2], shape[1], hspot, base])
    else:
        model.forestshape = np.array([lai, 0.001, 3.0, 3.0, hspot, 0.0])
    model.set_vegetation_input(kind, *(shape if kind == "row" else (0.5, 0.5, 1.0)), row_azimuth, "parametric")
    model.set_structural_input(shapes)
    model.set_spectral_input(1.0, 1.0)
    model._collect_optical_geometry = True
    view = min(max(abs(float(vza)), 1.0e-3), 89.999)
    model.set_angular_input(np.array([view]), np.array([saa if abs(float(vza)) <= 1.0e-8 else vaa]), float(sza), float(saa))
    model.calculate_component_direct_emissivity(0, 2)
    beam = (1.0 - diffuse_fraction) * float(sunlit)
    f = float(diffuse_fraction)
    materials = (float(rho), float(tau), float(soil))

    first = float(np.dot(model._optical_first_coefficients, materials))
    terms = TerrainExchange.scatter_terms(model, 0)
    sky_ground = 1.0 - terms["ws"] / terms["pPlaneV"] if terms["pPlaneV"] > 0 else 0.0
    beta = float(np.mean(np.arctan(shapes[:, 0] / shapes[:, 1])))
    cosine = math.cos(beta)
    facet_lai, facet_shape = _facet_canopy(kind, lai, shape, cosine)
    canopy = 0.0
    for lit_weight, weight, facet, v, s, _, va, sa, ra in model._optical_canopy_parts:
        part_lai, part_shape = (facet_lai, facet_shape) if facet else (lai, shape)
        if beam and lit_weight > 0.0:
            canopy += beam * lit_weight * float(_beam(kind, part_lai, part_shape, *materials, v, va, s, sa, ra)[0].sum())
        if f and weight > 0.0:
            canopy += f * weight * (terms["eu"] if facet else sky_ground) * \
                float(_sky(kind, part_lai, part_shape, *materials, v, va, ra)[0].sum())

    facet_beam, facet_albedo = _facet_sources(kind, facet_lai, facet_shape, *materials, math.degrees(beta),
                                              float(sza), float(saa), float(row_azimuth), model.facet_azimuths)
    ground_beam = _beam(kind, lai, shape, *materials, 0.0, 0.0, min(float(sza), 89.0), float(saa), float(row_azimuth))[1]
    ground_albedo = _sky(kind, lai, shape, *materials, 0.0, 0.0, float(row_azimuth))[1]
    x_ff = terms["p"] * terms["i0v"] * facet_albedo
    x_fg = terms["ws"] * ground_albedo
    x_gf = terms["ed"] * terms["i0v"] * facet_albedo
    terrain = beam * ((x_ff * terms["s_ww"] + x_fg * terms["s_ws"]) * facet_beam +
                      x_gf * terms["ground_sunlit"] * ground_beam)
    terrain += f * ((x_ff + x_fg) * terms["eu"] * facet_albedo + x_gf * sky_ground * ground_albedo)
    return float(beam * first + canopy + terrain)
