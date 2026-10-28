"""Street + four-wall Lambertian radiosity for urban thermal scattering.

All scattering orders are solved on the reciprocal view factors of
``optical_urban.wall_exchange``; each node emits with its own area-mean sunlit
fraction. Equal-height roofs see only the sky, so they emit directly and take no
part in the exchange. A street canopy enters the street node as an area mixture
weighted by its hemispherical leaf coverage; leaf transmittance is neglected (TIR).
The direct (first-order) term is the unchanged urban geometric-optical model.
"""
import math

import numpy as np

from .optical_urban import obstruction, wall_exchange

_WALL_OFFSETS = np.array([0.0, 90.0, 180.0, 270.0])


def _shape(model):
    shapes = np.asarray(model.shapes, dtype=float)
    if shapes.ndim != 2 or shapes.shape[0] != 1:
        raise ValueError("Urban radiosity currently requires one building class.")
    length, width, height, density, bearing = shapes[0][:5]
    return float(length), float(width), float(height), float(density), float(bearing) % 360.0


def _sunlit_fractions(shape, sza, saa):
    """Area-mean sunlit fractions of the street and of the four walls."""
    if sza >= 90.0:
        return 0.0, np.zeros(4)
    depth = obstruction(shape, sza, saa)
    street = math.exp(-depth)
    # A sun-facing wall point at height z is blocked with depth (1 - z/H)·depth.
    wall = 1.0 if depth < 1.0e-12 else -math.expm1(-depth) / depth
    facing = np.cos(np.deg2rad(saa - (shape[4] + _WALL_OFFSETS))) > 0.0
    return street, np.where(facing, wall, 0.0)


def _wall_view_weights(shape, vaa):
    length, width, _, _, bearing = shape
    weights = np.array([length, width, length, width]) * np.maximum(
        np.cos(np.deg2rad(vaa - (bearing + _WALL_OFFSETS))), 0.0)
    total = float(weights.sum())
    return weights / total if total > 0.0 else np.full(4, 0.25)


def _street_canopy(model):
    """Flux-weighted hemispherical leaf coverage of the street and its sunlit leaf share."""
    nz, na = model.coverage_grid
    zenith, azimuth = np.meshgrid((np.arange(nz) + 0.5) * 90.0 / nz, np.arange(na) * 360.0 / na)
    weights = (np.cos(np.deg2rad(zenith)) * np.sin(np.deg2rad(zenith))).ravel()
    fractions = model.canopy_proportions(zenith.ravel(), azimuth.ravel(), model.sza, model.saa,
                                         row_geometry="legacy")
    leaf = float(np.average(fractions[2] + fractions[3], weights=weights))
    sunlit = float(np.average(fractions[2], weights=weights)) / leaf if leaf > 0.0 else 0.0
    return leaf, sunlit


def scattering_operator(model):
    """Reflected radiance leaving the street and each wall, per component temperature (5 x n)."""
    shape = _shape(model)
    vegetated = hasattr(model, "Eveg") and float(model.forestshape[0]) > 0.0
    street_sun, wall_sun = _sunlit_fractions(shape, float(model.sza), float(model.saa))
    factors, _, _ = wall_exchange(shape, 0.0)
    source = np.zeros((5, 8 if hasattr(model, "Eveg") else 6))
    street_e = float(model.Estreat)
    source[0, 2:4] = street_e * np.array([street_sun, 1.0 - street_sun])
    if vegetated:
        cover, leaf_sun = _street_canopy(model)
        source[0, 2:4] *= 1.0 - cover
        lit = street_sun * leaf_sun
        source[0, 6:8] = cover * float(model.Eveg) * np.array([lit, 1.0 - lit])
        street_e = (1.0 - cover) * street_e + cover * float(model.Eveg)
    for k in range(4):
        source[1 + k, 0:2] = float(model.Ewall) * np.array([wall_sun[k], 1.0 - wall_sun[k]])
    reflectance = 1.0 - np.r_[street_e, np.full(4, float(model.Ewall))]
    radiosity = np.linalg.solve(np.eye(5) - reflectance[:, None] * factors, source)
    return shape, reflectance[:, None] * (factors @ radiosity)


def effective_component_emissivity(model):
    """Direct + radiosity component weights per view, in the model's component order.

    Radiance = sum_c w_c B(T_c). With all temperatures equal, sum_c w_c is the
    directional emissivity under a black sky.
    """
    shape, reflected = scattering_operator(model)
    vegetated = hasattr(model, "Eveg")
    emissivity = [model.Ewall, model.Ewall, model.Estreat, model.Estreat, model.Eroof, model.Eroof]
    if vegetated:
        emissivity += [model.Eveg, model.Eveg]
    emissivity = np.asarray(emissivity, dtype=float)
    rows = []
    for k in range(model.n_angle):
        direct = np.asarray(model.calculate_direct_emissivity(k, 2), dtype=float)
        visible = np.divide(direct, emissivity, out=np.zeros_like(direct), where=emissivity > 0.0)
        street = visible[2] + visible[3] + (visible[6] + visible[7] if vegetated else 0.0)
        walls = (visible[0] + visible[1]) * _wall_view_weights(shape, float(model.vaa_[k]))
        components = direct + street * reflected[0] + walls @ reflected[1:]
        if vegetated:
            canopy = np.zeros_like(components)
            canopy[2], canopy[3], canopy[6], canopy[7] = model._canopy_scatter(k)
            budget = max(0.0, 1.0 - float(components.sum()))
            if canopy.sum() > budget:
                canopy *= budget / canopy.sum()
            components = components + canopy
        rows.append(components)
    return np.asarray(rows)
