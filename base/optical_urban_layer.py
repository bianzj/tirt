"""Homogeneous/row canopy boundary coupled to urban wall/street radiosity.

The existing plane canopy kernels supply directional reflection. Vegetation
occupies the specified street fraction and acts as a planar street boundary;
this closure does not resolve individual plant/building intersections.
"""
import math
from functools import lru_cache

import numpy as np


def parameters(settings):
    from . import canopy_closure
    kind = str(settings.get('vegetation_model', 'hom')).lower()
    kind = {'homogeneous': 'hom', 'raw': 'row', 'forest': 'crown'}.get(kind, kind)
    closure = canopy_closure.closure(settings)
    crown = None
    if kind == 'crown':
        _, _, _, shape, _, base = canopy_closure.structure(settings)
        crown = (*shape, base)
    return (kind, float(settings.get('lai', 0)), float(settings.get('hspot', .1)),
            float(settings.get('row_width', .5)), float(settings.get('row_blank', .5)),
            float(settings.get('row_height', 1.)), float(settings.get('row_azimuth', 0)),
            str(settings.get('optical_row_method', 'conditioned_corridor')), closure, crown)


def _standard_shape(p):
    kind, crown = p[0], p[9]
    if kind == 'crown':
        return crown[:3], crown[3]
    return (p[3], p[4], p[5]) if kind == 'row' else None, None


@lru_cache(maxsize=4096)
def response(p, solar, solar_azimuth, view, view_azimuth, rho, tau, soil, diffuse):
    kind, lai, hotspot, width, gap, height, bearing, method, closure, _ = p
    if closure != 'legacy':
        from . import canopy_closure
        shape, base = _standard_shape(p)
        return canopy_closure.reflectance(kind, lai, hotspot, shape, rho, tau, soil, view, view_azimuth,
                                          solar, solar_azimuth, bearing, base, diffuse, True, level=closure)
    from .radiation import directional_reflectance
    local = dict(surface_model='plane', vegetation_model=kind, surface_type=0,
                 distribution=2 if kind == 'row' else 1, lai=lai, hspot=hotspot,
                 row_width=width, row_blank=gap, row_height=height,
                 row_azimuth=bearing, optical_row_method=method, canopy_closure=closure)
    return directional_reflectance(local, solar, solar_azimuth, view, view_azimuth,
                                   rho, tau, soil, diffuse_fraction=diffuse,
                                   sunlit=True, use_hapke=False)


@lru_cache(maxsize=1024)
def albedo(p, solar, solar_azimuth, rho, tau, soil, diffuse):
    if p[8] != 'legacy':
        from . import canopy_closure
        shape, _ = _standard_shape(p)
        if diffuse:
            return canopy_closure.sky_albedo(p[0], p[1], shape, rho, tau, soil, p[6])
        return canopy_closure.beam_albedo(p[0], p[1], shape, rho, tau, soil, solar, solar_azimuth, p[6])
    # Angular quadrature of the analytical canopy response, without a spatial
    # crown grid. Cache it across view directions and time-series iterations.
    from .optical_crown import quadrature
    cosine, weights, azimuths = quadrature(12, 24)
    return float(sum(2*mu*weight/len(azimuths)*sum(
        response(p, solar, solar_azimuth, math.degrees(math.acos(mu)), float(azimuth)+p[6],
                 rho, tau, soil, diffuse) for azimuth in azimuths)
        for mu, weight in zip(cosine, weights)))


def reflectance(settings, view, view_azimuth, solar, solar_azimuth,
                rho_roof, rho_wall, rho_street, rho_leaf, tau_leaf):
    from .optical_urban import building_shape, obstruction, joint_obstruction, wall_exchange
    shape = building_shape(settings)
    length, width, height, density, bearing = shape
    footprint = length*width*density
    ground = 1-footprint
    # The vegetation model itself describes the whole street surface: homogeneous,
    # row and crown gaps are already part of its directional response.
    cover = 1.0
    diffuse = float(settings.get('shortwave_diffuse_fraction', 0.))
    sunlit = bool(settings.get('shortwave_sunlit_surface', True))
    if not 0 <= diffuse <= 1:
        raise ValueError('Urban diffuse fraction must be in [0, 1].')
    if not 0 <= solar < 90 or not 0 <= view < 90:
        raise ValueError('Urban optical transport requires above-horizon directions.')
    p = parameters(settings)
    beam = (1-cover)*rho_street+cover*response(
        p, solar, solar_azimuth, view, view_azimuth, rho_leaf, tau_leaf, rho_street, 0.)
    # Retain the selected plane kernel's diffuse closure, including the
    # effective-depth approximation used by the existing row model.
    sky_response = (1-cover)*rho_street+cover*response(
        p, solar, solar_azimuth, view, view_azimuth, rho_leaf, tau_leaf, rho_street, 1.)
    if footprint == 0 or height == 0:
        return float((1-diffuse)*sunlit*beam+diffuse*sky_response)
    beam_albedo = (1-cover)*rho_street+cover*albedo(
        p, solar, solar_azimuth, rho_leaf, tau_leaf, rho_street, 0.)
    sky_albedo = (1-cover)*rho_street+cover*albedo(
        p, solar, solar_azimuth, rho_leaf, tau_leaf, rho_street, 1.)
    view_depth = obstruction(shape, view, view_azimuth)
    solar_depth = obstruction(shape, solar, solar_azimuth)
    relative = abs((view_azimuth-solar_azimuth+180) % 360-180)
    joint = joint_obstruction(view, solar, relative, view_depth, solar_depth)
    normals = np.array([bearing, bearing+90, bearing+180, bearing+270])
    visible = np.array([length, width, length, width])*np.maximum(
        np.cos(np.deg2rad(view_azimuth-normals)), 0.)
    visible /= max(float(visible.sum()), 1e-300)
    incidence = np.maximum(np.cos(np.deg2rad(solar_azimuth-normals)), 0.)*math.tan(math.radians(solar))
    factors, sky, walls = wall_exchange(shape, 0.)
    wall_energy = ground*(-math.expm1(-solar_depth))
    projected = float(walls@incidence)
    wall_sources = wall_energy*incidence/projected if projected > 0 else np.zeros(4)
    reflectances = np.r_[sky_albedo, np.full(4, rho_wall)]
    beam_sources = np.r_[beam_albedo*math.exp(-solar_depth), rho_wall*wall_sources]*sunlit
    sources = (1-diffuse)*beam_sources+diffuse*reflectances*sky
    radiosity = np.linalg.solve(np.eye(5)-reflectances[:, None]*factors, sources)
    higher = factors@radiosity
    wall_mutual = math.exp(-joint_obstruction(view, solar, relative,
        footprint+.5*view_depth, footprint+.5*solar_depth)+footprint+.5*view_depth)
    wall_response = ((1-diffuse)*sunlit*rho_wall*wall_mutual*incidence
                     +rho_wall*(diffuse*sky[1:]+higher[1:]))
    street_visible = ground*math.exp(-view_depth)
    wall_visible = max(0., 1-footprint-street_visible)
    street_response = ((1-diffuse)*sunlit*ground*math.exp(-joint)*beam
                       +street_visible*sky_response*(diffuse*sky[0]+higher[0]))
    return float(footprint*rho_roof*((1-diffuse)*sunlit+diffuse)
                 +wall_visible*(visible@wall_response)+street_response)
