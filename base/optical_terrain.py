"""Optical forest transport on statistical cone facets and their flat ground."""
import math
from functools import lru_cache
from types import SimpleNamespace

import numpy as np

from .optical_crown import crown_reflectance, quadrature, transport
from .solar import slope_direction
from .terrain_rt import _shape_values, _projected_gap, _mountain_visibility
from .terrain_exchange import TerrainExchange


@lru_cache(maxsize=128)
def exchange_geometry(shapes):
    array = np.asarray(shapes, dtype=float)
    model = SimpleNamespace(shapes=array, Es=1., n_shape=len(array), n_dim=3, n_part=10,
                            vza_=np.asarray([.001]), vaa_=np.asarray([0.]), sza=0., saa=0.)
    terms = TerrainExchange.scatter_terms(model, 0)
    ground_hit = terms['ws']/terms['pPlaneV'] if terms['pPlaneV'] > 0 else 0.
    matrix = np.asarray([[terms['p'], terms['ed']], [ground_hit, 0.]], dtype=float)
    # Optical radiosity requires flux view factors, rather than the legacy
    # view-weighted thermal terms. Ground cosine-weighted occlusion and
    # surface-area reciprocity define both directions of this exchange.
    occupied = sum(density*math.pi*radius**2 for height, radius, density in shapes)
    facet_area = sum(density*math.pi*radius*math.hypot(height, radius)
                     for height, radius, density in shapes)
    if facet_area > 0:
        cosine, weights, _ = quadrature(64)
        ground_hit = sum(2*mu*weight*(1-_projected_gap(shapes, math.degrees(math.acos(mu)), 0., 0.)[0]/(1-occupied))
                         for mu, weight in zip(cosine, weights))
        matrix[1,0] = ground_hit
        matrix[0,1] = ground_hit*(1-occupied)/facet_area
    if np.any(matrix < -1e-12) or np.any(matrix.sum(axis=1) > 1+1e-10):
        raise ValueError('Non-conservative terrain view factors.')
    return matrix


def forest_albedo(settings, solar, rho, tau, soil_rho):
    shape = tuple(float(settings[key]) for key in ('lai', 'stand_density', 'height', 'width', 'hspot'))
    if shape[0] <= 1e-8 or shape[1] == 0:
        return float(soil_rho)
    mode = settings.get('optical_crown_geometry', 'volume')
    source = settings.get('optical_crown_source', 'sunlit_volume' if mode == 'volume' else 'diffusion_profile')
    scales = None
    if mode == 'clumping_dual':
        scales = (float(settings.get('optical_hotspot_inner', shape[-1])),
                  float(settings.get('optical_hotspot_outer', shape[3]/shape[2])))
    _, budget = transport(shape, float(solar), float(rho), float(tau), float(soil_rho),
                          int(settings.get('optical_crown_quadrature', 12)), mode, source, scales)
    return budget['reflected']


def diffuse_albedo(settings, rho, tau, soil_rho):
    cosine, weights, _ = quadrature(int(settings.get('optical_crown_quadrature', 12)))
    return sum(2*value*weight*forest_albedo(settings, math.degrees(math.acos(value)), rho, tau, soil_rho)
               for value, weight in zip(cosine, weights))


def terrain_crown_reflectance(settings, sza, saa, vza, vaa, rho, tau, soil_rho,
                              diffuse_fraction=0., sunlit=True):
    if not 0 <= diffuse_fraction <= 1:
        raise ValueError('Optical diffuse fraction must be in [0, 1].')
    shapes = _shape_values(settings)
    if all(density == 0 for _, _, density in shapes):
        return crown_reflectance(dict(settings, surface_model='plane'), sza, vza,
                                  abs((vaa-saa+180.) % 360.-180.), rho, tau, soil_rho,
                                  diffuse_fraction, sunlit)
    if sum(density*math.pi*radius*radius for _, radius, density in shapes) >= 1:
        raise ValueError('Cone base fractions must sum to less than one.')
    relative = abs((vaa-saa+180.) % 360.-180.)
    plane, plane_sun, _ = _projected_gap(shapes, vza, sza, relative)
    mountain_raw, _, details = _mountain_visibility(shapes, vza, sza, relative)
    factors = exchange_geometry(shapes)
    sky_flux = 1-factors.sum(axis=1)
    local = dict(settings, surface_model='plane')
    # Only this terrain adapter opts into the beam-conditioned first source.
    # Plane/single-slope forests and all thermal transport retain their mode.
    primary = settings.get('optical_terrain_primary', 'beam_profile')
    if primary not in ('beam_profile', 'uniform'):
        raise ValueError('Unknown optical terrain primary source: '+str(primary))
    if primary == 'beam_profile' and local.get('optical_crown_source', 'diffusion_profile') == 'diffusion_profile':
        local['optical_crown_source'] = 'beam_diffusion_profile'
    plane_fraction = plane_sun / plane if plane > 0 else 0.
    direct = crown_reflectance(local, sza, vza, relative, rho, tau, soil_rho,
                               0., sunlit) * plane_fraction
    result = plane * (1-diffuse_fraction) * direct
    if diffuse_fraction:
        result += plane * diffuse_fraction * sky_flux[1] * crown_reflectance(
            local, sza, vza, relative, rho, tau, soil_rho, 1., sunlit)
    scale = (1-plane) / mountain_raw if mountain_raw > 0 else 0.
    facet_direct_flux, facet_area = 0., 0.
    facet_sky_view = 0.
    facet_sky_albedo = 0.
    for alpha, area, solar_area in details:
        slope = math.degrees(alpha)
        cosine_slope = math.cos(math.radians(slope))
        facet = dict(local, lai=local['lai']*cosine_slope,
                     stand_density=local['stand_density']*cosine_slope)
        azimuths = np.arange(0., 360., 5.)
        cos_view = (cosine_slope*math.cos(math.radians(vza)) +
                    math.sin(math.radians(slope))*math.sin(math.radians(vza))*
                    np.cos(np.radians(vaa-azimuths)))
        visible = cos_view > 1e-9
        weights = cos_view[visible]/sum(cos_view[visible])
        illumination = min(1., max(0., solar_area/area)) if area > 0 else 0.
        facet_result = 0.
        mean_direct = 0.
        for aspect in azimuths:
            local_sun, _, lit = slope_direction(sza, saa, slope, aspect)
            if sunlit and lit:
                incident_ratio = math.cos(math.radians(local_sun))/math.cos(math.radians(sza))
                mean_direct += incident_ratio*forest_albedo(facet, min(local_sun, 89.9), rho, tau, soil_rho)/len(azimuths)
        physical_area = area/cosine_slope
        facet_direct_flux += physical_area*illumination*mean_direct
        facet_area += physical_area
        facet_sky_albedo += physical_area*diffuse_albedo(facet, rho, tau, soil_rho)
        for weight, aspect in zip(weights, azimuths[visible]):
            local_view, view_azimuth, _ = slope_direction(vza, vaa, slope, aspect)
            local_sun, sun_azimuth, lit = slope_direction(sza, saa, slope, aspect)
            relative_local = abs((view_azimuth-sun_azimuth+180.) % 360.-180.)
            if sunlit and lit:
                incident_ratio = math.cos(math.radians(local_sun))/math.cos(math.radians(sza))
                facet_result += weight*(1-diffuse_fraction)*illumination*incident_ratio*crown_reflectance(
                    facet, min(local_sun, 89.9), min(local_view, 89.9), relative_local,
                    rho, tau, soil_rho)
            sky_view = crown_reflectance(facet, min(local_sun, 89.9), min(local_view, 89.9),
                                         relative_local, rho, tau, soil_rho, 1., sunlit)
            facet_sky_view += scale*area*weight*sky_view
            facet_result += weight*diffuse_fraction*sky_flux[0]*sky_view
        result += scale*area*facet_result
    if facet_area > 0:
        albedos = np.asarray([facet_sky_albedo/facet_area, diffuse_albedo(local, rho, tau, soil_rho)])
        direct_flux = np.asarray([facet_direct_flux/facet_area,
                                  plane_fraction*forest_albedo(local, sza, rho, tau, soil_rho) if sunlit else 0.])
        sources = (1-diffuse_fraction)*direct_flux+diffuse_fraction*albedos*sky_flux
        outgoing = np.linalg.solve(np.eye(2)-albedos[:, None]*factors, sources)
        incident = factors@outgoing
        ground_sky_view = crown_reflectance(local, sza, vza, relative, rho, tau, soil_rho, 1., sunlit)
        result += facet_sky_view*incident[0]+plane*ground_sky_view*incident[1]
    return float(result)
