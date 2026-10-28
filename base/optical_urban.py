"""Sunlit urban facets and conservative wall/street optical radiosity."""
import math
from functools import lru_cache

import numpy as np

from .optical_crown import (_single, primary, exchange_matrix, geometry, quadrature,
                            crown_reflectance)


def building_shape(settings):
    values = tuple(float(settings.get(key, default)) for key, default in (
        ('urban_building_length', 10.), ('urban_building_width', 10.),
        ('urban_building_height', 20.), ('urban_density', .0015)))
    azimuth = float(settings.get('urban_building_azimuth', settings.get('building_azimuth', 0.))) % 360.
    if not np.isfinite((*values, azimuth)).all() or min(values) < 0:
        raise ValueError('Urban dimensions and density must be finite and nonnegative.')
    if values[0]*values[1]*values[3] >= 1:
        raise ValueError('Building footprint must be less than one.')
    return (*values, azimuth)


def obstruction(shape, zenith, azimuth, origin=0.):
    length, width, height, density, bearing = shape
    relative = math.radians(azimuth-bearing)
    cross_section = length*abs(math.cos(relative))+width*abs(math.sin(relative))
    return density*max(height-origin, 0.)*math.tan(math.radians(zenith))*cross_section


def joint_obstruction(view, solar, relative, view_depth, solar_depth):
    tangent_view, tangent_solar = math.tan(math.radians(view)), math.tan(math.radians(solar))
    denominator = tangent_view+tangent_solar
    distance = math.sqrt(max(0., tangent_view*tangent_view+tangent_solar*tangent_solar-
                              2*tangent_view*tangent_solar*math.cos(math.radians(relative))))
    overlap = distance/denominator if denominator > 1e-12 else 0.
    return max(view_depth, solar_depth)+min(view_depth, solar_depth)*overlap


@lru_cache(maxsize=128)
def wall_exchange(shape, source_height):
    length, width, height, density, bearing = shape
    ground_area = 1-length*width*density
    wall_areas = density*height*np.array([length, width, length, width])
    cosine, weights, azimuths = quadrature(20, 72)
    ground_sky = sum(2*value*weight*sum(math.exp(-obstruction(shape, math.degrees(math.acos(value)),
                    azimuth+bearing, source_height)) for azimuth in azimuths)/len(azimuths)
                    for value, weight in zip(cosine, weights))
    ground_to_wall = (1-ground_sky)*wall_areas/wall_areas.sum()
    wall_to_ground = ground_area*ground_to_wall/wall_areas
    normals = np.deg2rad(np.array([bearing, bearing+90, bearing+180, bearing+270]))
    directions = (np.arange(720)+.5)*2*np.pi/720+math.radians(bearing)
    incoming = np.maximum(np.cos(directions[:, None]-normals), 0.)
    outgoing = np.maximum(-np.cos(directions[:, None]-normals), 0.)
    kernel = (incoming.T@outgoing)*wall_areas[:, None]*wall_areas[None, :]
    kernel = (kernel+kernel.T)/2
    target = wall_areas*(1-2*wall_to_ground)
    if np.min(target) < -1e-10:
        raise ValueError('Invalid urban view-factor reciprocity.')
    scale = np.ones(4)
    for _ in range(400):
        totals = scale*(kernel@scale)
        if np.max(abs(totals-target)) < 1e-11:
            break
        scale *= np.sqrt(np.maximum(target, 0.)/np.maximum(totals, 1e-300))
    mutual = scale[:, None]*kernel*scale[None, :]/wall_areas[:, None]
    matrix = np.zeros((5, 5))
    matrix[0, 1:] = ground_to_wall
    matrix[1:, 0] = wall_to_ground
    matrix[1:, 1:] = mutual
    sky = 1-matrix.sum(axis=1)
    if np.min(matrix) < -1e-10 or np.min(sky) < -1e-10:
        raise ValueError('Non-conservative urban radiosity geometry.')
    for array in (matrix, sky, wall_areas):
        array.setflags(write=False)
    return matrix, sky, wall_areas


def forest_parameters(settings):
    shape = tuple(float(settings[key]) for key in ('lai', 'stand_density', 'height', 'width', 'hspot'))
    mode = settings.get('optical_crown_geometry', 'clumping_dual')
    source = settings.get('optical_crown_source', 'sunlit_volume' if mode == 'volume' else 'diffusion_profile')
    scales = ((float(settings.get('optical_hotspot_inner', shape[-1])),
               float(settings.get('optical_hotspot_outer', shape[3]/shape[2])))
              if mode == 'clumping_dual' else None)
    return shape, mode, source, scales, int(settings.get('optical_crown_quadrature', 12))


@lru_cache(maxsize=1024)
def forest_origins(shape, solar, rho, tau, soil_rho, mode, source, scales, order):
    first, escaped, secondary, leaf_scale = primary(shape, solar, rho, tau, soil_rho,
                                                   order, mode, source, scales)
    matrix, exits = exchange_matrix(shape, rho, tau, soil_rho, order, source)
    ground_remaining = max(0., first[2]-escaped[1])
    backward, forward = (2*rho+tau)/3, (rho+2*tau)/3
    ground_secondary = np.array([forward*ground_remaining, backward*ground_remaining, 0.])
    leaf_higher = np.linalg.solve(np.eye(3)-matrix, secondary-ground_secondary)
    ground_higher = np.linalg.solve(np.eye(3)-matrix, ground_secondary)
    fluxes = np.array([escaped[0]+exits@leaf_higher, escaped[1]+exits@ground_higher])
    for array in (leaf_higher, ground_higher, fluxes):
        array.setflags(write=False)
    return leaf_higher, ground_higher, fluxes, leaf_scale


def forest_components(parameters, solar, view, relative, rho, tau, soil_rho,
                      illumination_leaf=1., illumination_soil=1.):
    shape, mode, source, scales, order = parameters
    leaf_higher, ground_higher, _, leaf_scale = forest_origins(
        shape, solar, rho, tau, soil_rho, mode, source, scales, order)
    single = _single(shape, solar, view, relative, rho, tau, soil_rho, mode, scales)
    single[0] *= leaf_scale
    higher = illumination_leaf*leaf_higher+illumination_soil*ground_higher
    from .gap import gap_probability_crown_analytical
    gap = float(gap_probability_crown_analytical(*shape[:4], view))
    ratio = 1.
    if source == 'diffusion_profile':
        from .optical_crown_profile import source_profile, escape_ratios
        depth, attenuation, _, _ = source_profile(shape, rho, tau, order)
        ratio = escape_ratios(view, depth, attenuation)[0]
    return single, np.array([higher[0]*(1-gap)/shape[0]*ratio, higher[2]*gap])


@lru_cache(maxsize=1024)
def forest_sky_components(parameters, view, rho, tau, soil_rho):
    cosine, weights, azimuths = quadrature(parameters[-1])
    result = np.zeros(2)
    for value, weight in zip(cosine, weights):
        solar = float(math.degrees(math.acos(value)))
        for azimuth in azimuths:
            single, multiple = forest_components(parameters, solar, view, float(azimuth), rho, tau, soil_rho)
            result += 2*value*weight/len(azimuths)*(single+multiple)
    return result


def urban_reflectance(settings, view, view_azimuth, solar, solar_azimuth, rho_roof, rho_wall,
                       rho_street, rho_leaf, tau_leaf):
    kind = str(settings.get('vegetation_model', 'bare')).strip().lower()
    from . import canopy_closure
    standard_canopy = (canopy_closure.closure(settings) != 'legacy' and float(settings.get('lai', 0.)) > 1e-8 and
                       (kind in ('hom', 'homogeneous', 'row', 'raw') or
                        (kind in ('crown', 'forest') and float(settings.get('stand_density', 0.)) > 0)))
    if standard_canopy or (kind in ('hom', 'homogeneous', 'row', 'raw') and float(settings.get('lai', 0.)) > 1e-8):
        from .optical_urban_layer import reflectance
        return reflectance(settings, view, view_azimuth, solar, solar_azimuth,
                           rho_roof, rho_wall, rho_street, rho_leaf, tau_leaf)
    shape = building_shape(settings)
    length, width, height, density, bearing = shape
    footprint = length*width*density
    ground_area = 1-footprint
    forest_fraction = (1.0 if kind in ('crown', 'forest')
                       and float(settings.get('lai', 0.)) > 1e-8
                       and float(settings.get('stand_density', 0.)) > 0 else 0.0)
    diffuse = float(settings.get('shortwave_diffuse_fraction', 0.))
    sunlit = bool(settings.get('shortwave_sunlit_surface', True))
    if not 0 <= diffuse <= 1:
        raise ValueError('Optical diffuse fraction must be in [0, 1].')
    if footprint == 0 or height == 0:
        if forest_fraction:
            return ((1-forest_fraction)*rho_street*((1-diffuse)*sunlit+diffuse)+
                    forest_fraction*crown_reflectance(dict(settings, surface_model='plane'), solar, view,
                    abs((view_azimuth-solar_azimuth+180.) % 360.-180.), rho_leaf, tau_leaf, rho_street, diffuse, sunlit))
        return rho_street*((1-diffuse)*sunlit+diffuse)
    relative = abs((view_azimuth-solar_azimuth+180.) % 360.-180.)
    view_ground = obstruction(shape, view, view_azimuth)
    solar_ground = obstruction(shape, solar, solar_azimuth)
    joint_ground = joint_obstruction(view, solar, relative, view_ground, solar_ground)
    street_visible = ground_area*math.exp(-view_ground)
    normals = np.array([bearing, bearing+90, bearing+180, bearing+270])
    view_weights = np.array([length, width, length, width])*np.maximum(np.cos(np.deg2rad(view_azimuth-normals)), 0.)
    view_weights /= max(view_weights.sum(), 1e-300)
    incidence = np.maximum(np.cos(np.deg2rad(solar_azimuth-normals)), 0.)*math.tan(math.radians(solar))
    crown_height = 0.
    leaf_visible, soil_visible = 0., street_visible
    forest_beam, forest_sky = np.zeros(2), np.zeros(2)
    forest_flux = rho_street*math.exp(-solar_ground)
    forest_albedo = rho_street
    incident_forest = math.exp(-solar_ground)
    if forest_fraction:
        parameters = forest_parameters(settings)
        crown_height = float(settings.get('crown_base_height', 3.))+float(settings['height'])
        if crown_height+float(settings['height']) >= height:
            raise ValueError('This urban closure requires tree crowns below equal-height rooftops.')
        view_leaf = obstruction(shape, view, view_azimuth, crown_height)
        solar_leaf = obstruction(shape, solar, solar_azimuth, crown_height)
        joint_leaf = joint_obstruction(view, solar, relative, view_leaf, solar_leaf)
        illumination = np.exp(-np.array([solar_leaf, solar_ground]))
        single, multiple = forest_components(parameters, solar, view, relative, rho_leaf, tau_leaf,
                                              rho_street, *illumination)
        sky_parts = forest_sky_components(parameters, view, rho_leaf, tau_leaf, rho_street)
        fluxes = forest_origins(parameters[0], solar, rho_leaf, tau_leaf, rho_street,
                                parameters[1], parameters[2], parameters[3], parameters[4])[2]
        forest_beam = ground_area*(single*np.exp(-np.array([joint_leaf, joint_ground]))+
                                  multiple*np.exp(-np.array([view_leaf, view_ground])))
        forest_sky = ground_area*sky_parts*np.exp(-np.array([view_leaf, view_ground]))
        forest_flux = float(fluxes@illumination)
        from .gap import gap_probability_crown_analytical
        gap = float(gap_probability_crown_analytical(*parameters[0][:4], view))
        solar_gap = float(gap_probability_crown_analytical(*parameters[0][:4], solar))
        leaf_visible = forest_fraction*ground_area*math.exp(-view_leaf)*(1-gap)
        soil_visible = street_visible*(1-forest_fraction+forest_fraction*gap)
        incident_forest = float(np.array([1-solar_gap, solar_gap])@illumination)
        from .optical_terrain import diffuse_albedo
        forest_albedo = diffuse_albedo(settings, rho_leaf, tau_leaf, rho_street)
        crown_height *= forest_fraction*(1-geometry(parameters[0], parameters[-1])[1])
    wall_visible = max(0., 1-footprint-soil_visible-leaf_visible)
    wall_mutual = math.exp(-joint_obstruction(view, solar, relative,
                    footprint+.5*view_ground, footprint+.5*solar_ground)+footprint+.5*view_ground)
    factors, sky, wall_areas = wall_exchange(shape, crown_height)
    albedos = np.r_[(1-forest_fraction)*rho_street+forest_fraction*forest_albedo, np.full(4, rho_wall)]
    incoming_ground = ((1-forest_fraction)*math.exp(-solar_ground)+forest_fraction*incident_forest)
    wall_energy = ground_area*(1-incoming_ground)
    projected_walls = float(wall_areas@incidence)
    wall_sources = (wall_energy*incidence/projected_walls if projected_walls > 0 else np.zeros(4))
    bare_wall_energy = ground_area*(-math.expm1(-solar_ground))
    tree_shade = wall_energy/bare_wall_energy if bare_wall_energy > 0 else 1.
    sources_beam = np.r_[(1-forest_fraction)*rho_street*math.exp(-solar_ground)+forest_fraction*forest_flux,
                         rho_wall*wall_sources]*sunlit
    sources = (1-diffuse)*sources_beam+diffuse*albedos*sky
    radiosity = np.linalg.solve(np.eye(5)-albedos[:, None]*factors, sources)
    incoming_multiple = factors@radiosity
    wall_response = (1-diffuse)*sunlit*rho_wall*wall_mutual*incidence*tree_shade+rho_wall*(diffuse*sky[1:]+incoming_multiple[1:])
    ground_response = ((1-diffuse)*sunlit*((1-forest_fraction)*ground_area*rho_street*math.exp(-joint_ground)+forest_fraction*forest_beam.sum())+
                       (diffuse*sky[0]+incoming_multiple[0])*((1-forest_fraction)*street_visible*rho_street+forest_fraction*forest_sky.sum()))
    return float(footprint*rho_roof*((1-diffuse)*sunlit+diffuse)+wall_visible*(view_weights@wall_response)+ground_response)
