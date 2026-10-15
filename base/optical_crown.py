"""Optical-only four-component synthesis for statistical ellipsoidal forests.

First scattering defaults to crown-volume joint sun/view geometry. The
optical_crown_geometry='clumping' option reuses the thermal model's clumping
and hotspot equations, with no spatial crown integration by default.
The 'clumping_dual' option separates crown-scale and leaf-scale hotspots
using closed-form projected-area and optical-depth averages.
Higher orders use a conservative three-channel exchange closure: upward
foliage, downward foliage, and upward Lambertian ground. The diffusion_profile
source retains the diffuse transport eigenmode and directional source depth;
beam_diffusion_profile also retains the direct source's downward depth profile.
uniform retains the earlier spatially constant source. This is
an approximation, not a reconstruction of a particular finite tree layout.
No thermal source, temperature, emissivity, or thermal dispatch is changed.
"""
import math
from functools import lru_cache

import numpy as np

from .gap import gap_probability_crown_analytical, gap_probability_crown_inside_voxel, gap_probability_crown_outside_voxel
from .physicsF import ellipsoid_grid_crown_FRT
from .proportion import proportion_bidirectional_crown_one, proportion_bidirectional_crown_voxel_one
from .vnir import spherical_coefficients
from .optical_hotspot import dual_hotspot_fractions


def quadrature(order=12, azimuths=24):
    x, w = np.polynomial.legendre.leggauss(order)
    mu, weight = (x+1)/2, w/2
    return mu, weight, (np.arange(azimuths)+.5)*360/azimuths


@lru_cache(maxsize=8192)
def fractions(lai, density, hcr, rcr, hspot, sza, vza, raa, geometry_mode='volume', hotspot_scales=None):
    """Visible sunlit/shaded soil and foliage, in the thermal component order."""
    if geometry_mode not in ('volume', 'clumping', 'clumping_dual'):
        raise ValueError('Optical crown geometry must be volume, clumping or clumping_dual.')
    if geometry_mode == 'clumping_dual':
        return dual_hotspot_fractions(lai,density,hcr,rcr,hspot,sza,vza,raa,hotspot_scales)
    gv = float(gap_probability_crown_analytical(lai,density,hcr,rcr,vza))
    # Exact coincident-ray limit also avoids the inherited circle-overlap
    # arccos roundoff at native solar angles such as 30.00000000000001 deg.
    separation = ((math.tan(math.radians(vza))-math.tan(math.radians(sza)))**2
                  + 4*math.tan(math.radians(vza))*math.tan(math.radians(sza))*math.sin(math.radians(raa)/2)**2)
    if separation < 1e-14:
        return gv, 0., 1-gv, 0.
    if geometry_mode == 'clumping':
        result = proportion_bidirectional_crown_one(
            lai, density, max(hspot, 1e-12), hcr, rcr,
            np.array([vza]), np.array([sza]), np.array([raa]))
    else:
        result = proportion_bidirectional_crown_voxel_one(
            lai, density, max(hspot, 1e-12), 2*hcr, hcr, rcr,
            np.array([vza]), np.array([sza]), np.array([raa]), np.array([0.]))
    values = np.array(result, dtype=float).ravel()
    if not np.isfinite(values).all() or abs(values.sum()-1) > 1e-8:
        raise ValueError(f'Invalid optical crown component probabilities: sza={sza}, vza={vza}, raa={raa}, values={values}.')
    # The inherited volume approximation can overshoot joint probabilities
    # at grazing views. Enforce probability bounds in this optical adapter;
    # do not alter the shared function or its thermal results.
    gs = float(gap_probability_crown_analytical(lai,density,hcr,rcr,sza))
    soil_sun = float(np.clip(values[0], max(0.,gv+gs-1), min(gv,gs)))
    leaf_sun = float(np.clip(values[2], 0., 1-gv))
    return soil_sun, gv-soil_sun, leaf_sun, 1-gv-leaf_sun


def _single(shape, sza, vza, raa, rho, tau, soil_rho, geometry_mode='volume', hotspot_scales=None):
    # The statistical forest is azimuthally symmetric. Use the same unsigned
    # separation as the public API, including in hemisphere integrations;
    # otherwise the finite crown cubature introduces a spurious handedness.
    raa = abs((float(raa)+180.)%360.-180.)
    p = fractions(*shape, float(sza), float(vza), float(raa), geometry_mode, hotspot_scales)
    _, ko, sob, sof = spherical_coefficients(float(sza), float(vza), float(raa))
    return np.array([(sob*rho+sof*tau)/ko*p[2], soil_rho*p[0]])


@lru_cache(maxsize=256)
def geometry(shape, order=12):
    lai, density, hcr, rcr, _ = shape
    mu, weights, _ = quadrature(order)
    gap = gap_probability_crown_analytical(lai, density, hcr, rcr, np.rad2deg(np.arccos(mu)))
    ground_escape = float(np.dot(2*mu*weights, gap))
    # Isotropic leaf normals imply G=1/2. This identity makes directional
    # escape consistent with projected foliage cover and flux reciprocity.
    leaf_escape = (1-ground_escape)/lai
    if not 0 < leaf_escape <= 1+1e-8:
        raise ValueError('Invalid leaf escape probability.')
    return min(leaf_escape, 1.0), ground_escape


@lru_cache(maxsize=512)
def sunlit_down_escape(shape, sza, order=12):
    """Ground-reaching probability from sunlit leaf volumes (hemisphere mean).

    Within the downward hemisphere the first source is assumed isotropic;
    its total upward/downward energy still distinguishes rho from tau.
    """
    lai, density, hcr, rcr, _ = shape
    vox, _ = ellipsoid_grid_crown_FRT()
    x, y, z = vox[0]*rcr, vox[1]*rcr, vox[2]*hcr
    weights = vox[3]/sum(vox[3])
    lad = lai/(density*4*np.pi/3*rcr*rcr*hcr)

    def transmission(zz, zenith, azimuth):
        inside, length, _, _ = gap_probability_crown_inside_voxel(x, y, zz, lad, hcr, rcr, zenith, azimuth)
        th, ph = np.deg2rad(zenith), np.deg2rad(azimuth)
        outside = gap_probability_crown_outside_voxel(density,
            x+length*np.sin(th)*np.cos(ph), y+length*np.sin(th)*np.sin(ph),
            zz+length*np.cos(th)+hcr, lad, 2*hcr, hcr, rcr, zenith, azimuth)[0]
        return inside*outside

    illuminated = weights*transmission(z, sza, 0.)
    mu, w, az = quadrature(order, 12)
    down = np.zeros_like(z)
    for m, weight in zip(mu, w):
        for a in az:
            # Reflect z about the crown midplane to trace toward the ground.
            down += weight/len(az)*transmission(-z, np.rad2deg(np.arccos(m)), a)
    value = float(np.dot(illuminated, down)/sum(illuminated))
    if not 0 <= value <= 1:
        raise ValueError('Invalid downward first-scattering escape.')
    return value


@lru_cache(maxsize=512)
def primary(shape, sza, rho, tau, rs, order=12, geometry_mode='volume', source_mode='sunlit_volume', hotspot_scales=None):
    """First-scattered energy and the source of second and higher orders."""
    lai, density, hcr, rcr, _ = shape
    gs = float(gap_probability_crown_analytical(lai, density, hcr, rcr, sza))
    ks = .5/math.cos(math.radians(sza))
    beta = (ks+1/3)/(2*ks)
    first = np.array([(1-gs)*(beta*rho+(1-beta)*tau),
                      (1-gs)*((1-beta)*rho+beta*tau), gs*rs])
    mu, w, az = quadrature(order)
    escaped = np.zeros(2)
    for m, weight in zip(mu, w):
        vza = float(np.rad2deg(np.arccos(m)))
        for angle in az:
            escaped += 2*m*weight/len(az)*_single(shape, sza, vza, angle, rho, tau, rs, geometry_mode, hotspot_scales)
    # Escape of the directly lit sources is integrated from their actual
    # directional first scattering, preserving the crown shadow correlation.
    # A finite legacy crown cubature can overestimate directional foliage
    # radiance under nearly horizontal illumination. Normalize that angular
    # source only when its integral exceeds the energy available to it.
    # This is an energy bound, independent of any reference reflectance.
    leaf_scale = min(1., first[0]/escaped[0]) if escaped[0] > 0 else 1.
    escaped[0] *= leaf_scale
    remaining_up = first[0]-escaped[0]
    remaining_ground = first[2]-escaped[1]
    if min(remaining_up, remaining_ground) < -1e-7:
        raise ValueError('First-scattering directional integral exceeds its source energy.')
    remaining_up, remaining_ground = max(remaining_up, 0.), max(remaining_ground, 0.)
    if source_mode == 'sunlit_volume':
        downward_escape = sunlit_down_escape(shape, float(sza), order)
    elif source_mode == 'beam_diffusion_profile':
        from .optical_crown_profile import beam_down_escape
        downward_escape = beam_down_escape(shape, float(sza), order)
    elif source_mode in ('uniform', 'diffusion_profile'):
        # Same mean leaf escape as the higher-order closure. This removes
        # the spatial integral, but approximates the illuminated leaf source
        # by a uniform one and therefore loses its depth dependence.
        downward_escape = geometry(shape, order)[0]
    else:
        raise ValueError('Unknown optical crown source: '+str(source_mode))
    to_ground = first[1]*downward_escape
    remaining_down = first[1]-to_ground
    back, forward = (2*rho+tau)/3, (rho+2*tau)/3
    incident_up = remaining_up+remaining_ground
    secondary = np.array([forward*incident_up+back*remaining_down,
                          back*incident_up+forward*remaining_down, rs*to_ground])
    return first, escaped, secondary, leaf_scale


def exchange_matrix(shape, rho, tau, rs, order=12, source_mode='uniform'):
    leaf_escape, ground_escape = geometry(shape, order)
    upward_escape = downward_escape = leaf_escape
    if source_mode in ('diffusion_profile', 'beam_diffusion_profile'):
        from .optical_crown_profile import source_profile
        _, _, upward_escape, downward_escape = source_profile(shape, rho, tau, order)
    back, forward = (2*rho+tau)/3, (rho+2*tau)/3
    matrix = np.array([[forward*(1-upward_escape), back*(1-downward_escape), forward*(1-ground_escape)],
                       [back*(1-upward_escape), forward*(1-downward_escape), back*(1-ground_escape)],
                       [0., rs*downward_escape, 0.]])
    exits = np.array([upward_escape, 0., ground_escape])
    # Every emitted packet is escaped, absorbed, or scattered again exactly once.
    if np.any(matrix.sum(axis=0)+exits > 1+1e-10):
        raise ValueError('Non-conservative crown exchange operator.')
    return matrix, exits


@lru_cache(maxsize=512)
def transport(shape, sza, rho, tau, rs, order=12, geometry_mode='volume', source_mode='sunlit_volume', hotspot_scales=None):
    first, escaped, secondary, leaf_scale = primary(shape, sza, rho, tau, rs, order, geometry_mode, source_mode, hotspot_scales)
    matrix, exits = exchange_matrix(shape, rho, tau, rs, order, source_mode)
    higher = np.linalg.solve(np.eye(3)-matrix, secondary)
    absorption_first = 1-float(first.sum())
    absorption_second = float(first.sum()-escaped.sum()-secondary.sum())
    absorption_higher = float(np.dot(1-matrix.sum(axis=0)-exits, higher))
    reflected = float(escaped.sum()+np.dot(exits, higher))
    residual = reflected+absorption_first+absorption_second+absorption_higher-1
    if higher.min() < -1e-10 or min(absorption_first, absorption_second, absorption_higher) < -1e-8 or abs(residual) > 1e-9:
        raise ValueError('Crown radiative energy balance failed.')
    upward_escape = downward_escape = geometry(shape, order)[0]
    if source_mode in ('diffusion_profile', 'beam_diffusion_profile'):
        from .optical_crown_profile import source_profile
        _, _, upward_escape, downward_escape = source_profile(shape, rho, tau, order)
    return higher, dict(reflected=reflected, absorbed=1-reflected, balance_residual=residual,
                        first_escape=float(escaped.sum()), multiple_escape=float(exits@higher),
                        first_leaf_energy_scale=leaf_scale,
                        leaf_escape=geometry(shape,order)[0], ground_escape=geometry(shape,order)[1],
                        leaf_escape_up=upward_escape, leaf_escape_down=downward_escape, source_mode=source_mode)


def _beam(shape, sza, vza, raa, rho, tau, rs, order, geometry_mode='volume', source_mode='sunlit_volume', hotspot_scales=None):
    first = _single(shape, sza, vza, raa, rho, tau, rs, geometry_mode, hotspot_scales)
    higher, budget = transport(shape, sza, rho, tau, rs, order, geometry_mode, source_mode, hotspot_scales)
    first[0] *= budget['first_leaf_energy_scale']
    lai, density, hcr, rcr, _ = shape
    gv = float(gap_probability_crown_analytical(lai, density, hcr, rcr, vza))
    # Uniform diffuse foliage source: normalized directional escape; ground
    # is Lambertian. These integrate to the same exits used by the matrix.
    profile_ratio = 1.
    if source_mode in ('diffusion_profile', 'beam_diffusion_profile'):
        from .optical_crown_profile import source_profile, escape_ratios
        depth, attenuation, _, _ = source_profile(shape, rho, tau, order)
        profile_ratio = escape_ratios(vza, depth, attenuation)[0]
    multiple = higher[0]*(1-gv)/lai*profile_ratio + higher[2]*gv
    return float(first.sum()+multiple)


@lru_cache(maxsize=2048)
def _sky(shape, vza, rho, tau, rs, order, geometry_mode='volume', source_mode='sunlit_volume', hotspot_scales=None):
    mu, w, az = quadrature(order)
    value = 0.
    # Isotropic sky is a flux-weighted hemisphere of collimated sources.
    for m, weight in zip(mu, w):
        sza = float(np.rad2deg(np.arccos(m)))
        for angle in az:
            value += 2*m*weight/len(az)*_beam(shape, sza, vza, angle, rho, tau, rs, order, geometry_mode, source_mode, hotspot_scales)
    return value


def crown_reflectance(settings, sza, vza, raa, rho, tau, rs,
                      diffuse_fraction=0., sunlit=True):
    f = float(diffuse_fraction)
    if not np.isfinite([f,rho,tau,rs,sza,vza,raa]).all():
        raise ValueError('Optical coefficients and angles must be finite.')
    if not 0 <= sza < 90 or not 0 <= vza < 90:
        raise ValueError('Crown transport requires sun and view above the horizon.')
    if not 0 <= f <= 1 or min(rho,tau,rs) < 0 or rho+tau > 1+1e-12 or rs > 1:
        raise ValueError('Invalid optical illumination or material coefficients.')
    geometry_mode = str(settings.get('optical_crown_geometry', 'volume')).strip().lower()
    source_mode = str(settings.get('optical_crown_source',
                      'sunlit_volume' if geometry_mode == 'volume' else 'diffusion_profile')).strip().lower()
    if geometry_mode not in ('volume', 'clumping', 'clumping_dual'):
        raise ValueError('Optical crown geometry must be volume, clumping or clumping_dual.')
    if source_mode not in ('sunlit_volume', 'uniform', 'diffusion_profile', 'beam_diffusion_profile'):
        raise ValueError('Unknown optical crown source: '+str(source_mode))
    if geometry_mode == 'clumping_dual' and source_mode == 'sunlit_volume':
        raise ValueError('clumping_dual requires a parameterized source to avoid crown spatial integration.')
    lai = float(settings['lai'])
    if lai <= 1e-8:
        return rs*((1-f)*bool(sunlit)+f)
    shape = (lai,float(settings['stand_density']),float(settings['height']),
             float(settings['width']),float(settings['hspot']))
    if not np.isfinite(shape).all() or min(shape[:4]) <= 0 or shape[4] < 0:
        raise ValueError('Invalid ellipsoidal crown geometry.')
    hotspot_scales = None
    if geometry_mode == 'clumping_dual':
        hotspot_scales = (float(settings.get('optical_hotspot_inner',shape[4])),
                          float(settings.get('optical_hotspot_outer',shape[3]/shape[2])))
        if not np.isfinite(hotspot_scales).all() or min(hotspot_scales) < 0:
            raise ValueError('Optical inner/outer hotspot widths must be finite and nonnegative.')
    order = int(settings.get('optical_crown_quadrature',12))
    if not 4 <= order <= 32:
        raise ValueError('Optical crown quadrature must be between 4 and 32.')
    result = 0.
    if sunlit and f < 1:
        result += (1-f)*_beam(shape,float(sza),float(vza),float(raa),rho,tau,rs,order,geometry_mode,source_mode,hotspot_scales)
    if f > 0:
        result += f*_sky(shape,float(vza),rho,tau,rs,order,geometry_mode,source_mode,hotspot_scales)
    return float(result)
