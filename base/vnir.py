"""Optical observation transport, independent of the thermal correction.

Spherical leaves, Lambertian soil, isotropic sky. The two-stream boundary
problem retains canopy/soil multiple scattering; direct single scattering
uses sun/view joint transmission. Row/crown inputs retain their geometric
four-component fractions with an effective-slab diffuse approximation.
"""
import math
from functools import lru_cache

import numpy as np
from scipy.linalg import expm, solve_banded

from .rt_transfer import volscatt

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)


@lru_cache(maxsize=4096)
def spherical_coefficients(sza, vza, raa):
    cs, cv = math.cos(math.radians(sza)), math.cos(math.radians(vza))
    ks, ko = 0.5 / max(cs, 1e-6), 0.5 / max(cv, 1e-6)
    # Uniform leaf-normal cosine is the spherical inclination distribution.
    sob = sof = 0.0
    for mu, weight in zip((_NODES + 1) / 2, _WEIGHTS / 2):
        _, _, fr, ft = volscatt(sza, vza, raa, math.degrees(math.acos(mu)))
        sob += weight * math.pi * fr / (cs * cv)
        sof += weight * math.pi * ft / (cs * cv)
    return ks, ko, sob, sof


def joint_transmission(lai, ks, ko, hspot, sza, vza, raa):
    """Return bottom joint gap and its depth integral (units of LAI).

    Exponential path overlap, with stable exact-hotspot and zero-size limits.
    hspot=0 explicitly disables correlation, including at coincident angles.
    """
    ts, tv = math.tan(math.radians(sza)), math.tan(math.radians(vza))
    distance = math.sqrt(max((ts-tv)**2 + 4*ts*tv*math.sin(math.radians(raa)/2)**2, 0))
    if hspot <= 0:
        k = ks + ko
        return math.exp(-k*lai), -math.expm1(-k*lai) / k
    alpha = distance / hspot * 2 / (ks + ko)
    x = np.r_[(_NODES + 1)/2, 1.0]
    overlap = x if alpha < 1e-10 else -np.expm1(-alpha*x)/alpha
    joint = np.exp(-lai*(ks+ko)*x + lai*math.sqrt(ks*ko)*overlap)
    return float(joint[-1]), float(lai * np.dot(_WEIGHTS/2, joint[:-1]))


@lru_cache(maxsize=2048)
def diffuse_field(lai, ks, rho, tau, rs, direct, diffuse, layers=128):
    """Solve down/up fluxes with D(0)=sky, U(L)=rs*(D(L)+beam(L)).

    Exact layer transfer matrices and a banded boundary solve remain stable
    for black leaves and conservative scattering (rho+tau=1, rs=1).
    """
    dz = lai / layers
    bf = 1/3
    back = (1+bf)/2*rho + (1-bf)/2*tau
    forward = (1-bf)/2*rho + (1+bf)/2*tau
    sb = (ks+bf)/2*rho + (ks-bf)/2*tau
    sf = (ks-bf)/2*rho + (ks+bf)/2*tau
    matrix = np.array([[-(1-forward), back, sf],
                       [-back, 1-forward, -sb], [0, 0, -ks]])
    transfer = expm(matrix*dz)
    r = -transfer[1,0]/transfer[1,1]
    t = 1/transfer[1,1]
    qu = -transfer[1,2]/transfer[1,1]
    qd = transfer[0,2] - transfer[0,1]*transfer[1,2]/transfer[1,1]
    size = 2*(layers+1)
    ab, rhs = np.zeros((7,size)), np.zeros(size)
    def put(i,j,value):
        ab[3+i-j,j] = value
    put(0,0,1)
    rhs[0] = diffuse
    edges = np.arange(layers+1)*dz
    beam = direct*np.exp(-ks*edges)
    for i in range(layers):
        down, up = 2*i+2, 2*i+1
        put(down,2*i+2,1); put(down,2*i,-t); put(down,2*i+3,-r)
        rhs[down] = qd*beam[i]
        put(up,2*i+1,1); put(up,2*i,-r); put(up,2*i+3,-t)
        rhs[up] = qu*beam[i]
    put(size-1,size-1,1); put(size-1,size-2,-rs)
    rhs[-1] = rs*beam[-1]
    field = solve_banded((3,3),ab,rhs).reshape(-1,2)
    # Exact midpoint fields avoid linear-interpolation errors in thick layers.
    mid = (expm(matrix*dz/2) @ np.vstack((field[:-1].T,beam[:-1])))[:2]
    for array in (edges, field, mid):
        array.setflags(write=False)
    return edges[:-1]+dz/2, mid[0], mid[1], float(field[-1,0]), float(field[0,1])


def canopy_reflectance(lai, hspot, rho, tau, rs, sza, vza, raa,
                       fractions, solar_gap, homogeneous=True,
                       diffuse_fraction=0.0, sunlit=True, soil_direct=None):
    """BRF for direct illumination; HDRF when a diffuse fraction is supplied.

    No [0,1] clipping: directional reflectance factors can exceed unity.
    Shaded components receive the solved diffuse field, not a new material rho.
    """
    f = float(diffuse_fraction)
    if not 0 <= f <= 1:
        raise ValueError('Optical diffuse fraction must be in [0,1].')
    if min(rho,tau,rs) < 0 or rho+tau > 1+1e-12 or rs > 1:
        raise ValueError('Invalid optical material coefficients.')
    direct = (1-f) if sunlit else 0.0
    rd = rs if soil_direct is None else soil_direct
    if lai <= 1e-10:
        return direct*rd + f*rs
    ks, ko, sob, sof = spherical_coefficients(float(sza),float(vza),float(raa))
    soil_sun, soil_shade, leaf_sun, leaf_shade = fractions
    view_gap = float(soil_sun+soil_shade)
    if homogeneous:
        joint, path = joint_transmission(lai,ks,ko,hspot,sza,vza,raa)
        single = direct*((sob*rho+sof*tau)*path + rd*joint)
        depth, view_k = lai, ko
    else:
        # Existing GO fractions retain row/crown shadows. Diffuse transport
        # uses the slab with the same solar gap; this is not full 3-D RT.
        depth = max(-math.log(max(float(solar_gap),1e-12))/ks, 1e-8)
        view_k = -math.log(max(view_gap,1e-12))/depth
        single = direct*((sob*rho+sof*tau)/ko*leaf_sun + rd*soil_sun)
    z, down, up, ground_down, _ = diffuse_field(depth,ks,rho,tau,rs,direct,f)
    scale = view_k/ko
    vb = ((ko+1/3)/2*rho + (ko-1/3)/2*tau)*scale
    vf = ((ko-1/3)/2*rho + (ko+1/3)/2*tau)*scale
    scattered = float(np.sum((vb*down+vf*up)*np.exp(-view_k*z))*depth/len(z))
    # All visible soil receives diffuse light; only jointly sunlit soil gets
    # the direct term above. The upwelling field includes soil interreflection.
    scattered += rs*ground_down*view_gap
    return max(float(single+scattered),0.0)
