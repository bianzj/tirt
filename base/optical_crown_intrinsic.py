"""Opt-in crown-intrinsic source depth for the finite-terrain optical adapter.

The uniform-source escape still includes statistical crown-group gaps.
The exponential source varies over the central crown chord, D=2*h*LAD,
instead of a depth inferred from the gap-diluted whole-forest escape.
This is a central-chord source approximation, not a spherical volume solution.
No shared optical or thermal profile function is patched or replaced.
"""
import math
from functools import lru_cache
import numpy as np
from . import optical_crown as crown
from .gap import gap_probability_crown_analytical as gap
from .optical_crown_profile import escape_ratios, beam_down_escape as old_down_escape


def intrinsic_depth(shape):
    lai,density,hcr,rcr,_=shape
    lad=lai/(density*4*math.pi/3*hcr*rcr*rcr)
    return 2*hcr*lad


@lru_cache(maxsize=512)
def source_profile(shape,rho,tau,order):
    depth=intrinsic_depth(shape)
    backward,forward=(2*rho+tau)/3,(rho+2*tau)/3
    attenuation=math.sqrt(max(0.,(1-forward)**2-backward**2))
    cosine,weights,_=crown.quadrature(order)
    upward_escape=downward_escape=0.
    for mu,weight in zip(cosine,weights):
        view=math.degrees(math.acos(mu))
        up,down=escape_ratios(view,depth,attenuation)
        uniform=2*mu*weight*(1-float(gap(*shape[:4],view)))/shape[0]
        upward_escape+=uniform*up
        downward_escape+=uniform*down
    if not (0<=upward_escape<=1+1e-10 and 0<=downward_escape<=1+1e-10):
        raise ValueError('Intrinsic crown source has invalid escaped flux.')
    return depth,attenuation,float(upward_escape),float(downward_escape)


@lru_cache(maxsize=512)
def beam_down_escape(shape,solar,order):
    depth=intrinsic_depth(shape)
    attenuation=.5/math.cos(math.radians(solar))
    cosine,weights,_=crown.quadrature(order)
    value=0.
    for mu,weight in zip(cosine,weights):
        view=math.degrees(math.acos(mu))
        uniform=2*mu*weight*(1-float(gap(*shape[:4],view)))/shape[0]
        value+=uniform*escape_ratios(view,depth,attenuation)[1]
    if not 0<=value<=1:
        raise ValueError('Intrinsic beam source has invalid downward escape.')
    return float(value)


@lru_cache(maxsize=512)
def primary(shape,solar,rho,tau,soil,order,mode,scales):
    first,escaped,secondary,scale=crown.primary(shape,solar,rho,tau,soil,order,mode,
                                               'beam_diffusion_profile',scales)
    # GO first escape and the first emitted energy remain exact copies of
    # the previous adapter. Only the downward source-to-collision handover
    # changes with its crown-intrinsic depth.
    delta=first[1]*(beam_down_escape(shape,solar,order)-old_down_escape(shape,solar,order))
    backward,forward=(2*rho+tau)/3,(rho+2*tau)/3
    secondary=secondary+delta*np.array([-backward,-forward,soil])
    if secondary.min()<-1e-10:
        raise ValueError('Intrinsic source produces negative collided energy.')
    return first,escaped,secondary,scale


def exchange_matrix(shape,rho,tau,soil,order):
    _,_,eu,ed=source_profile(shape,rho,tau,order)
    eg=crown.geometry(shape,order)[1]
    back,forward=(2*rho+tau)/3,(rho+2*tau)/3
    matrix=np.array([[forward*(1-eu),back*(1-ed),forward*(1-eg)],
                     [back*(1-eu),forward*(1-ed),back*(1-eg)],
                     [0.,soil*ed,0.]])
    exits=np.array([eu,0.,eg])
    if np.any(matrix.sum(axis=0)+exits>1+1e-10):
        raise ValueError('Intrinsic crown exchange is not conservative.')
    return matrix,exits


@lru_cache(maxsize=512)
def transport(shape,solar,rho,tau,soil,order,mode,scales):
    first,escaped,secondary,scale=primary(shape,solar,rho,tau,soil,order,mode,scales)
    matrix,exits=exchange_matrix(shape,rho,tau,soil,order)
    higher=np.linalg.solve(np.eye(3)-matrix,secondary)
    absorbed=1-first.sum()+first.sum()-escaped.sum()-secondary.sum()+\
             (1-matrix.sum(axis=0)-exits)@higher
    reflected=float(escaped.sum()+exits@higher)
    residual=float(reflected+absorbed-1)
    if higher.min()<-1e-10 or absorbed<-1e-8 or abs(residual)>1e-9:
        raise ValueError('Intrinsic crown energy balance failed.')
    return higher,dict(reflected=reflected,absorbed=float(absorbed),balance_residual=residual,
                       first_leaf_energy_scale=scale)


def beam(shape,solar,view,relative,rho,tau,soil,order,mode,scales):
    first=crown._single(shape,solar,view,relative,rho,tau,soil,mode,scales)
    higher,budget=transport(shape,solar,rho,tau,soil,order,mode,scales)
    first[0]*=budget['first_leaf_energy_scale']
    depth,attenuation,_,_=source_profile(shape,rho,tau,order)
    gv=float(gap(*shape[:4],view))
    return float(first.sum()+higher[0]*(1-gv)/shape[0]*escape_ratios(view,depth,attenuation)[0]+higher[2]*gv)


@lru_cache(maxsize=2048)
def sky(shape,view,rho,tau,soil,order,mode,scales):
    mu,w,az=crown.quadrature(order)
    return sum(2*m*weight/len(az)*beam(shape,float(np.degrees(np.arccos(m))),view,float(a),
               rho,tau,soil,order,mode,scales) for m,weight in zip(mu,w) for a in az)


def reflectance(settings,solar,view,relative,rho,tau,soil,f=0.,sunlit=True):
    if float(settings['lai'])<=1e-8:
        return float(soil*((1-f)*bool(sunlit)+f))
    from .optical_terrain_ports import parameters
    shape,order,mode,scales=parameters(settings)
    return ((1-f)*beam(shape,solar,view,relative,rho,tau,soil,order,mode,scales) if sunlit else 0.)+\
           (f*sky(shape,view,rho,tau,soil,order,mode,scales) if f else 0.)


def albedo(settings,solar,rho,tau,soil):
    from .optical_terrain_ports import parameters
    shape,order,mode,scales=parameters(settings)
    return transport(shape,solar,rho,tau,soil,order,mode,scales)[1]['reflected']


def diffuse_albedo(settings,rho,tau,soil):
    mu,w,_=crown.quadrature(int(settings.get('optical_crown_quadrature',12)))
    return sum(2*m*weight*albedo(settings,float(np.degrees(np.arccos(m))),rho,tau,soil)
               for m,weight in zip(mu,w))
