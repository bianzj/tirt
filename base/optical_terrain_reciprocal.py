"""Optical trial using the shared four-state reciprocal forest geometry.

First scattering retains separate rho/tau spherical-leaf phase coefficients
and dual-hotspot visibility. Higher orders use an isotropic mean leaf albedo
rho+tau. Thermal emissivities are not substituted for optical reflectance.
"""
import math
from functools import lru_cache
import numpy as np
from .forest_exchange import exchange_operator,view_geometry,outer_width,higher_orders
from .terrain_rt import _shape_values


def parameters(settings):
    shape=tuple(float(settings[k]) for k in ('lai','stand_density','height','width','hspot'))
    if settings.get('optical_crown_geometry','clumping_dual')!='clumping_dual':
        raise ValueError('Reciprocal optical terrain requires clumping_dual.')
    bottom=float(settings.get('crown_bottom',shape[2]))
    inner=float(settings.get('optical_hotspot_inner',shape[4]))
    outer=outer_width(shape[3],shape[2],bottom,settings.get('optical_hotspot_outer'))
    shapes=_shape_values(settings)
    if len(shapes)!=1 or min(shapes[0][:2])<=0 or shapes[0][2]<=0:
        raise ValueError('Reciprocal optical trial requires one populated cone class.')
    if min(shape[:2])<0 or not np.isfinite(shape).all() or inner<0 or not math.isfinite(inner):
        raise ValueError('Invalid forest structure or hotspot scale.')
    return shapes,shape,bottom,inner,outer


@lru_cache(maxsize=1024)
def state(shapes,shape,bottom,inner,outer,solar,rho,tau,soil):
    matrix,areas=exchange_operator(shapes,*shape[:4])
    q_s=np.asarray(view_geometry(shapes,shape,bottom,inner,outer,solar,solar,0.)[0])
    albedo=np.array([soil,rho+tau,soil,rho+tau])
    source=np.divide(albedo*q_s,areas,out=np.zeros(4),where=areas>0)
    higher=higher_orders(matrix,albedo,source)
    total=source+higher
    escaped=float(areas@((1-matrix.sum(axis=1))*total))
    absorbed=float((1-albedo)@q_s+areas@((1-albedo)*(matrix@total)))
    if min(escaped,absorbed)<-1e-10 or abs(escaped+absorbed-1)>1e-9:
        raise ValueError('Reciprocal optical state flux budget failed.')
    return higher,dict(incident=1.,escaped=escaped,absorbed=absorbed,
        residual=escaped+absorbed-1,state_first_escape=float(areas@((1-matrix.sum(axis=1))*source)))


def reflectance(settings,sza,saa,vza,vaa,rho,tau,soil,diffuse_fraction=0.,sunlit=True,return_budget=False):
    scalars=[sza,saa,vza,vaa,rho,tau,soil,diffuse_fraction]
    if (not np.isfinite(scalars).all() or not 0<=sza<90 or not 0<=vza<90 or
            min(rho,tau,soil,diffuse_fraction)<0 or rho+tau>1+1e-12 or max(soil,diffuse_fraction)>1):
        raise ValueError('Invalid optical direction or material coefficients.')
    shapes,shape,bottom,inner,outer=parameters(settings)
    relative=abs((vaa-saa+180)%360-180)
    q,phase,_=view_geometry(shapes,shape,bottom,inner,outer,float(sza),float(vza),float(relative))
    q=np.asarray(q)
    first=float(np.dot(phase,[rho,tau,soil]))
    higher,budget=state(shapes,shape,bottom,inner,outer,float(sza),float(rho),float(tau),float(soil))
    beam=first+float(q@higher) if sunlit else 0.
    f=diffuse_fraction
    sky=0.
    sky_budget=dict(incident=0.,escaped=0.,absorbed=0.)
    if f:
        matrix,areas=exchange_operator(shapes,*shape[:4])
        albedo=np.array([soil,rho+tau,soil,rho+tau])
        source=albedo*(1-matrix.sum(axis=1))
        total=source+higher_orders(matrix,albedo,source)
        sky=float(q@total)
        boundary=1-matrix.sum(axis=1)
        sky_budget=dict(incident=float(areas@boundary),escaped=float(areas@(boundary*total)),
            absorbed=float(areas@((1-albedo)*(boundary+matrix@total))))
    result=(1-f)*beam+f*sky
    if not return_budget:
        return result
    budget={key:(1-f)*float(sunlit)*budget[key]+f*sky_budget[key] for key in ('incident','escaped','absorbed')}
    budget['residual']=budget['escaped']+budget['absorbed']-budget['incident']
    budget=dict(budget,directional_first=first,directional_higher=float(q@higher),
        flux_budget_definition='state radiosity; angular first kernel is a separate statistical moment',
        higher_orders_leaf_phase='isotropic_rho_plus_tau',geometry_shared_with_thermal=True)
    return result,budget
