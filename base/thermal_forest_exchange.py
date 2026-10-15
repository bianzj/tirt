"""Thermal sources on shared reciprocal whole-crown geometry."""
import math
import numpy as np
from .forest_exchange import (canopy_transmission, canopy_block, exchange_operator,
    source_sunlit, outer_width, higher_orders)

def scatter_components(model):
    lai,density,radius,height,hspot,bottom=map(float,model.forestshape[:6])
    shapes=tuple(tuple(map(float,row)) for row in model.shapes)
    matrix,areas=exchange_operator(shapes,lai,density,height,radius)
    q=np.asarray(model._reciprocal_first_hit_states,dtype=float)
    eps=np.array([model.Es,model.Em,model.Es,model.Em],dtype=float)
    reflectance=1-eps
    inner=hspot if model.thermal_hotspot_inner is None else float(model.thermal_hotspot_inner)
    outer=outer_width(radius,height,bottom,model.thermal_hotspot_outer)
    beta=math.atan2(shapes[0][0],shapes[0][1]);cosine=math.cos(beta)
    sun_slope=source_sunlit(lai*cosine,density*cosine,height,radius,hspot,bottom,inner,outer,
        float(model.sza),math.degrees(beta))
    sun_flat=source_sunlit(lai,density,height,radius,hspot,bottom,inner,outer,float(model.sza),0.)
    sources=np.zeros((4,4))
    for i,(soil_sun,leaf_sun) in enumerate((sun_slope,sun_flat)):
        sources[2*i,:2]=model.Es*np.array([soil_sun,1-soil_sun])
        sources[2*i+1,2:]=model.Em*np.array([leaf_sun,1-leaf_sun])
    transport=reflectance[:,None]*matrix
    scattered=q@higher_orders(matrix,reflectance,sources)
    if np.any(scattered<-1e-12) or not np.isfinite(scattered).all():
        raise ValueError('Invalid reciprocal thermal scattering.')
    model.thermal_exchange_budget=dict(exchange=matrix.copy(),areas=areas.copy(),
        first_hit=q.copy(),source=sources.copy(),transport=transport.copy(),
        scattered=scattered.copy())
    return scattered
