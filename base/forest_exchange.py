"""Shared reciprocal whole-crown geometry for optical and thermal transport.

States: slope soil/leaves, flat soil/leaves. Local canopy returns and escape
to other terrain patches share one probability budget. No crown cubature,
reference BRF, temperature fitting or Monte Carlo is used by this closure.
"""
import math
from functools import lru_cache

import numpy as np

from .gap import gap_probability_crown_analytical
from .thermal_crown import component_fractions
from .solar import slope_direction
from .optical_terrain import exchange_geometry


def outer_width(radius,height,bottom,override=None):
    """One crown-centre height rule for both spectral adapters."""
    if not np.isfinite([radius,height,bottom]).all() or min(radius,height)<=0 or bottom<0:
        raise ValueError('Positive crown semiaxes and nonnegative crown bottom required.')
    width=radius/(height+bottom) if override is None else float(override)
    if not math.isfinite(width) or width<0:
        raise ValueError('Hotspot width must be finite and nonnegative.')
    return width


def higher_orders(matrix,albedo,source):
    """Shared closure; material coefficients and sources stay separate."""
    matrix=np.asarray(matrix,dtype=float)
    albedo=np.asarray(albedo,dtype=float)
    source=np.asarray(source,dtype=float)
    if (matrix.shape!=(4,4) or albedo.shape!=(4,) or source.shape[0]!=4 or
            not np.isfinite(matrix).all() or not np.isfinite(albedo).all() or
            not np.isfinite(source).all() or min(matrix.min(),albedo.min(),source.min())<0 or
            albedo.max()>1 or matrix.sum(axis=1).max()>1+1e-10):
        raise ValueError('Invalid whole-crown exchange coefficients.')
    transport=albedo[:,None]*matrix
    return np.linalg.solve(np.eye(4)-transport,transport@source)


@lru_cache(maxsize=4096)
def view_geometry(shapes,shape,bottom,inner,outer,solar,view,relative):
    """Extract common first-hit probabilities and optical phase kernels.

    The black geometry probe does not evaluate temperature or thermal scatter.
    """
    from .terrain_forest import TerrainForest
    lai,density,height,radius,hspot=shape
    model=TerrainForest()
    model.set_hotspot_geometry('clumping_dual',inner,outer)
    model.set_radiative_transport('reciprocal_exchange')
    model.forestshape=np.array([lai,density,radius,height,hspot,bottom])
    model.set_structural_input(shapes)
    model.set_spectral_input(1.,1.)
    model._collect_optical_geometry=True
    model.set_angular_input(np.array([max(view,.001)]),np.array([relative if view else 0.]),solar,0.)
    components=model.calculate_component_direct_emissivity(0,2)
    return (tuple(model._reciprocal_first_hit_states),tuple(model._optical_first_coefficients),
            tuple(map(float,components)))


@lru_cache(maxsize=256)
def canopy_transmission(lai,density,height,radius):
    if lai<=0 or density<=0:
        return 1.
    x,w=np.polynomial.legendre.leggauss(8)
    mu=(x+1)/2;weights=w/2
    gap=gap_probability_crown_analytical(lai,density,height,radius,np.degrees(np.arccos(mu)))
    return float(np.dot(2*mu*weights,gap))


def canopy_block(lai,transmission):
    """Area reciprocity: A_leaf/A_ground=2*LAI for two-sided leaf emission."""
    if lai<=0 or transmission>=1.:
        return np.zeros((2,2)),np.array([1.,0.]),np.array([1.,0.])
    escape=(1-transmission)/(2*lai)
    if not 0<=escape<=.5+1e-10:
        raise ValueError('Invalid reciprocal canopy escape probability.')
    local=np.array([[0.,1-transmission],[escape,1-2*escape]])
    return local,np.array([transmission,escape]),np.array([transmission,1-transmission])


@lru_cache(maxsize=128)
def exchange_operator(shapes,lai,density,height,radius):
    """Allocate local and macro routes once, with area reciprocity."""
    if len(shapes)!=1:
        raise ValueError('Reciprocal thermal forest currently requires one cone class.')
    h,r,n=shapes[0]
    cosine=r/math.hypot(h,r)
    coverage=n*math.pi*r*r
    areas=np.array([coverage/cosine,1-coverage])
    local_lai=np.array([lai*cosine,lai])
    transmission=[canopy_transmission(local_lai[i],density*(cosine if i==0 else 1),height,radius)
                  for i in range(2)]
    blocks=[canopy_block(local_lai[i],transmission[i]) for i in range(2)]
    macro=exchange_geometry(shapes)
    matrix=np.zeros((4,4))
    for i,(local,escape,_) in enumerate(blocks):
        matrix[2*i:2*i+2,2*i:2*i+2]=local
        for j,(_,_,entry) in enumerate(blocks):
            matrix[2*i:2*i+2,2*j:2*j+2]+=macro[i,j]*np.outer(escape,entry)
    state_areas=np.repeat(areas,2)*np.ravel(np.column_stack((np.ones(2),2*local_lai)))
    flux=state_areas[:,None]*matrix
    if (matrix.min()<-1e-12 or matrix.sum(axis=1).max()>1+1e-10 or
            not np.allclose(flux,flux.T,atol=1e-10,rtol=1e-10)):
        raise ValueError('Thermal exchange violates probability or area reciprocity.')
    return matrix,state_areas


@lru_cache(maxsize=512)
def source_sunlit(lai,density,height,radius,hspot,bottom,inner,outer,sza,slope):
    """Small fixed angular average for source populations, independent of view."""
    x,w=np.polynomial.legendre.leggauss(4)
    mu=(x+1)/2;weights=w/2
    azimuth=np.arange(0.,360.,60.)
    views,relative=np.meshgrid(np.degrees(np.arccos(mu)),np.minimum(azimuth,360-azimuth),indexing='ij')
    wt=np.repeat(2*mu*weights,len(azimuth))/len(azimuth)
    settings=dict(lai=lai,stand_density=density,height=height,width=radius,hspot=hspot,
        thermal_crown_geometry='clumping_dual',thermal_hotspot_inner=inner,thermal_hotspot_outer=outer)
    mean=np.zeros(4)
    aspects=np.arange(0.,360.,15.) if slope else [0.]
    for aspect in aspects:
        sun,_,lit=slope_direction(sza,0.,slope,aspect) if slope else (sza,0.,True)
        p=component_fractions(settings,min(sun,89.9),views.ravel(),relative.ravel(),sunlit=lit)
        mean+=p@wt/len(aspects)
    soil=mean[0]/max(float(mean[:2].sum()),1e-15)
    leaf=mean[2]/max(float(mean[2:].sum()),1e-15)
    return float(soil),float(leaf)
