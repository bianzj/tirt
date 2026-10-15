"""Two-scale PARAMETERIZED optical crown geometry; no spatial quadrature.

Other crowns are a Poisson process of porous projected disks with the same
areas and single-direction gaps as gap_probability_crown_analytical. Their
overlap supplies the outer hotspot. A separate leaf-scale correlation
supplies the joint transmission inside a shared crown.

For foliage, use a closed-form uniform optical-depth average in the host
crown and half of the stand optical depth for other crowns. This is a new
conditional-mean closure, not an exact average of spatially varying joint
transmissions. The equivalent disks and half-depth closure are approximations.
Shared thermal geometry/hotspot routines are deliberately not changed.
"""
import math

from .gap import gap_probability_crown_analytical


def attenuation_mean(depth):
    """Integral of exp(-depth*u), u=0..1, evaluated without quadrature."""
    return -math.expm1(-depth)/depth if depth > 1e-8 else 1-depth/2+depth*depth/6


def inner_hotspot(depth_v, depth_s, separation, width, correction):
    """Joint optical depth inside a shared crown; width is dimensionless.

    Geometric-mean overlap follows the legacy path hotspot convention.
    Bounds enforce max(depth_v,depth_s) <= joint <= depth_v+depth_s.
    """
    correlation = (1. if separation == 0 else 0.) if width == 0 else attenuation_mean(separation*correction/width)
    overlap = min(depth_v, depth_s, math.sqrt(depth_v*depth_s)*correlation)
    return depth_v+depth_s-overlap


def disk_overlap(area_v, area_s, distance):
    """Exact overlap of equal-area surrogate disks, including tangent limits."""
    rv,rs=math.sqrt(area_v/math.pi),math.sqrt(area_s/math.pi)
    if distance >= rv+rs:return 0.
    if distance <= abs(rv-rs):return min(area_v,area_s)
    a=max(-1.,min(1.,(distance*distance+rv*rv-rs*rs)/(2*distance*rv)))
    b=max(-1.,min(1.,(distance*distance+rs*rs-rv*rv)/(2*distance*rs)))
    radicand=(-distance+rv+rs)*(distance+rv-rs)*(distance-rv+rs)*(distance+rv+rs)
    value=rv*rv*math.acos(a)+rs*rs*math.acos(b)-.5*math.sqrt(max(0.,radicand))
    return max(0.,min(min(area_v,area_s),value))


def outer_hotspot(density, area_v, area_s, depth_v, depth_s, depth_joint,
                  separation, radius, width):
    """Joint stand optical depth from crown overlap and porous interception.

    Width controls projected crown displacement; the geometric default
    radius/vertical_semiaxis gives displacement=vertical_semiaxis*separation.
    Internal transmission enters the shared area once, avoiding double counting.
    """
    distance = (0. if separation == 0 else math.inf) if width == 0 else radius*separation/width
    common=disk_overlap(area_v,area_s,distance)
    iv,is_,ivs=-math.expm1(-depth_v),-math.expm1(-depth_s),-math.expm1(-depth_joint)
    dv,ds=density*area_v*iv,density*area_s*is_
    joint=density*((area_v-common)*iv+(area_s-common)*is_+common*ivs)
    return max(dv,ds,min(dv+ds,joint)),dv,ds


def dual_hotspot_fractions(lai,density,hcr,rcr,hspot,sza,vza,raa,scales=None):
    """Sunlit/shaded soil and leaves, with unchanged visible ground/leaf totals.

    scales=(inner_width,outer_width). Both are dimensionless angular scales;
    defaults are the input hspot and the crown aspect ratio rcr/hcr. No
    reference reflectance or fitted coefficient enters these defaults.
    """
    inner,outer=(hspot,rcr/hcr) if scales is None else scales
    if not all(math.isfinite(x) and x >= 0 for x in (inner,outer)):
        raise ValueError('Optical inner/outer hotspot widths must be finite and nonnegative.')
    sv,vv,az=map(math.radians,(sza,vza,raa))
    ms,mv=math.cos(sv),math.cos(vv)
    ts,tv=math.tan(sv),math.tan(vv)
    separation=math.sqrt(max(0.,(tv-ts)**2+4*tv*ts*math.sin(az/2)**2))
    av=math.pi*rcr*math.sqrt(rcr*rcr+(hcr*tv)**2)
    ass=math.pi*rcr*math.sqrt(rcr*rcr+(hcr*ts)**2)
    kv,ks=.5*lai/(density*av*mv),.5*lai/(density*ass*ms)
    correction=2/(.5/mv+.5/ms)
    joint_inner=inner_hotspot(kv,ks,separation,inner,correction)
    joint_outer,dv,_=outer_hotspot(density,av,ass,kv,ks,joint_inner,separation,rcr,outer)
    gv=float(gap_probability_crown_analytical(lai,density,hcr,rcr,vza))
    gs=float(gap_probability_crown_analytical(lai,density,hcr,rcr,sza))
    soil_sun=max(max(0.,gv+gs-1),min(gv,gs,math.exp(-joint_outer)))
    # Analytic average through the host crown, conditional on leaf visibility.
    inside=attenuation_mean(joint_inner)/attenuation_mean(kv)
    # Half-stand-depth approximation for the other crowns above a typical leaf.
    outside=math.exp(-.5*max(0.,joint_outer-dv))
    leaf_sun=(1-gv)*max(0.,min(1.,inside*outside))
    return soil_sun,gv-soil_sun,leaf_sun,1-gv-leaf_sun
