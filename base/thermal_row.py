"""Thermal row geometry conditioned on the joint empty corridor.

Only the four direct-emission proportions change. Material emissivities,
component temperatures and the established thermal scattering remain intact.
"""
import math
import numpy as np

from .gap import gap_probability_row_analytical
from .optical_row import clear_corridor
from .vnir import joint_transmission


def conditioned_proportions(lai,hspot,width,blank,height,vza,vaa,sza,saa,row_azimuth,legacy):
    """Preserve row view marginals while correcting sun/shade correlation.

    The empty phase stays sunlit soil. In the occupied phase, infer optical
    depth and view extinction from the conditional solar and view gaps.
    Joint transmission gives soil sun fraction and its path integral gives
    visible sunlit foliage. Along rows this is the exact empty+slab limit.
    """
    if width<=0 or blank<0 or height<=0 or lai<0 or hspot<0:
        raise ValueError('Invalid thermal row geometry.')
    v,a,s,b,az = [x.reshape(-1) for x in np.broadcast_arrays(
        np.asarray(vza,float),np.asarray(vaa,float),np.asarray(sza,float),
        np.asarray(saa,float),np.asarray(row_azimuth,float))]
    if not np.isfinite([lai,hspot,width,blank,height]).all() or not np.isfinite([v,a,s,b,az]).all():
        raise ValueError('Thermal row inputs must be finite.')
    result = np.array(legacy,dtype=float,copy=True).reshape(4,-1)
    if result.shape[1]!=len(v):
        raise ValueError('Thermal row proportions do not match angles.')
    if blank==0 or lai==0:
        return result
    # Match the existing geometry's one-sided horizon convention.
    v = np.clip(v,0.,89.999); s = np.clip(s,0.,89.999)
    solar = gap_probability_row_analytical(lai,width,blank,height,s,abs(b-az)%180.)
    for j in range(len(v)):
        clear = clear_corridor(width,blank,height,s[j],b[j],v[j],a[j],az[j])
        if clear<=1e-14:
            continue
        cover = 1-clear
        view_gap = float(result[0,j]+result[1,j])
        gs = float(np.clip((solar[j]-clear)/cover,1e-12,1.))
        gv = float(np.clip((view_gap-clear)/cover,1e-12,1.))
        ks = .5/math.cos(math.radians(s[j]))
        depth = max(-math.log(gs)/ks,1e-10)
        ko = -math.log(gv)/depth
        raa = abs((a[j]-b[j]+180.)%360.-180.)
        joint,path = joint_transmission(depth,ks,ko,hspot,float(s[j]),float(v[j]),float(raa))
        soil_sun = float(np.clip(clear+cover*joint,0.,view_gap))
        leaf_sun = float(np.clip(cover*ko*path,0.,1-view_gap))
        conditioned = np.array([soil_sun,view_gap-soil_sun,leaf_sun,1-view_gap-leaf_sun])
        # Same geometrical corridor survival used in the optical adapter.
        # This returns the old closure exactly when the corridor is closed.
        weight = min(1.,clear*(width+blank)/blank)
        result[:,j] = (1-weight)*result[:,j]+weight*conditioned
    return result
