import numpy as np
import scipy.integrate as sci
from gap import *
from hotspot import *
from utils import *
from scatter import *
import matplotlib.pyplot as plt
import quadpy
from emissivity import *


'''
冠层辐射
unit：四组分
individual: 六组分
endmember：八组分
layer：标识，标识其中的上层植被部分是分层计算的
'''


def radiance_direct(Ecom,Rcom):
    number_component = np.size(Rcom)
    number_angle = np.shape(Ecom[0])[0]
    Rcomnew = np.zeros([number_component,number_angle])
    for k in range(number_component):
        Rcomtemp = Rcom[k]
        Rcomtemp = np.transpose(np.tile(Rcomtemp,[number_angle,1]))
        Rcomnew[k,:] = np.sum(np.asarray(Ecom[k])*Rcomtemp,axis=0)
    # Rcomnew = np.asarray([Ecom[k,:]*Rcom[k] for k in range(number_component)])
    Rcomnew = np.sum(Rcomnew,axis=0)
    return Rcomnew

def radiance_scatter(Ecom,Rcom):
    Ecom =np.asarray(Ecom)
    Rcom = np.asarray(Rcom)
    number_component = np.size(Rcom)
    Rcomnew = np.asarray([Ecom[k,:]*Rcom[k] for k in range(number_component)])
    Rcomnew = np.sum(Rcomnew,axis=0)
    return Rcomnew




def radiance_scatter_hom_endmember_sail(Pcom,
                             lai_crown, lai_trunk, lai_unveg,
                             ec, et, ev, es, Tss,Tsh,Tcs,Tch,Tvs,Tvh,Tts,Tth,
                             vza0, sza0, vaa0,
                             saa0=0, alg_crown = 54, alg_unveg = 54):
    Pcom = np.asarray(Pcom)
    Pss = Pcom[0,:]
    Psh = Pcom[1,:]
    Pvs = Pcom[4,:]
    Pvh = Pcom[5,:]
    Pcs = Pcom[2,:]
    Pch = Pcom[3,:]
    Pts = Pcom[6,:]
    Pth = Pcom[7,:]
    rad = multiple_scattering_hom_sail_rad_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                                     lai_crown, lai_trunk, lai_unveg,
                                                     ec, et, ev, es, Tss, Tsh, Tcs, Tch, Tvs, Tvh, Tts, Tth,
                                                     vza0, sza0, vaa0)
    return rad

def radiance_scatter_crown_endmember_sail(Pcom,
                             lai_crown, lai_trunk, lai_unveg,std_crown,hcr_crown,rcr_crown,
                             ec, et, ev, es, Tss,Tsh,Tcs,Tch,Tvs,Tvh,Tts,Tth,
                             vza0, sza0, vaa0,
                             saa0=0, alg_crown = 54, alg_unveg = 54):

    Pcom = np.asarray(Pcom)
    Pss = Pcom[0,:]
    Psh = Pcom[1,:]
    Pvs = Pcom[4,:]
    Pvh = Pcom[5,:]
    Pcs = Pcom[2,:]
    Pch = Pcom[3,:]
    Pts = Pcom[6,:]
    Pth = Pcom[7,:]
    rad = multiple_scattering_crown_sail_rad_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                                       lai_crown, lai_trunk, lai_unveg, std_crown, hcr_crown, rcr_crown,
                                                       ec, et, ev, es, Tss, Tsh, Tcs, Tch, Tvs, Tvh, Tts, Tth,
                                                       vza0, sza0, vaa0)

    return rad

