import numpy as np
import scipy.integrate as sci
from rt.gap import *
from rt.hotspot import *
from rt.scatter import *
import matplotlib.pyplot as plt
from rt.emissivity import *


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

