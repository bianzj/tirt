import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt

rd = np.pi/180.0
def hapke_function(angle,K,albedo):
    emis = 1-albedo
    term = 2*np.cos(angle*rd)/K
    emis_dir = np.sqrt(emis) *(1+term)/(1+term*np.sqrt(emis))
    return emis_dir

def P_distribution(angle,sigmas = 1.0):
    sigmas2 = sigmas * sigmas
    tan2 = np.tan(angle*rd) * np.tan(angle*rd)
    return 1.0/(2.0*np.pi*sigmas2) * np.exp(-tan2/(2.0*sigmas2))

def Hemisphere(targetza,targetaa):
    za = np.arange(0,90)
    aa = np.arange(0,360)
    za_,aa_ = np.meshgrid(za,aa)
    P_ = P_distribution(za_)
    cosangle = np.cos(za_*rd)*np.cos(targetza*rd)+\
               np.sin(za_*rd)*np.sin(targetza*rd)*np.cos((aa_-targetaa)*rd)
    eweight = 0
    esum = 0
    angle = np.arccos(cosangle)
    ind = angle > 0
    esum = np.sum(cosangle[ind]*P_[ind])
    eweight = np.sum(P_[ind])
    result = esum/eweight/np.cos(targetza*rd)
    result = np.sum(cosangle[ind]*1.0*rd*1.0*rd*P_[ind])
    result = result/2.0/np.pi/eweight
    return result

phi = 0.35
K = -np.log(1-1.209*np.power(phi,2/3.0))/(1.209*np.power(phi,2.0/3.0))
K = 1.0
vza = np.arange(0,90)
emis_dir = hapke_function(vza,K,0.05)
# plt.plot(vza,emis_dir)
# plt.show()



n = 1.0
Ra = 0.5
es = 0.9
R = 1.0/(1+(1.25*Ra)*(1.25*Ra)*np.pi*np.pi*n*n)
er = 1.0/(1.0+(1.0/es-1.0)*R)

print(es,er)

sigma = 1.0
l = 1.0
sigmas = np.sqrt(sigma)/l

vza = np.arange(0,90)
emis = np.asarray([Hemisphere(targetza,0) for targetza in vza])


Pd = P_distribution(vza,0.5)
plt.plot(vza,emis)
plt.show()


np.sum(np.cos(vza*rd))
