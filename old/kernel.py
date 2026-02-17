
from sklearn.linear_model import LinearRegression
from scipy.optimize import root,fsolve
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import fmin
import numpy as np

#########################################
#### constant
rd = np.pi/180.0

############################################
####  组分差异核
############################################


def kernel_LSF_analytical(vza):
    '''
    一种LSF核的参数化方案
    :param vza: 观测天顶角
    :return: LSF核值
    '''
    ### 这个参数化方案是JL.Roujean 提出的
    cthetv = np.cos(np.deg2rad(vza))
    Klsf = 0.069*cthetv*cthetv-0.215*cthetv+1.176
    return Klsf

def kernel_LSF(vza):
    '''
    LSF核的原始计算方法
    :param vza: 观测天顶角
    :return: LSF核
    '''
    ### 这个方案是Su. 提出的，对hapke的简化，考虑了上层和下层间的差异
    cthetv = np.cos(np.deg2rad(vza))
    Klsf = ((1 + 2 * cthetv) / (np.sqrt(0.96) + 2 * 0.96 * cthetv) - 0.25 * cthetv
            / (1 + 2 * cthetv) + 0.15 * (1 - np.exp(-0.75 / cthetv)))
    return Klsf

def kernel_Vin(vza):
    '''
    Vinnikov 提出的简单的考虑，只是cos函数
    :param vza:
    :return: Vin 核值
    '''
    Kvin = 1 - np.cos(np.deg2rad(vza))
    return Kvin

############################################################
##### 热点核
############################################################

def kernel_LiDense(vza,sza,raa, h = 2, r = 1, b = 1):
    '''

    :param vza:
    :param sza:
    :param raa:
    :return:
    '''
    ind = raa > 180
    raa[ind] = 360 - raa[ind]

    ### 椭球到球体的转换
    vza = np.arctan(b/r * np.tan(np.deg2rad(vza))) * 180 / np.pi
    sza = np.arctan(b / r * np.tan(np.deg2rad(sza))) * 180 / np.pi

    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    cphivs = np.cos(np.deg2rad(raa))
    sphivs = np.sin(np.deg2rad(raa))
    tgths = np.tan(np.deg2rad(sza))
    tgthv = np.tan(np.deg2rad(vza))

    calpha =cthets * cthetv + sthetv * sthets * cphivs
    DD = np.power(tgths, 2) + np.power(tgthv, 2) - 2 * tgths * tgthv * cphivs
    DD[DD < 0] = 0
    D = np.sqrt(DD)

    cost = h / b  *np.sqrt(D * D + np.power(tgthv * tgths * sphivs, 2)) / (1.0 / cthets + 1.0 / cthetv)
    cost[cost < -1] = -0.999
    cost[cost > 1] = 0.999

    t = np.arccos(cost)
    O = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * (1.0 / cthets + 1.0 / cthetv)
    Kli = (1 + calpha) * (1.0 / cthets * 1.0 / cthetv) / (1.0 / cthetv + 1 / cthets - O) - 2
    ### 热点消除
    Kli[sza >= 75] = 0
    return Kli


def kernel_LiSparse(vza, sza, raa, h = 2, r = 1, b = 1):
    '''

    :param vza:
    :param sza:
    :param raa:
    :return:
    '''
    ind = raa > 180
    raa[ind] = 360 - raa[ind]

    ### 椭球到球体的转换
    vza = np.arctan(b/r * np.tan(np.deg2rad(vza))) * 180 / np.pi
    sza = np.arctan(b / r * np.tan(np.deg2rad(sza))) * 180 / np.pi

    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    cphivs = np.cos(np.deg2rad(raa))
    sphivs = np.sin(np.deg2rad(raa))
    tgths = np.tan(np.deg2rad(sza))
    tgthv = np.tan(np.deg2rad(vza))

    DD = np.power(tgths, 2) + np.power(tgthv, 2) - 2 * tgths * tgthv * cphivs
    DD[DD < 0] = 0
    D = np.sqrt(DD)
    cost = h / b * np.sqrt(D * D + np.power(tgthv * tgths * sphivs, 2)) / (1.0 / cthets + 1.0 / cthetv)
    cost[cost < -1] = -0.999
    cost[cost > 1] = 0.999
    calpha = cthetv * cthets + sthetv * sthets * cphivs
    t = np.arccos(cost)
    O = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * (1.0 / cthets + 1.0 / cthetv)
    Kli = O - 1.0/cthetv - 1.0/cthets + 0.5*(1+calpha) *(1.0/cthets)*(1.0/cthetv)

    ### 热点消除
    Kli[sza >= 75] = 0

    return Kli

def kernel_RLf(vza,sza,raa):


    tgths = np.tan(np.deg2rad(sza))
    tgthv = np.tan(np.deg2rad(vza))
    cphivs = np.cos(np.deg2rad(raa))
    f = np.sqrt((np.power(tgthv, 2) + np.power(tgths, 2) - 2 * tgths * tgthv * cphivs))
    return f

####################################################
#### 多次散射项
####################################################



def kernel_RossThick(vza, sza, raa):
    '''
    Ross 简化模型，用于多次散射项计算
    :param sza:
    :param vza:
    :param raa:
    :return:
    '''
    ind = raa > 180
    raa[ind] = 360 - raa[ind]
    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    cphivs = np.cos(np.deg2rad(raa))


    calpha = cthetv * cthets + sthetv * sthets * cphivs
    alpha = np.arccos(calpha)
    salpha = np.sin(alpha)
    Kross = ((np.pi/2.0 - alpha)*calpha + salpha)/(cthetv+cthets)-np.pi/4.0


    return Kross




def kernel_RossThin():
    pass









#####################################################
#### 核驱动模型组合
####################################################

def kernel_LSFLiDense(vza,sza,raa):
    '''
    LSF组分温度差 和 Li 热点核
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''
    Klsf = kernel_LSF(vza)
    Kli = kernel_LiDense(vza,sza,raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Klsf,Kli,Kiso




def kernel_LSFLiDenseRossThick(vza,sza,raa):
    '''
    LSF组分温度差 和 Li 热点核
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''
    Klsf = kernel_LSF(vza)
    Kli = kernel_LiDense(vza,sza,raa)
    Kross = kernel_RossThick(vza,sza,raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Klsf,Kli,Kross,Kiso



def kernel_LSFLiSparse(vza,sza,raa):
    '''
    LSF 组分差异核 LiSparse 热点核
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''
    Klsf = kernel_LSF(vza)
    Kli = kernel_LiSparse(vza,sza,raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Klsf,Kli,Kiso

def kernel_RossThickLiSparse(vza,sza,raa):
    '''
    体散射核核几何光学核，分别为Ross和Li
    :param vza:
    :param sza:
    :param raa:
    :return: 体散射核核几何光学核
    '''

    Kvol = kernel_RossThick(vza,sza, raa)
    Kgeo = kernel_LiSparse(vza, sza, raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Kvol,Kgeo,Kiso

def kernel_RossThickLiDense(vza,sza,raa):
    '''
    体散射核核几何光学核，分别为Ross和Li
    :param vza:
    :param sza:
    :param raa:
    :return: 体散射核核几何光学核
    '''

    Kvol = kernel_RossThick(vza, sza, raa)
    Kgeo = kernel_LiDense(vza, sza, raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Kvol,Kgeo,Kiso

def kernel_VinLiDense(vza,sza,raa):
    '''
    体散射核核几何光学核，分别为Ross和Li
    :param vza:
    :param sza:
    :param raa:
    :return: 体散射核核几何光学核
    '''

    Kcom = kernel_Vin(vza)
    Kgeo = kernel_LiDense(vza, sza, raa)
    number = np.size(vza)
    Kiso = np.ones(number)
    return Kcom,Kgeo,Kiso
