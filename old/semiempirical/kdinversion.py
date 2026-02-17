
'''
核驱动模拟拟合
'''
from scipy.optimize import root,fsolve
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import fmin
from old.kernel import *

'''
进行拟合或反演的时候可以采用线性或非线性的方式
可以是多角度或者多像元的方式进行
一般情况：行代表了不同的角度，列代表不同的像元
xxx 为矩阵
xxx0 为参考矩阵
这里包括2中:
1种是直接，没有后缀
1种是差异，后缀diff
其中直接的多用于可见光、差异的多用于热红外

Kvol:  体散射核
Kgeo: 几何光学核
Kcom: 组分差异核
Kiso: 各向同性核

'''


#########################################
### constant

### 夜晚阈值，太阳天顶角大于该值认为是夜晚，没有热点效应
sza_night_threshold = 75.0
rd = np.pi/180.0

################################
###  核驱动模拟非线性拟合用的函数
################################

### (Ross/ LSF /Vin) * (LiDense / LiSparse / RL)

def fun_RossThickLiDense(x,y,vza,sza,raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kvol,Kgeo,Kiso = kernel_RossThickLiDense(vza,sza,raa)
    sim = a * Kvol + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

def fun_RossThickLiSparse(x,y,vza,sza,raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kvol,Kgeo,Kiso = kernel_RossThickLiSparse(vza,sza,raa)
    sim = a * Kvol + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

def fun_RossThickRL(x,y,vza,sza,raa):
    '''
    反演VinRLE方法的b和k值，其中a值是通过Vin方法直接反演得到的
    :param x: 待拟合的核系数
    :param lst:  温度
    :param a: 系数a
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 相对方位角
    :param rad: 下行长波辐射、短波辐射、下行辐射
    :return:
    '''
    a = x[0]
    b = x[1]
    c = x[2]
    k = x[3]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    Kvol = kernel_RossThick(vza,sza,raa)
    Kgeo = np.exp(-k*f)

    sim = a * Kvol + b  * Kgeo + c
    mea = y
    res = sim - mea
    return np.sum(res*res)

def fun_LSFLiSparse(x, y, vza, sza, raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kcom, Kgeo, Kiso = kernel_LSFLiSparse(vza, sza, raa)
    sim = a * Kcom + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

def fun_LSFLiDense(x, y, vza, sza, raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kcom, Kgeo, Kiso = kernel_LSFLiDense(vza, sza, raa)
    sim = a * Kcom + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

def fun_LSFLiDenseRossThick(x, y, vza, sza, raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]

    Kcom, Kgeo, Kvol, Kiso = kernel_LSFLiDenseRossThick(vza, sza, raa)
    sim = a * Kcom + b * Kgeo + c * Kiso + d * Kvol
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)


def fun_LSFRL(x,y,vza,sza,raa):
    '''
    反演VinRLE方法的b和k值，其中a值是通过Vin方法直接反演得到的
    :param x: 待拟合的核系数
    :param lst:  温度
    :param a: 系数a
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 相对方位角
    :param rad: 下行长波辐射、短波辐射、下行辐射
    :return:
    '''
    a = x[0]
    b = x[1]
    c = x[2]
    k = x[3]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    phi = kernel_LSF(vza)


    sim = a * phi + b  * np.exp(-k*f) + c
    mea = y
    res = sim - mea
    return np.sum(res*res)

def fun_VinRL(x,y,vza,sza,raa):
    '''
    反演VinRLE方法的b和k值，其中a值是通过Vin方法直接反演得到的
    :param x: 待拟合的核系数
    :param lst:  温度
    :param a: 系数a
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 相对方位角
    :param rad: 下行长波辐射、短波辐射、下行辐射
    :return:
    '''
    a = x[0]
    b = x[1]
    c = x[2]
    k = x[3]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    phi = kernel_Vin(vza)

    sim = a * phi + b  * np.exp(-k*f) + c
    mea = y
    res = sim - mea
    return np.sum(res*res)

def fun_VinLiDense(x,y,vza,sza,raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kvol,Kgeo,Kiso = kernel_VinLiDense(vza,sza,raa)
    sim = a * Kvol + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

def fun_VinLiSparse(x,y,vza,sza,raa):
    '''
    用于可见光和近红外波段观测的拟合，体散射核核几何光学核
    :param x:
    :param refl:
    :param vza:
    :param sza:
    :param raa:
    :return:
    '''

    a = x[0]
    b = x[1]
    c = x[2]

    Kvol,Kgeo,Kiso = kernel_RossThickLiSparse(vza,sza,raa)
    sim = a * Kvol + b * Kgeo + c * Kiso
    mea = y
    ind = sza < sza_night_threshold
    res = sim[ind] - mea[ind]

    return np.sum(res * res)

### (Ross/ LSF) * (LiDense / LiSparse / RL)
### RL 及其变体

def fun_Vin_dif(x,y,vza):
    '''
     Vinnikov 的直接拟合，多用于夜晚的Vin系数模拟,这里
     lst和vza分别是要处理的观测和输入；
     lst0和vza0是从观测中选中的某个观测或输入；
    :param x:
    :param lst:
    :param vza:
    :param vza0:
    :return:
    '''

    y0 = y[0,:]
    vza0 = vza[0,:]
    y = y[1:,:]
    vza = vza[1:,:]
    phi = kernel_Vin(vza)
    phi0 = kernel_Vin(vza0)
    res = (phi*y0 - phi0*y)*x

    return np.sum(res**2)

def fun_LSF_dif(x,y,vza):
    '''
     Vinnikov 的直接拟合，多用于夜晚的Vin系数模拟,这里
     lst和vza分别是要处理的观测和输入；
     lst0和vza0是从观测中选中的某个观测或输入；
    :param x:
    :param lst:
    :param vza:
    :param vza0:
    :return:
    '''

    y0 = y[0,:]
    vza0 = vza[0,:]
    y = y[1:,:]
    vza = vza[1:,:]
    phi = kernel_LSF(vza)
    phi0 = kernel_LSF(vza0)
    res = (phi*y0 - phi0*y)*x

    return np.sum(res**2)

def fun_VinRLE_bk_dif(x,y,a,vza,sza,raa,var):
    '''
    反演VinRLE方法的b和k值，其中a值是通过Vin方法直接反演得到的
    :param x: 待拟合的核系数
    :param lst:  温度
    :param a: 系数a
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 相对方位角
    :param rad: 下行长波辐射、短波辐射、下行辐射
    :return:
    '''

    b = x[0]
    k = x[1]


    y0 = y[:, 0]
    y = y[:, 1:]
    vza0 = vza[:, 0]
    vza = vza[:, 1:]
    sza = sza[:,1:]
    raa = raa[:,1:]
    sza0 = sza[:,0]
    raa0 = raa[:,0]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    fn = np.tan(np.deg2rad(sza))
    phi = kernel_Vin(vza)
    f0 = kernel_RLf(vza0,sza0,raa0)
    # fn0 = np.tan(np.deg2rad(sza0))
    phi0 =kernel_Vin(vza0)


    sim = a * (phi*y0 - phi0*y0) + b * var *(np.exp(-k*f)-np.exp(-k*f0))/(1-np.exp(-k*fn))
    mea = y0 -y
    res = sim - mea
    return np.sum(res*res)

def fun_VinRLE_abk_dif(x,y,vza,sza,raa,var):

    a = x[0]
    b = x[1]
    k = x[2]

    y0 = y[:, 0]
    y = y[:, 1:]
    vza0 = vza[:, 0]
    vza = vza[:, 1:]
    sza = sza[:, 1:]
    raa = raa[:, 1:]
    sza0 = sza[:, 0]
    raa0 = raa[:, 0]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza, sza, raa)
    fn = np.tan(np.deg2rad(sza))
    phi = kernel_Vin(vza)
    f0 = kernel_RLf(vza0, sza0, raa0)
    # fn0 = np.tan(np.deg2rad(sza0))
    phi0 = kernel_Vin(vza0)

    sim = a * (phi * y0 - phi0 * y0) + b * var * (np.exp(-k * f) - np.exp(-k * f0)) / (1 - np.exp(-k * fn))
    mea = y0 - y
    ind = sza > 80
    sim[ind] = a * (phi * y0 - phi0 * y0)
    res = sim - mea

    return np.sum(res*res)

def fun_LSFRLE_bk_dif(x,y,a,vza,sza,raa,var):

    b = x[0]
    k = x[1]


    y0 = y[:, 0]
    y = y[:, 1:]
    vza0 = vza[:, 0]
    vza = vza[:, 1:]
    sza = sza[:,1:]
    raa = raa[:,1:]
    sza0 = sza[:,0]
    raa0 = raa[:,0]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    fn = np.tan(np.deg2rad(sza))
    phi = kernel_LSF_analytical(vza)
    f0 = kernel_RLf(vza0,sza0,raa0)
    # fn0 = np.tan(np.deg2rad(sza0))
    phi0 =kernel_LSF_analytical(vza0)


    sim = a * (phi*y0 - phi0*y0) + b * var *(np.exp(-k*f)-np.exp(-k*f0))/(1-np.exp(-k*fn))
    mea = y0 -y
    res = sim - mea
    return np.sum(res*res)

def fun_LSFRLE_abk_dif(x,y,vza,sza,raa,var):

    a = x[0]
    b = x[1]
    k = x[2]


    y0 = y[:, 0]
    y = y[:, 1:]
    vza0 = vza[:, 0]
    vza = vza[:, 1:]
    sza = sza[:,1:]
    raa = raa[:,1:]
    sza0 = sza[:,0]
    raa0 = raa[:,0]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    fn = np.tan(np.deg2rad(sza))
    phi = kernel_LSF_analytical(vza)
    f0 = kernel_RLf(vza0,sza0,raa0)
    # fn0 = np.tan(np.deg2rad(sza0))
    phi0 =kernel_LSF_analytical(vza0)


    sim = a * (phi*y0 - phi0*y0) + b * var *(np.exp(-k*f)-np.exp(-k*f0))/(1-np.exp(-k*fn))
    mea = y0 -y
    res = sim - mea
    return np.sum(res*res)

def fun_RossThickRLE_abk_dif(x,y,vza,sza,raa,var):

    a = x[0]
    b = x[1]
    k = x[2]


    y0 = y[:, 0]
    y = y[:, 1:]
    vza0 = vza[:, 0]
    vza = vza[:, 1:]
    sza = sza[:,1:]
    raa = raa[:,1:]
    sza0 = sza[:,0]
    raa0 = raa[:,0]

    ### xxx0是作为参考，xxx是变量
    f = kernel_RLf(vza,sza,raa)
    fn = np.tan(np.deg2rad(sza))
    phi = kernel_RossThick(vza,sza,raa)
    f0 = kernel_RLf(vza0,sza0,raa0)
    # fn0 = np.tan(np.deg2rad(sza0))
    phi0 =kernel_RossThick(vza0,sza0,raa0)


    sim = a * (phi*y0 - phi0*y0) + b * var *(np.exp(-k*f)-np.exp(-k*f0))/(1-np.exp(-k*fn))
    mea = y0 -y
    res = sim - mea
    return np.sum(res*res)

#############################################
### fitting USING FMIN
###############################################

#### DIRECT VNIR
def fitting_RossThickLiSparse(y,vza,sza,raa,coeffa = 0.5,coeffb = 0.5, coeffc = 300):
    res = fmin(fun_RossThickLiSparse,x0 = [coeffa,coeffb,coeffc],args=(y,vza,sza,raa),disp=0)
    return res

def fitting_RossThickLiDense(y, vza, sza, raa,coeffa = 0.5,coeffb = 0.5, coeffc = 300):
    res = fmin(fun_RossThickLiDense, x0=[coeffa, coeffb, coeffc], args=(y, vza, sza, raa), disp=0)
    return res

def fitting_RossThickRL(y, vza, sza, raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300, coeffk = 0.5):
    res = fmin(fun_RossThickRL, x0 = [coeffa, coeffb, coeffc,coeffk], args=(y, vza, sza, raa), disp=0)
    return res

def fitting_LSFLiSparse(y,vza,sza,raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300):
    res = fmin(fun_LSFLiSparse,x0 = [coeffa, coeffb, coeffc],args=(y,vza,sza,raa),disp=0)
    return res

def fitting_LSFLiDenseRossThick(y,vza,sza,raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300,coeffd = 0.5):
    res = fmin(fun_LSFLiDenseRossThick,x0 = [coeffa, coeffb, coeffc, coeffd],args=(y,vza,sza,raa),disp=0)
    return res

def fitting_LSFLiDense(y, vza, sza, raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300):
    res = fmin(fun_LSFLiDense, x0=[coeffa, coeffb, coeffc], args=(y, vza, sza, raa), disp=0)
    return res

def fitting_LSFRL(y,vza,sza,raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300, coeffk = 0.5):
    res = fmin(fun_LSFRL, x0=[coeffa, coeffb, coeffc,coeffk], args=(y, vza, sza, raa), disp=0)
    return  res

def fitting_VinLiDense(y,vza,sza,raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300):
    res = fmin(fun_VinLiSparse,x0 = [coeffa, coeffb, coeffc],args=(y,vza,sza,raa),disp=0)
    return res

def fitting_VinRL(y,vza,sza,raa,coeffa = 0.5, coeffb = 0.5, coeffc = 300, coeffk = 0.5):
    res = fmin(fun_VinRL,x0 = [coeffa, coeffb, coeffc,coeffk],args=(y,vza,sza,raa),disp=0)
    return  res

### DIFFERENCE TIR

def fitting_Vin_dif(y,vza,coeffa = 0.05):
    res = fmin(fun_Vin_dif, x0 = [coeffa],args = (y,vza), disp = 0)
    return res

def fitting_LSF_dif(y,vza,coeffa = 0.05):
    res = fmin(fun_LSF_dif, x0 = [coeffa], args = (y,vza), disp = 0)
    return res

def fitting_VinRLE_bk_dif(y, a, vza,sza,raa, var,coeffb = 0.5, coeffk = 0.5):
    res = fmin(fun_VinRLE_bk_dif, x0 = [coeffb,coeffk],args = (y,a,vza,sza,raa,var),disp = 0)
    return res

def fitting_VinRLE_abk_dif(y, vza,sza,raa,var,coeffa = 0.5, coeffb = 0.5, coeffk = 0.5):
    res = fmin(fun_VinRLE_abk_dif, x0 = [coeffa,coeffb,coeffk], args = (y,vza,sza,raa,var),disp = 0)
    return res

def fitting_LSFRLE_bk_dif(y, a, vza,sza,raa, var,coeffb = 0.5, coeffk = 0.5):
    res = fmin(fun_LSFRLE_bk_dif, x0 = [coeffb,coeffk],args = (y,a,vza,sza,raa,var),disp = 0)
    return res

def fitting_LSFRLE_abk_dif(y, vza,sza,raa,var,coeffa = 0.5, coeffb = 0.5, coeffk = 0.5):
    res = fmin(fun_LSFRLE_abk_dif, x0 = [coeffa,coeffb,coeffk], args = (y,vza,sza,raa,var),disp = 0)
    return res

def fitting_RossThickRLE_abk_dif(y, vza,sza,raa,var, coeffa = 0.5, coeffb = 0.5, coeffk = 0.5):
    res = fmin(fun_RossThickRLE_abk_dif, x0 = [coeffa,coeffb,coeffk], args = (y,vza,sza,raa,var),disp = 0)
    return res


###################################################
#### fitting USING MINIMIZE
###################################################

def fitting_VinLiSparse_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_VinLiSparse,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x


def fitting_VinLiDense_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_VinLiDense,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x


def fitting_RossThickLiDense_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_RossThickLiDense,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x

def fitting_RossThickLiSparse_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_RossThickLiSparse,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x

def fitting_LSFLiDense_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 275,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_LSFLiDense,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x

def fitting_LSFLiSparse_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 275,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350) ):
    x0 = np.asarray([coeffa,coeffb,coeffc])
    bnds = (bnda,bndb,bndc)
    res = minimize(fun_LSFLiSparse,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x


def fitting_LSFRL_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,coeffk = 0.5,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350),bndk = (-100,100) ):
    x0 = np.asarray([coeffa,coeffb,coeffc,coeffk])
    bnds = (bnda,bndb,bndc,bndk)
    res = minimize(fun_LSFRL,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x

def fitting_RossThickRL_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,coeffk = 0.5,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350),bndk = (-100,100) ):
    x0 = np.asarray([coeffa,coeffb,coeffc,coeffk])
    bnds = (bnda,bndb,bndc,bndk)
    res = minimize(fun_RossThickRL,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x

def fitting_VinRL_min(y,vza,sza,raa, method = 'powell',
                        coeffa = 0.5,coeffb = 0.5,coeffc = 300,coeffk = 0.5,
                        bnda =(-100,100),bndb = (-100,100),bndc = (0,350),bndk = (-100,100) ):
    x0 = np.asarray([coeffa,coeffb,coeffc,coeffk])
    bnds = (bnda,bndb,bndc,bndk)
    res = minimize(fun_VinRL,x0,
                   args=(y,vza,sza,raa),
                   method = method,
                   bounds = bnds)

    return res.x
