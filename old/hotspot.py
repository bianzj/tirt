import numpy
import numpy as np
import scipy.integrate as sci
from gap import *

'''
热点/双向透过率
layer:均质场景观测角度决定；
pl：  异质场景路径长度决定；
area: 异质场景几何投影+路径长度决定；

'''

def hotspot_layer(lai, hspot, vza, sza, vsa, CIs=1, CIv = 1, G = 0.5):
    '''
    热点函数，用于土壤热点计算
    :param lai: 叶面积指数
    :param hspot: 热点系数
    :param vza:  观测天顶角
    :param sza:  太阳天顶角
    :param vsa:  观测-太阳方位角
    :param CIi:  聚集指数
    :param CIv:  聚集指数
    :param G:  投影系数
    :return: 双向透过率
    '''

    cthetv = np.cos(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    sthetv = np.sin(np.deg2rad(vza))
    cphivs = np.cos(np.deg2rad(vsa))

    ### correction 是SAIL 采纳的系数
    kv = G/cthetv
    ks = G/cthets
    correction = 2.0 / (kv + ks)

    H = 1.0
    d = H * hspot / correction
    delta = np.power(cthets, -2) + np.power(cthetv, -2) - 2 * (cthets * cthetv + sthets * sthetv * cphivs) / (cthets * cthetv)
    delta[delta<0] = 0
    # delta = 1.0 * np.sqrt(delta)
    if type(delta) == numpy.float64:
        alpha = H * np.sqrt(delta) / d
        w = (1 - np.exp(- alpha)) / alpha
    else:
        w = delta*1.0
        ind0 = np.where(delta < 0.0001)
        ind1 = np.where(delta > 0.0001)
        if np.size(ind0)>0:
            w[ind0] = 1
        if np.size(ind1)>0:
            alpha = H * np.sqrt(delta[ind1]) / d[ind1]
            w[ind1] = (1 - np.exp(- alpha)) / alpha

    pls = CIs * G / cthets
    plv = CIv * G / cthetv
    overlapping = np.sqrt(plv * pls)*w
    gapvs = np.exp(-( pls + plv - overlapping) * lai)
    return gapvs

def hotspot_analytical(lai, hspot, vza, sza, vsa, CIi = 1, CIv = 1, G = 0.5):
    '''
    多层植被热点的累计策略，假设为上下两层，上层可视假设全部被光照，下层可视认为层间独立
    :param lai:  叶面积指数
    :param hspot:  热点系数
    :param vza:  观测天顶角
    :param sza:  太阳天顶角
    :param vsa:  冠岑-太阳方位角
    :param CIi:  聚集指数
    :param CIv:  聚集指数
    :param G:  投影系数
    :return:
    '''
    cthetv = np.cos(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    sthetv = np.sin(np.deg2rad(vza))
    cphivs = np.cos(np.deg2rad(vsa))

    bv = np.exp(-G * lai * CIv / cthetv)
    bi = np.exp(-G * lai * CIi / cthets)

    ### 经验系数
    index = 0.82 - 0.04 * np.floor(lai / 0.5)

    kv = G / cthetv
    ki = G / cthets
    correction = 2.0 / (kv + ki)

    H = 1.0
    # d = H * hspot * np.ones(np.size(correction))
    d = H * hspot / correction

    delta = 1.0 * np.sqrt(
        np.power(cthets, -2) + np.power(cthetv, -2) - 2 * (cthets * cthetv + sthets * sthetv * cphivs) / (
                    cthets * cthetv))
    bv_upper = 1 - (1 - bv) * index
    bi_upper = 1 - (1 - bi) * index
    hv_upper = -np.log(bv_upper) * cthetv / G / lai  ### 在观测方向的上层比例
    hi_upper = -np.log(bi_upper) * cthets / G / lai  ### 在太阳方向的上层比例
    h_upper = np.sqrt(hv_upper * hi_upper)  # 上层高度比例
    # h_upper = (hv_upper + hi_upper)*0.5
    # h_upper = hv_upper  ### 替代，这个可能会跟SAIL更接近
    h_bottom = 1 - h_upper  # 下层高度比例
    H_upper = H * h_upper  ### 上层高度
    H_bottom = H * h_bottom  ### 下层高度
    lai_upper = lai * h_upper  ### 上层的LAI
    lai_bottom = lai * h_bottom  ### 下层的LAI
    ### w_upper = (1 - np.exp(- H_upper * delta / d)) * d / (H_upper * delta)
    w_upper = delta * 1.0
    if type(w_upper) == numpy.float64:
        if delta < 0.0001:
            w_upper = 1
        else:
            w_upper = (1 - np.exp(- H * delta / d)) * d / (H * delta)
    else:
        ind0 = np.where(delta < 0.0001)
        ind1 = np.where(delta > 0.0001)
        if np.size(ind0) > 0:
            w_upper[ind0] = 1
        if np.size(ind1) > 0:
            ### 即使分层不能改变高度和半径的比例
            # alpha = H_upper[ind1] * delta[ind1] / d[ind1]
            alpha = H * delta[ind1] / d[ind1]
            w_upper[ind1] = (1 - np.exp(- alpha)) / alpha

    overlapping = G * np.sqrt(1.0 * CIi * CIv / (cthets * cthetv)) * w_upper
    ### 计算上层植被的光照可视背景
    sunlit_fraction_background_upper = np.exp(-(CIi * G / cthets + CIv * G / cthetv - overlapping) * lai_upper)
    ### w_bottom = (1 - np.exp(-H_bottom * delta / d)) * d / (H_bottom * delta)
    w_bottom = delta * 1.0
    if type(w_bottom) == numpy.float64:
        if delta < 0.0001:
            w_bottom = 1
        else:
            w_bottom = (1 - np.exp(- H * delta / d)) * d / (H * delta)
    else:
        ind0 = np.where(delta < 0.0001)
        ind1 = np.where(delta > 0.0001)
        if np.size(ind0) > 0:
            w_bottom[ind0] = 1
        if np.size(ind1) > 0:
            ### 即使分层不能改变高度和半径的比例
            # alpha = H_bottom[ind1] * delta[ind1] / d[ind1]
            alpha = H * delta[ind1] / d[ind1]
            w_bottom[ind1] = (1 - np.exp(- alpha)) / alpha
    ### 计算下层植被的光照比例
    sunlit_fraction_vegetation_bottom = 1 - np.exp(
        -G * np.sqrt(1.0 * CIv * CIi / (cthetv * cthets)) * w_bottom * lai_bottom)

    ### 上层的透过率，1-bv_upper 上层的不透过率，假设造成这些不透过率的植被就是光照植被
    bv_upper = np.exp(-G * lai_upper * CIv / cthetv)
    gapvs = (1 - bv_upper + sunlit_fraction_vegetation_bottom * sunlit_fraction_background_upper)

    return gapvs

def hotspot_numerical(lai, hspot, vza, sza, vsa, CI = 1, G = 0.5):
    '''
    多层植被热点的累计策略，每一层都通过hotspot_layer计算
    这里是采用了SAIL的方案，并不是按照物理距离分层，而是按照透过率聚集分层
    :param lai: 叶面积指数
    :param hspot: 热点系数
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param vsa: 观测-太阳方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''


    uv = np.cos(np.deg2rad(vza))
    ui = np.cos(np.deg2rad(sza))
    si = np.sin(np.deg2rad(sza))
    sv = np.sin(np.deg2rad(vza))
    cospsi = np.cos(np.deg2rad(vsa))
    delta = 1.0 * np.sqrt(np.power(ui, -2) + np.power(uv, -2) - 2 * (ui * uv + si * sv * cospsi) / (ui * uv))

    bv = np.exp(-G*CI*lai/uv)
    ki = G / ui
    kv = G / uv
    alpha = (delta / hspot) * 2.0 /(kv + ki)

    if type(alpha) == numpy.float64:
        if alpha < 0.0001: sunlit_fraction = 1
        else:
            x1 = 0
            y1 = 0
            f1 = 1.0
            fhot = lai * CI * np.sqrt(ki * kv)
            fint = (1.0 - np.exp(-alpha)) * 0.05
            sumint = 0
            for i in range(1,20):
                if i < 20:
                    x2 = - np.log(1.0 - i * fint) / alpha
                else:
                    x2 = 1
                y2 = -(ki + kv) * lai * CI * x2 + fhot * (1 - np.exp(-alpha * x2)) / alpha
                f2 = np.exp(y2)
                sumint = sumint + (f2 - f1) * (x2 - x1) / (y2 - y1)
                x1 = x2 * 1.0
                y1 = y2 * 1.0
                f1 = f2 * 1.0
            sunlit_fraction = kv * CI * lai * sumint / (1 - bv)
    else:
        sunlit_fraction = alpha * 1.0
        ind0 = np.where(alpha < 0.0001)
        ind1 = np.where(alpha > 0.0001)
        if np.size(ind0) > 0:
            sunlit_fraction[ind0] = 1.0
        if np.size(ind1) > 1:
            x1 = 0
            y1 = 0
            f1 = 1.0
            fhot = lai * CI * np.sqrt(ki * kv)
            fint = (1.0 - np.exp(-alpha)) * 0.05
            sumint = 0
            for i in range(1,21):
                if i < 20:
                    x2 = - np.log(1.0 - i * fint[ind1]) / alpha[ind1]
                else:
                    x2 = 1
                y2 = -(ki[ind1] + kv[ind1]) * lai * CI * x2 + fhot[ind1] * (1 - np.exp(-alpha[ind1] * x2)) / alpha[ind1]
                f2 = np.exp(y2)
                sumint = sumint + (f2 - f1) * (x2 - x1) / (y2 - y1)
                x1 = x2 * 1.0
                y1 = y2 * 1.0
                f1 = f2 * 1.0
            sunlit_fraction[ind1] = kv[ind1] * CI * lai * sumint / (1 - bv[ind1])

    return sunlit_fraction

def hotspot_layer_function(x):
    '''
    单层的热点函数，与hotspot_layer 相同，差别在于部分系数被放在全局变量
    :param x:  距离系数 0- 1
    :return:
    '''
    global lai0,G0,hspot0,sza0,vza0,raa0
    lai = lai0 *1.0
    G = G0 *1.0
    hspot = hspot0 * 1.0
    sza = sza0 * 1.0
    vza = vza0 * 1.0
    raa = raa0 * 1.0
    eps = 0.001
    rd = np.pi/180.0
    sthetv = np.sin(vza * rd)
    cthetv = np.cos(vza * rd)
    sthets = np.sin(sza * rd)
    cthets = np.cos(sza * rd)
    cospsi = np.cos(raa * rd)
    tants = sthets / cthets
    tanto = sthetv / cthetv
    dso = np.sqrt(tants * tants + tanto * tanto - 2.0 * tants * tanto * cospsi)
    ko = G / cthetv
    ks = G / cthets

    if dso < eps:
        pso = np.exp((ks + ko) * lai * (x) - np.sqrt(ko * ks) * lai * (x))
    else:
        alf = (dso / hspot) * 2 * (ks + ko)
        pso = np.exp((ks + ko) * lai * (x) + np.sqrt(ko * ks) * lai / alf * (1 - np.exp(x * alf)))
    return pso

def hotspot_integrated(lai, hspot, vza, sza, raa,  G = 0.5):
    '''
	通过积分的方式解决层间问题
    :param lai: 叶面积指数
    :param hspot: 热点系数
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 观测-太阳方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    nz = 20
    xl = np.linspace(0,-1,nz)
    dx = np.float(1.0/nz)
    gapvs = np.zeros(nz)
    global lai0,G0,hspot0,sza0,vza0,raa0
    lai0 = lai*1.0
    G0 = G*1.0
    hspot0 = hspot*1.0
    sza0 = sza*1.0
    vza0 = vza*1.0
    raa0 = raa*1.0

    for j in range(nz):
        t1 = xl[j]-dx
        t2 = xl[j]
        gapvs[j] = sci.quad(hotspot_layer_function, t1, t2)[0]
        gapvs[j] = gapvs[j]/np.float(dx)

    return gapvs

def hotspot_vegetation_volume(lai, sza, CI = 1.0, G = 0.5):
    '''
    植被体系的整体的光照比例
    :param lai: 叶面积指数
    :param sza: 太阳天顶角
    :param CI:  聚集指数
    :param G:  投影系数
    :return:  光照比例
    '''
    sthets = np.cos(np.deg2rad(sza))
    ### 植被层的透过率
    gap_probability_illuminate = np.exp(-G * lai * CI / sthets[0])
    ### 1-透过率：植被阻挡概率； 植被的叶面积指数情况，计算的只是光照比例
    sunlit_fraction_volume = (1 - gap_probability_illuminate) / (G * lai * CI) * sthets[0]
    return sunlit_fraction_volume

def hotspot_hom_voxel(lai, hspot, vza, sza, vsa,CIv = 1,CIs = 1, G = 0.5):
    '''
	通过体元的方式计算热点
    :param lai: 叶面积指数
    :param number_voxel: 体元数目
    :param hspot: 热点系数
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param raa: 观测-太阳方位角
    :param CI:  聚集指数
    :param G:   投影系数
    :return:
    '''
    ### 体素场景的双向透过率，因为是层所以整体的和单层的计算方法一致
    gapvs = hotspot_layer(lai, hspot, vza, sza, vsa,CIv,CIs,G)

    return gapvs

def hotspot_path_length(pl_v, pl_s,density,hspot, vza, sza, vsa, G = 0.5):
    '''
	基于路径长度的热点函数，可以用于并不是均质层状分布的场景
	!!! CI聚集指数这里没有用，这是用于计算路径长度用的
	:param pl_v:  观测方向路径长度
	:param pl_s:  太阳方向路径长度
	:param density:  树冠体密度
	:param hspot:  热点系数
	:param vza:  观测天顶角
	:param sza:  太阳天顶角
	:param vsa:  观测-太阳方位角
	:param CI:  聚集指数
	:param G:   投影系数
	:return:
    '''


    eps = 0.0001
    number = np.size(pl_v)
    uv = np.cos(np.deg2rad(vza))
    ui = np.cos(np.deg2rad(sza))
    si = np.sin(np.deg2rad(sza))
    sv = np.sin(np.deg2rad(vza))
    cospsi = np.cos(np.deg2rad(vsa))
    cosalpha = si*sv*cospsi + uv*ui
    # delta = np.power(pl_s-pl_v,2)+(1-cosalpha)*2*pl_s*pl_v
    delta = np.power(pl_v,2) + np.power(pl_s,2) - 2*cosalpha* pl_v * pl_s
    delta[delta<0] = 0
    delta = np.sqrt(delta)

    '''校正，并且计算重叠概率'''
    '''bg = sigma * h / d'''
    # correct = np.zeros(number)
    # w = np.zeros(number)
    # ind = pl > eps
    # if np.sum(ind) > 0:
    # 	correct[ind] = np.sqrt(pl[ind]) / hspot
    # ind = correct < eps
    # if np.sum(ind) > 0:
    # 	w[ind] = 1.0 - correct[ind] * 0.5
    # ind = correct > eps
    # if np.sum(ind) > 0:
    # 	w[ind] = (1 - np.exp(-correct[ind])) / correct[ind]

    kv = G / uv
    ki = G / ui
    correction = 2.0 / (kv + ki)
    sigma = (delta / hspot)* correction
    w = np.ones(number)
    ind = sigma > eps
    if np.sum(ind) > 0:
        w[ind] = (1-np.exp(-sigma[ind]))/sigma[ind]
    ind = sigma < eps
    if np.sum(ind) > 0:
        w[ind] = 1.0 - sigma[ind] * 0.5
    ### 计算投影重合比例
    overlapping = np.zeros(number)
    ind = pl_v * pl_s > 0
    if np.sum(ind)>0:
        overlapping[ind] = np.sqrt(pl_v[ind] * pl_s[ind]) * w[ind]
    pl_hotspot = (pl_v + pl_s-overlapping)
    gap_hotspot = np.exp(-pl_hotspot * density * G)

    return gap_hotspot,pl_hotspot

def hotspot_projection_area(std, density, hspot, xv, yv, zv, plv, upAreav, hcrv_up, interv,
                            xs, ys, zs, pls, upAreas, hcrs_up, inters, vza, sza, vsa, G = 0.5, CI = 1.0):
    '''
    树冠的热点函数,通过投影面积的方式,相比于路径长度，考虑了树冠形状的影响
    :param std: 森林树密度
    :param density:  树冠的体密度
    :param hspot:  热点因子
    :param xv:  树冠表面x
    :param yv: 树冠表面y
    :param zv:  树冠表面z
    :param plv:  路径长度 观测方向
    :param upAreav:  观测方向 投影面积
    :param hcrv_up:  等效半径长度
    :param interv:   阻挡因子 1- gap
    :param xs:   树冠表面x
    :param ys:   树冠表面y
    :param zs:   树冠表面z
    :param pls:  路径长度 太阳平面
    :param upAreas:  太阳方向 投影面积
    :param hcrs_up:  等效半径长度
    :param inters:  阻挡一直
    :param vza:  观测天顶角
    :param sza:  太阳天顶角
    :param vsa:  观测-太阳方位角
    :param G:    投影系数
    :param CI:   聚集指数
    :return:
    '''
    eps = 0.0001
    sthetv = np.sin(np.deg2rad(vza))
    cthetv = np.cos(np.deg2rad(vza))
    sthets = np.sin(np.deg2rad(sza))
    cthets = np.cos(np.deg2rad(sza))
    sphivs = np.sin(np.deg2rad(vsa))
    cphivs = np.cos(np.deg2rad(vsa))
    calph = sthets * sthetv * cphivs + cthetv * cthets
    # ti = si / ui
    # tv = sv / uv


    #####################################
    #### 计算观测和太阳方向两个投影圆圈的关系
    #####################################
    x7 = xv + hcrv_up * cphivs
    y7 = yv + hcrv_up * sphivs
    x8 = xs + hcrs_up
    y8 = ys * 1.0
    ### 两个投影圆圈中心之间的距离
    sl8 = np.sqrt((x7 - x8) * (x7 - x8) + (y7 - y8) * (y7 - y8))

    ### 确定大圆圈和小圆圈机器半径
    largeArea = upAreav * 1.0
    ind = upAreav < upAreas
    if np.sum(ind) > 0: largeArea[ind] = upAreas[ind]
    smallArea = upAreav * 1.0
    ind = upAreav > upAreas
    if np.sum(ind) > 0: smallArea[ind] = upAreas[ind]
    large_rcr = np.sqrt(largeArea / np.pi)
    small_rcr = np.sqrt(smallArea / np.pi)
    number = np.size(large_rcr)
    scomm = np.zeros(number)
    ### 小圆在大圆里边
    ind = sl8 <= (large_rcr - small_rcr)
    if np.sum(ind) > 0:
        scomm[ind] = smallArea[ind]
    ind = sl8 < eps
    if np.sum(ind) > 0:
        scomm[ind] = smallArea[ind]
    ### 大圆和小圆之间有重叠区
    ind = (sl8 < (large_rcr + small_rcr)) * (sl8 > (large_rcr - small_rcr))
    if np.sum(ind) > 0:
        ppp = (large_rcr[ind] + small_rcr[ind] + sl8[ind]) * .5
        angl4 = np.arccos((large_rcr[ind] * large_rcr[ind] + sl8[ind] * sl8[ind] - small_rcr[ind] * small_rcr[ind]) / (
                    2.0 * large_rcr[ind] * sl8[ind]))
        angl8 = np.arccos(
            (large_rcr[ind] * large_rcr[ind] + small_rcr[ind] * small_rcr[ind] - sl8[ind] * sl8[ind]) / (
                        2.0 * large_rcr[ind] * small_rcr[ind]))
        angl5 = (np.pi - angl4 - angl8)
        ss4 = angl4 * large_rcr[ind] * large_rcr[ind]
        ss5 = angl5 * small_rcr[ind] * small_rcr[ind]
        ss6 = np.sqrt(ppp * (ppp - large_rcr[ind]) * (ppp - small_rcr[ind]) * (ppp - sl8[ind]))
        scomm[ind] = ss4 + ss5 - 2.0 * ss6


    #####################################
    #### 计算太阳和观测球体重叠区域的联合透过率
    #####################################

    rlls = pls
    rllv = plv
    delta = (rllv - rlls) * (rllv - rlls) + (1.0 - calph) * 2.0 * rlls * rllv
    # delta = np.power(plv,2) + np.power(pls,2) - 2*calph* plv * pls
    delta[delta<0] = 0
    delta = np.sqrt(delta)

    kv = G / cthetv
    ki = G / cthets
    correction = 2.0 / (kv + ki)
    sigma = (delta / hspot)* correction
    w = np.ones(number)
    ind = sigma > eps
    if np.sum(ind) > 0:
        w[ind] = (1-np.exp(-sigma[ind]))/sigma[ind]
    ind = sigma < eps
    if np.sum(ind) > 0:
        w[ind] = 1- sigma[ind] * 0.5

    overlapping = np.sqrt(plv*pls) * w
    pl_hotspot = (pls + plv - overlapping)
    gap_hotspot = np.exp(-pl_hotspot*density * G)

    CICI = 1 - CI
    ### intercept 等效不透过率、等效遮挡概率
    if np.abs(CICI)> 0.001:	intercept = -np.log(1.0-(1.0-gap_hotspot)*CICI)/CICI
    else: intercept = 1- gap_hotspot

    ### sumu 是很多不同类型植被的集合，这里只有1中植被类型，所以就直接是结果
    sumu = 0
    uu = std * (interv * (upAreav - scomm) + inters * (upAreas - scomm) + intercept * scomm)
    sumu = sumu + uu

    gapvs = np.exp(-sumu)
    return gapvs

def hotspot_projection_area_under(std, density, hspot, xv, yv, zv, plv, upAreav, hcrv_up, interv,
                                  xs, ys, zs, pls, upAreas, hcrs_up, inters,
                                  hc, hcr, dbh_trunk, stemProjv, stemProjs,
                                  vza, sza, vsa, G = 0.5, CI = 1.0):
    '''
    树冠的热点函数，路径长度，树冠形状，外加树干影响，这里就是FRT里边的 spooj.f
    :param std: 森林树密度
    :param density:  树冠的体密度
    :param hspot:  热点因子
    :param xv:  树冠表面x
    :param yv: 树冠表面y
    :param zv:  树冠表面z
    :param plv:  路径长度 观测方向
    :param upAreav:  观测方向 投影面积
    :param hcrv_up:  等效半径长度
    :param interv:   阻挡因子 1- gap
    :param xs:   树冠表面x
    :param ys:   树冠表面y
    :param zs:   树冠表面z
    :param pls:  路径长度 太阳平面
    :param upAreas:  太阳方向 投影面积
    :param hcrs_up:  等效半径长度
    :param inters:  阻挡一直
    :param vza:  观测天顶角
    :param sza:  太阳天顶角
    :param vsa:  观测-太阳方位角
    :param G:    投影系数
    :param CI:   聚集指数
    :return:
    '''
    number_voxel = np.size(vza)
    sthetv = np.sin(np.deg2rad(vza))
    cthetv = np.cos(np.deg2rad(vza))
    sthets = np.sin(np.deg2rad(sza))
    cthets = np.cos(np.deg2rad(sza))
    sphivs = np.sin(np.deg2rad(vsa))
    cphivs = np.cos(np.deg2rad(vsa))
    calph = sthets * sthetv * cphivs + cthetv * cthets
    tgths = sthets / cthets
    tgthv = sthetv / cthetv
    eps = 0.0001

    #########################################################################
    ##################  树冠之间的overlaping
    #########################################################################

    ###  the distance between projection centres
    ### s方向是2，v方向是1
    x7 = xv + hcrv_up * cphivs
    y7 = yv + hcrv_up * sphivs
    x8 = xs + hcrs_up
    y8 = ys * 1.0
    sl8 = np.sqrt((x7 - x8) * (x7 - x8) + (y7 - y8) * (y7 - y8))

    largeArea = upAreav * 1.0
    ind = upAreav < upAreas
    if np.sum(ind) > 0: largeArea[ind] = upAreas[ind]
    smallArea = upAreav * 1.0
    ind = upAreav > upAreas
    if np.sum(ind) > 0: smallArea[ind] = upAreas[ind]

    large_rcr = np.sqrt(largeArea / np.pi)
    small_rcr = np.sqrt(smallArea / np.pi)
    number = np.size(large_rcr)
    scomm = np.zeros(number)
    ind = sl8 <= (large_rcr - small_rcr)
    if np.sum(ind) > 0:
        scomm[ind] = smallArea[ind]
    ind = sl8 < eps
    if np.sum(ind) > 0:
        scomm[ind] = smallArea[ind]

    ind = (sl8 < (large_rcr + small_rcr)) * (sl8 > (large_rcr - small_rcr))
    if np.sum(ind) > 0:
        ppp = (large_rcr[ind] + small_rcr[ind] + sl8[ind]) * .5
        angl4 = np.arccos((large_rcr[ind] * large_rcr[ind] + sl8[ind] * sl8[ind] - small_rcr[ind] * small_rcr[ind]) / (
                    2.0 * large_rcr[ind] * sl8[ind]))
        angl8 = np.arccos(
            (large_rcr[ind] * large_rcr[ind] + small_rcr[ind] * small_rcr[ind] - sl8[ind] * sl8[ind]) / (
                        2.0 * large_rcr[ind] * small_rcr[ind]))
        angl5 = (np.pi - angl4 - angl8)
        ss4 = angl4 * large_rcr[ind] * large_rcr[ind]
        ss5 = angl5 * small_rcr[ind] * small_rcr[ind]
        ss6 = np.sqrt(ppp * (ppp - large_rcr[ind]) * (ppp - small_rcr[ind]) * (ppp - sl8[ind]))
        scomm[ind] = ss4 + ss5 - 2.0 * ss6


    ###! external hot-spot correction
    rlls = pls
    rllv = plv
    delta = (rllv - rlls) * (rllv - rlls) + (1.0 - calph) * 2.0 * rlls * rllv
    delta[delta<0] = 0
    delta = np.sqrt(delta)
    kv = G / cthetv
    ki = G / cthets
    correction = 2.0 / (kv + ki)
    alpha = (delta / hspot)* correction
    w = np.ones(number)
    ind = alpha > eps
    if np.sum(ind) > 0:
        w[ind] = (1-np.exp(-alpha[ind]))/alpha[ind]
    overlapping = np.sqrt(plv*pls) * w
    pl_hotspot = (pls + plv - overlapping)
    gap_hotspot = np.exp(-pl_hotspot*density * G)

    CICI = 1 - CI
    if np.abs(CICI)> 0.001:	inter = -np.log(1.0-(1.0-gap_hotspot)*CICI)/CICI
    else: inter = 1- gap_hotspot


    ###############################################################################
    ####  树冠和树干之间的重叠 面积
    ###############################################################################

    slty1 = np.zeros(number)
    slty2 = np.zeros(number)
    scty = np.zeros(number)

    hcb = hc - hcr * 2
    ind = zv < hcb
    if np.sum(ind) > 0: slty1[ind] = (hcb - zv[ind]) * tgthv
    ind = zs < hcb
    if np.sum(ind) > 0: slty2[ind] = (hcb - zs[ind]) * tgths

    ind1 = (stemProjv > 0) * (stemProjs > 0.0)
    ind2 = (cphivs > eps) * (np.abs(yv - ys) < eps)

    small_slty = slty1 * 1.0
    ind = slty1 > slty2
    if np.sum(ind) > 0: small_slty[ind] = slty2[ind]
    ind = small_slty > eps
    tstphi = np.zeros(number)
    tstphi[ind] = dbh_trunk / small_slty[ind]

    x5 = xs + slty2
    x6 = xv + slty1 * cphivs
    x56min = x5 * 1.0
    ind = x56min > x6
    if np.sum(ind) > 0: x56min[ind] = x6[ind]
    x12max = xv * 1.0
    ind = x12max < xs
    if np.sum(ind) > 0: x12max[ind] = xs[ind]

    ind3 = (np.abs(xv - xs) < eps) * (sphivs > tstphi)
    indtemp =  (sphivs > eps)
    ind4 = (np.abs(xv - xs) < eps) * (x56min > x12max)
    ind = ind1 * ind2 * ind3 * indtemp
    if  np.sum(ind) > 0: scty[ind] = dbh_trunk * dbh_trunk * 0.50 / sphivs[ind]
    ind = ind1 * ind2 * ind4
    if np.sum(ind) > 0: scty[ind] = (x56min[ind] - x12max[ind]) * dbh_trunk

    ind5 = (np.abs(xv - xs) > eps) * (sphivs > eps)
    ind6 = (np.abs(xv - xs) > eps) * (x56min > x12max)
    ind = ind1 * ind2 * ind6
    if np.sum(ind) > 0: scty[ind] = (x56min[ind] - x12max[ind]) * dbh_trunk

    styx = stemProjv * 1.0
    ind = styx > stemProjs
    if np.sum(ind) > 0: styx[ind] = stemProjs[ind]
    ind = styx > scty
    if np.sum(ind) > 0: styx[ind] = scty[ind]
    styx = stemProjv + stemProjs - styx


    sumu = 0
    uu = std * (interv * (upAreav - scomm) + inters * (upAreas - scomm) + inter * scomm + styx)
    sumu = sumu + uu

    bidirectional_gap = np.exp(-sumu)
    return bidirectional_gap



