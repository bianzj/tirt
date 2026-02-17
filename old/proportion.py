import numpy as np
import scipy.integrate as sci
from gap import *
from hotspot import *
from utils import *
# from scatter import *
import matplotlib.pyplot as plt
import quadpy

'''
冠层内组分的可视比例
proportion_directional: 表示单向的；
proportion_bidirectional: 表示双向的；
voxel: 表示体素的进体素的出；
one： 表示参数化方法直接算组分；
voxel_one: 表示体素方法直接积分到组分；
voxel_layer: 表示体素方法积分到了等效层；
_unit: 表示四组分,默认四组分，直接隐去；
_individual: 表示六组分，四组分外加树干；
_endmember:表示八组分，四组分外加树干和低层植被，端元可以直接是一组分、两组分、三组分。。。；

'''

####################################################################
###### Analytical Solution
####################################################################


def proprotion_directional_crown_one(lai,std,hcr,rcr,xza, G = 0.5):
    '''
    计算单向的组分可视比例，观测方向或者太阳方向
    :param lai:
    :param std:
    :param hcr:
    :param rcr:
    :param xza:
    :param G:
    :return:
    '''
    Psoil = gap_probability_crown_analytical(lai, std, hcr, rcr, xza, G)
    Pleaf = 1- Psoil
    return Psoil,Pleaf

def proportion_directional_hom_one(lai,xza,CI = 1.0, G = 0.5):
    Psoil = gap_probability_hom_analytical(lai,xza)
    Pleaf = 1-Psoil
    return Psoil,Pleaf

def proportion_directional_row_one(lai,row_width, row_blank, row_height, xza, xaa, raa):
    xra = np.abs(xaa - raa)
    Psoil = gap_probability_row_analytical(lai,row_width,row_blank,row_height,xza,xra)
    Pleaf = 1 - Psoil
    return Psoil, Pleaf

def proportion_bidirectional_crown_one(lai, std, hspot, hcr, rcr, vza, sza, vsa, G = 0.5):
    '''
    参数化方法，计算森林冠层四组分的可视比例
    :param lai:
    :param std:
    :param hspot:
    :param hcr:
    :param rcr:
    :param vza:
    :param sza:
    :param vsa:
    :param G:
    :return:
    '''

    bv = gap_probability_crown_analytical(lai, std, hcr, rcr, vza)
    bi = gap_probability_crown_analytical(lai, std, hcr, rcr, sza)
    Psoil = bv
    Pleaf = 1 - Psoil
    uv = np.cos(np.deg2rad(vza))
    ui = np.cos(np.deg2rad(sza))
    CIv = -np.log(bv) * uv / (lai * G)
    CIi = -np.log(bi) * ui / (lai * G)
    Psoil_sunlit = hotspot_layer(lai, hspot, vza, sza, vsa, CIi, CIv)
    Pleaf_sunlit = hotspot_analytical(lai, hspot, vza, sza, vsa, CIi, CIv)
    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded = Pleaf - Pleaf_sunlit

    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

def proportion_bidirectional_hom_one(lai,hspot,vza,sza,vsa):
    Psoil = gap_probability_hom_analytical(lai, vza)
    Pleaf = 1 - Psoil
    Psoil_sunlit = hotspot_layer(lai, hspot, vza, sza, vsa)
    Pleaf_sunlit = hotspot_analytical(lai, hspot, vza, sza, vsa)
    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded = Pleaf - Pleaf_sunlit
    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

def proportion_bidirectional_row_one(lai,hspot,row_width,row_blank,row_height,vza,vaa,sza,saa,raa,G = 0.5):

    vsa = np.abs(saa - vaa)
    sra = np.abs(saa - raa) % 180
    vra = np.abs(vaa - raa) % 180
    vsa[vsa > 180] = 360 - vsa[vsa > 180]

    bv = gap_probability_row_analytical(lai, row_width, row_blank, row_height, vza, vra)
    bi = gap_probability_row_analytical(lai, row_width, row_blank, row_height, sza, sra)
    Psoil = bv
    Pleaf = 1 - Psoil
    cthetv = np.cos(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    CIv = -np.log(bv) * cthetv / (lai * G)
    CIi = -np.log(bi) * cthets / (lai * G)
    Psoil_sunlit = hotspot_layer(lai, hspot, vza, sza, vsa, CIi, CIv)
    Pleaf_sunlit = hotspot_analytical(lai, hspot, vza, sza, vsa, CIi, CIv)
    ####################################################################
    ##### '''Interesting !!!'''
    ####################################################################

    saa_ = saa - raa
    vaa_ = vaa - raa
    saa_[saa_ < 0] = 360 + saa_[saa_ < 0]
    vaa_[vaa_ < 0] = 360 + vaa_[vaa_ < 0]
    saa_[saa_ > 360] = saa_[saa_ > 360] - 360
    vaa_[vaa_ > 360] = vaa_[vaa_ > 360] - 360
    shadow_s = row_height * np.tan(np.deg2rad(sza)) * np.sin(np.deg2rad(sra))
    shadow_v = row_height * np.tan(np.deg2rad(vza)) * np.sin(np.deg2rad(vra))
    shadow_v[shadow_v > row_blank] = row_blank
    shadow_s[shadow_s > row_blank] = row_blank
    '''同向'''
    ind = (row_blank > shadow_v) * (row_blank > shadow_s) * ((saa_ - 180) * (vaa_ - 180) > 0)
    correction = np.ones(np.size(vza))
    if np.sum(ind) > 0:
        shadow = shadow_v * 1.0
        shadow[shadow < shadow_s] = shadow_s[shadow < shadow_s]
        correction[ind] = (row_width + row_blank) / (shadow[ind] + row_width)
    '''不同向'''
    ind = ((saa_ - 180) * (vaa_ - 180) < 0) * ((shadow_s + shadow_v) < row_blank)
    if np.sum(ind) > 0:
        shadow = shadow_s + shadow_v
        correction[ind] = (row_width + row_blank) / (shadow[ind] + row_width)
    '''顺垄'''
    ind = (shadow_s * shadow_v == 0) * (row_blank > shadow_v) * (row_blank > shadow_s)
    if np.sum(ind) > 0:
        shadow = shadow_v * 1.0
        shadow[shadow < shadow_s] = shadow_s[shadow < shadow_s]
        correction[ind] = (row_width + row_blank) / (shadow[ind] + row_width)
    testi = (shadow_s + row_width) / (row_blank + row_width)
    testv = (shadow_v + row_width) / (row_width + row_blank)
    Ksoil_sunlit = Psoil_sunlit / Psoil
    Ksoil_sunlit = Ksoil_sunlit / correction + (1 - 1.0 / correction)
    Psoil_sunlit = Psoil * Ksoil_sunlit

    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded = Pleaf - Pleaf_sunlit
    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

####################################################################
##### Voxel-based solution
####################################################################
#### directional

def proportion_directional_crown_voxel(lai, std, hcr, rcr, xza0, xaa0 = 0, G = 0.5, ifleaf = 1):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param std:  森林树密度
    :param hc:   冠层高度
    :param hcr:   垂直方向的树冠半径
    :param rcr:   水平方向的树冠半径
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param vaa0:   观测方位角
    :param sza0:   太阳天顶角
    :param saa0:   太阳方位角
    :param G:      投影系数[
    :return:  组分比例
    '''


    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    density = lai / (std * vol_crown)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    voxels = scheme.points
    vo = scheme.weights *vol_crown
    number_voxel = np.size(vo)
    number_angle = np.size(xza0)
    number =  number_angle * number_voxel


    xza = np.repeat(xza0, number_voxel)
    xaa = np.repeat(xaa0, number_voxel)

    x = np.tile(voxels[0,:] * rcr,number_angle)
    y = np.tile(voxels[1,:] * rcr,number_angle)
    z = np.tile(voxels[2,:] * hcr,number_angle)
    vo = np.tile(vo, number_angle)
    x_soil = np.asarray([0])
    y_soil = np.asarray([0])
    z_soil = np.asarray([0])

    cthetx = np.cos(np.deg2rad(xza))
    sthetx = np.sin(np.deg2rad(xza))
    tgthv = np.tan(np.deg2rad(xza))
    cphiv = np.cos(np.deg2rad(xaa))
    sphiv = np.sin(np.deg2rad(xaa))

    hc = hcr * 2
    ### 观测方向的透过率 = 观测方向树冠内透过率 * 观测方向树冠间/树冠外的透过率
    gapx_inside_up,plv_inside_up,gapv_inside_down,plv_inside_down = gap_probability_crown_inside_voxel(x,y,z,density,hcr,rcr,xza,xaa)
    xv_up = x + plv_inside_up * sthetx * cphiv
    yv_up = y + plv_inside_up * sthetx * sphiv
    zv_up = z + plv_inside_up * cthetx + hcr
    gapv_outside_up,plv_outside_up,upAreav,hcrv_up,interv = \
        gap_probability_crown_outside_voxel(std, xv_up, yv_up, zv_up, density, hc, hcr,rcr, xza, xaa)
    gapv = gapx_inside_up * gapv_outside_up

    base_voxel = std * vo * density * G / cthetx
    Pleaf_voxel =  gapv * base_voxel
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Pleaf = np.sum(Pleaf_voxel, axis=0)

    gapv_soil, plv_outside_soil, upAreav_soil, hcrv_soil, interv_soil = \
        gap_probability_crown_outside_voxel(std, x_soil, y_soil, z_soil, density, hc, hcr, rcr, xza0, xaa0)
    Psoil = gapv_soil

    if ifleaf == 0:
        ### 比例偏差通过土壤平差
        Psoil = 1 - Pleaf

    else:
        ### 比例偏差通过叶片平差
        Pleafnew = 1-Psoil
        Kcorrect = Pleafnew / Pleaf
        Pleaf_voxel = Pleaf_voxel * Kcorrect

    return Psoil,Pleaf_voxel

def proportion_directional_hom_voxel(lai0, vza0,  CI = 1.0, G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    number_voxel = 10
    number_angle = np.size(vza0)
    vza = np.repeat(vza0, number_voxel)
    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))

    dlai = lai0 / number_voxel
    step = 1.0 / number_voxel
    lai = np.linspace(step * 0.5, 1 - step * 0.5, number_voxel) * lai0
    lai = np.tile(lai,number_angle)

    Gleaf_voxel = gap_probability_hom_voxel(lai, vza)
    Pleaf_voxel = Gleaf_voxel * dlai * G * CI / cthetv
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Psoil = gap_probability_hom_analytical(lai0, vza0)


    return Psoil,Pleaf_voxel

def proportion_directional_row_voxel(lai, row_width, row_blank, row_height, vza0, vaa0, raa0,
						nvw = 5, nvh = 10, nvb = 5,  G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    rw = row_width
    rb = row_blank
    rh = row_height
    rs = rw + rb
    nvw = nvw
    nvh = nvh
    nvb = nvb
    nvs = nvw + nvb
    vra0 = np.abs(vaa0 - raa0) % 180
    laie = lai * rs / rw
    density = laie / rh

    scale_w = rw / nvw
    scale_h = rh / nvh
    scale_s = rb / nvb

    length_layer_width = np.linspace(0.5, nvw - 0.5, nvw) * scale_w
    length_layer_blank = np.linspace(0.5, nvb - 0.5, nvb) * scale_s + rw
    height_layer = np.linspace(0.5, nvh - 0.5, nvh) * scale_h
    length_voxel_vegetation = np.tile(length_layer_width, (nvh))
    height_voxel_vegetation = np.asarray(np.tile(np.transpose(np.asmatrix(height_layer)), (1, nvw)))


    ###'''植被和背景体素延展'''
    length_voxel_vegetation = np.reshape(length_voxel_vegetation, -1)
    height_voxel_vegetation = np.reshape(height_voxel_vegetation, -1)
    length_voxel_background = np.hstack([length_layer_width, length_layer_blank])
    height_voxel_background = np.repeat(rh, nvs)
    # length_voxel = np.hstack([length_voxel_vegetation, length_voxel_background])
    # height_voxel = np.hstack([height_voxel_vegetation, height_voxel_background])
    number_voxel_vegetation = np.size(length_voxel_vegetation)
    number_voxel_background = np.size(length_voxel_background)
    number_angle = np.size(vza0)
    length_voxel_background = np.tile(length_voxel_background, number_angle)
    height_voxel_background = np.tile(height_voxel_background, number_angle)
    length_voxel_vegetation = np.tile(length_voxel_vegetation, number_angle)
    height_voxel_vegetation = np.tile(height_voxel_vegetation, number_angle)

    ### 根据数目进行角度的展开
    vza_vegetation = np.repeat(vza0, number_voxel_vegetation)
    vra_vegetation = np.repeat(vra0, number_voxel_vegetation)
    vza_background = np.repeat(vza0, number_voxel_background)
    vra_background = np.repeat(vra0, number_voxel_background)

    ### 透过率和路径长度
    Gvleaf_voxel,plv_leaf = gap_probability_row_voxel(density,row_width,row_blank,row_height,
                length_voxel_vegetation,height_voxel_vegetation,vza_vegetation,vra_vegetation,0,0)
    cthetv = np.cos(np.deg2rad(vza_vegetation))
    ### 体素的LAI = (density * rh / nvh * rs / nvs)
    Pvleaf_voxel = Gvleaf_voxel * density * G /cthetv * (rs/ nvs) *(rh / nvh)
    Pvsoil_voxel,plv_soil = gap_probability_row_voxel(density,row_width,row_blank,row_height,
               length_voxel_background, height_voxel_background,vza_background,vra_background,0,0)
    ### 土壤的集成，因为每个的面积是相同的，所以就直接取数值平均，理应是按照距离进行加权平均
    Pvsoil_voxel = Pvsoil_voxel * (1.0/nvs)

    Pvleaf_voxel = np.transpose(np.reshape(Pvleaf_voxel, [number_angle, number_voxel_vegetation]))
    Pvsoil_voxel = np.transpose(np.reshape(Pvsoil_voxel, [number_angle, number_voxel_background]))

    return Pvsoil_voxel,Pvleaf_voxel

def proportion_directional_crown_voxel_one(lai, std, hcr, rcr, xza0, xaa0 = 0, G = 0.5, ifleaf = 1):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param std:  森林树密度
    :param hc:   冠层高度
    :param hcr:   垂直方向的树冠半径
    :param rcr:   水平方向的树冠半径
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param vaa0:   观测方位角
    :param sza0:   太阳天顶角
    :param saa0:   太阳方位角
    :param G:      投影系数[
    :return:  组分比例
    '''

    ### 因为 xaa0 的默认值设定，需要对其进行大小控制
    if np.size(xaa0) != np.size(xza0):
        xaa0 = np.resize(xaa0,np.size(xza0))
    Psoil,Pleaf_voxel = proportion_directional_crown_voxel(lai,std,hcr,rcr,xza0,xaa0,G,ifleaf)
    ### 将体素进行集成
    Pleaf = np.sum(Pleaf_voxel,axis=0)

    return Psoil,Pleaf

def proportion_directional_hom_voxel_one(lai0, vza0,  CI = 1.0, G = 0.5):
    '''
    冠层场景的组分观测比例，将体素的进行组装集成
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    Psoil, Pleaf_voxel = proportion_directional_hom_voxel(lai0,vza0,CI, G)
    Pleaf = np.sum(Pleaf_voxel,axis = 0)


    return Psoil,Pleaf

def proportion_directional_row_voxel_one(lai,row_width, row_blank, row_height, vza0, vaa0, raa0,
						nvw = 5, nvh = 10, nvb = 5,  G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''
    ### 这里是体素的集成比例
    Psoil_voxel,Pleaf_voxel = proportion_directional_row_voxel(lai,row_width,row_blank,row_height,vza0,vaa0,raa0,nvw,nvh,nvb,G)
    Psoil = np.sum(Psoil_voxel,axis=0)
    Pleaf = np.sum(Pleaf_voxel,axis=0)

    return Psoil,Pleaf


def proportion_directional_crown_voxel_layer(lai, std, hcr, rcr, xza0,number_layer = 5, xaa0 = 0, G = 0.5, ifleaf = 1):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param std:  森林树密度
    :param hc:   冠层高度
    :param hcr:   垂直方向的树冠半径
    :param rcr:   水平方向的树冠半径
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param vaa0:   观测方位角
    :param sza0:   太阳天顶角
    :param saa0:   太阳方位角
    :param G:      投影系数[
    :return:  组分比例
    '''
    Psoil,Pleaf_voxel = proportion_directional_crown_voxel(lai,std,hcr,rcr,xza0,xaa0,G,ifleaf)
    ### 计算设置层的观测比例，按照垂直高度进行的分层
    number_angle = np.size(xza0)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    vz = scheme.points[2,:]
    number_voxel = np.size(vz)
    Pleaf_layer = np.zeros([number_layer,number_angle])
    dlh = 1 * 2 / number_layer
    index_layer = (vz+1) // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer,:] = np.sum(Pleaf_voxel[ind,:], axis = 0)
    return Psoil,Pleaf_layer

def proportion_directional_hom_voxel_layer(lai0, xza0, number_layer = 5, CI = 1.0, G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    ### 假设分为5层，然后计算结果，按照距离的分层
    Psoil, Pleaf_voxel = proportion_directional_hom_voxel(lai0, xza0, CI, G)
    number_voxel = np.size(Pleaf_voxel)
    number_angle = np.size(xza0)
    vz = np.linspace(0.5,number_voxel-0.5,number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    Pleaf_layer = np.zeros([number_layer,number_angle])
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer,:] = np.sum(Pleaf_voxel[ind,:],axis = 0)


    return Psoil,Pleaf_layer

def proportion_directional_row_voxel_layer(lai, row_width, row_blank, row_height, xza0, xaa0, raa0, number_layer = 5,
                                           nvw = 5, nvh = 10, nvb = 5, G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    Psoil_voxel,Pleaf_voxel = proportion_directional_row_voxel(lai, row_width, row_blank, row_height, xza0, xaa0, raa0, nvw, nvh, nvb, G)
    Psoil = np.sum(Psoil_voxel,axis=0)

    ### 垄行做的的分层，这里是按照垂直高度的分层
    number_voxel = nvh
    number_angle = np.size(xza0)
    Pleaf_layer = np.zeros([number_layer,number_angle])
    vz = np.linspace(0.5,number_voxel-0.5,number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer,:] = np.sum(Pleaf_voxel[ind,:],axis = 0)

    return Psoil,Pleaf_layer

###### bidirectional

def proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0=0, G = 0.5, ifleaf = 1):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param std:  森林树密度
    :param hc:   冠层高度
    :param hcr:   垂直方向的树冠半径
    :param rcr:   水平方向的树冠半径
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param vaa0:   观测方位角
    :param sza0:   太阳天顶角
    :param saa0:   太阳方位角
    :param G:      投影系数
    :return:  组分比例
    '''

    if np.size(saa0) != np.size(sza0):
        saa0 = np.resize(saa0,np.size(sza0))

    hcc = hc - hcr
    hcb = hc - hcr * 2
    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    density = lai / (std * vol_crown)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    voxels = scheme.points
    ### scheme.weights 的累计是1，乘以体积，则是体素的真实的体积
    vo = scheme.weights *vol_crown
    number_voxel = np.size(vo)
    number_angle = np.size(vza0)
    number =  number_angle * number_voxel

    vsa0 = np.abs(saa0 - vaa0)
    vza = np.repeat(vza0,number_voxel)
    vaa = np.repeat(vaa0,number_voxel)
    sza = np.repeat(sza0,number_voxel)
    saa = np.repeat(saa0,number_voxel)
    vsa = np.repeat(vsa0,number_voxel)
    ### 冠层的体素坐标和体积
    x = np.tile(voxels[0,:] * rcr,number_angle)
    y = np.tile(voxels[1,:] * rcr,number_angle)
    z = np.tile(voxels[2,:] * hcr,number_angle)
    v = np.tile(vo, number_angle)
    ### 土壤的体素的位置
    x_soil = np.asarray([0])
    y_soil = np.asarray([0])
    z_soil = np.asarray([0])

    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    sphiv = np.sin(np.deg2rad(vaa))
    cphiv = np.cos(np.deg2rad(vaa))
    tgthv = np.tan(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    sphis = np.sin(np.deg2rad(saa))
    cphis = np.cos(np.deg2rad(saa))
    tgths = np.tan(np.deg2rad(sza))


    ### 观测方向的透过率 = 观测方向树冠内透过率 * 观测方向树冠间/树冠外的透过率
    gapv_inside_up,plv_inside_up,gapv_inside_down,plv_inside_down = gap_probability_crown_inside_voxel(x,y,z,density,hcr,rcr,vza,vaa)
    xv_up = x + plv_inside_up * sthetv * cphiv
    yv_up = y + plv_inside_up * sthetv * sphiv
    zv_up = z + plv_inside_up * cthetv + hcr
    gapv_outside_up,plv_outside_up,upAreav,hcrv_up,interv = \
        gap_probability_crown_outside_voxel(std, xv_up, yv_up, zv_up, density, hc, hcr,rcr, vza, vaa)
    gapv = gapv_inside_up * gapv_outside_up

    ### 太阳方向的透过率 = 太阳方向树冠内透过率 * 太阳方向树冠间/树冠外的透过率
    gaps_inside_up,pls_inside_up,gaps_inside_down,pls_inside_down = gap_probability_crown_inside_voxel(x,y,z,density,hcr,rcr,sza,saa)
    xs_up = x + pls_inside_up * sthets * cphis
    ys_up = y + pls_inside_up * sthets * sphis
    zs_up = z + pls_inside_up * cthets + hcr
    gaps_outside_up,pls_outside_up,upAreas,hcrs_up,inters = \
        gap_probability_crown_outside_voxel(std, xs_up, ys_up, zs_up, density, hc, hcr,rcr, sza, saa)
    gaps = gaps_inside_up * gaps_outside_up

    ### 树冠内热点的集成
    gapvs_inside_up,pl_hotspot = hotspot_path_length(plv_inside_up,pls_inside_up,density,hspot,vza,sza,vsa)
    gapvs_outside_up = hotspot_projection_area(std, density, hspot,
                                               xv_up, yv_up, zv_up, plv_outside_up, upAreav, hcrv_up, interv,
                                               xs_up, ys_up, zs_up, pls_outside_up, upAreas, hcrs_up, inters,
                                               vza, sza, vsa)
    gapvs = gapvs_inside_up * gapvs_outside_up

    ### 体积密度和遮挡
    base_voxel = std * v * density * G / cthetv
    Pleaf_voxel =  gapv * base_voxel
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Pleaf = np.sum(Pleaf_voxel, axis=0)
    Pleaf_sunlit_voxel = gapvs * base_voxel
    Pleaf_sunlit_voxel = np.transpose(np.reshape(Pleaf_sunlit_voxel, [number_angle, number_voxel]))
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)

    ### 这里采用树冠的结果，而不是under，没有用到x和y坐标信息，即不考虑自身树冠的遮挡
    gaps_soil, pls_outside_soil, upAreas_soil, hcrs_soil, inters_soil = \
        gap_probability_crown_outside_voxel(std, x_soil, y_soil, z_soil, density, hc, hcr, rcr, sza0, saa0)
    gapv_soil, plv_outside_soil, upAreav_soil, hcrv_soil, interv_soil = \
        gap_probability_crown_outside_voxel(std, x_soil, y_soil, z_soil, density, hc, hcr, rcr, vza0, vaa0)
    gapvs_soil = hotspot_projection_area(std, density, hspot, x_soil, y_soil, z_soil, plv_outside_soil, upAreav_soil, hcrv_soil, interv_soil,
                                         x_soil, y_soil, z_soil, pls_outside_soil, upAreas_soil, hcrs_soil, inters_soil, vza0, sza0, vsa0)


    ### 土壤的可视比例
    Psoil = gapv_soil
    Psoil_sunlit = gapvs_soil

    if ifleaf == 0:
        ### 比例偏差通过土壤平差
        Psoil = 1 - Pleaf
        Psoil_sunlit = gapvs_soil/gapv_soil * Psoil
    else:
        ### 比例偏差通过叶片平差
        Pleafnew = 1-Psoil
        Kcorrect = Pleafnew / Pleaf
        Pleaf_voxel  = Pleaf_voxel * Kcorrect
        Pleaf_sunlit_voxel = Pleaf_sunlit_voxel * Kcorrect

    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded_voxel = Pleaf_voxel - Pleaf_sunlit_voxel

    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_voxel,Pleaf_shaded_voxel

def proportion_bidirectional_hom_voxel(lai0,hspot, vza0, sza0, vsa0,  CI = 1.0, G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例 光照土壤，阴影土壤，光照植被，阴影植被
    '''
    number_voxel = 20
    number_angle = np.size(vza0)
    vza = np.repeat(vza0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    cthetv = np.cos(np.deg2rad(vza))
    # sthetv = np.sin(np.deg2rad(vza))
    # tgthv = np.tan(np.deg2rad(vza))
    # cthets = np.cos(np.deg2rad(sza))
    # sthets = np.sin(np.deg2rad(sza))
    # tgths = np.tan(np.deg2rad(sza))

    dlai = lai0 / number_voxel
    step = 1.0 / number_voxel
    lai = np.linspace(step * 0.5, 1 - step * 0.5, number_voxel) * lai0
    lai = np.tile(lai,number_angle)

    Gleaf_voxel = gap_probability_hom_voxel(lai, vza)
    Pleaf_voxel = Gleaf_voxel * dlai * G * CI / cthetv
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Pleaf = np.sum(Pleaf_voxel, axis = 0)

    Gleaf_sunlit_voxel = hotspot_hom_voxel(lai, hspot, vza, sza, vsa)
    Pleaf_sunlit_voxel = Gleaf_sunlit_voxel * dlai * G * CI / cthetv
    Pleaf_sunlit_voxel = np.transpose(np.reshape(Pleaf_sunlit_voxel, [number_angle, number_voxel]))
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis = 0)

    Psoil = gap_probability_hom_analytical(lai0, vza0)
    Psoil_sunlit = hotspot_layer(lai0, hspot, vza0, sza0, vsa0)
    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded_voxel = Pleaf_voxel - Pleaf_sunlit_voxel

    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_voxel,Pleaf_shaded_voxel

def proportion_bidirectional_row_voxel(lai, hspot,row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
						nvw = 5, nvh = 10, nvb = 5,  G = 0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    rw = row_width
    rb = row_blank
    rh = row_height
    rs = rw + rb
    nvw = nvw
    nvh = nvh
    nvb = nvb
    nvs = nvw + nvb
    ### 观测和太阳方位角差
    vsa0 = np.abs(saa0 - vaa0)
    ### 太阳和垄行方位角差
    sra0 = np.abs(saa0 - raa0) % 180
    ### 观测和垄行方方位角差
    vra0 = np.abs(vaa0 - raa0) % 180
    laie = lai * rs / rw
    ### 体密度
    density = laie / rh

    scale_w = rw / nvw
    scale_h = rh / nvh
    scale_s = rb / nvb

    length_layer_width = np.linspace(0.5, nvw - 0.5, nvw) * scale_w
    length_layer_blank = np.linspace(0.5, nvb - 0.5, nvb) * scale_s + rw
    height_layer = np.linspace(0.5, nvh - 0.5, nvh) * scale_h
    length_voxel_vegetation = np.tile(length_layer_width, (nvh))
    height_voxel_vegetation = np.asarray(np.tile(np.transpose(np.asmatrix(height_layer)), (1, nvw)))


    '''植被和背景体素集成'''
    length_voxel_vegetation = np.reshape(length_voxel_vegetation, -1)
    height_voxel_vegetation = np.reshape(height_voxel_vegetation, -1)
    length_voxel_background = np.hstack([length_layer_width, length_layer_blank])
    height_voxel_background = np.repeat(rh, nvs)
    # length_voxel = np.hstack([length_voxel_vegetation, length_voxel_background])
    # height_voxel = np.hstack([height_voxel_vegetation, height_voxel_background])
    number_voxel_vegetation = np.size(length_voxel_vegetation)
    number_voxel_background = np.size(length_voxel_background)
    number_angle = np.size(vza0)


    length_voxel_background = np.tile(length_voxel_background, number_angle)
    height_voxel_background = np.tile(height_voxel_background, number_angle)
    length_voxel_vegetation = np.tile(length_voxel_vegetation, number_angle)
    height_voxel_vegetation = np.tile(height_voxel_vegetation, number_angle)

    vza_vegetation = np.repeat(vza0, number_voxel_vegetation)
    sza_vegetation = np.repeat(sza0, number_voxel_vegetation)
    vra_vegetation = np.repeat(vra0, number_voxel_vegetation)
    sra_vegetation = np.repeat(sra0, number_voxel_vegetation)
    vsa_vegetation = np.repeat(vsa0, number_voxel_vegetation)

    vza_background = np.repeat(vza0, number_voxel_background)
    sza_background = np.repeat(sza0, number_voxel_background)
    vra_background = np.repeat(vra0, number_voxel_background)
    sra_background = np.repeat(sra0, number_voxel_background)
    vsa_background = np.repeat(vsa0, number_voxel_background)

    ### 观测方向和太阳方向的透过率，太阳方向的时候ifverse = 1
    Gvleaf_voxel,plv_leaf = gap_probability_row_voxel(density,row_width,row_blank,row_height,
                length_voxel_vegetation,height_voxel_vegetation,vza_vegetation,vra_vegetation,vsa_vegetation)
    Gsleaf_voxel,pls_leaf = gap_probability_row_voxel(density,row_width,row_blank,row_height,
                length_voxel_vegetation, height_voxel_vegetation,sza_vegetation, sra_vegetation, vsa_vegetation,1)
    ### 观测和太阳方向的热点
    Gvleaf_sunlit_voxel,plvs_leaf = hotspot_path_length(plv_leaf,pls_leaf,density,hspot,vza_vegetation,sza_vegetation,vsa_vegetation)
    cthetv = np.cos(np.deg2rad(vza_vegetation))
    #########################################################################################
    ##### ?????????????????????????????????????????????????????????????????????????????
    ########################################################################################
    ### 从透过率到组分比例
    # Pvleaf_voxel = Gvleaf_voxel * density * G /cthetv * (rw / nvw) *(rh/nvh)
    # Pvleaf_sunlit_voxel = Gvleaf_sunlit_voxel * density * G / cthetv * (rw/nvw) *(rh/nvh)
    Pvleaf_voxel = Gvleaf_voxel * density * G /cthetv * (1.0/nvs)*(rh/nvh)
    Pvleaf_sunlit_voxel = Gvleaf_sunlit_voxel * density * G / cthetv *(1.0/nvs)*(rh/nvh)

    Pvsoil_voxel,plv_soil = gap_probability_row_voxel(density,row_width,row_blank,row_height,
               length_voxel_background, height_voxel_background,vza_background,vra_background,vsa_background)
    Pssoil_voxel,pls_soil = gap_probability_row_voxel(density,row_width,row_blank,row_height,
                length_voxel_background, height_voxel_background,sza_background,sra_background,vsa_background,1)
    Pvsoil_sunlit_voxel,plvs_soil = hotspot_path_length(plv_soil, pls_soil, density, hspot, vza_background, sza_background,vsa_background)


    Pvleaf_sunlit_voxel = np.transpose(np.reshape(Pvleaf_sunlit_voxel, [number_angle, number_voxel_vegetation]))
    Pvleaf_sunlit = np.sum(Pvleaf_sunlit_voxel, axis=0)
    Pvleaf_voxel = np.transpose(np.reshape(Pvleaf_voxel, [number_angle, number_voxel_vegetation]))
    Pvleaf = np.sum(Pvleaf_voxel, axis=0)
    Pvsoil_voxel = np.transpose(np.reshape(Pvsoil_voxel, [number_angle, number_voxel_background]))
    Pvsoil = np.sum(Pvsoil_voxel, axis=0)/number_voxel_background
    Pvsoil_sunlit_voxel = np.transpose(np.reshape(Pvsoil_sunlit_voxel, [number_angle, number_voxel_background]))
    Pvsoil_sunlit = np.sum(Pvsoil_sunlit_voxel, axis=0)/number_voxel_background
    Pvsoil_shaded = Pvsoil - Pvsoil_sunlit
    Pvleaf_shaded_voxel = Pvleaf_voxel - Pvleaf_sunlit_voxel

    return Pvsoil_sunlit,Pvsoil_shaded,Pvleaf_sunlit_voxel,Pvleaf_shaded_voxel

def proportion_bidirectional_crown_voxel_one(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0=0, G = 0.5, ifleaf = 1):
    '''
    冠层组分比例的集成，光照/阴影组分的比例
    :param lai:
    :param std:
    :param hspot:
    :param hc:
    :param hcr:
    :param rcr:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param saa0:
    :param G:
    :param ifleaf:
    :return:
    '''
    Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0, G, ifleaf)

    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis = 0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis = 0)
    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

def proportion_bidirectional_hom_voxel_one(lai0,hspot, vza0, sza0, vsa0, CI = 1.0, G = 0.5):
    '''
    均质场景，光照/阴影组分的比例
    :param lai0:
    :param hspot:
    :param vza0:
    :param sza0:
    :param vsa0:
    :param CI:
    :param G:
    :return:
    '''
    Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_hom_voxel(lai0, hspot, vza0, sza0, vsa0,  CI, G)
    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis = 0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis = 0)
    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

def proportion_bidirectional_row_voxel_one(lai, hspot,row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
						nvw = 5, nvh = 10, nvb = 5,  G = 0.5):
    '''
    垄行场景，光照/阴影组分的比例
    :param lai:
    :param hspot:
    :param row_width:
    :param row_blank:
    :param row_height:
    :param vza0:
    :param vaa0:
    :param sza0:
    :param saa0:
    :param raa0:
    :param nvw:
    :param nvh:
    :param nvb:
    :param G:
    :return:
    '''
    Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_voxel, Pleaf_shaded_voxel  = \
        proportion_bidirectional_row_voxel(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                       nvw, nvh, nvb, G)
    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis=0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)

    return Psoil_sunlit,Psoil_shaded, Pleaf_sunlit, Pleaf_shaded

def proportion_bidirectional_crown_voxel_layer(
        lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, number_layer = 5, saa0=0, G=0.5, ifleaf=1):
    '''
    组分比例的分层结果，通过垂直高度数值计算
    :param lai:
    :param std:
    :param hspot:
    :param hc:
    :param hcr:
    :param rcr:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param number_layer:
    :param saa0:
    :param G:
    :param ifleaf:
    :return:
    '''


    Psoil_sunlit,Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0, G, ifleaf)

    number_angle = np.size(vza0)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    vz = scheme.points[2,:]
    number_voxel = np.size(vz)
    Pleaf_shaded_layer = np.zeros([number_layer,number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer,number_angle])
    dlh = 1 * 2 / number_layer
    index_layer = (vz+1) // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer,:] = np.sum(Pleaf_shaded_voxel[ind,:],axis = 0)
            Pleaf_sunlit_layer[klayer,:] = np.sum(Pleaf_sunlit_voxel[ind,:],axis = 0)

    return Psoil_sunlit,Psoil_shaded, Pleaf_sunlit_layer,Pleaf_shaded_layer


def proportion_bidirectional_hom_voxel_layer(lai0, hspot, vza0, sza0, vsa0, number_layer = 5, CI=1.0, G=0.5):
    '''
    组分比例的分层结果，通过垂直高度数值计算
    :param lai0:
    :param hspot:
    :param vza0:
    :param sza0:
    :param vsa0:
    :param number_layer:
    :param CI:
    :param G:
    :return:
    '''
    Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_hom_voxel(lai0, hspot, vza0, sza0, vsa0, CI, G)


    [number_voxel,number_angle] = np.shape(Pleaf_sunlit_voxel)
    # number_angle = np.size(vza0)
    vz = np.linspace(0.5,number_voxel-0.5,number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    Pleaf_shaded_layer = np.zeros([number_layer,number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer,number_angle])
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer,:] = np.sum(Pleaf_shaded_voxel[ind,:], axis = 0)
            Pleaf_sunlit_layer[klayer,:] = np.sum(Pleaf_sunlit_voxel[ind,:], axis = 0)

    return Psoil_sunlit,Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer


def proportion_bidirectional_row_voxel_layer(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                             number_layer = 5, nvw=5, nvh=10, nvb=5, G=0.5):
    '''
       组分比例的分层结果，通过垂直高度数值计算
    :param lai:
    :param hspot:
    :param row_width:
    :param row_blank:
    :param row_height:
    :param vza0:
    :param vaa0:
    :param sza0:
    :param saa0:
    :param raa0:
    :param number_layer:
    :param nvw:
    :param nvh:
    :param nvb:
    :param G:
    :return:
    '''
    Psoil_sunlit_voxel, Psoil_shaded_voxel, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_row_voxel(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0, nvw, nvh, nvb, G)

    number_voxel = nvh
    number_angle = np.size(vza0)

    Psoil_shaded = np.sum(Psoil_shaded_voxel, axis=0)
    Psoil_sunlit = np.sum(Psoil_sunlit_voxel, axis=0)
    Pleaf_shaded_layer = np.zeros([number_layer,number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer,number_angle])
    ### 分层划分，取整
    vz = np.linspace(0.5,number_voxel-0.5,number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer,:] = np.sum(Pleaf_shaded_voxel[ind,:],axis = 0)
            Pleaf_sunlit_layer[klayer,:] = np.sum(Pleaf_sunlit_voxel[ind,:],axis = 0)

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer

def proportion_bidirectional_trunk_voxel_individual(lai,std,hspot,hc,hcr,rcr,hc_trunk,dbh_trunk,vza0,sza0,vaa0, saa0 = 0):
    hcc = hc - hcr
    hcb = hc - hcr * 2
    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    area_crown = np.pi * rcr * rcr
    density = lai / (std * vol_crown)
    number_trunk = 10
    scale_trunk = hc_trunk / number_trunk
    number_angle = np.size(vza0)
    zt = np.linspace(0.5, number_trunk - 0.5, number_trunk) * scale_trunk
    xt = np.zeros(number_trunk * number_angle)
    yt = np.zeros(number_trunk * number_angle)
    zt = np.tile(zt, number_angle)
    vsa0 = np.abs(saa0 - vaa0)
    vza = np.repeat(vza0, number_trunk)
    vaa = np.repeat(vaa0, number_trunk)
    sza = np.repeat(sza0, number_trunk)
    saa = np.repeat(saa0, number_trunk)
    vsa = np.repeat(vsa0, number_trunk)

    cthetv = np.cos(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    sphiv = np.sin(np.deg2rad(vaa))
    cphiv = np.cos(np.deg2rad(vaa))
    tgthv = np.tan(np.deg2rad(vza))
    cthets = np.cos(np.deg2rad(sza))
    sthets = np.sin(np.deg2rad(sza))
    sphis = np.sin(np.deg2rad(saa))
    cphis = np.cos(np.deg2rad(saa))
    tgths = np.tan(np.deg2rad(sza))

    #### 观测方向
    gapv_inside_up, plv_inside_up = gap_probability_under_inside_voxel(xt, yt, zt, density, hcr, rcr, vza, vaa)
    xv_up = xt + plv_inside_up * sthetv * cphiv
    yv_up = yt + plv_inside_up * sthetv * sphiv
    zv_up = zt + plv_inside_up * cthetv + hc - hcr
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv, stemProjv = \
        gap_probability_under_outside_voxel(std, xv_up, yv_up, zv_up, density, hc, hcr, rcr, dbh_trunk, hc_trunk, vza,
                                            vaa)
    gapv = gapv_inside_up * gapv_outside_up

    #### 光照方向
    gaps_inside_up, pls_inside_up = gap_probability_under_inside_voxel(xt, yt, zt, density, hcr, rcr, sza, saa)
    xs_up = xt + pls_inside_up * sthets * cphis
    ys_up = yt + pls_inside_up * sthets * sphis
    zs_up = zt + pls_inside_up * cthets + hc - hcr
    gaps_outside_up, pls_outside_up, upAreas, hcrs_up, inters, stemProjs = \
        gap_probability_under_outside_voxel(std, xs_up, ys_up, zs_up, density, hc, hcr, rcr, dbh_trunk, hc_trunk, sza,
                                            saa)
    gaps = gaps_inside_up * gaps_outside_up

    ### 树冠内热点的集成
    gapvs_inside_up, pl_hotspot = hotspot_path_length(plv_inside_up, pls_inside_up, density, hspot, vza, sza, vsa)
    gapvs_outside_up = hotspot_projection_area_under(std, density, hspot,
                                                     xv_up, yv_up, zv_up, plv_outside_up, upAreav, hcrv_up, interv,
                                                     xs_up, ys_up, zs_up, pls_outside_up, upAreas, hcrs_up, inters,
                                                     hc, hcr, dbh_trunk, stemProjv, stemProjs,
                                                     vza, sza, vsa)
    gapvs = gapvs_inside_up * gapvs_outside_up
    phi = (np.pi - np.deg2rad(vsa)) / np.pi
    Ptrunk_voxel = gapv * dbh_trunk * std * tgthv * scale_trunk
    Ptrunk_sunlit_voxel = gapvs * dbh_trunk * std * phi * tgthv * scale_trunk
    Ptrunk_voxel = np.transpose(np.reshape(Ptrunk_voxel, [number_angle, number_trunk]))
    # Ptrunk = np.sum(Ptrunk_voxel, axis=0)
    Ptrunk_sunlit_voxel = np.transpose(np.reshape(Ptrunk_sunlit_voxel, [number_angle, number_trunk]))
    # Ptrunk_sunlit = np.sum(Ptrunk_sunlit_voxel, axis=0)

    #############################################################################################################
    ##### 背景场
    zb = np.tile(np.asarray([0]),number_angle)
    xb = np.tile(np.asarray([0]),number_angle)
    yb = np.tile(np.asarray([0]),number_angle)

    #### 观测方向
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv, stemProjv = \
        gap_probability_under_outside_voxel(std, xb, yb, zb, density, hc, hcr, rcr,dbh_trunk,hc_trunk, vza0, vaa0)
    Psoil = gapv_outside_up

    #### 光照方向
    gaps_outside_up, pls_outside_up, upAreas, hcrs_up, inters, stemProjs = \
        gap_probability_under_outside_voxel(std, xb, yb, zb, density, hc, hcr, rcr,dbh_trunk,hc_trunk,  sza0, saa0)
    gaps = gaps_outside_up

    ### 树冠内热点的集成
    Psoil_sunlit = hotspot_projection_area_under(std, density, hspot,
                                               xb, yb, zb, plv_outside_up, upAreav, hcrv_up, interv,
                                               xb, yb, zb, pls_outside_up, upAreas, hcrs_up, inters,
                                               hc, hcr, dbh_trunk, stemProjv, stemProjs,
                                               vza0, sza0, vsa0)

    Psoil_shaded = Psoil - Psoil_sunlit
    Ptrunk_shaded_voxel = Ptrunk_voxel - Ptrunk_sunlit_voxel

    return Psoil_sunlit, Psoil_shaded, Ptrunk_sunlit_voxel, Ptrunk_shaded_voxel

def proportion_bidirectional_trunk_voxel_one_individual(lai,std,hspot,hc,hcr,rcr,hc_trunk,dbh_trunk,vza0,sza0,vaa0, saa0 = 0):
    Psoil_sunlit,Psoil_shaded,Ptrunk_sunlit_voxel, Ptrunk_shaded_voxel = \
        proportion_bidirectional_trunk_voxel_individual(lai, std, hspot, hc, hcr, rcr, hc_trunk, dbh_trunk, vza0, sza0,vaa0,saa0)
    Ptrunk_shaded = np.sum(Ptrunk_shaded_voxel,axis = 0)
    Ptrunk_sunlit = np.sum(Ptrunk_sunlit_voxel,axis = 0)
    return Psoil_sunlit,Psoil_shaded,Ptrunk_sunlit,Ptrunk_shaded


def proportion_bidirectional_trunk():
    pass


###############################
### 树冠+树干 的个体 + 土壤
###############################

def proportion_bidirectional_crown_voxel_one_individual(lai, std, hspot, hc, hcr, rcr, hc_trunk, dbh_trunk, vza0, sza0, vaa0, iftrunk = 0, saa0=0, G = 0.5, ifleaf = 1):
    '''
    树冠+树干+土壤场景的组分可视比例
    :param lai:
    :param std:
    :param hspot:
    :param hc:
    :param hcr:
    :param rcr:
    :param hc_trunk:
    :param dbh_trunk:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param saa0:
    :param iftrunk:
    :param G:
    :param ifleaf:
    :return:
    '''

    Psoil_sunlit,Psoil_shaded, Pleaf_sunlit,Pleaf_shaded = \
        proportion_bidirectional_crown_voxel_one(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0)




    number_angle = np.size(vza0)
    Ptrunk_shaded = np.zeros(number_angle)
    Ptrunk_sunlit = np.zeros(number_angle)
    if iftrunk==1:
        ###############################################################################################################
        Psoil_sunlit,Psoil_shaded,Ptrunk_sunlit,Ptrunk_shaded = \
            proportion_bidirectional_trunk_voxel_one_individual(lai,std,hspot,hc,hcr,rcr,hc_trunk,dbh_trunk,vza0,sza0,vaa0,saa0)


    # Pcomsum = Psoil_sunlit + Psoil_shaded + Ptrunk_sunlit + Ptrunk_shaded
    # plt.plot(Pcomsum)
    # plt.show()

    #############################################
    ##### 这里有一点偏差，需要后期定位
    ##############################################


    Pleaf = Pleaf_shaded + Pleaf_sunlit
    Ptrunk = Ptrunk_sunlit + Ptrunk_shaded
    Psoil = Psoil_sunlit + Psoil_shaded

    Ksoil_sunlit = Psoil_sunlit/Psoil
    Ksoil_shaded = Psoil_shaded/Psoil
    Psoil = 1- Pleaf - Ptrunk
    Psoil_shaded = Psoil * Ksoil_shaded
    Psoil_sunlit = Psoil * Ksoil_sunlit

    return Psoil_sunlit,Psoil_shaded, Pleaf_sunlit,Pleaf_shaded,Ptrunk_sunlit,Ptrunk_shaded

def proportion_bidirectional_crown_voxel_layer_individual(lai, std, hspot, hc, hcr,
                                                          rcr, hc_trunk, dbh_trunk,
                                                          vza0, sza0, vaa0, number_layer, iftrunk = 0, saa0=0, G = 0.5, ifleaf = 1):
    '''
    树冠+树干+土壤场景的组分可视比例
    :param lai:
    :param std:
    :param hspot:
    :param hc:
    :param hcr:
    :param rcr:
    :param hc_trunk:
    :param dbh_trunk:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param saa0:
    :param iftrunk:
    :param G:
    :param ifleaf:
    :return:
    '''

    Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_layer,Pleaf_shaded_layer =  \
        proportion_bidirectional_crown_voxel_layer(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, number_layer, saa0)

    #

    Ptrunk_shaded = 0
    Ptrunk_sunlit = 0
    if iftrunk == 1:
        ###############################################################################################################
        Psoil_sunlit, Psoil_shaded, Ptrunk_sunlit, Ptrunk_shaded = \
            proportion_bidirectional_trunk_voxel_one_individual(lai, std, hspot, hc, hcr, rcr, hc_trunk, dbh_trunk,
                                                                vza0, sza0, vaa0, saa0)
    Pleaf_shaded = np.sum(Pleaf_shaded_layer,axis=0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_layer,axis = 0)
    Pleaf = Pleaf_shaded + Pleaf_sunlit
    Ptrunk = Ptrunk_sunlit + Ptrunk_shaded
    Psoil = Psoil_sunlit + Psoil_shaded
    Ksoil_sunlit = Psoil_sunlit/Psoil
    Ksoil_shaded = Psoil_shaded/Psoil

    Psoil = 1- Pleaf - Ptrunk
    Psoil_shaded = Psoil * Ksoil_shaded
    Psoil_sunlit = Psoil * Ksoil_sunlit


    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer, Ptrunk_sunlit, Ptrunk_shaded


def proportion_bidirectional_hom_voxel_one_individual(lai0, hspot, hc_trunk, dbh_trunk, std_trunk, vza0, sza0, vaa0, iftrunk=0, saa0=0, number_voxel = 10, CI = 1.0, G = 0.5):
    '''
    计算均质冠层+树干+土壤类型场景的组分可视比例
    :param lai0:
    :param hspot:
    :param hc_trunk:
    :param dbh_trunk:
    :param std_trunk:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param saa0:
    :param iftrunk:
    :param number_voxel:
    :param CI:
    :param G:
    :return:
    '''

    vsa0 = np.abs(vaa0 - saa0)
    Psoil_sunlit,Psoil_shaded, Pleaf_sunlit,Pleaf_shaded =  proportion_bidirectional_hom_voxel_one(lai0, hspot, vza0, sza0, vsa0)

    Psoil = Psoil_sunlit + Psoil_shaded
    Ptrunk_shaded = 0
    Ptrunk_sunlit = 0
    Ptrunk = 0
    if iftrunk == 1:
        Gtrunk = gap_probability_trunk_analytical(hc_trunk,dbh_trunk,std_trunk,vza0)
        phi = (np.pi - np.deg2rad(vsa0)) / np.pi
        Ptrunk = (1-Gtrunk) * Psoil
        Ptrunk_sunlit = (Ptrunk * phi) * Psoil_sunlit
        Ptrunk_shaded = Ptrunk - Ptrunk_sunlit
        Psoil = Psoil * Gtrunk
        Psoil_sunlit = Psoil_sunlit * Gtrunk
        Psoil_shaded = Psoil - Psoil_sunlit

    # Pleaf =  Pleaf
    # Pleaf_sunlit =  Pleaf_sunlit
    # Psoil =  Psoil
    # Psoil_sunlit =  Psoil_sunlit

    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded, Ptrunk_sunlit,Ptrunk_shaded

def proportion_bidirectional_hom_voxel_layer_individual(lai0, hspot, hc_trunk, dbh_trunk, std_trunk, vza0, sza0, vaa0,number_layer = 5, iftrunk=0, saa0=0,  CI = 1.0, G = 0.5):
    '''
    计算均质冠层+树干+土壤类型场景的组分可视比例
    :param lai0:
    :param hspot:
    :param hc_trunk:
    :param dbh_trunk:
    :param std_trunk:
    :param vza0:
    :param sza0:
    :param vaa0:
    :param saa0:
    :param iftrunk:
    :param number_voxel:
    :param CI:
    :param G:
    :return:
    '''

    vsa0 = np.abs(vaa0 - saa0)
    Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_layer,Pleaf_shaded_layer =  proportion_bidirectional_hom_voxel_layer(lai0, hspot, vza0, sza0, vsa0,number_layer)

    Ptrunk = 0
    Psoil = Psoil_sunlit + Psoil_shaded

    Ptrunk_sunlit = 0
    if iftrunk == 1:
        Gtrunk = gap_probability_trunk_analytical(hc_trunk,dbh_trunk,std_trunk,vza0)
        phi = (np.pi - np.deg2rad(vsa0)) / np.pi
        Ptrunk = (1-Gtrunk) * Psoil
        Ptrunk_sunlit = (Ptrunk * phi) * Psoil_sunlit
        Ptrunk_shaded = Ptrunk - Ptrunk_sunlit
        Psoil = Psoil * Gtrunk
        Psoil_sunlit = Psoil_sunlit * Gtrunk
        Psoil_shaded = Psoil - Psoil_sunlit

    # Pleaf =  Pleaf
    # Pleaf_sunlit =  Pleaf_sunlit
    # Psoil =  Psoil
    # Psoil_sunlit =  Psoil_sunlit

    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit_layer,Pleaf_shaded_layer,Ptrunk_sunlit,Ptrunk_shaded


###############################
### 树冠+树干+低层植被 的个体
###############################


def proportion_bidirectional_hom_voxel_one_endmember(lai_crown, hspot_crown, hc_trunk, dbh_trunk, std_trunk,
                                                     lai_unveg, hspot_unveg, vza, sza, vaa, saa=0, G_crown = 0.5, G_unveg = 0.5):
    '''
    端元上的组分可视化比例，包括上层的均质场景+下层的均质场景
    :param lai_crown:
    :param hspot_crown:
    :param hc_trunk:
    :param dbh_trunk:
    :param std_trunk:
    :param lai_unveg:
    :param hspot_unveg:
    :param vza:
    :param sza:
    :param vaa:
    :param iftrunk:
    :param saa:
    :param G_crown:
    :param G_unveg:
    :return:
    '''

    iftrunk = 0
    if dbh_trunk * hc_trunk > 0.1: iftrunk = 1
    Pbackground_sunlit,Pbackground_shaded,Pcrown_sunlit,Pcrown_shaded,Ptrunk_sunlit,Ptrunk_shaded = \
        proportion_bidirectional_hom_voxel_one_individual(lai_crown, hspot_crown, hc_trunk, dbh_trunk, std_trunk, vza, sza, vaa, iftrunk, saa)

    Pbackground = Pbackground_shaded + Pbackground_sunlit
    vsa = np.abs(vaa - saa)
    Psoil_sunlit,Psoil_shaded,Punveg_sunlit,Punveg_shaded = \
        proportion_bidirectional_hom_voxel_one(lai_unveg, hspot_unveg, vza, sza, vsa)
    Psoil = Psoil_sunlit + Psoil_shaded
    Punveg = Punveg_shaded + Punveg_sunlit
    Psoil_sunlit = Pbackground_sunlit * Psoil_sunlit
    Psoil_shaded = Pbackground * Psoil - Psoil_sunlit
    Punveg_sunlit = Pbackground_sunlit * Punveg_sunlit
    Punveg_shaded = Pbackground * Punveg - Punveg_sunlit

    return  Psoil_sunlit,Psoil_shaded,Pcrown_sunlit,Pcrown_shaded,Punveg_sunlit,Punveg_shaded,Ptrunk_sunlit,Ptrunk_shaded


def proportion_bidirectional_hom_voxel_layer_endmember(lai_crown, hspot_crown, hc_trunk, dbh_trunk, std_trunk,
                                                       lai_unveg, hspot_unveg, vza, sza, vaa, number_layer = 5, saa=0, G_crown = 0.5, G_unveg = 0.5):
    '''
    端元上的组分可视化比例，包括上层的均质场景+下层的均质场景
    :param lai_crown:
    :param hspot_crown:
    :param hc_trunk:
    :param dbh_trunk:
    :param std_trunk:
    :param lai_unveg:
    :param hspot_unveg:
    :param vza:
    :param sza:
    :param vaa:
    :param iftrunk:
    :param saa:
    :param G_crown:
    :param G_unveg:
    :return:
    '''

    iftrunk = 0
    if dbh_trunk * hc_trunk > 0.1: iftrunk = 1
    Pbackground_sunlit, Pbackground_shaded, Pcrown_sunlit_layer, Pcrown_shaded_layer, Ptrunk_sunlit, Ptrunk_shaded = \
        proportion_bidirectional_hom_voxel_layer_individual(lai_crown, hspot_crown, hc_trunk, dbh_trunk, std_trunk, vza,
                                                          sza, vaa, number_layer, iftrunk, saa)

    Pbackground = Pbackground_shaded + Pbackground_sunlit
    vsa = np.abs(vaa - saa)
    Psoil_sunlit,Psoil_shaded,Punveg_sunlit,Punveg_shaded = \
        proportion_bidirectional_hom_voxel_one(lai_unveg, hspot_unveg, vza, sza, vsa)
    Psoil = Psoil_sunlit + Psoil_shaded
    Punveg = Punveg_shaded + Punveg_sunlit
    Psoil_sunlit = Pbackground_sunlit * Psoil_sunlit
    Psoil_shaded = Pbackground * Psoil - Psoil_sunlit
    Punveg_sunlit = Pbackground_sunlit * Punveg_sunlit
    Punveg_shaded = Pbackground * Punveg - Punveg_sunlit

    return Psoil_sunlit, Psoil_shaded, Pcrown_sunlit_layer, Pcrown_shaded_layer, Punveg_sunlit, Punveg_shaded, Ptrunk_sunlit, Ptrunk_shaded


def proportion_bidirectional_crown_voxel_one_endmember(lai_crown, std_crown, hspot_crown, hc_crown, hcr_crown, rcr_crown, hc_trunk, dbh_trunk,
                                                       lai_unveg, hspot_unveg, vza, sza, vaa, saa=0, G_crown = 0.5, G_unveg = 0.5):
    '''
    端元上的组分可视化比例，包括上层的树冠场景+下层的均质场景
    :param lai_crown:
    :param std_crown:
    :param hspot_crown:
    :param hc_crown:
    :param hcr_crown:
    :param rcr_crown:
    :param hc_trunk:
    :param dbh_trunk:
    :param lai_unveg:
    :param hspot_unveg:
    :param vza:
    :param sza:
    :param vaa:
    :param iftrunk:
    :param saa:
    :param G_crown:
    :param G_unveg:
    :return:
    '''

    iftrunk = 0
    if dbh_trunk * hc_trunk > 0.1: iftrunk = 1
    Pbackground_sunlit,Pbackground_shaded,Pcrown_sunlit,Pcrown_shaded,Ptrunk_sunlit,Ptrunk_shaded = \
        proportion_bidirectional_crown_voxel_one_individual(lai_crown, std_crown, hspot_crown, hc_crown, hcr_crown, rcr_crown, hc_trunk, dbh_trunk,
                                                            vza, sza, vaa, iftrunk, saa)



    Pbackground = Pbackground_shaded + Pbackground_sunlit
    vsa = np.abs(vaa - saa)
    Psoil_sunlit,Psoil_shaded,Punveg_sunlit,Punveg_shaded = \
        proportion_bidirectional_hom_voxel_one(lai_unveg, hspot_unveg, vza, sza, vsa)
    Psoil = Psoil_sunlit + Psoil_shaded
    Punveg = Punveg_shaded + Punveg_sunlit
    Psoil_sunlit = Pbackground_sunlit * Psoil_sunlit
    Psoil_shaded = Pbackground * Psoil - Psoil_sunlit
    Punveg_sunlit = Pbackground_sunlit * Punveg_sunlit
    Punveg_shaded = Pbackground * Punveg - Punveg_sunlit

    return  Psoil_sunlit,Psoil_shaded,Pcrown_sunlit,Pcrown_shaded,Punveg_sunlit,Punveg_shaded,Ptrunk_sunlit,Ptrunk_shaded


def proportion_bidirectional_crown_voxel_layer_endmember(lai_crown, std_crown, hspot_crown, hc_crown, hcr_crown,
                                                         rcr_crown, hc_trunk, dbh_trunk,
                                                         lai_unveg, hspot_unveg, vza, sza, vaa, number_layer = 5, saa=0,
                                                         G_crown = 0.5, G_unveg = 0.5):
    '''
    端元上的组分可视化比例，包括上层的树冠场景+下层的均质场景
    :param lai_crown:
    :param std_crown:
    :param hspot_crown:
    :param hc_crown:
    :param hcr_crown:
    :param rcr_crown:
    :param hc_trunk:
    :param dbh_trunk:
    :param lai_unveg:
    :param hspot_unveg:
    :param vza:
    :param sza:
    :param vaa:
    :param iftrunk:
    :param saa:
    :param G_crown:
    :param G_unveg:
    :return:
    '''

    iftrunk = 0
    if dbh_trunk * hc_trunk > 0.1: iftrunk = 1
    Pbackground_sunlit, Pbackground_shaded, Pcrown_sunlit_layer, Pcrown_shaded_layer, Ptrunk_sunlit, Ptrunk_shaded = \
        proportion_bidirectional_crown_voxel_layer_individual(lai_crown, std_crown, hspot_crown, hc_crown, hcr_crown,
                                                            rcr_crown, hc_trunk, dbh_trunk,
                                                            vza, sza, vaa, number_layer, iftrunk, saa)

    Pbackground = Pbackground_shaded + Pbackground_sunlit
    vsa = np.abs(vaa - saa)
    Psoil_sunlit,Psoil_shaded,Punveg_sunlit,Punveg_shaded = \
        proportion_bidirectional_hom_voxel_one(lai_unveg, hspot_unveg, vza, sza, vsa)
    Psoil = Psoil_sunlit + Psoil_shaded
    Punveg = Punveg_shaded + Punveg_sunlit
    Psoil_sunlit = Pbackground_sunlit * Psoil_sunlit
    Psoil_shaded = Pbackground * Psoil - Psoil_sunlit
    Punveg_sunlit = Pbackground_sunlit * Punveg_sunlit
    Punveg_shaded = Pbackground * Punveg - Punveg_sunlit


    return Psoil_sunlit, Psoil_shaded, Pcrown_sunlit_layer, Pcrown_shaded_layer, Punveg_sunlit, Punveg_shaded, Ptrunk_sunlit, Ptrunk_shaded

