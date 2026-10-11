from functools import lru_cache

import numpy as np
import scipy.integrate as sci
from .gap import *
from .hotspot import *
from .path import *
# from scatter import *
import matplotlib.pyplot as plt

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



def _spheroid_interval(p, d, hcr, rcr):
    """Distances t >= 0 where rays p + t*d are inside the spheroid x^2/rcr^2 + y^2/rcr^2 + z^2/hcr^2 <= 1."""
    scale = np.array([1.0 / rcr, 1.0 / rcr, 1.0 / hcr])
    ps, ds = p * scale, d * scale
    a = ds @ ds
    b = ps @ ds
    disc = b * b - a * (np.sum(ps * ps, axis=-1) - 1.0)
    root = np.sqrt(np.maximum(disc, 0.0))
    t1 = np.maximum((-b - root) / a, 0.0)
    t2 = np.maximum((-b + root) / a, t1)
    hit = disc > 0.0
    return np.where(hit, t1, 0.0), np.where(hit, t2, 0.0)


def _crown_depths(p, view, sun, extinction, kappa, hcr, rcr):
    """Optical depths through one crown of the view ray, the sun ray and both together, rays leaving offsets p."""
    (dv, cv), (ds, cs) = view, sun
    t1v, t2v = _spheroid_interval(p, dv, hcr, rcr)
    t1s, t2s = _spheroid_interval(p, ds, hcr, rcr)
    tau_v, tau_s = extinction * (t2v - t1v), extinction * (t2s - t1s)
    # Leaf-scale (Kuusk) correlation where both rays are inside the crown, decaying with height above the ray origin.
    low = np.maximum(t1v * cv, t1s * cs)
    high = np.maximum(np.minimum(t2v * cv, t2s * cs), low)
    shared = (np.exp(-kappa * low) - np.exp(-kappa * high)) / kappa if kappa > 1.0e-12 else high - low
    overlap = np.minimum(extinction / np.sqrt(cv * cs) * shared, np.minimum(tau_v, tau_s))
    return tau_v, tau_s, tau_v + tau_s - overlap


def _crown_offsets(z, centre_height, hcr, rcr, tan_x, direction, n_cross=24, n_max=96):
    """Per height z: crown-centre offsets (as p = ray origin - crown centre) that can intersect a ray leaving z, and the cell area."""
    h_low, h_high = np.maximum(centre_height - hcr - z, 0.0), centre_height + hcr - z
    s_low, span = h_low * tan_x - rcr, (h_high - h_low) * tan_x + 2.0 * rcr
    n_long = int(np.clip(round(n_cross * span.max() / (2.0 * rcr)), n_cross, n_max))
    S = (s_low[:, None] + (np.arange(n_long) + 0.5) / n_long * span[:, None])[:, :, None]
    N = rcr * ((np.arange(n_cross) + 0.5) * 2.0 / n_cross - 1.0)[None, None, :]
    S, N = np.broadcast_arrays(S, N)
    Z = np.broadcast_to((z - centre_height)[:, None, None], S.shape)
    p = np.stack([-(S * direction[0] - N * direction[1]), -(S * direction[1] + N * direction[0]), Z], axis=-1)
    return p, span / n_long * 2.0 * rcr / n_cross


_CROWN_LEVELS, _CROWN_LEVEL_WEIGHTS = np.polynomial.legendre.leggauss(8)
_CROWN_RINGS = np.sqrt((np.arange(5) + 0.5) / 5)
_CROWN_PHI = (np.arange(12) + 0.5) * 2.0 * np.pi / 12


@lru_cache(maxsize=8192)
def _crown_fractions(lai, std, hspot, hcr, rcr, hcb, vza, sza, vsa, G):
    if lai <= 0.0 or std <= 0.0:
        return 1.0, 0.0, 0.0, 0.0
    rd = np.pi / 180.0
    centre = hcb + hcr
    extinction = G * lai / (std * 4.0 / 3.0 * np.pi * rcr * rcr * hcr)
    cv, cs = np.cos(vza * rd), np.cos(sza * rd)
    tan_v, tan_s = np.tan(vza * rd), np.tan(sza * rd)
    ev, es = np.array([np.cos(vsa * rd), np.sin(vsa * rd)]), np.array([1.0, 0.0])
    view = (np.array([np.sin(vza * rd) * ev[0], np.sin(vza * rd) * ev[1], cv]), cv)
    sun = (np.array([np.sin(sza * rd), 0.0, cs]), cs)
    # hspot = leaf size / crown depth; SAIL scaling of the correlation length, as in hotspot_layer
    kappa = np.hypot(*(tan_v * ev - tan_s * es)) * 2.0 / (G / cv + G / cs) / max(hspot * 2.0 * hcr, 1.0e-9)

    # Integrals over Poisson crown positions of (1 - transmittance), for the view ray and the view+sun pair leaving each height:
    # the ground first, then the Gauss levels of a crown.
    heights = np.concatenate([[0.0], centre + _CROWN_LEVELS * hcr])
    pv, area_v = _crown_offsets(heights, centre, hcr, rcr, tan_v, ev)
    ps, area_s = _crown_offsets(heights, centre, hcr, rcr, tan_s, es)
    # 1 - Tv - Ts + Tvs is non-zero only for crowns hit by both rays, so it is summed over the shorter strip only
    if tan_v <= tan_s:
        tau_v, tau_s, tau_vs = _crown_depths(pv, view, sun, extinction, kappa, hcr, rcr)
        j_v = np.sum(-np.expm1(-tau_v), axis=(1, 2)) * area_v
        t1, t2 = _spheroid_interval(ps, sun[0], hcr, rcr)
        j_s = np.sum(-np.expm1(-extinction * (t2 - t1)), axis=(1, 2)) * area_s
        area = area_v
    else:
        tau_v, tau_s, tau_vs = _crown_depths(ps, view, sun, extinction, kappa, hcr, rcr)
        j_s = np.sum(-np.expm1(-tau_s), axis=(1, 2)) * area_s
        t1, t2 = _spheroid_interval(pv, view[0], hcr, rcr)
        j_v = np.sum(-np.expm1(-extinction * (t2 - t1)), axis=(1, 2)) * area_v
        area = area_s
    j_vs = j_v + j_s - np.sum(-np.expm1(-tau_v) - np.expm1(-tau_s) + np.expm1(-tau_vs), axis=(1, 2)) * area

    soil = float(gap_probability_crown_analytical(lai, std, hcr, rcr, vza, G))
    soil_sunlit = soil * float(np.exp(-std * (j_vs[0] - j_v[0])))

    # Leaves: own crown transmittance times the Poisson gap of the other crowns, integrated over the crown volume.
    radius = rcr * np.sqrt(np.maximum(1.0 - _CROWN_LEVELS ** 2, 0.0))
    ring = radius[:, None] * _CROWN_RINGS[None, :]
    x = (ring[:, :, None] * np.cos(_CROWN_PHI)).reshape(radius.size, -1)
    y = (ring[:, :, None] * np.sin(_CROWN_PHI)).reshape(radius.size, -1)
    own = np.stack([x, y, np.broadcast_to((_CROWN_LEVELS * hcr)[:, None], x.shape)], axis=-1)
    own_v, _, own_vs = _crown_depths(own, view, sun, extinction, kappa, hcr, rcr)
    volume = np.pi * radius ** 2 * _CROWN_LEVEL_WEIGHTS * hcr
    with np.errstate(divide="ignore"):
        log_v = np.log(volume * np.mean(np.exp(-own_v), axis=1)) - std * j_v[1:]
        log_vs = np.log(volume * np.mean(np.exp(-own_vs), axis=1)) - std * j_vs[1:]
    shift = np.max(log_v)
    sunlit_share = min(float(np.sum(np.exp(log_vs - shift)) / np.sum(np.exp(log_v - shift))), 1.0)
    leaf = 1.0 - soil
    return soil_sunlit, soil - soil_sunlit, leaf * sunlit_share, leaf * (1.0 - sunlit_share)


def proportion_bidirectional_crown_one(lai, std, hspot, hcr, rcr, vza, sza, vsa, G = 0.5, hcb = None):
    '''
    离散树冠四组分比例：树冠为泊松分布的椭球（冠底高度 hcb，默认 hcr），冠内为浑浊介质。
    土壤和叶片的观测/太阳双向透过率对树冠位置精确积分，包含树冠尺度的阴影重叠和叶片尺度的热点。
    '''
    hcb = hcr if hcb is None else max(float(hcb), 0.0)
    vza, sza, vsa = np.broadcast_arrays(*(np.asarray(value, dtype=float) for value in (vza, sza, vsa)))
    results = np.zeros((4,) + vza.shape)
    for index in np.ndindex(vza.shape):
        relative = abs((float(vsa[index]) + 180.0) % 360.0 - 180.0)
        results[(slice(None),) + index] = _crown_fractions(float(lai), float(std), float(hspot), float(hcr), float(rcr), hcb,
                                                           float(vza[index]), float(sza[index]), relative, float(G))
    return results[0], results[1], results[2], results[3]

def proportion_bidirectional_hom_one(lai,hspot,vza,sza,vsa):
    Psoil = gap_probability_hom_analytical(lai, vza)
    Pleaf = 1 - Psoil
    Psoil_sunlit = hotspot_layer(lai, hspot, vza, sza, vsa)
    Pleaf_sunlit = hotspot_analytical(lai, hspot, vza, sza, vsa)
    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded = Pleaf - Pleaf_sunlit
    return Psoil_sunlit,Psoil_shaded,Pleaf_sunlit,Pleaf_shaded

def _row_inside(x, period, width):
    """Cumulative length inside the rows, C(x) = integral of 1[(t mod period) < width] over [0, x]."""
    return np.floor(x / period) * width + np.minimum(np.mod(x, period), width)


def _row_path_inside(x, slope, rise, period, width):
    """Height extent inside the rows of a ray that rises `rise` while moving `slope` per unit height from cross-row x."""
    moving = np.abs(slope) > 1.0e-12
    safe = np.where(moving, np.abs(slope), 1.0)
    along = np.abs(_row_inside(x + slope * rise, period, width) - _row_inside(x, period, width)) / safe
    return np.where(moving, along, rise * (np.mod(x, period) < width))


def proportion_bidirectional_row_one(lai,hspot,row_width,row_blank,row_height,vza,vaa,sza,saa,raa,G = 0.5):
    '''
    垄行四组分比例：在垂直垄向的周期横截面内，对观测和太阳视线穿过垄带的光学厚度做精确积分。
    垄内按浑浊介质处理（消光系数 G*LAI*(w+b)/(w*H)）；叶片尺度的热点重叠按 hotspot_layer 的重叠权重叠加在两条路径的共享部分上。
    '''
    vza, vaa, sza, saa, raa = (np.atleast_1d(np.asarray(value, dtype=float)) for value in (vza, vaa, sza, saa, raa))
    period = row_width + row_blank
    if 1.0 * row_blank / period < 0.03:
        vsa = np.abs(saa - vaa)
        vsa = np.where(vsa > 180, 360 - vsa, vsa)
        return proportion_bidirectional_hom_one(lai, hspot, vza, sza, vsa)

    rd = np.pi / 180.0
    extinction = G * lai * period / row_width / row_height
    x0 = (np.arange(96) + 0.5) / 96 * period
    z = (np.arange(48) + 0.5) / 48 * row_height
    dz = row_height / 48
    X0, Z = np.meshgrid(x0, z, indexing="ij")

    results = np.zeros((4, vza.size))
    for index in range(vza.size):
        cv, cs = np.cos(vza[index] * rd), np.cos(sza[index] * rd)
        tan_v, tan_s = np.tan(vza[index] * rd), np.tan(sza[index] * rd)
        slope_v = tan_v * np.sin((vaa[index] - raa[index]) * rd)
        slope_s = tan_s * np.sin((saa[index] - raa[index]) * rd)

        # leaf-scale hotspot weight, same as hotspot_layer (relative azimuth view-sun)
        delta = tan_s ** 2 + tan_v ** 2 - 2 * tan_s * tan_v * np.cos((saa[index] - vaa[index]) * rd)
        correction = 2.0 / (G / cv + G / cs)
        alpha = np.sqrt(max(delta, 0.0)) * correction / hspot
        w = float(np.where(alpha > 1.0e-12, -np.expm1(-alpha) / np.where(alpha > 1.0e-12, alpha, 1.0), 1.0))

        tau_v0 = extinction / cv * _row_path_inside(x0, slope_v, row_height, period, row_width)
        tau_s0 = extinction / cs * _row_path_inside(x0, slope_s, row_height, period, row_width)
        soil_visible = float(np.mean(np.exp(-tau_v0)))
        soil_sunlit = float(np.mean(np.exp(-tau_v0 - tau_s0 + np.sqrt(tau_v0 * tau_s0) * w)))

        position = X0 + slope_v * Z
        inside = np.mod(position, period) < row_width
        rise = row_height - Z
        tau_v = extinction / cv * _row_path_inside(position, slope_v, rise, period, row_width)
        tau_s = extinction / cs * _row_path_inside(position, slope_s, rise, period, row_width)
        leaf_sunlit = float(np.mean(np.sum(extinction / cv * inside * np.exp(-tau_v - tau_s + np.sqrt(tau_v * tau_s) * w) * dz, axis=1)))

        leaf_visible = 1.0 - soil_visible
        leaf_sunlit = min(max(leaf_sunlit, 0.0), leaf_visible)
        results[:, index] = (soil_sunlit, soil_visible - soil_sunlit, leaf_sunlit, leaf_visible - leaf_sunlit)
    return results[0], results[1], results[2], results[3]


####################################################################
##### Voxel-based solution
####################################################################
#### directional


####################################################################
##### Voxel-based solution
####################################################################
#### directional

def proportion_directional_crown_voxel(lai, std, hcr, rcr, xza0, xaa0=0, G=0.5, ifleaf=1):
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
    # scheme = quadpy.sn.stroud_1967_7_c(3)
    # voxels = scheme.points
    # vo = scheme.weights *vol_crown

    voxels, number_voxel = ellipsoid_grid_crown_FRT()
    vo = voxels[3, :] / np.sum(voxels[3, :]) * vol_crown

    number_voxel = np.size(vo)
    number_angle = np.size(xza0)
    number = number_angle * number_voxel

    xza = np.repeat(xza0, number_voxel)
    xaa = np.repeat(xaa0, number_voxel)

    x = np.tile(voxels[0, :] * rcr, number_angle)
    y = np.tile(voxels[1, :] * rcr, number_angle)
    z = np.tile(voxels[2, :] * hcr, number_angle)
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
    gapx_inside_up, plv_inside_up, gapv_inside_down, plv_inside_down = gap_probability_crown_inside_voxel(x, y, z,
                                                                                                          density, hcr,
                                                                                                          rcr, xza, xaa)
    xv_up = x + plv_inside_up * sthetx * cphiv
    yv_up = y + plv_inside_up * sthetx * sphiv
    zv_up = z + plv_inside_up * cthetx + hcr
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv = \
        gap_probability_crown_outside_voxel(std, xv_up, yv_up, zv_up, density, hc, hcr, rcr, xza, xaa)
    gapv = gapx_inside_up * gapv_outside_up

    base_voxel = std * vo * density * G / cthetx
    Pleaf_voxel = gapv * base_voxel
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
        Pleafnew = 1 - Psoil
        Kcorrect = Pleafnew / Pleaf
        Pleaf_voxel = Pleaf_voxel * Kcorrect

    return Psoil, Pleaf_voxel


def proportion_directional_hom_voxel(lai0, vza0, CI=1.0, G=0.5):
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
    lai = np.tile(lai, number_angle)

    Gleaf_voxel = gap_probability_hom_voxel(lai, vza)
    Pleaf_voxel = Gleaf_voxel * dlai * G * CI / cthetv
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Psoil = gap_probability_hom_analytical(lai0, vza0)

    return Psoil, Pleaf_voxel


def proportion_directional_row_voxel(lai, row_width, row_blank, row_height, vza0, vaa0, raa0,
                                     nvw=5, nvh=10, nvb=5, G=0.5):
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
    Gvleaf_voxel, plv_leaf = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_vegetation, height_voxel_vegetation, vza_vegetation,
                                                       vra_vegetation, 0, 0)
    cthetv = np.cos(np.deg2rad(vza_vegetation))
    ### 体素的LAI = (density * rh / nvh * rs / nvs)
    Pvleaf_voxel = Gvleaf_voxel * density * G / cthetv * (rs / nvs) * (rh / nvh)
    Pvsoil_voxel, plv_soil = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_background, height_voxel_background, vza_background,
                                                       vra_background, 0, 0)
    ### 土壤的集成，因为每个的面积是相同的，所以就直接取数值平均，理应是按照距离进行加权平均
    Pvsoil_voxel = Pvsoil_voxel * (1.0 / nvs)

    Pvleaf_voxel = np.transpose(np.reshape(Pvleaf_voxel, [number_angle, number_voxel_vegetation]))
    Pvsoil_voxel = np.transpose(np.reshape(Pvsoil_voxel, [number_angle, number_voxel_background]))

    return Pvsoil_voxel, Pvleaf_voxel


def proportion_directional_crown_voxel_one(lai, std, hcr, rcr, xza0, xaa0=0, G=0.5, ifleaf=1):
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
        xaa0 = np.resize(xaa0, np.size(xza0))
    Psoil, Pleaf_voxel = proportion_directional_crown_voxel(lai, std, hcr, rcr, xza0, xaa0, G, ifleaf)
    ### 将体素进行集成
    Pleaf = np.sum(Pleaf_voxel, axis=0)

    return Psoil, Pleaf


def proportion_directional_hom_voxel_one(lai0, vza0, CI=1.0, G=0.5):
    '''
    冠层场景的组分观测比例，将体素的进行组装集成
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param vza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    Psoil, Pleaf_voxel = proportion_directional_hom_voxel(lai0, vza0, CI, G)
    Pleaf = np.sum(Pleaf_voxel, axis=0)

    return Psoil, Pleaf


def proportion_directional_row_voxel_one(lai, row_width, row_blank, row_height, vza0, vaa0, raa0,
                                         nvw=5, nvh=10, nvb=5, G=0.5):
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
    Psoil_voxel, Pleaf_voxel = proportion_directional_row_voxel(lai, row_width, row_blank, row_height, vza0, vaa0, raa0,
                                                                nvw, nvh, nvb, G)
    Psoil = np.sum(Psoil_voxel, axis=0)
    Pleaf = np.sum(Pleaf_voxel, axis=0)

    return Psoil, Pleaf


def proportion_directional_crown_voxel_layer(lai, std, hcr, rcr, xza0, number_layer=5, xaa0=0, G=0.5, ifleaf=1):
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
    Psoil, Pleaf_voxel = proportion_directional_crown_voxel(lai, std, hcr, rcr, xza0, xaa0, G, ifleaf)
    ### 计算设置层的观测比例，按照垂直高度进行的分层
    number_angle = np.size(xza0)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    vz = scheme.points[2, :]
    number_voxel = np.size(vz)
    Pleaf_layer = np.zeros([number_layer, number_angle])
    dlh = 1 * 2 / number_layer
    index_layer = (vz + 1) // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer, :] = np.sum(Pleaf_voxel[ind, :], axis=0)
    return Psoil, Pleaf_layer


def proportion_directional_hom_voxel_layer(lai0, xza0, number_layer=5, CI=1.0, G=0.5):
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
    vz = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    Pleaf_layer = np.zeros([number_layer, number_angle])
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer, :] = np.sum(Pleaf_voxel[ind, :], axis=0)

    return Psoil, Pleaf_layer


def proportion_directional_row_voxel_layer(lai, row_width, row_blank, row_height, xza0, xaa0, raa0, number_layer=5,
                                           nvw=5, nvh=10, nvb=5, G=0.5):
    '''
    冠层场景的组分观测比例
    :param lai:  叶面积指数
    :param hspot:  热点因子
    :param xza0:   观测天顶角
    :param sza0:   太阳天顶角
    :param G:      投影系数
    :return:  组分比例
    '''

    Psoil_voxel, Pleaf_voxel = proportion_directional_row_voxel(lai, row_width, row_blank, row_height, xza0, xaa0, raa0,
                                                                nvw, nvh, nvb, G)
    Psoil = np.sum(Psoil_voxel, axis=0)

    ### 垄行做的的分层，这里是按照垂直高度的分层
    number_voxel = nvh
    number_angle = np.size(xza0)
    Pleaf_layer = np.zeros([number_layer, number_angle])
    vz = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_layer[klayer, :] = np.sum(Pleaf_voxel[ind, :], axis=0)

    return Psoil, Pleaf_layer


###### bidirectional

def proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0=0, G=0.5, ifleaf=1):
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
        saa0 = np.resize(saa0, np.size(sza0))

    hcc = hc - hcr
    hcb = hc - hcr * 2
    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    density = lai / (std * vol_crown)
    # scheme = quadpy.sn.stroud_1967_7_c(3)
    # voxels = scheme.points

    voxels, number_voxel = ellipsoid_grid_crown_FRT()
    vo = voxels[3, :] / np.sum(voxels[3, :]) * vol_crown

    ### scheme.weights 的累计是1，乘以体积，则是体素的真实的体积
    # vo = scheme.weights *vol_crown
    number_voxel = np.size(vo)
    number_angle = np.size(vza0)
    number = number_angle * number_voxel

    vsa0 = np.abs(saa0 - vaa0)
    vza = np.repeat(vza0, number_voxel)
    vaa = np.repeat(vaa0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    saa = np.repeat(saa0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    ### 冠层的体素坐标和体积
    x = np.tile(voxels[0, :] * rcr, number_angle)
    y = np.tile(voxels[1, :] * rcr, number_angle)
    z = np.tile(voxels[2, :] * hcr, number_angle)
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
    gapv_inside_up, plv_inside_up, gapv_inside_down, plv_inside_down = gap_probability_crown_inside_voxel(x, y, z,
                                                                                                          density, hcr,
                                                                                                          rcr, vza, vaa)
    xv_up = x + plv_inside_up * sthetv * cphiv
    yv_up = y + plv_inside_up * sthetv * sphiv
    zv_up = z + plv_inside_up * cthetv + hcr
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv = \
        gap_probability_crown_outside_voxel(std, xv_up, yv_up, zv_up, density, hc, hcr, rcr, vza, vaa)
    gapv = gapv_inside_up * gapv_outside_up

    ### 太阳方向的透过率 = 太阳方向树冠内透过率 * 太阳方向树冠间/树冠外的透过率
    gaps_inside_up, pls_inside_up, gaps_inside_down, pls_inside_down = gap_probability_crown_inside_voxel(x, y, z,
                                                                                                          density, hcr,
                                                                                                          rcr, sza, saa)
    xs_up = x + pls_inside_up * sthets * cphis
    ys_up = y + pls_inside_up * sthets * sphis
    zs_up = z + pls_inside_up * cthets + hcr
    gaps_outside_up, pls_outside_up, upAreas, hcrs_up, inters = \
        gap_probability_crown_outside_voxel(std, xs_up, ys_up, zs_up, density, hc, hcr, rcr, sza, saa)
    gaps = gaps_inside_up * gaps_outside_up

    ### 树冠内热点的集成
    gapvs_inside_up, pl_hotspot = hotspot_path_length(plv_inside_up, pls_inside_up, density, hspot, vza, sza, vsa)
    gapvs_outside_up = hotspot_projection_area(std, density, hspot,
                                               xv_up, yv_up, zv_up, plv_outside_up, upAreav, hcrv_up, interv,
                                               xs_up, ys_up, zs_up, pls_outside_up, upAreas, hcrs_up, inters,
                                               vza, sza, vsa)
    gapvs = gapvs_inside_up * gapvs_outside_up

    ### 体积密度和遮挡
    base_voxel = std * v * density * G / cthetv
    Pleaf_voxel = gapv * base_voxel
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
    gapvs_soil = hotspot_projection_area(std, density, hspot, x_soil, y_soil, z_soil, plv_outside_soil, upAreav_soil,
                                         hcrv_soil, interv_soil,
                                         x_soil, y_soil, z_soil, pls_outside_soil, upAreas_soil, hcrs_soil, inters_soil,
                                         vza0, sza0, vsa0)

    ### 土壤的可视比例
    Psoil = gapv_soil
    Psoil_sunlit = gapvs_soil

    if ifleaf == 0:
        ### 比例偏差通过土壤平差
        Psoil = 1 - Pleaf
        Psoil_sunlit = gapvs_soil / gapv_soil * Psoil
    else:
        ### 比例偏差通过叶片平差
        Pleafnew = 1 - Psoil
        Kcorrect = Pleafnew / Pleaf
        Pleaf_voxel = Pleaf_voxel * Kcorrect
        Pleaf_sunlit_voxel = Pleaf_sunlit_voxel * Kcorrect

    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded_voxel = Pleaf_voxel - Pleaf_sunlit_voxel

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel


def proportion_bidirectional_hom_voxel(lai0, hspot, vza0, sza0, vsa0, CI=1.0, G=0.5):
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
    lai = np.tile(lai, number_angle)

    Gleaf_voxel = gap_probability_hom_voxel(lai, vza)
    Pleaf_voxel = Gleaf_voxel * dlai * G * CI / cthetv
    Pleaf_voxel = np.transpose(np.reshape(Pleaf_voxel, [number_angle, number_voxel]))
    Pleaf = np.sum(Pleaf_voxel, axis=0)

    Gleaf_sunlit_voxel = hotspot_hom_voxel(lai, hspot, vza, sza, vsa)
    Pleaf_sunlit_voxel = Gleaf_sunlit_voxel * dlai * G * CI / cthetv
    Pleaf_sunlit_voxel = np.transpose(np.reshape(Pleaf_sunlit_voxel, [number_angle, number_voxel]))
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)

    Psoil = gap_probability_hom_analytical(lai0, vza0)
    Psoil_sunlit = hotspot_layer(lai0, hspot, vza0, sza0, vsa0)
    Psoil_shaded = Psoil - Psoil_sunlit
    Pleaf_shaded_voxel = Pleaf_voxel - Pleaf_sunlit_voxel

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel


def proportion_bidirectional_row_voxel(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                       nvw=5, nvh=10, nvb=5, G=0.5):
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
    Gvleaf_voxel, plv_leaf = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_vegetation, height_voxel_vegetation, vza_vegetation,
                                                       vra_vegetation, vsa_vegetation)
    Gsleaf_voxel, pls_leaf = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_vegetation, height_voxel_vegetation, sza_vegetation,
                                                       sra_vegetation, vsa_vegetation, 1)
    ### 观测和太阳方向的热点
    Gvleaf_sunlit_voxel, plvs_leaf = hotspot_path_length(plv_leaf, pls_leaf, density, hspot, vza_vegetation,
                                                         sza_vegetation, vsa_vegetation)
    cthetv = np.cos(np.deg2rad(vza_vegetation))
    #########################################################################################
    ##### ?????????????????????????????????????????????????????????????????????????????
    ########################################################################################
    ### 从透过率到组分比例
    # Pvleaf_voxel = Gvleaf_voxel * density * G /cthetv * (rw / nvw) *(rh/nvh)
    # Pvleaf_sunlit_voxel = Gvleaf_sunlit_voxel * density * G / cthetv * (rw/nvw) *(rh/nvh)
    Pvleaf_voxel = Gvleaf_voxel * density * G / cthetv * (1.0 / nvs) * (rh / nvh)
    Pvleaf_sunlit_voxel = Gvleaf_sunlit_voxel * density * G / cthetv * (1.0 / nvs) * (rh / nvh)

    Pvsoil_voxel, plv_soil = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_background, height_voxel_background, vza_background,
                                                       vra_background, vsa_background)
    Pssoil_voxel, pls_soil = gap_probability_row_voxel(density, row_width, row_blank, row_height,
                                                       length_voxel_background, height_voxel_background, sza_background,
                                                       sra_background, vsa_background, 1)
    Pvsoil_sunlit_voxel, plvs_soil = hotspot_path_length(plv_soil, pls_soil, density, hspot, vza_background,
                                                         sza_background, vsa_background)

    Pvleaf_sunlit_voxel = np.transpose(np.reshape(Pvleaf_sunlit_voxel, [number_angle, number_voxel_vegetation]))
    Pvleaf_sunlit = np.sum(Pvleaf_sunlit_voxel, axis=0)
    Pvleaf_voxel = np.transpose(np.reshape(Pvleaf_voxel, [number_angle, number_voxel_vegetation]))
    Pvleaf = np.sum(Pvleaf_voxel, axis=0)
    Pvsoil_voxel = np.transpose(np.reshape(Pvsoil_voxel, [number_angle, number_voxel_background]))
    Pvsoil = np.sum(Pvsoil_voxel, axis=0) / number_voxel_background
    Pvsoil_sunlit_voxel = np.transpose(np.reshape(Pvsoil_sunlit_voxel, [number_angle, number_voxel_background]))
    Pvsoil_sunlit = np.sum(Pvsoil_sunlit_voxel, axis=0) / number_voxel_background
    Pvsoil_shaded = Pvsoil - Pvsoil_sunlit
    Pvleaf_shaded_voxel = Pvleaf_voxel - Pvleaf_sunlit_voxel

    return Pvsoil_sunlit, Pvsoil_shaded, Pvleaf_sunlit_voxel, Pvleaf_shaded_voxel


def proportion_bidirectional_crown_voxel_one(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0=0, G=0.5, ifleaf=1):
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
    Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0, G, ifleaf)

    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis=0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)
    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit, Pleaf_shaded


def proportion_bidirectional_hom_voxel_one(lai0, hspot, vza0, sza0, vsa0, CI=1.0, G=0.5):
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
    Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_hom_voxel(lai0, hspot, vza0, sza0, vsa0, CI, G)
    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis=0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)
    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit, Pleaf_shaded


def proportion_bidirectional_row_voxel_one(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                           nvw=5, nvh=10, nvb=5, G=0.5):
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
    Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_row_voxel(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                           nvw, nvh, nvb, G)
    Pleaf_shaded = np.sum(Pleaf_shaded_voxel, axis=0)
    Pleaf_sunlit = np.sum(Pleaf_sunlit_voxel, axis=0)

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit, Pleaf_shaded


def proportion_bidirectional_crown_voxel_layer(
        lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, number_layer=5, saa0=0, G=0.5, ifleaf=1):
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

    Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_voxel, Pleaf_shaded_voxel = \
        proportion_bidirectional_crown_voxel(lai, std, hspot, hc, hcr, rcr, vza0, sza0, vaa0, saa0, G, ifleaf)

    number_angle = np.size(vza0)
    scheme = quadpy.sn.stroud_1967_7_c(3)
    vz = scheme.points[2, :]
    number_voxel = np.size(vz)
    Pleaf_shaded_layer = np.zeros([number_layer, number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer, number_angle])
    dlh = 1 * 2 / number_layer
    index_layer = (vz + 1) // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer, :] = np.sum(Pleaf_shaded_voxel[ind, :], axis=0)
            Pleaf_sunlit_layer[klayer, :] = np.sum(Pleaf_sunlit_voxel[ind, :], axis=0)

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer

def proportion_bidirectional_hom_voxel_layer(lai0, hspot, vza0, sza0, vsa0, number_layer=5, CI=1.0, G=0.5):
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

    [number_voxel, number_angle] = np.shape(Pleaf_sunlit_voxel)
    # number_angle = np.size(vza0)
    vz = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    Pleaf_shaded_layer = np.zeros([number_layer, number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer, number_angle])
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer, :] = np.sum(Pleaf_shaded_voxel[ind, :], axis=0)
            Pleaf_sunlit_layer[klayer, :] = np.sum(Pleaf_sunlit_voxel[ind, :], axis=0)

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer
def proportion_bidirectional_row_voxel_layer(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                             number_layer=5, nvw=5, nvh=10, nvb=5, G=0.5):
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
        proportion_bidirectional_row_voxel(lai, hspot, row_width, row_blank, row_height, vza0, vaa0, sza0, saa0, raa0,
                                           nvw, nvh, nvb, G)

    number_voxel = nvh
    number_angle = np.size(vza0)

    Psoil_shaded = np.sum(Psoil_shaded_voxel, axis=0)
    Psoil_sunlit = np.sum(Psoil_sunlit_voxel, axis=0)
    Pleaf_shaded_layer = np.zeros([number_layer, number_angle])
    Pleaf_sunlit_layer = np.zeros([number_layer, number_angle])
    ### 分层划分，取整
    vz = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    dlh = np.int(number_voxel / number_layer)
    index_layer = vz // dlh
    for klayer in range(number_layer):
        ind = klayer == index_layer
        if np.sum(ind) > 0:
            Pleaf_shaded_layer[klayer, :] = np.sum(Pleaf_shaded_voxel[ind, :], axis=0)
            Pleaf_sunlit_layer[klayer, :] = np.sum(Pleaf_sunlit_voxel[ind, :], axis=0)

    return Psoil_sunlit, Psoil_shaded, Pleaf_sunlit_layer, Pleaf_shaded_layer
