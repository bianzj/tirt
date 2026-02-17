import numpy as np
import scipy.integrate as sci
from path import *
from utils import *

'''
gap probability 透过率
inside: 表示穿过了相干冠层
outside: 表示穿过了不相干冠层
'''


def leaf_inclination_distribution_function(ala):
    '''
    叶倾角分布函数，通过度数ala计算2.5,7.5.。。87.5度的分布概率，总概率的和为1
    :param ala: 角度，度数
    :return: 中间度数，出现概率
    '''

    ### 计算角度步长的中间代表步长

    n = 18
    tx2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]
    tx1 = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
    tx2 = np.asarray(tx2)
    tx1 = np.asarray(tx1)
    x = (tx2 + tx1) / 2.0
    tl1 = tx1 * (np.pi / 180.0)
    tl2 = tx2 * (np.pi / 180.0)
    excent = np.exp(-1.6184e-5 * np.power(ala, 3.) + 2.1145e-3 * np.power(ala, 2.0) - 1.2390e-1 * ala + 3.2491)

    freq = np.zeros(n)
    for i in range(n):
        x1 = excent / (np.sqrt(1.0 + excent * excent * np.tan(tl1[i]) * np.tan(tl1[i])))
        x2 = excent / (np.sqrt(1.0 + excent * excent * np.tan(tl2[i]) * np.tan(tl2[i])))
        if (excent == 1):
            freq[i] = abs(np.cos(tl1[i]) - np.cos(tl2[i]))
        else:
            alpha = excent / np.sqrt(abs(1.0 - excent * excent))
            alpha2 = alpha * alpha
            x12 = x1 * x1
            x22 = x2 * x2
            if (excent > 1):
                alpx1 = np.sqrt(alpha2 + x12)
                alpx2 = np.sqrt(alpha2 + x22)
                dum = x1 * alpx1 + alpha2 * np.log(x1 + alpx1)
                freq[i] = abs(dum - (x2 * alpx2 + alpha2 * np.log(x2 + alpx2)))
            else:
                almx1 = np.sqrt(alpha2 - x12)
                almx2 = np.sqrt(alpha2 - x22)
                dum = x1 * almx1 + alpha2 * np.arcsin(x1 / alpha)
                freq[i] = abs(dum - (x2 * almx2 + alpha2 * np.arcsin(x2 / alpha)))
    sum0 = np.sum(freq)
    freq0 = freq / sum0
    return x, freq0


#############################################################################
##### Analytical Solution
#############################################################################

#### 均质结构

def gap_probability_hom_analytical(lai, xza, CI=1.0, G=0.5):
    '''
    计算均质场景的透过率
    :param lai: 叶面积指数
    :param xza: 观测/太阳天顶角
    :param CI: 聚集指数
    :param G: 投影比例
    :return: 透过率
    '''
    cthetx = np.cos(np.deg2rad(xza))
    return np.exp(-lai * G * CI / cthetx)


def gap_probability_hom_hemisphere_analytical(lai, CI=1.0):
    '''
    均质场景的半球平均透过率，一种参数化估算方案
    :param lai: 叶面积指数
    :param CI:  聚集指数
    :return: 半球平均透过率
    '''
    ### 等效结果，通过拟合得到, 见 FR97 et al. FR02 et al.
    ### coeff 相当于投影比例 G
    coeff = 0.825
    return np.exp(-coeff * lai * CI)


#### 树冠结构

def gap_probability_crown_analytical(lai, std, hcr, rcr, xza, G=0.5):
    '''
    树冠的透过率
    :param lai: 叶面积指数
    :param std: 森林树密度
    :param hcr: 垂直方向的半径
    :param rcr: 水平方向的半径
    :param xza: 天顶角
    :param G: 投影系数
    :return: 透过率
    '''
    #### FRT
    tgthx = np.tan(np.deg2rad(xza))
    cthetx = np.cos(np.deg2rad(xza))
    ### 计算椭球树冠在地面的投影面积，即需要计算透过率的部分
    areav = np.sqrt(rcr * rcr + np.power(hcr * tgthx, 2)) * np.pi * rcr
    ### lai /(areav*std) 植被部分的等效LAI
    bv_in = np.exp(-lai / (areav * std) * G / cthetx)
    ### 面积*面积的透过率，这里面积就是考虑了几何形状，然后以该几何形状作为个体；而不是通过每个小面元作为个体。
    ### 原理一样，但是处理的基本单元不同，每个叶片还是叶片的集合*叶片的等效不透过率
    bv = np.exp(-areav * (1 - bv_in) * std)

    #### Li Ju Cai et al. 不太推荐使用，有偏差可能
    # area0 = np.pi * np.power(rcr, 2)
    # slai = lai / area0 / std
    # tv = np.tan(np.deg2rad(vza))
    # areav = np.sqrt(rcr * rcr + np.power(hcr * tv, 2)) * np.pi * rcr
    ### 椭球体到球体的转换，之后可以直接用
    # theta = np.arctan(hcr / rcr * np.tan(np.deg2rad(vza)))
    # uv = np.cos(theta)
    # bv_betweeen = np.exp(-std * area0/uv)
    # bv_in = np.exp(-slai/uv * G)
    # bv = bv_betweeen + (1-bv_betweeen) * bv_in

    return bv


#### 垄行结构

def gap_probability_row_analytical(lai, row_width, row_blank, row_height, xza, xra, CI=1.0, G=0.5):
    '''
    垄行植被类型的透过率，参数化解决方案,是多个区域的平均
    :param lai: 叶面积指数
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return: 透过率
    '''
    row_spacing = row_width + row_blank
    laie = lai * row_spacing / row_width
    density = laie / row_height

    ### homgeneous vegetation
    ### 阈值判断，当小于0.03时候，认为垄行情况退化为均质
    if 1.0 * row_blank / row_spacing < 0.03:
        uv = np.cos(np.deg2rad(xza))
        bv = np.exp(-G * lai * CI / uv)
        return bv

    ### row-planted crop for angle = 0
    ### 观测方位角方向和垄行方向相同的时候，就是情况0
    ### 观测天顶角方向是0时，也是情况0
    number_angle = np.size(xza)
    bv = np.zeros(number_angle)
    ind = (xra % 180 == 0)
    if np.sum(ind):
        bv[ind] = gap_probability_row_zero(row_width, row_blank, row_height, density, xza[ind])
    ind = (xza == 0)
    if np.sum(ind):
        bv[ind] = gap_probability_row_zero(row_width, row_blank, row_height, density, xza[ind])

    ### 其他情况分别处理 ind0是排除以上情况类型；ind是ind0和各个情况
    ind0 = (xra % 180 != 0) * (xza != 0)
    if np.sum(ind0) > 0:
        row_width_ = np.zeros(number_angle)
        row_blank_ = np.zeros(number_angle)
        row_height_ = np.zeros(number_angle)
        number = np.zeros(number_angle)

        sr = np.sin(np.deg2rad(xra[ind0] % 180))
        row_width_[ind0] = row_width / sr
        row_blank_[ind0] = row_blank / sr
        row_spacing_ = row_width_ + row_blank_
        row_height_[ind0] = row_height * np.tan(np.deg2rad(xza[ind0]))
        number[ind0] = np.floor(row_height_[ind0] / row_spacing_[ind0])

        ### 详见 论文
        ### part 1
        ind = ind0 * ((number * row_spacing_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_width_)) * (row_width_ <= row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part14(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                 xza[ind])
        ### part 2
        ind = ind0 * ((number * row_spacing_ + row_width_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_blank_)) * (row_width_ <= row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part2(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                xza[ind])
        ### part 3
        ind = ind0 * ((number * row_spacing_ + row_blank_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_spacing_)) * (row_width_ <= row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part36(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                 xza[ind])
        ### part 4
        ind = ind0 * ((number * row_spacing_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_blank_)) * (row_width_ > row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part14(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                 xza[ind])
        ### part 5
        ind = ind0 * ((number * row_spacing_ + row_blank_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_width_)) * (row_width_ > row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part5(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                xza[ind])
        ### part 6
        ind = ind0 * ((number * row_spacing_ + row_width_) <= row_height_) * (
                    row_height_ <= (number * row_spacing_ + row_spacing_)) * (row_width_ > row_blank_)
        if np.sum(ind) > 0:
            bv[ind] = gap_probability_row_part36(row_width_[ind], row_blank_[ind], row_height, number[ind], density,
                                                 xza[ind])

    return bv


def gap_probability_row_zero(row_width, row_blank, row_height, density, xza, CI=1.0, G=0.5):
    '''
    垂直方向、天顶角0、平行垄行方向
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    row_spacing = row_width + row_blank
    cthetx = np.cos(np.deg2rad(xza))
    bv = (row_width * np.exp(-G * row_height * density * CI / cthetx) + row_blank) / row_spacing
    return bv


def gap_probability_row_part14(row_width_, row_blank_, row_height, number, density, xza, G=0.5):
    '''
    条带1/4
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    row_spacing_ = row_width_ + row_blank_
    uv = np.cos(np.deg2rad(xza))
    tv = np.tan(np.deg2rad(xza))
    sv = np.sin(np.deg2rad(xza))
    ctg = 1.0 / tv
    p1 = (number * row_spacing_ + row_width_ - row_height * tv - 2 * sv / G / density) * \
         np.exp(-G * density * (row_height - number * row_blank_ * ctg) / uv)
    p2 = (number * row_spacing_ + row_blank_ - row_height * tv + 2 * sv / G / density) * \
         np.exp(-number * G * density * row_width_ * ctg / uv)
    return (p1 + p2) / row_spacing_


def gap_probability_row_part2(row_width_, row_blank_, row_height, number, density, xza, G=0.5):
    '''
    条带2
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    row_spacing_ = row_width_ + row_blank_
    uv = np.cos(np.deg2rad(xza))
    tv = np.tan(np.deg2rad(xza))
    sv = np.sin(np.deg2rad(xza))
    ctg = 1.0 / tv
    p1 = (-number * row_spacing_ - row_width_ + row_height * tv - 2 * sv / G / density) * \
         np.exp(-(number + 1) * G * density * row_width_ * ctg / uv)
    p2 = (number * row_spacing_ + row_blank_ - row_height * tv + 2 * sv / G / density) * \
         np.exp(-number * G * density * row_width_ * ctg / uv)
    return (p1 + p2) / row_spacing_


def gap_probability_row_part36(row_width_, row_blank_, row_height, number, density, xza, G=0.5):
    '''
    条带3/6
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    row_spacing_ = row_width_ + row_blank_
    uv = np.cos(np.deg2rad(xza))
    tv = np.tan(np.deg2rad(xza))
    sv = np.sin(np.deg2rad(xza))
    ctg = 1.0 / tv
    p1 = (-number * row_spacing_ - row_width_ + row_height * tv - 2 * sv / G / density) * \
         np.exp(-(number + 1) * G * density * row_width_ * ctg / uv)
    p2 = (-number * row_spacing_ - row_blank_ + row_height * tv + 2 * sv / G / density) * \
         np.exp(-G * density * (row_height - (number + 1) * row_blank_ * ctg) / uv)
    return (p1 + p2) / row_spacing_


def gap_probability_row_part5(row_width_, row_blank_, row_height, number, density, xza, G=0.5):
    '''
    条带5
    :param row_width: 行植被宽度
    :param row_blank: 行空白宽度
    :param row_height: 行高度
    :param xza: 太阳/观测天顶角
    :param xra: 太阳/观测-垄行相对方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return:
    '''
    row_spacing_ = row_width_ + row_blank_
    uv = np.cos(np.deg2rad(xza))
    tv = np.tan(np.deg2rad(xza))
    sv = np.sin(np.deg2rad(xza))
    ctg = 1.0 / tv
    p1 = (-number * row_spacing_ - row_blank_ + row_height * tv + 2 * sv / G / density) * \
         np.exp(- G * density * (row_height - (number + 1) * row_blank_ * ctg) / uv)
    p2 = (number * row_spacing_ + row_width_ - row_height * tv - 2 * sv / G / density) * \
         np.exp(- G * density * (row_height - number * row_blank_ * ctg) / uv)
    return (p1 + p2) / row_spacing_


### 树干结构

def gap_probability_trunk_analytical(hc_trunk, dbh_trunk, std_trunk, xza):
    '''
    参数化的方法计算树干的透过率，假设树干是圆柱
    :param hc_trunk:  树干高度
    :param dbh_trunk:  树干宽度
    :param std_trunk:  树干密度
    :param xza:  观测/太阳天顶角
    :return:  透过率
    '''

    ### 树干在观测上的投影面积比例
    tgthx = np.tan(np.deg2rad(xza))
    ### 研究区树干的投影面积
    Strunk = hc_trunk * dbh_trunk * tgthx * std_trunk

    ### 与树冠的投影相同，只是树干的透过率是0
    return np.exp(-Strunk)


##################################################################################
###### Voxel-based Solution
##################################################################################

#### 均质场景

def gap_probability_hom_voxel(lai, xza, CI=1.0, G=0.5):
    '''
    均质场景的体元的透过率
    :param lai: 叶面积指数
    :param number_voxel: 体元层数
    :param xza:  天顶角
    :param CI:  聚集指数
    :param G:  投影系数
    :return:  所有体元的透过率
    '''

    ### 整体透过率和矩阵体素的透过率相同
    cthetx = np.cos(np.deg2rad(xza))
    gap = np.exp(-lai * CI * G / cthetx)

    return gap


#### 树冠场景

def gap_probability_crown_inside_voxel(xi, yj, zk, density, hcr, rcr, xza, xaa, G=0.5):
    '''
    计算树冠内各个体元传过树冠内部分，到达树冠表面的透过率，向上和向下分别的结果
    :param xi: 坐标
    :param yj: 坐标
    :param zk: 坐标
    :param density: 树冠体密度
    :param hcr: 垂直方向半径
    :param rcr:  水平方向半径
    :param xza:  天顶角
    :param xaa:  方位角
    :param G:  投影系数
    :return:  向上透过率，向上的路径长度，向下透过率，向下路径长度(负的)
    '''
    pl_inside_up, pl_inside_down = path_length_crown_inside_voxel(xi, yj, zk, xza, xaa, hcr, rcr)
    ### pl_inside_up 正的； pl*density = lai/cthetx
    gap_inside_up = np.exp(-pl_inside_up * density * G)
    ### pl_inside_down 负值; pl*density = lai/cthetx
    gap_inside_down = np.exp(pl_inside_down * density * G)
    return gap_inside_up, pl_inside_up, gap_inside_down, pl_inside_down


def gap_probability_crown_outside_voxel(
        std, xi_up, yj_up, zk_up, density, hc, hcr, rcr, xza, xaa, CI=1.0, G=0.5):
    '''
    从各个体元从树冠表面，穿过其他冠层，到传感器/冠层顶的透过率
    :param std: 森林树密度
    :param xi_up: 坐标
    :param yj_up: 坐标
    :param zk_up: 坐标
    :param density: 树冠内体密度
    :param hc: 树高
    :param hcr: 垂直方向半径
    :param rcr: 水平方向半径
    :param xza: 天顶角
    :param xaa: 方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return: 向上方向的透过率，向上方向的路径长度，上层植被投影面积，上层植被被假设成树冠的投影的等效半径，上层树冠的等效拦截概率
    '''
    ### 高度以上冠层的等效路径长度和投影面积
    pl_outside_up, upArea = path_length_crown_outside_voxel(xi_up, yj_up, zk_up, hc, hcr, rcr, xza, xaa)
    gap_outside_up = np.exp(-pl_outside_up * density * G)
    ### 没有直接用透过率，而是假设以上植被为新的树冠，计算了树冠的平均投影和平均的比例
    ### 计算等效阻挡概率（不透过率概率），CICI是考虑聚集指数情况
    ### CI = 1，表明完全的随机状态，小于1是聚集
    CICI = 1 - CI
    if np.abs(CICI) > 0.001:
        inter = -np.log(1.0 - (1.0 - gap_outside_up) * CICI) / CICI
    else:
        inter = 1 - gap_outside_up
    gap_outside_up = np.exp(-std * inter * upArea)

    ### hcr_up 上层冠层等效成树冠，然后在地面投影的等效边长（而不是树冠本身的等效边长，因为乘tgthx）
    temp = zk_up * 1.0
    hcb = hc - hcr * 2
    indd = (temp < hcb)
    tgthx = np.tan(np.deg2rad(xza))
    if np.sum(indd) > 0: temp[indd] = hcb
    hcr_up = ((hc + temp) * 0.5 - zk_up) * tgthx

    return gap_outside_up, pl_outside_up, upArea, hcr_up, inter


def gap_probability_crown_outside_voxel_corrected(
        std, xi_up, yj_up, zk_up, density, hc, hcr, rcr, za, aa, CI=1.0,G=0.5):
    '''
    从树冠表面穿过冠层到传感器/冠层顶的透过率,修正方法，树冠和树冠不能重叠
    :param std: 森林树密度
    :param xi_up: 坐标
    :param yj_up: 坐标
    :param zk_up: 坐标
    :param density: 树冠内体密度
    :param hc: 树高
    :param hcr: 垂直方向半径
    :param rcr: 水平方向半径
    :param za: 天顶角
    :param aa: 方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return: 透过率
    '''
    pl_outside_up, upArea = path_length_crown_outside_voxel_corrected(xi_up, yj_up, zk_up, hc, hcr, rcr, za, aa)
    gap_outside_up = np.exp(-pl_outside_up * density * G)
    CICI = 1 - CI
    if np.abs(CICI) > 0.001:
        inter = -np.log(1.0 - (1.0 - gap_outside_up) * CICI) / CICI
    else:
        inter = 1 - gap_outside_up
    gap_outside_up = np.exp(-std * inter * upArea)

    temp = zk_up * 1.0
    hcb = hc - hcr * 2
    indd = (temp < hcb)
    tgthx = np.tan(np.deg2rad(za))
    if np.sum(indd) > 0: temp[indd] = hcb
    hcr_up = ((hc + temp) * 0.5 - zk_up) * tgthx

    return gap_outside_up, pl_outside_up, upArea, hcr_up, inter


#### 垄行场景

def gap_probability_row_voxel(
        density, row_width, row_blank, row_height, length, height, xza, raa, sva,ifverse=0,  G=0.5):
    '''
    垄行场景背景组分体元的透过率
    :param lai: 植被指数
    :param row_width: 垄宽度
    :param row_blank: 垄空白宽度
    :param row_height:  垄高度
    :param xza:  天顶角
    :param xaa:  方位角
    :param sva:  相对方位角
    :param nvw:  体元数目
    :param nvh:  体元数目
    :param nvb:  体元数目
    :param CI: 聚集指数
    :param G: 投影系数
    :return: 透过率
    '''

    ### 计算各个体素离开植被冠层的路径长度，ifverse是正向/逆向的问题
    pl_voxel = path_length_row_voxel(length, height, row_width, row_blank, row_height, xza, raa, sva, ifverse)
    Gvoxel = np.exp(-pl_voxel * density * G)
    return Gvoxel, pl_voxel


###################################################################################
##### Under for 4SAIL, MGP and thermalFRT
###################################################################################

def gap_probability_under_inside_voxel(xi, yj, zk, density, hcr, rcr, xza, xaa, G = 0.5):
    '''
    计算树冠及下方体素到达树冠表面的透过率,穿过自己所属的树冠的透过率及其路径长度，如果是没有穿过，那透过率就是1，路径长度是0
    :param xi: 坐标
    :param yj: 坐标
    :param zk: 坐标
    :param density: 树冠体密度
    :param hcr: 垂直方向半径
    :param rcr:  水平方向半径
    :param xza:  天顶角
    :param xaa:  方位角
    :param G:  投影系数
    :return:  透过率
    '''
    ### 树冠往上和往下方的路径长度
    pl_inside_up, pl_inside_down = path_length_crown_inside_voxel(xi, yj, zk, xza, xaa, hcr, rcr)
    ### 确认所计算高度要低于树冠底层的高度，才是np.abs(pl_inside_up - pl_inside_down)
    pl = np.abs(pl_inside_up - pl_inside_down)
    gap_inside_up = np.exp(-pl * density * G)
    return gap_inside_up, pl


def gap_probability_under_outside_voxel(std, xi_up, yj_up, zk_up, density, hc, hcr, rcr, dbh_trunk, hc_trunk, xza, xaa,
                                        CI=1.0, G=0.5):
    '''
    从某一层到达树冠顶部的透过率，这里除了树冠的路径，还包括了树干的阻挡，但是这里是没有穿自身冠层的
    under, 表示了位置属性，而并不是特指树干类型
    :param std: 森林树密度
    :param xi_up: 坐标
    :param yj_up: 坐标
    :param zk_up: 坐标
    :param density: 树冠内体密度
    :param hc: 树高
    :param hcr: 垂直方向半径
    :param rcr: 水平方向半径
    :param xza: 天顶角
    :param xaa: 方位角
    :param CI: 聚集指数
    :param G: 投影系数
    :return: 透过率
    '''

    ### 树干投影比例
    tgthx = np.tan(np.deg2rad(xza))
    ### 树冠透过率
    pl_outside_up, upArea = path_length_crown_outside_voxel(xi_up, yj_up, zk_up, hc, hcr, rcr, xza, xaa)
    gap_outside_up = np.exp(-pl_outside_up * density * G)

    ### 没有直接用透过率，而是假设以上植被为新的树冠，计算了树冠的平均投影和平均的比例
    CICI = 1 - CI
    if np.abs(CICI) > 0.001:
        inter = -np.log(1.0 - (1.0 - gap_outside_up) * CICI) / CICI
    else:
        inter = 1 - gap_outside_up

    ### 树干的投影面积，树干圆柱的校正，而非长方形
    hcc = hc - hcr
    stemArea = projection_area_stem(zk_up, hcc, dbh_trunk, hc_trunk)
    stemProj = (dbh_trunk / 2.0) * (dbh_trunk / 2.0) * np.pi / 2.0 * np.ones(np.size(xza))
    ind = xza > 0
    if np.sum(ind) > 0:
        stemProj[ind] = stemProj[ind] + stemArea[ind] * tgthx[ind]

    ### 透过率情况，包括了树冠的等效不透过率+树干的透过率
    gap_outside_up = np.exp(-std * (inter * upArea + stemProj))

    ### 投影的等效半径
    temp = zk_up * 1.0
    hcb = hc - hcr * 2
    indd = (temp < hcb)
    if np.sum(indd) > 0: temp[indd] = hcb
    hcr_up = ((hc + temp) * 0.5 - zk_up) * tgthx

    return gap_outside_up, pl_outside_up, upArea, hcr_up, inter, stemProj

