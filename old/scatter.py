import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utils import *
from gap import *
from hotspot import *
import quadpy
from proportion import *

''' 
多次散射项
首先实现变量扩展的方法；
然后才开始矩阵方法转化；
'''

def multiple_scattering_analytical(lai,vza,refl_soil,refl_leaf):
    '''
    参数化的方案计算多次散射
    :param lai: 叶面积指数
    :param vza: 观测天顶角
    :param refl_soil:  土壤反射率
    :param refl_leaf:  叶片反射率
    :return:  冠层的多次散射项
    '''

    bv = gap_probability_hom_analytical(lai, vza)
    M = gap_probability_hom_hemisphere_analytical(lai)
    alpha = np.asarray([0.2855375, 0.2885375, 0.2964427, 0.3003953, 0.3083004, 0.3201581, 0.3399209, 0.3715415, 0.4189723, 1])
    vza_index_bottom = np.asarray(vza / 10, dtype=np.int)
    vza_index_top = vza_index_bottom + 1

    ratio = (vza % 10) / 10.0
    alphanew = alpha[vza_index_bottom] * (1-ratio) + alpha[vza_index_top] * ratio
    multiple_scattering_emissivity = bv * (1-M)* refl_soil  + (1-alphanew)*(1-bv*M) *(1-bv)*refl_leaf
    multiple_scattering_emissivity = multiple_scattering_emissivity *(1-refl_leaf)

    return multiple_scattering_emissivity


def multiple_scattering_analytical_sunlit(lai,vza,sza,refl_soil,refl_leaf):
    '''
	参数化的方案计算多次散射,区分光照和阴影
    :param lai: 叶面积指数
    :param vza: 观测天顶角
    :param sza: 太阳天顶角
    :param Kleaf: 比例系数
    :param refl_soil: 土壤反射率
    :param refl_leaf: 叶片反射率
    :return:
    '''
    bv = gap_probability_hom_analytical(lai, vza)
    M = gap_probability_hom_hemisphere_analytical(lai)
    alpha = np.asarray([0.2885375, 0.2885375, 0.2964427, 0.3003953, 0.3083004, 0.3201581, 0.3399209, 0.3715415, 0.4189723, 1])
    vza_index = np.asarray(vza / 10, dtype=np.int)
    Vsunlit = hotspot_vegetation_volume(lai, sza)
    Vshaded = 1 - Vsunlit
    multiple_scattering_emissivity_sunlit = bv * (1-M)* refl_soil*Vsunlit + (1-alpha[vza_index])*(1-bv*M) *(1-bv)*refl_leaf*Vsunlit
    multiple_scattering_emissivity_sunlit = multiple_scattering_emissivity_sunlit *(1-refl_leaf)
    multiple_scattering_emissivity_shaded = bv * (1-M)* refl_soil*Vshaded + (1-alpha[vza_index])*(1-bv*M) *(1-bv)*refl_leaf*Vshaded
    multiple_scattering_emissivity_shaded = multiple_scattering_emissivity_shaded *(1-refl_leaf)

    return multiple_scattering_emissivity_sunlit,multiple_scattering_emissivity_shaded


def multiple_scattering_voxel(lai,number_vegetation,vza,refl_soil,refl_leaf,number_background = 1,CI = 1.0,G = 0.5):

    '''
	均质场景体元的多次散射项，生成了大矩阵，表明了任意两个体元的相互影响，透过率，反射率，透过率，这是均质场景的多次散射计算法方法
	这种方法跟其他的比，很特别
	E * M * R * B
	:param lai:  叶面积指数
	:param number_vegetation:
	:param vza:
	:param refl_soil:
	:param refl_leaf:
	:param number_background:
	:param CI:
	:param G:
	:return:
	layer 是层的标识
	matrix 是矩阵的标识
	gap 是透过率
	lai 是层的物质密度
	计算的是每个体素的接收项
	E_matrix 体素的发射项
	R_matrix 体素的反射项
	M_matrix 体素到待求算体素的路程
	'''

    ### 植被层 number_voxel + 土壤层 1
    number_voxel = number_vegetation + number_background
    number_angle = np.size(vza)

    dlai = lai / number_vegetation


    lai_layer = np.hstack([np.linspace(0.5,number_vegetation-0.5,number_vegetation),number_vegetation])

    lai_matrix_up = np.tile(np.transpose(np.asmatrix(lai_layer)), (1,number_voxel))
    lai_matrix_down = number_voxel - lai_matrix_up
    uv_angle = np.cos(np.deg2rad(vza))
    uv_matrix_angle = np.tile(np.transpose(uv_angle),(number_voxel,1))
    lai_matrix_angle = np.tile(np.transpose(np.asmatrix(lai_layer)),(1,number_angle))
    B_matrix = np.exp(-G*CI*lai_matrix_angle*dlai/uv_matrix_angle)


    lai_matrix_scatter = np.tile(lai_layer,(number_voxel,1))
    lai_matrix_scatter = np.abs(lai_matrix_scatter-lai_matrix_up)
    m_matrix_scatter = gap_probability_hom_hemisphere_analytical(lai_matrix_scatter * dlai)
    m_layer = gap_probability_hom_hemisphere_analytical(dlai)

    emis_soil = 1- refl_soil
    emis_leaf = 1- refl_leaf

    E_layer = np.hstack([np.repeat(emis_leaf,number_vegetation) * dlai * 0.825, emis_soil*1.0])
    E_matrix = np.tile(np.transpose(np.asmatrix(E_layer)),(1,number_voxel))
    R_layer = np.hstack([np.repeat(refl_leaf,number_vegetation) * dlai * G , refl_soil*1.0])
    R_matrix = np.tile((np.asmatrix(R_layer)),(number_voxel,1))

    R_matrix = np.triu(R_matrix,1) + np.diag(np.diag(R_matrix))
    R_matrix[-1] = 0
    M_matrix = (m_matrix_scatter)


    result1 = np.multiply(E_matrix,M_matrix)
    result2 = np.multiply(result1,R_matrix)
    result3 = np.matmul(result2,B_matrix)
    emis_angle = np.sum(result3, axis=0)

    return np.squeeze(np.array(emis_angle))

def multiple_scattering_voxel_sunlit(lai,number_vegetation,vza,sza,refl_soil,refl_leaf,number_background = 1,CI = 1.0,G = 0.5):
    '''
    均质场景体元的多次散射项，生成了大矩阵，表明了任意两个体元的相互影响，发射项i * 半球透过率i-》j * 反射率j * 透过率j -》 sensor
	E * M * R * B
    :param lai: 叶面积指数
    :param number_vegetation: 植被层数
    :param vza:  观测天顶角
    :param sza:  太阳天顶角
    :param refl_soil:  土壤反射率
    :param refl_leaf:  叶片反射率
    :param number_background:  土壤层数
    :param CI:  聚集指数
    :param G:   投影系数
    :return:   多次散射项，区分光照和阴影
    '''
    ### 植被层 number_voxel + 土壤层 1
    number_voxel = number_vegetation + number_background
    number_angle = np.size(vza)

    dlai = lai / number_vegetation
    lai_layer = np.hstack([np.linspace(0.5,number_vegetation-0.5,number_vegetation),number_vegetation])
    lai_matrix_up = np.tile(np.transpose(np.asmatrix(lai_layer)), (1,number_voxel))
    lai_matrix_down = number_voxel - lai_matrix_up
    uv_angle = np.cos(np.deg2rad(vza))
    uv_matrix_angle = np.tile(np.transpose(uv_angle),(number_voxel,1))
    lai_matrix_angle = np.tile(np.transpose(np.asmatrix(lai_layer)),(1,number_angle))
    B_matrix = np.exp(-G*CI*lai_matrix_angle*dlai/uv_matrix_angle) / uv_matrix_angle


    lai_matrix_scatter = np.tile(lai_layer,(number_voxel,1))
    lai_matrix_scatter = np.abs(lai_matrix_scatter-lai_matrix_up)
    m_matrix_scatter = gap_probability_hom_hemisphere_analytical(lai_matrix_scatter * dlai)
    m_layer = gap_probability_hom_hemisphere_analytical(dlai)

    emis_soil = 1- refl_soil
    emis_leaf = 1- refl_leaf

    E_layer = np.hstack([np.repeat(emis_leaf,number_vegetation) * dlai * 0.825, emis_soil*1.0])
    E_matrix = np.tile(np.transpose(np.asmatrix(E_layer)),(1,number_voxel))
    R_layer = np.hstack([np.repeat(refl_leaf,number_vegetation) * dlai * G , refl_soil*1.0])
    R_matrix = np.tile((np.asmatrix(R_layer)),(number_voxel,1))
    # R_matrix[-1,:] = 0

    R_matrix = np.triu(R_matrix,1) + np.diag(np.diag(R_matrix))
    R_matrix[-1] = 0
    M_matrix = (m_matrix_scatter)

    Vsunlit = hotspot_vegetation_volume(lai_matrix_up, sza)
    Vshaded = 1- Vsunlit

    E_matrix_sunlit = np.multiply(E_matrix,Vsunlit)
    E_matrix_shaded = np.multiply(E_matrix,Vshaded)

    result1_sunlit = np.multiply(E_matrix_sunlit,M_matrix)
    result1_shaded = np.multiply(E_matrix_shaded,M_matrix)
    result2_sunlit = np.multiply(result1_sunlit,R_matrix)
    result2_shaded = np.multiply(result1_shaded,R_matrix)

    result3_sunlit = np.matmul(result2_sunlit,B_matrix)
    result3_shaded = np.matmul(result2_shaded,B_matrix)

    emis_angle_sunlit = np.sum(result3_sunlit, axis=0)
    emis_angle_shaded = np.sum(result3_shaded, axis=0)

    return np.squeeze(np.array(emis_angle_sunlit)),np.squeeze(np.array(emis_angle_shaded))

def multiple_scattering_hom_spectral_invariance_endmember(lai_crown, lai_unveg, hc_trunk, dbh_trunk, std_trunk,
                                                          ec, et, ev, es,
                                                          vza0, sza0, vaa0, saa0 = 0, G_crown = 0.5, G_unveg = 0.5):


    eps = 0.0001
    number_voxel = 10
    number_angle = np.size(vza0)
    number_hemisphere = 18
    if2 = 2.0
    au = 0.5
    ad = 0.5

    hza0 = np.linspace(0.5, number_hemisphere - 0.5, number_hemisphere) * np.pi / 2.0 / number_hemisphere
    hz0 = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    fweight0 = np.sin(hza0)*np.cos(hza0) * np.pi / 2.0 / number_hemisphere
    tgthh0 = np.tan(hza0)
    ctheth0 = np.cos(hza0)

    if lai_unveg < eps: lai_unveg = eps
    if lai_crown < eps: lai_crown = eps

    #### expand
    vsa0 = np.abs(vaa0 - saa0)
    vza = np.repeat(vza0,number_voxel)
    sza = np.repeat(sza0,number_voxel)
    vsa = np.repeat(vsa0,number_voxel)
    hz = np.tile(hz0, number_angle)
    cthetv = np.cos(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))
    cthetv0 = np.cos(np.deg2rad(vza0))
    tgthv0 = np.tan(np.deg2rad(vza0))
    tgths0 = np.tan(np.deg2rad(sza0))
    cthets0 = np.cos(np.deg2rad(sza0))

    '''
    光谱不变理论：
    从观测方向出发发射光子，出现2种情况：2）光子碰到植被；1）光子透过植被；
    1）光子透过植被，表明这部分是植被下层的贡献；
    2）光子碰到植被，表明这部分是植被这层的贡献，除了碰撞的植被的发射项，光子继续反射，出现3种情况包括：
        a) 从上层逃逸 eu；
        b) 从下层逃逸 ed；
        c) 碰到植被的其他部分 p； 
    经过这层植被作为跳板进入到传感器的贡献： 上层复合贡献*eu + 下层复合贡献 * ed + 本层贡献 * p
    
    这里通过对其分层或者分体元的方式分别进行计算，然后进行累计，每层或体元计算如下：
    i. 选择层或者体元计算dlai 和 dS
    ii. 计算观测方向的透过率fv；
    iii. 计算上半球的平均透过率 mu；
    iv. 计算下半球的平均透过率 md；
    v. 计算往上方向的逃逸概率 fv * dS * au * mu   到达该层概率；在该层碰撞到概率；往上半球方向概率；从该层到上层逃逸概率
    vi. 计算往下方向的逃逸概率 fv * dS * ad * md 同上
    vii. 计算该层的吸收概率 p = 1 - eu - ed 
    '''
    ### 观测方向的拦截概率
    Strunk = dbh_trunk * hc_trunk * tgthv0 * std_trunk
    ind = vza0 < 0.001
    if np.sum(ind) > 0: Strunk[ind] =(dbh_trunk/2.0)*(dbh_trunk/2.0)*np.pi*std_trunk
    ftrunk = np.exp(-Strunk)  # gap frequency of trunk in a direction
    fcrown = np.exp(-lai_crown * G_crown / cthetv0)  # gap frequency of vegetation in a direction
    fveg = np.exp(-lai_unveg * G_unveg / cthetv0)  # gap frequency of vegetation in a direction
    Strunk_s = dbh_trunk * hc_trunk * tgths0 * std_trunk
    ind = sza0 < 0.001
    if np.sum(ind) > 0: Strunk_s[ind] = (dbh_trunk / 2.0) * (dbh_trunk / 2.0) * np.pi * std_trunk
    ftrunk_s = np.exp(-Strunk_s)  # gap frequency of trunk in a direction
    fcrown_s = np.exp(-lai_crown * G_crown / cthets0)  # gap frequency of vegetation in a direction
    fveg_s = np.exp(-lai_unveg * G_unveg / cthets0)  # gap frequency of vegetation in a direction
    i0t = 1 - ftrunk  # interception probability
    i0c = 1 - fcrown  # interception probability
    i0v = 1 - fveg

    ### 半球方向的平均拦截概率
    tempT = dbh_trunk * hc_trunk * tgthh0 * std_trunk  # projection of tree trunk in hemisphere space
    tempV = lai_unveg * G_unveg / ctheth0  # projection of tree crown in hemisphere space
    tempC = lai_crown * G_crown / ctheth0  # projection of tree crown in hemisphere space
    i0ttemp = 1 - np.exp(-tempT)  # interception probability of trunk
    i0tplus = if2 * np.sum(i0ttemp * fweight0)  # average interception probability
    i0ctemp = 1 - np.exp(-tempC)  # interception probability of vegetation
    i0cplus = if2 * np.sum(i0ctemp * fweight0)  # average interception probability
    i0vtemp = 1 - np.exp(-tempV)  # interception probability of vegetation
    i0vplus = if2 * np.sum(i0vtemp * fweight0)  # average interception probability


    ####################################################
    #### Crown
    dlai = lai_crown * G_crown  / number_voxel
    dS =  dlai / cthetv
    fv = np.exp(- hz * dS)

    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0,number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)
    cthethtemp = np.cos(hzatemp)
    futemp = np.exp(-hztemp * dlai / cthethtemp)
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai / cthethtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    eu = np.transpose(np.reshape(mu * dS * au * fv,[number_angle,number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv,[number_angle,number_voxel]))
    edc = np.sum(ed,axis=0)/i0c
    euc = np.sum(eu,axis=0)/i0c
    pc = 1- edc - euc
    ####################################################
    ####### Trunk
    dlai = dbh_trunk * hc_trunk * std_trunk / number_voxel
    dS = dlai * tgthv
    fv = np.exp(- hz * dS )
    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0,number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)

    tgthhtemp = np.tan(hzatemp)
    futemp = np.exp(-hztemp * dlai * tgthhtemp )
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai * tgthhtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    eu = np.transpose(np.reshape(mu * dS * au * fv,[number_angle,number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv,[number_angle,number_voxel]))
    edt = np.sum(ed,axis=0)/i0t
    eut = np.sum(eu,axis=0)/i0t
    pt = 1- edt - eut

    ####################################################
    #### Veg
    dlai = lai_unveg * G_unveg  / number_voxel
    dS =  dlai / cthetv
    fv = np.exp(- hz * dS)

    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0,number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)
    cthethtemp = np.cos(hzatemp)
    futemp = np.exp(-hztemp * dlai / cthethtemp)
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai / cthethtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    eu = np.transpose(np.reshape(mu * dS * au * fv,[number_angle,number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv,[number_angle,number_voxel]))
    edv = np.sum(ed,axis=0)/i0v
    euv = np.sum(eu,axis=0)/i0v
    pv = 1- edv - euv

    '''
    组分间的多次项，这里只考虑了单次散射
    组分的发射项： 组分发射率 * 组分的半球拦截概率
    组分的反射项： 组分反射率 * 来自上层反射率贡献 （i0x * eux ）
                  组分反射率 * 来自下层反射率贡献  (i0x * edx )
                  组分反射率 * 来自中层反射率贡献  (i0x * px)
     ！！！ 多次散射可以通过雅克比方法或者高斯赛德尔方法求解
    '''

    ####################################################################
    ### 光照和阴影树干
    ###################################################################
    ### 树干到树干部分 发射率et，反射率(1-et) 吸收项i0t * pt，上层透过率  1-i0c
    emit0 = et * (1-et) * (i0t * pt) * (1-i0c)
    ### 树干到下层植被
    emit1 = et * i0tplus * (i0v * euv) * (1-ev) * (1-i0t) * (1-i0c)
    ### 树干到土壤部分, 发射率 et, 向下方向的贡献 i0t * ed, 植被的半球透过率， 土壤反射，下层植被的方向透过率，树干的方向透过率，树冠的方向透过率
    emit2 = et * i0tplus * (1-i0vplus)* (1-es) * (1-i0v) *(1-i0t) *(1-i0c)
    ### 树干到树冠部分  发射率， 往上方向贡献概率， 树冠的拦截概率，树冠的反射率
    ### 可以被忽略
    emit3 = et * i0tplus * (i0c * edc) * (1 - ec)
    emit = emit0 + emit1 + emit2 + emit3
    vsraa = (180.0 - vsa0) / 180.0  # sunlit part becasue of raa
    vhraa = 1 - vsraa
    fts = (1 - np.exp(-dbh_trunk * hc_trunk * std_trunk * tgths0)) / (dbh_trunk * hc_trunk * std_trunk * tgths0)
    emits = emit * vhraa * fts * fcrown_s  # shaded back->trunk->sensor
    emith = emit - emits  # sunlit back-> trunk -> sensor


    ######################################################################
    #### 光照和阴影树冠
    ######################################################################
    ### 树冠到树冠 树冠发射率， 树冠本层的贡献率，树冠的反射率
    emic0 = ec * (i0c * pc) *(1-ec)
    ### 树冠到树干 树冠发射率，树冠下层的贡献率，树干的拦截概率，树干的反射率，树冠的方向透过率
    emic1 = ec * i0cplus * (i0t * eut) * (1-et) * (1-i0c)
    ### 树冠到下层植被,树冠发射率， 树冠下层的贡献率，树干的透过率，下层植被的半球拦截概率，下层植被的反射率，树干方向透过率，树冠方向透过率
    emic2 = ec * i0cplus * (1-i0tplus) * (i0v * euv) * (1-ev) * (1-i0t) * (1-i0c)
    ### 树冠到土壤 树冠发射率，树冠贡献比例，树干的半球透过率，下层植被的半球透过率，土壤的反射率，下层植被的方向透过率，树冠的方向透过率，树干的方向透过率
    emic3 = ec * i0cplus * (1-i0tplus) * (1-i0vplus) * (1-es) * (1-i0v) * (1-i0t) * (1-i0c)
    emic = emic0 + emic1 + emic2 + emic3
    fss = (1-fcrown_s)/(lai_crown*G_crown/cthets0)
    emics = emic * fss
    emich = emic - emics


    ######################################################################
    #### 光照和阴影下层植被
    ######################################################################
    ### 下层植被到下层植被,下层植被发射率，这层内贡献，下层植被反射率，树干的方向透过率，冠层的方向透过率
    emiv0 = ev * (i0v * pv) *(1-ev) * (1-i0t) * (1-i0c)
    ### 下层植被到土壤， 下层植被发射率，下层植被贡献，土壤透过率，树干的方向透过率，下层植被的方向透过率，树冠的方向透过率
    emiv1 = ev * i0vplus * (1-es) * (1-i0t) * (1-i0c) *(1-i0v)
    ### 下层植被到树干,发射，贡献，反射率，方向透过率
    ### 可以被忽略
    emiv2 = ev * i0vplus * (i0t * edt) * (1-et) * (1-i0c)
    ### 下层植被到上层植被，发射率，发射项，树冠半球透过，树冠下层贡献，树冠反射率
    ### 可以被忽略
    emiv3 = ev * i0vplus * (1-i0tplus) * (i0c * edc) * (1-ec)

    emiv = emiv0 + emiv1 + emiv2 + emiv3
    fss = (1-fveg_s) / (lai_unveg*G_unveg/cthets0) *fcrown_s*ftrunk_s
    emivs = fss * emiv
    emivh = emiv - emivs

    #####################################################################
    ##### 光照和阴影土壤
    #####################################################################
    ### 土壤到下层植被，发射率，下层植被到下层，下层植被的反射率，树干的方向透过率，树冠的方向透过率
    ### 可以被忽略
    emis1 = es * (i0v * edv) * (1-ev) * (1-i0t) * (1-i0c)
    ### 土壤到树干,发射率，下层植被半球透过率，树冠下层贡献率，树冠的方向透过率
    ### 可以被忽略
    emis2 = es * (1-i0vplus) * (i0t * edt) * (1-et) * (1-i0c)
    ### 土壤到冠层,发射率，下层植被半球透过率，树干的半球透过率，树冠的下层贡献率，树冠的反射率
    ### 可以被忽略
    emis3 = es * (1-i0vplus) * (1-i0tplus) * (i0c * edc) * (1-ec)

    emis = emis1 + emis2 + emis3
    fss = fcrown_s * ftrunk_s * fveg_s
    emiss = fss * emis
    emish = emis - emiss

    memis = emiv + emit + emis + emic


    return emiss,emish,emics,emich,emivs,emivh,emits,emith,memis

def multiple_scattering_crown_spectral_invariance_endmember(lai_crown, std_crown, hc, hcr, rcr, lai_unveg, hc_trunk, dbh_trunk,
                                                            ec, et, ev, es,
                                                            vza0, sza0, vaa0, saa0 = 0, G_crown = 0.5, G_unveg = 0.5):
    '''
    光谱不变理论计算组分有效发射率，可以得到各个组分的结果，是解析方法
    '''

    eps = 0.0001
    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    scheme = quadpy.sn.stroud_1967_7_c(3)
    voxels = scheme.points
    vo = scheme.weights * vol_crown
    density = lai_crown / (std_crown * vol_crown)
    number_voxel = np.size(vo)
    number_angle = np.size(vza0)
    number_hemisphere = 18
    std_trunk = std_crown
    if2 = 2.0
    au = 0.5
    ad = 0.5

    hza0 = np.linspace(0.5, number_hemisphere - 0.5, number_hemisphere) * np.pi / 2.0 / number_hemisphere
    haa0 = np.repeat(0,number_hemisphere)*np.pi / 2.0 / number_hemisphere
    hz0 = np.linspace(0.5, number_voxel - 0.5, number_voxel)

    fweight0 = np.sin(hza0) * np.cos(hza0) * np.pi / 2.0 / number_hemisphere
    tgthh0 = np.tan(hza0)
    ctheth0 = np.cos(hza0)

    if lai_unveg < eps: lai_unveg = eps
    if lai_crown < eps: lai_crown = eps

    x0 = (voxels[0, :] * rcr)
    y0 = (voxels[1, :] * rcr)
    z0 = (voxels[2, :] * hcr)

    #### expand
    vsa0 = np.abs(vaa0 - saa0)
    vza = np.repeat(vza0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    x = np.tile(x0,number_angle)
    y = np.tile(y0,number_angle)
    z = np.tile(z0,number_angle)

    hz = np.tile(hz0, number_angle)
    cthetv = np.cos(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    sphivs = np.sin(np.deg2rad(vsa))
    cphivs = np.cos(np.deg2rad(vsa))

    cthetv0 = np.cos(np.deg2rad(vza0))
    tgthv0 = np.tan(np.deg2rad(vza0))
    tgths0 = np.tan(np.deg2rad(sza0))
    cthets0 = np.cos(np.deg2rad(sza0))

    '''
    光谱不变理论：
    从观测方向出发发射光子，出现2种情况：2）光子碰到植被；1）光子透过植被；
    1）光子透过植被，表明这部分是植被下层的贡献；
    2）光子碰到植被，表明这部分是植被这层的贡献，除了碰撞的植被的发射项，光子继续反射，出现3种情况包括：
        a) 从上层逃逸 eu；
        b) 从下层逃逸 ed；
        c) 碰到植被的其他部分 p； 
    经过这层植被作为跳板进入到传感器的贡献： 上层复合贡献*eu + 下层复合贡献 * ed + 本层贡献 * p

    这里通过对其分层或者分体元的方式分别进行计算，然后进行累计，每层或体元计算如下：
    i. 选择层或者体元计算dlai 和 dS
    ii. 计算观测方向的透过率fv；
    iii. 计算上半球的平均透过率 mu；
    iv. 计算下半球的平均透过率 md；
    v. 计算往上方向的逃逸概率 fv * dS * au * mu   到达该层概率；在该层碰撞到概率；往上半球方向概率；从该层到上层逃逸概率
    vi. 计算往下方向的逃逸概率 fv * dS * ad * md 同上
    vii. 计算该层的吸收概率 p = 1 - eu - ed 
    
    
    temp： 表示是半球的所有/逐个角度
    plus:  表示半球所有角度的积分
    
    number_angle number_voxel number_hemi
    先number_Angle*number_voxel,number_hemi
    再number_angle,number_voxel
    '''
    ### 观测方向的拦截概率
    Strunk = dbh_trunk * hc_trunk * tgthv0 * std_trunk
    ind = Strunk < 0.001
    if np.sum(ind) > 0: Strunk[ind] = (dbh_trunk / 2.0) * (dbh_trunk / 2.0) * np.pi * std_trunk
    ftrunk = np.exp(-Strunk)
    fcrown,temp = proportion_directional_crown_voxel_one(lai_crown,std_crown,hcr,rcr,vza0)
    # fcrown = 0.5
    fveg = np.exp(-lai_unveg * G_unveg / cthetv0)  # gap frequency of vegetation in a direction
    Strunk_s = dbh_trunk * hc_trunk * tgths0 * std_trunk
    ind = Strunk_s < 0.001
    if np.sum(ind) > 0: Strunk_s[ind] = (dbh_trunk / 2.0) * (dbh_trunk / 2.0) * np.pi * std_trunk
    ftrunk_s = np.exp(-Strunk_s)  # gap frequency of trunk in a direction
    fcrown_s = np.exp(-lai_crown * G_crown / cthets0)  # gap frequency of vegetation in a direction
    fveg_s = np.exp(-lai_unveg * G_unveg / cthets0)  # gap frequency of vegetation in a direction
    i0t = 1 - ftrunk  # interception probability
    i0c = 1 - fcrown  # interception probability
    i0v = 1 - fveg

    ### 半球方向的平均拦截概率
    tempT = dbh_trunk * hc_trunk * tgthh0 * std_trunk  # projection of tree trunk in hemisphere space
    tempV = lai_unveg * G_unveg / ctheth0  # projection of tree crown in hemisphere space
    fcrown0,temp = proportion_directional_crown_voxel_one(lai_crown,std_crown,hcr,rcr,hza0)
    i0ctemp = 1 - fcrown0
    i0ttemp = 1 - np.exp(-tempT)  # interception probability of trunk
    i0tplus = if2 * np.sum(i0ttemp * fweight0)  # average interception probability
    i0cplus = if2 * np.sum(i0ctemp * fweight0)  # average interception probability
    i0vtemp = 1 - np.exp(-tempV)  # interception probability of vegetation
    i0vplus = if2 * np.sum(i0vtemp * fweight0)  # average interception probability

    ####################################################
    #### Crown

    dlai0 = vo * density * G_crown
    dlai = np.tile(dlai0, number_angle)
    dS = dlai / cthetv

    ### 观测方向的透过率 = 观测方向树冠内透过率 * 观测方向树冠间/树冠外的透过率
    gapv_inside_up, plv_inside_up,gapv_inside_down,plv_inside_down = gap_probability_crown_inside_voxel(x, y, z, density, hcr, rcr, vza, vsa)
    xv_up = x + plv_inside_up * sthetv * cphivs
    yv_up = y + plv_inside_up * sthetv * sphivs
    zv_up = z + plv_inside_up * cthetv + hcr
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv = \
        gap_probability_crown_outside_voxel(std_crown, xv_up, yv_up, zv_up, density, hc, hcr, rcr, vza, vsa)
    fv = gapv_inside_up * gapv_outside_up

    ### 半球方向的问题
    hzatemp = np.tile(hza0, number_angle * number_voxel)
    haatemp = np.tile(haa0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    cthethtemp = np.cos(hzatemp)
    sthethtemp = np.sin(hzatemp)
    cphivstemp = np.cos(haatemp)
    sphivstemp = np.sin(haatemp)
    ztemp = np.repeat(z, number_hemisphere)
    ytemp = np.repeat(y, number_hemisphere)
    xtemp = np.repeat(x, number_hemisphere)

    gapvtemp_inside_up, plvtemp_inside_up,gapvtemp_inside_down, plvtemp_inside_down = \
        gap_probability_crown_inside_voxel(xtemp, ytemp, ztemp, density, hcr, rcr, hzatemp, haatemp)
    ### 往上
    xvtemp_up = xtemp + plvtemp_inside_up * sthethtemp * cphivstemp
    yvtemp_up = ytemp + plvtemp_inside_up * sthethtemp * sphivstemp
    zvtemp_up = ztemp + plvtemp_inside_up * cthethtemp + hcr
    gapvtemp_outside_up, plvtemp_outside_up, upAreavtemp, hcrvtemp_up, intervtemp_up = \
        gap_probability_crown_outside_voxel(std_crown, xvtemp_up, yvtemp_up, zvtemp_up, density, hc, hcr, rcr, hzatemp, haatemp)
    futemp = gapvtemp_inside_up * gapvtemp_outside_up
    ### 往下
    xvtemp_down = xtemp + plvtemp_inside_down * sthethtemp * cphivstemp
    yvtemp_down = ytemp + plvtemp_inside_down * sthethtemp * sphivstemp
    zvtemp_down = ztemp + plvtemp_inside_down * cthethtemp + hcr
    gapvtemp_outside_down, plvtemp_outside_down, downAreavtemp, hcrvtemp_down, intervtemp_down = \
        gap_probability_crown_outside_voxel(std_crown, xvtemp_down, yvtemp_down, zvtemp_down, density, hc, hcr, rcr, hzatemp,haatemp)
    fdtemp = gapvtemp_inside_down * gapvtemp_outside_down

    mu = if2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel * number_angle, number_hemisphere])),
                    axis=0)
    md = if2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel * number_angle, number_hemisphere])),
                    axis=0)
    eu = np.transpose(np.reshape(mu * dS * au * fv, [number_angle, number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv, [number_angle, number_voxel]))
    edc = np.sum(ed, axis=0) / i0c
    euc = np.sum(eu, axis=0) / i0c
    pc = 1 - edc - euc




    ####################################################
    ####### Trunk
    dlai = dbh_trunk * hc_trunk * std_trunk / number_voxel
    dS = dlai * tgthv
    fv = np.exp(- hz * dS)
    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)
    tgthhtemp = np.tan(hzatemp)
    futemp = np.exp(-hztemp * dlai * tgthhtemp)
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai * tgthhtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    eu = np.transpose(np.reshape(mu * dS * au * fv,[number_angle,number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv,[number_angle,number_voxel]))
    edt = np.sum(ed,axis=0)/i0t
    eut = np.sum(eu,axis=0)/i0t
    pt = 1 - edt - eut

    ####################################################
    #### Veg
    dlai = lai_unveg * G_unveg / number_voxel
    dS = dlai / cthetv
    fv = np.exp(- hz * dS)

    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)
    cthethtemp = np.cos(hzatemp)
    futemp = np.exp(-hztemp * dlai / cthethtemp)
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai / cthethtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel*number_angle,number_hemisphere])),axis= 0)
    eu = np.transpose(np.reshape(mu * dS * au * fv,[number_angle,number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv,[number_angle,number_voxel]))
    edv = np.sum(ed,axis=0)/i0v
    euv = np.sum(eu,axis=0)/i0v
    pv = 1 - edv - euv

    ####################################################################
    ### 光照和阴影树干
    ###################################################################
    ### 树干到树干部分 发射率et，反射率(1-et) 吸收项i0t * pt，上层透过率  1-i0c
    emit0 = et * (1 - et) * (i0t * pt) * (1 - i0c)
    ### 树干到下层植被
    emit1 = et * i0tplus * (i0v * euv) * (1 - ev) * (1 - i0t) * (1 - i0c)
    ### 树干到土壤部分, 发射率 et, 向下方向的贡献 i0t * ed, 植被的半球透过率， 土壤反射，下层植被的方向透过率，树干的方向透过率，树冠的方向透过率
    emit2 = et * i0tplus * (1 - i0vplus) * (1 - es) * (1 - i0v) * (1 - i0t) * (1 - i0c)
    ### 树干到树冠部分  发射率， 往上方向贡献概率， 树冠的拦截概率，树冠的反射率
    ### 可以被忽略
    emit3 = et * i0tplus * (i0c * edc) * (1 - ec)
    emit = emit0 + emit1 + emit2 + emit3
    vsraa = (180.0 - vsa0) / 180.0  # sunlit part becasue of raa
    vhraa = 1 - vsraa
    fts = (1 - np.exp(-dbh_trunk * hc_trunk * std_trunk * tgths0)) / (dbh_trunk * hc_trunk * std_trunk * tgths0)
    emits = emit * vhraa * fts * fcrown_s  # shaded back->trunk->sensor
    emith = emit - emits  # sunlit back-> trunk -> sensor

    ######################################################################
    #### 光照和阴影树冠
    ######################################################################
    ### 树冠到树冠 树冠发射率， 树冠本层的贡献率，树冠的反射率
    emic0 = ec * (i0c * pc) * (1 - ec)
    ### 树冠到树干 树冠发射率，树冠下层的贡献率，树干的拦截概率，树干的反射率，树冠的方向透过率
    emic1 = ec * i0cplus * (i0t * eut) * (1 - et) * (1 - i0c)
    ### 树冠到下层植被,树冠发射率， 树冠下层的贡献率，树干的透过率，下层植被的半球拦截概率，下层植被的反射率，树干方向透过率，树冠方向透过率
    emic2 = ec * i0cplus * (1 - i0tplus) * (i0v * euv) * (1 - ev) * (1 - i0t) * (1 - i0c)
    ### 树冠到土壤 树冠发射率，树冠贡献比例，树干的半球透过率，下层植被的半球透过率，土壤的反射率，下层植被的方向透过率，树冠的方向透过率，树干的方向透过率
    emic3 = ec * i0cplus * (1 - i0tplus) * (1 - i0vplus) * (1 - es) * (1 - i0v) * (1 - i0t) * (1 - i0c)
    emic = emic0 + emic1 + emic2 + emic3
    fss = (1 - fcrown_s) / (lai_crown * G_crown / cthets0)
    emics = emic * fss
    emich = emic - emics

    ######################################################################
    #### 光照和阴影下层植被
    ######################################################################
    ### 下层植被到下层植被,下层植被发射率，这层内贡献，下层植被反射率，树干的方向透过率，冠层的方向透过率
    emiv0 = ev * (i0v * pv) * (1 - ev) * (1 - i0t) * (1 - i0c)
    ### 下层植被到土壤， 下层植被发射率，下层植被贡献，土壤透过率，树干的方向透过率，下层植被的方向透过率，树冠的方向透过率
    emiv1 = ev * i0vplus * (1 - es) * (1 - i0t) * (1 - i0c) * (1 - i0v)
    ### 下层植被到树干,发射，贡献，反射率，方向透过率
    ### 可以被忽略
    emiv2 = ev * i0vplus * (i0t * edt) * (1 - et) * (1 - i0c)
    ### 下层植被到上层植被，发射率，发射项，树冠半球透过，树冠下层贡献，树冠反射率
    ### 可以被忽略
    emiv3 = ev * i0vplus * (1 - i0tplus) * (i0c * edc) * (1 - ec)

    emiv = emiv0 + emiv1 + emiv2 + emiv3
    fss = (1 - fveg_s) / (lai_unveg * G_unveg / cthets0) * fcrown_s * ftrunk_s
    emivs = fss * emiv
    emivh = emiv - emivs

    #####################################################################
    ##### 光照和阴影土壤
    #####################################################################
    ### 土壤到下层植被，发射率，下层植被到下层，下层植被的反射率，树干的方向透过率，树冠的方向透过率
    ### 可以被忽略
    emis1 = es * (i0v * edv) * (1 - ev) * (1 - i0t) * (1 - i0c)
    ### 土壤到树干,发射率，下层植被半球透过率，树冠下层贡献率，树冠的方向透过率
    ### 可以被忽略
    emis2 = es * (1 - i0vplus) * (i0t * edt) * (1 - et) * (1 - i0c)
    ### 土壤到冠层,发射率，下层植被半球透过率，树干的半球透过率，树冠的下层贡献率，树冠的反射率
    ### 可以被忽略
    emis3 = es * (1 - i0vplus) * (1 - i0tplus) * (i0c * edc) * (1 - ec)

    emis = emis1 + emis2 + emis3
    fss = fcrown_s * ftrunk_s * fveg_s
    emiss = fss * emis
    emish = emis - emiss

    memis = emiv + emit + emis + emic

    return emics, emich, emits, emith, emivs, emivh, emiss, emish, memis


def multiple_scattering_hom_spectral_invariance(lai_crown, ec,  es, vza0, sza0, vaa0, saa0=0, G_crown=0.5, ):
    eps = 0.0001
    number_voxel = 10
    number_angle = np.size(vza0)
    number_hemisphere = 18
    if2 = 2.0
    au = 0.5
    ad = 0.5

    hza0 = np.linspace(0.5, number_hemisphere - 0.5, number_hemisphere) * np.pi / 2.0 / number_hemisphere
    hz0 = np.linspace(0.5, number_voxel - 0.5, number_voxel)
    fweight0 = np.sin(hza0) * np.cos(hza0) * np.pi / 2.0 / number_hemisphere
    tgthh0 = np.tan(hza0)
    ctheth0 = np.cos(hza0)


    if lai_crown < eps: lai_crown = eps

    #### expand
    vsa0 = np.abs(vaa0 - saa0)
    vza = np.repeat(vza0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    hz = np.tile(hz0, number_angle)
    cthetv = np.cos(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))
    cthetv0 = np.cos(np.deg2rad(vza0))
    tgthv0 = np.tan(np.deg2rad(vza0))
    tgths0 = np.tan(np.deg2rad(sza0))
    cthets0 = np.cos(np.deg2rad(sza0))

    '''
    光谱不变理论：
    从观测方向出发发射光子，出现2种情况：2）光子碰到植被；1）光子透过植被；
    1）光子透过植被，表明这部分是植被下层的贡献；
    2）光子碰到植被，表明这部分是植被这层的贡献，除了碰撞的植被的发射项，光子继续反射，出现3种情况包括：
        a) 从上层逃逸 eu；
        b) 从下层逃逸 ed；
        c) 碰到植被的其他部分 p； 
    经过这层植被作为跳板进入到传感器的贡献： 上层复合贡献*eu + 下层复合贡献 * ed + 本层贡献 * p

    这里通过对其分层或者分体元的方式分别进行计算，然后进行累计，每层或体元计算如下：
    i. 选择层或者体元计算dlai 和 dS
    ii. 计算观测方向的透过率fv；
    iii. 计算上半球的平均透过率 mu；
    iv. 计算下半球的平均透过率 md；
    v. 计算往上方向的逃逸概率 fv * dS * au * mu   到达该层概率；在该层碰撞到概率；往上半球方向概率；从该层到上层逃逸概率
    vi. 计算往下方向的逃逸概率 fv * dS * ad * md 同上
    vii. 计算该层的吸收概率 p = 1 - eu - ed 
    '''
    ### 观测方向的拦截概率

    fcrown = np.exp(-lai_crown * G_crown / cthetv0)  # gap frequency of vegetation in a direction
    fcrown_s = np.exp(-lai_crown * G_crown / cthets0)  # gap frequency of vegetation in a direction
    i0c = 1 - fcrown  # interception probability


    ### 半球方向的平均拦截概率

    tempC = lai_crown * G_crown / ctheth0  # projection of tree crown in hemisphere space
    i0ctemp = 1 - np.exp(-tempC)  # interception probability of vegetation
    i0cplus = if2 * np.sum(i0ctemp * fweight0)  # average interception probability

    ####################################################
    #### Crown
    dlai = lai_crown * G_crown / number_voxel
    dS = dlai / cthetv
    fv = np.exp(- hz * dS)

    hzatemp = np.tile(hza0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    hztemp = np.repeat(hz, number_hemisphere)
    cthethtemp = np.cos(hzatemp)
    futemp = np.exp(-hztemp * dlai / cthethtemp)
    fdtemp = np.exp(-(number_voxel - hztemp) * dlai / cthethtemp)
    mu = 2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel * number_angle, number_hemisphere])),
                    axis=0)
    md = 2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel * number_angle, number_hemisphere])),
                    axis=0)
    eu = np.transpose(np.reshape(mu * dS * au * fv, [number_angle, number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv, [number_angle, number_voxel]))
    edc = np.sum(ed, axis=0) / i0c
    euc = np.sum(eu, axis=0) / i0c
    pc = 1 - edc - euc

    '''
    组分间的多次项，这里只考虑了单次散射
    组分的发射项： 组分发射率 * 组分的半球拦截概率
    组分的反射项： 组分反射率 * 来自上层反射率贡献 （i0x * eux ）
                  组分反射率 * 来自下层反射率贡献  (i0x * edx )
                  组分反射率 * 来自中层反射率贡献  (i0x * px)
     ！！！ 多次散射可以通过雅克比方法或者高斯赛德尔方法求解
    '''



    ######################################################################
    #### 光照和阴影树冠
    ######################################################################
    ### 树冠到树冠 树冠发射率， 树冠本层的贡献率，树冠的反射率
    emic0 = ec * (i0c * pc) * (1 - ec)
    ### 树冠到树干 树冠发射率，树冠下层的贡献率，树干的拦截概率，树干的反射率，树冠的方向透过率
    emic1 = ec * i0cplus * (1-es) * (1-i0c)

    emic = emic0 + emic1
    fss = (1 - fcrown_s) / (lai_crown * G_crown / cthets0)
    emics = emic * fss
    emich = emic - emics



    memis = emics + emich

    return emics, emich, memis


def multiple_scattering_crown_spectral_invariance(lai_crown, std_crown, hc, hcr, rcr,
                                                   ec,  es, vza0, sza0, vaa0, saa0=0, G_crown=0.5):
    '''
    光谱不变理论计算组分有效发射率，可以得到各个组分的结果，是解析方法
    '''

    eps = 0.0001
    vol_crown = 4.0 * np.pi / 3 * rcr * rcr * hcr
    scheme = quadpy.sn.stroud_1967_7_c(3)
    voxels = scheme.points
    vo = scheme.weights * vol_crown
    density = lai_crown / (std_crown * vol_crown)
    number_voxel = np.size(vo)
    number_angle = np.size(vza0)
    number_hemisphere = 18
    std_trunk = std_crown
    if2 = 2.0
    au = 0.5
    ad = 0.5

    hza0 = np.linspace(0.5, number_hemisphere - 0.5, number_hemisphere) / number_hemisphere * 90
    haa0 = np.repeat(0, number_hemisphere)  / number_hemisphere * 90
    hz0 = np.linspace(0.5, number_voxel - 0.5, number_voxel)

    stheth0 = np.sin(np.deg2rad(hza0))
    ctheth0 = np.cos(np.deg2rad(hza0))
    fweight0 = stheth0 * ctheth0 * np.pi / 2.0 / number_hemisphere



    if lai_crown < eps: lai_crown = eps

    x0 = (voxels[0, :] * rcr)
    y0 = (voxels[1, :] * rcr)
    z0 = (voxels[2, :] * hcr)

    #### expand
    vsa0 = np.abs(vaa0 - saa0)
    vza = np.repeat(vza0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    x = np.tile(x0, number_angle)
    y = np.tile(y0, number_angle)
    z = np.tile(z0, number_angle)

    cthetv = np.cos(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    sphivs = np.sin(np.deg2rad(vsa))
    cphivs = np.cos(np.deg2rad(vsa))

    cthetv0 = np.cos(np.deg2rad(vza0))
    tgthv0 = np.tan(np.deg2rad(vza0))
    tgths0 = np.tan(np.deg2rad(sza0))
    cthets0 = np.cos(np.deg2rad(sza0))

    '''
    1.光谱不变理论
    从观测方向出发发射光子，出现2种情况：2）光子碰到植被；1）光子透过植被；
    1）光子透过植被，表明这部分是植被下层的贡献；
    2）光子碰到植被，表明这部分是植被这层的贡献，除了碰撞的植被的发射项，光子继续反射，出现3种情况包括：
        a) 从上层逃逸 eu；
        b) 从下层逃逸 ed；
        c) 碰到植被的其他部分 p； 
    经过这层植被作为跳板进入到传感器的贡献： 上层复合贡献*eu + 下层复合贡献 * ed + 本层贡献 * p
    2. 体素离散和积分：
    i. 选择层或者体元计算dlai 和 dS
    ii. 计算观测方向的透过率fv；
    iii. 计算上半球的平均透过率 mu；
    iv. 计算下半球的平均透过率 md；
    v. 计算往上方向的逃逸概率 fv * dS * au * mu   到达该层概率；在该层碰撞到概率；往上半球方向概率；从该层到上层逃逸概率
    vi. 计算往下方向的逃逸概率 fv * dS * ad * md 同上
    vii. 计算该层的吸收概率 p = 1 - eu - ed 
    3. 半球离散
    temp： 表示是半球的所有/逐个角度
    plus:  表示半球所有角度的积分
    4. 数组折叠
    长数组拆分方法（number_angle*number_voxel*number_hemi）
    先number_Angle*number_voxel,number_hemi： 
    np.transpose(np.reshape(xxx, [number_voxel * number_angle, number_hemisphere])
    再number_angle,number_voxel
    np.transpose(np.reshape(xxx, [number_angle, number_voxel]))
    5. 数组延展
    从角度到体素再到半球是np.repeat,是单独重复；
    从体素到角度，或从半球到体素到角度是np.tile,是批重复
    '''
    ### 观测方向的拦截概率i0c, fcrown 和 fcrown_S 分别是透过率
    fcrown, temp = proportion_directional_crown_voxel_one(lai_crown, std_crown, hcr, rcr, vza0)
    i0c = 1 - fcrown
    fcrown_s, temp = proportion_directional_crown_voxel_one(lai_crown, std_crown, hcr, rcr, sza0)

    ### 半球方向的平均拦截概率
    fcrowntemp, temp = proportion_directional_crown_voxel_one(lai_crown, std_crown, hcr, rcr, hza0)
    i0ctemp = 1 - fcrowntemp
    i0cplus = if2 * np.sum(i0ctemp * fweight0)  # average interception probability

    ####################################################
    #### Crown
    ### dlai单个体素的体密度
    dlai0 = vo * density * G_crown
    dlai = np.tile(dlai0, number_angle)
    ### 体素在空间上的路径长度
    dS = dlai / cthetv * std_crown

    ### 观测方向的透过率 = 观测方向树冠内透过率 * 观测方向树冠间/树冠外的透过率
    gapv_inside_up, plv_inside_up, gapv_inside_down, plv_inside_down = \
        gap_probability_crown_inside_voxel(x, y, z, density, hcr, rcr, vza, vsa)
    xv_up = x + plv_inside_up * sthetv * cphivs
    yv_up = y + plv_inside_up * sthetv * sphivs
    zv_up = z + plv_inside_up * cthetv + hcr
    gapv_outside_up, plv_outside_up, upAreav, hcrv_up, interv = \
        gap_probability_crown_outside_voxel(std_crown, xv_up, yv_up, zv_up, density, hc, hcr, rcr, vza, vsa)
    fv = gapv_inside_up * gapv_outside_up

    ### 半球方向的问题

    hzatemp = np.tile(hza0, number_angle * number_voxel)
    haatemp = np.tile(haa0, number_angle * number_voxel)
    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    cthethtemp = np.cos(np.deg2rad(hzatemp))
    sthethtemp = np.sin(np.deg2rad(hzatemp))
    cphivstemp = np.cos(np.deg2rad(haatemp))
    sphivstemp = np.sin(np.deg2rad(haatemp))
    ztemp = np.repeat(z, number_hemisphere)
    ytemp = np.repeat(y, number_hemisphere)
    xtemp = np.repeat(x, number_hemisphere)

    gapvtemp_inside_up, plvtemp_inside_up, gapvtemp_inside_down, plvtemp_inside_down = \
        gap_probability_crown_inside_voxel(xtemp, ytemp, ztemp, density, hcr, rcr, hzatemp, haatemp)

    ### 往上
    xvtemp_up = xtemp + plvtemp_inside_up * sthethtemp * cphivstemp
    yvtemp_up = ytemp + plvtemp_inside_up * sthethtemp * sphivstemp
    zvtemp_up = ztemp + plvtemp_inside_up * cthethtemp + hcr
    gapvtemp_outside_up, plvtemp_outside_up, upAreavtemp, hcrvtemp_up, intervtemp_up = \
        gap_probability_crown_outside_voxel(std_crown, xvtemp_up, yvtemp_up, zvtemp_up, density, hc, hcr, rcr, hzatemp,
                                            haatemp)
    futemp = gapvtemp_inside_up * gapvtemp_outside_up
    ### 往下
    xvtemp_down = xtemp + plvtemp_inside_down * sthethtemp * cphivstemp
    yvtemp_down = ytemp + plvtemp_inside_down * sthethtemp * sphivstemp
    zvtemp_down = ztemp + plvtemp_inside_down * cthethtemp + hcr
    gapvtemp_outside_down, plvtemp_outside_down, downAreavtemp, hcrvtemp_down, intervtemp_down = \
        gap_probability_crown_outside_voxel(std_crown, xvtemp_down, yvtemp_down, zvtemp_down, density, hc, hcr, rcr,
                                            hzatemp, haatemp)
    fdtemp = gapvtemp_inside_down * gapvtemp_outside_down

    mu = if2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel * number_angle, number_hemisphere])),axis=0)
    md = if2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel * number_angle, number_hemisphere])),axis=0)
    eu = np.transpose(np.reshape(mu * dS * au * fv , [number_angle, number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv , [number_angle, number_voxel]))
    edc = np.sum(ed, axis=0) / i0c
    euc = np.sum(eu, axis=0) / i0c
    pc = 1 - edc - euc

    ######################################################################
    #### 光照和阴影树冠
    ######################################################################
    ### 树冠到树冠 树冠发射率， 树冠本层的贡献率，树冠的反射率
    emic0 = ec * (i0c * pc) * (1 - ec)
    ### 树冠到树干 树冠发射率，树冠下层的贡献率，树干的拦截概率，树干的反射率，树冠的方向透过率
    emic1 = ec * i0cplus * (1-es) * (1-i0c)
    emic = emic0 + emic1
    fss = (1 - fcrown_s) / (lai_crown * G_crown / cthets0)
    emics = emic * fss
    emich = emic - emics


    memis = emics + emich

    return emics,emich, memis


def multiple_scattering_row_spectral_invariance(lai_crown, row_width,row_blank,row_height,
                                                  ec, es, vza0, sza0, vaa0, saa0,raa0, G_crown=0.5):
    '''
    光谱不变理论计算组分有效发射率，可以得到各个组分的结果，是解析方法
    '''

    eps = 0.0001

    number_voxel_height = 50
    number_voxel_width = np.int(row_width/(row_height/number_voxel_height))
    number_voxel_blank = np.int(row_blank/(row_height/number_voxel_height))
    number_voxel = number_voxel_width * number_voxel_height
    number_angle = np.size(vza0)

    ### 与均质和树冠对称结构不同，除了观测天顶角0-90以5为间隔，还需要方位角的离散 0-90以10为间隔

    number_hza = 18
    number_haa = 9
    number_hemisphere = number_hza * number_haa
    if2 = 2.0
    au = 0.5
    ad = 0.5

    rw = row_width
    rb = row_blank
    rh = row_height
    rs = rw + rb
    nvw = number_voxel_width
    nvh = number_voxel_height
    nvb = number_voxel_blank
    nvs = nvw + nvb
    vsa0 = np.abs(saa0 - vaa0)
    sra0 = np.abs(saa0 - raa0) % 180
    vra0 = np.abs(vaa0 - raa0) % 180
    laie = lai_crown * rs / rw
    density = laie / rh

    scale_w = rw / nvw
    scale_h = rh / nvh


    length_layer_width = np.linspace(0.5, nvw - 0.5, nvw) * scale_w

    height_layer = np.linspace(0.5, nvh - 0.5, nvh) * scale_h
    length0 = np.tile(length_layer_width, (nvh))
    height0 = np.asarray(np.tile(np.transpose(np.asmatrix(height_layer)), (1, nvw)))
    height0 = np.reshape(height0,-1)
    length = np.tile(length0, number_angle)
    height = np.tile(height0, number_angle)


    ### 显示半球方位角的区分，然后是半球天顶角
    ### hza1 hza2 hza3 hza1 hza2 hza3...
    ### haa1 haa1 haa1 haa2 haa2 haa2...
    hza0 = np.tile(np.linspace(0.5, number_hza - 0.5, number_hza) / number_hza, number_haa) * 90
    haa0 = np.repeat(np.linspace(0.5, number_haa - 0.5, number_haa) / number_haa, number_hza) * 90
    raa0_h = np.resize(raa0,number_hemisphere)

    stheth0 = np.sin(np.deg2rad(hza0))
    ctheth0 = np.cos(np.deg2rad(hza0))
    fweight0 = stheth0 * ctheth0 * np.pi / 2.0 / number_hemisphere
    tgthh0 = np.tan(np.deg2rad(hza0))

    if lai_crown < eps: lai_crown = eps


    #### 将观测角度扩展到
    vsa0 = np.abs(vaa0 - saa0)
    vza = np.repeat(vza0, number_voxel)
    sza = np.repeat(sza0, number_voxel)
    vsa = np.repeat(vsa0, number_voxel)
    raa = np.repeat(raa0, number_voxel)

    cthetv = np.cos(np.deg2rad(vza))
    tgthv = np.tan(np.deg2rad(vza))
    sthetv = np.sin(np.deg2rad(vza))
    sphivs = np.sin(np.deg2rad(vsa))
    cphivs = np.cos(np.deg2rad(vsa))

    cthetv0 = np.cos(np.deg2rad(vza0))
    tgthv0 = np.tan(np.deg2rad(vza0))
    tgths0 = np.tan(np.deg2rad(sza0))
    cthets0 = np.cos(np.deg2rad(sza0))

    '''
    光谱不变理论：
    从观测方向出发发射光子，出现2种情况：2）光子碰到植被；1）光子透过植被；
    1）光子透过植被，表明这部分是植被下层的贡献；
    2）光子碰到植被，表明这部分是植被这层的贡献，除了碰撞的植被的发射项，光子继续反射，出现3种情况包括：
        a) 从上层逃逸 eu；
        b) 从下层逃逸 ed；
        c) 碰到植被的其他部分 p； 
    经过这层植被作为跳板进入到传感器的贡献： 上层复合贡献*eu + 下层复合贡献 * ed + 本层贡献 * p

    这里通过对其分层或者分体元的方式分别进行计算，然后进行累计，每层或体元计算如下：
    i. 选择层或者体元计算dlai 和 dS
    ii. 计算观测方向的透过率fv；
    iii. 计算上半球的平均透过率 mu；
    iv. 计算下半球的平均透过率 md；
    v. 计算往上方向的逃逸概率 fv * dS * au * mu   到达该层概率；在该层碰撞到概率；往上半球方向概率；从该层到上层逃逸概率
    vi. 计算往下方向的逃逸概率 fv * dS * ad * md 同上
    vii. 计算该层的吸收概率 p = 1 - eu - ed 


    temp： 表示是半球的所有/逐个角度
    plus:  表示半球所有角度的积分

    number_angle number_voxel number_hemi
    先number_Angle*number_voxel,number_hemi
    再number_angle,number_voxel
    '''


    ### 观测方向的拦截概率i0c, fcrown 和 fcrown_S 分别是观测和太阳方向的透过率
    fcrown, temp = proportion_directional_row_voxel_one(lai_crown,row_width,row_blank,row_height,vza0,vaa0,raa0,nvw,nvh,nvb)
    i0c = 1 - fcrown
    fcrown_s, temp = proportion_directional_row_voxel_one(lai_crown,row_width,row_blank,row_height,sza0,saa0,raa0,nvw,nvh,nvb)

    ### 半球方向的平均拦截概率
    fcrowntemp, temp = proportion_directional_row_voxel_one(lai_crown,row_width,row_blank,row_height,hza0,haa0,raa0_h,nvw,nvh,nvb)
    i0ctemp = 1 - fcrowntemp
    i0cplus = if2 * np.sum(i0ctemp * fweight0)  # average interception probability


    ####################################################
    #### Crown

    dlai0 =  density * G_crown * (rh/number_voxel_height)*(rw/number_voxel_width)
    dlai = np.tile(dlai0, number_angle*number_voxel)
    dS = dlai / cthetv
    fv,plv = gap_probability_row_voxel(density,row_width,row_blank,row_height,length,height,vza,raa,vsa)


    fweightemp = np.tile(fweight0, number_voxel * number_angle)
    vzatemp = np.repeat(vza, number_hemisphere)
    raatemp = np.repeat(raa, number_hemisphere)
    vsatemp = np.repeat(vsa, number_hemisphere)
    lengthtemp = np.repeat(length0, number_angle * number_hemisphere)
    heighttemp = np.repeat(height0, number_angle * number_hemisphere)

    futemp,plutemp = gap_probability_row_voxel(density,row_width,row_blank,row_height,lengthtemp,heighttemp,
                                         vzatemp,raatemp,vsatemp)
    fdtemp,pldtemp = gap_probability_row_voxel(density,row_width,row_blank,row_height,lengthtemp,row_height-heighttemp,
                                         vzatemp,raatemp,vsatemp)


    mu = if2 * np.sum(np.transpose(np.reshape(fweightemp * futemp, [number_voxel * number_angle, number_hemisphere])),axis=0)
    md = if2 * np.sum(np.transpose(np.reshape(fweightemp * fdtemp, [number_voxel * number_angle, number_hemisphere])),axis=0)
    eu = np.transpose(np.reshape(mu * dS * au * fv, [number_angle, number_voxel]))
    ed = np.transpose(np.reshape(md * dS * ad * fv, [number_angle, number_voxel]))
    edc = np.sum(ed, axis=0) / i0c
    euc = np.sum(eu, axis=0) / i0c
    pc = 1 - edc - euc

    ######################################################################
    #### 光照和阴影树冠
    ######################################################################
    ### 树冠到树冠 树冠发射率， 树冠本层的贡献率，树冠的反射率
    emic0 = ec * (i0c * pc) * (1 - ec)
    ### 树冠到树干 树冠发射率，树冠下层的贡献率，树干的拦截概率，树干的反射率，树冠的方向透过率
    emic1 = ec * i0cplus * (1 - es) * (1 - i0c)
    emic = emic0 + emic1
    fss = (1 - fcrown_s) / (lai_crown * G_crown / cthets0)
    emics = emic * fss
    emich = emic - emics

    memis = emics + emich

    return  emics, emich, memis


def volscatt(tts,tto,psi,ttl):
    rd = 3.1415926/180.0
    eps = 1e-6
    costs = np.cos(rd * tts)
    costo = np.cos(rd * tto)
    sints = np.sin(rd * tts)
    sinto = np.sin(rd * tto)
    cospsi = np.cos(rd * psi)
    psir = rd * psi
    costl = np.cos(rd * ttl)
    sintl = np.sin(rd * ttl)
    cs = costl * costs
    co = costl * costo
    ss = sintl * sints
    so = sintl * sinto
    number = np.shape(tto)

    cosbts = np.resize(5.0,number)
    cosbto = np.resize(5.0,number)

    ind = np.abs(ss) > eps
    if np.sum(ind) > 0: cosbts[ind] = -cs[ind]/ss[ind]
    ind = np.abs(so) > eps
    if np.sum(ind) > 0: cosbto[ind] = -co[ind]/so[ind]

    bts = np.resize(3.1415926,number)
    ds = cs
    ind = np.abs(cosbts) < 1.0
    if np.sum(ind) > 0:
        bts[ind] = np.arccos(cosbts[ind])
        ds[ind] = ss[ind]
    chi_s = 2.0 /np.pi * ((bts -np.pi * 0.5) * cs + np.sin(bts) * ss)

    bto = np.resize(3.1415926,number)
    doo = co
    ind = np.abs(cosbto) < 1.0
    if np.sum(ind) > 0:
        bto[ind] = np.arccos(cosbto[ind])
        doo[ind] = so[ind]
    chi_o = 2.0 / np.pi * ((bto - np.pi * 0.5) * co + np.sin(bto) * so)

    btran1 = abs(bts - bto)
    btran2 =np.pi - abs(bts + bto -np.pi)

    bt1 = psir * 1.0
    bt2 = btran1 *1.0
    bt3 = btran2 * 1.0
    ind = (psi > btran1) * (psi <= btran2)
    if np.sum(ind) > 0:
        bt1[ind] = btran1[ind] * 1.0
        bt2[ind] = psir[ind] * 1.0
        bt3[ind] = btran2[ind] * 1.0
    ind = (psi > btran1) * (psi > btran2)
    if np.sum(ind) > 0:
        bt1[ind] = btran1[ind]* 1.0
        bt3[ind] = psir[ind]* 1.0
        bt2[ind] = btran2[ind]* 1.0
    t1 = 2.0 * cs * co + ss * so * cospsi
    t2 = 0.0
    ind = bt2 >0
    if np.sum(ind): t2= np.sin(bt2) * (2.0 * ds * doo+ss * so * np.cos(bt1) * np.cos(bt3))
    denom = 2.0 *np.pi *np.pi
    frho = ((np.pi-bt2) * t1 + t2) / denom
    ftau = (-bt2 * t1 + t2) / denom
    frho[frho<0]=0.0
    ftau[ftau<0]=0.0
    return chi_o,chi_s,frho,ftau

def Jfunc1(k,l,t):
    dell=(k - l) * t
    Jout =dell * 1.0

    ind = Jout > 0.001
    if np.sum(ind)>0:
        Jout[ind] = (np.exp(-l * t) - np.exp(-k[ind] * t)) / (k[ind] - l)
    ind = Jout < 0.001
    if np.sum(ind)>0:
        Jout[ind] = 0.5 * t * (np.exp(-k[ind] * t) + np.exp(-l * t)) * (1 - dell[ind] * dell[ind] / 12.0)

    return Jout

def Jfunc2(k,l,t):
    Jout = (1-np.exp(-(k+l)*t))/(k+l)
    return Jout

def scattering_coefficient(lai, thm, rho, tau, sza, vza, raa, iftrunk=0):
    rd = np.pi / 180.0
    eps = 0.001
    sthetv = np.sin(vza * rd)
    cthetv = np.cos(vza * rd)
    sthets = np.sin(sza * rd)
    cthets = np.cos(sza * rd)
    sphi = np.sin(raa * rd)
    cphi = np.cos(raa * rd)
    calph = sthets * sthetv * cphi + cthetv * cthets
    tgthv = sthetv / cthetv
    tgths = sthets / cthets
    cts = cthets
    cto = cthetv
    tants = tgths
    tanto = tgthv
    ctscto = cts * cto
    cospsi = cphi
    sinpsi = sphi
    dso = np.sqrt(tants * tants + tanto * tanto -
                  2.0 * tants * tanto * cospsi)
    tto = vza
    tts = sza
    psi = raa

    ks = 0.0
    ko = 0.0
    bf = 0.0
    sob = 0.0
    sof = 0.0

    if (iftrunk == 0):
        na = 18
        lna = np.linspace(0, 17, 18)
        thetal = 2.5 + 5.0 * lna
        lidf = leaf_inclination_distribution_function(thm)
        for i in range(na):
            ttl = thetal[i]
            ctl = np.cos(ttl * rd)
            [chi_s, chi_o, frho, ftau] = \
                volscatt(tto, tts, psi, ttl)
            ksli = chi_s / cts
            koli = chi_o / cto
            sobli = frho * np.pi / ctscto
            sofli = ftau * np.pi / ctscto
            bfli = ctl * ctl
            ks = ks + ksli * lidf[i]
            ko = ko + koli * lidf[i]
            bf = bf + bfli * lidf[i]
            sob = sob + sobli * lidf[i]
            sof = sof + sofli * lidf[i]
        ks = 0.5 / cts
        ko = 0.5 / cto
    else:
        ttl = 89.9
        ctl = np.cos(ttl * rd)
        [chi_s, chi_o, frho, ftau] = \
            volscatt(tto, tts, psi, ttl)
        ks = chi_s / cts
        ko = chi_o / cto
        sob = frho * np.pi / ctscto
        sof = ftau * np.pi / ctscto
        bf = ctl * ctl

    sdb = 0.5 * (ks + bf)
    sdf = 0.5 * (ks - bf)
    dob = 0.5 * (ko + bf)
    dof = 0.5 * (ko - bf)
    ddb = 0.5 * (1.0 + bf)
    ddf = 0.5 * (1.0 - bf)

    sigb = ddb * rho + ddf * tau
    sigf = ddf * rho + ddb * tau
    att = 1.0 - sigf
    m2 = (att + sigb) * (att - sigb)

    ind = (m2 > 0)
    m2 = m2 * ind
    m = np.sqrt(m2)
    sb = sdb * rho + sdf * tau
    sf = sdf * rho + sdb * tau
    vb = dob * rho + dof * tau
    vf = dof * rho + dob * tau
    w = sob * rho + sof * tau
    e1 = np.exp(-m * lai)
    e2 = e1 * e1
    rinf = (att - m) / sigb
    rinf2 = rinf * rinf
    re = rinf * e1
    denom = 1.0 - rinf2 * e2
    J1ks = Jfunc1(ks, m, lai)
    J2ks = Jfunc2(ks, m, lai)
    J1ko = Jfunc1(ko, m, lai)
    J2ko = Jfunc2(ko, m, lai)
    Ps = (sf + sb * rinf) * J1ks
    Qs = (sf * rinf + sb) * J2ks
    Pv = (vf + vb * rinf) * J1ko
    Qv = (vf * rinf + vb) * J2ko

    pdd = rinf * (1.0 - e2) / denom
    tdd = (1.0 - rinf2) * e1 / denom
    tsd = (Ps - re * Qs) / denom
    psd = (Qs - re * Ps) / denom
    tdo = (Pv - re * Qv) / denom
    pdo = (Qv - re * Pv) / denom

    gammasdf = (1 + rinf) * (J1ks - re * J2ks) / denom
    gammasdb = (1 + rinf) * (-re * J1ks + J2ks) / denom

    return dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf


def multiple_scattering_crown_sail_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                             lai_crown, lai_trunk, lai_unveg, std_crown, hcr_crown, rcr_crown,
                                             ec, et, ev, es,
                                             vza0, sza0, vaa0,
                                             saa0=0, G_crown = 0.5, agl_crown = 54, agl_unveg = 54):
    '''通过4SAIL的方法计算多次散射项，通过JACOB方法进行了方程求解，是数值方法'''
    eps = 0.001
    number_crown = 10
    number_trunk = 10
    number_unveg = 10
    number_soil = 1
    number_layer = number_crown + number_trunk + number_unveg + number_soil

    vsa0 = np.abs(vaa0 - saa0)
    vza = np.tile(vza0,[number_layer,1])
    sza = np.tile(sza0,[number_layer,1])
    vsa = np.tile(vsa0,[number_layer,1])
    number_angle = np.size(vza0)

    Ems = np.zeros([number_layer,number_angle])  # thermal emission for each layer
    Eplui = np.zeros([number_layer,number_angle] ) # upward emission for each layer
    Emini = np.zeros([number_layer,number_angle])  # downward emission for each layer
    mi = np.zeros([number_layer,number_angle])
    LAIi = np.zeros([number_layer,number_angle])  # LAI for each layer, like lavd
    rinfi = np.zeros([number_layer,number_angle])
    fHc = np.zeros([number_layer,number_angle])   # probability for vegetation in each layer
    vbi = np.zeros([number_layer,number_angle])   # probability for backward hemisphere space
    vfi = np.zeros([number_layer,number_angle])    # probability for forward hemisphere space
    Po = np.zeros([number_layer,number_angle])  # visible proportions for each layer
    Rss = 10
    Rsh = 10
    Rcs = 10
    Rch = 10
    Rvs = 10
    Rvh = 10
    Rts = 10
    Rth = 10

    ks = Pss/(Pss+Psh)
    Emsi = Rss * ks + (1-ks) * Rsh
    kc = Pcs/(Pcs+Pch)
    Emci = Rcs * kc + (1-kc) * Rch
    kv = Pvs/(Pvs+Pvh)
    Emvi = Rvs * kv + (1-kv) * Rvh

    kt = np.zeros(np.size(kv))
    ind = Pth + Pts >0
    if np.sum(ind) >0 : kt[ind] = Pts[ind] / (Pts[ind] + Pth[ind])
    Emti = Rts * kt + (1-kt) * Rth

    Emsi = np.tile(Emsi,[number_soil,1])
    Emci = np.tile(Emci,[number_crown,1])
    Emvi = np.tile(Emvi,[number_unveg,1])
    Emti = np.tile(Emti, [number_trunk,1])

    layertemp = 0
    #################################################################
    #### Crown
    #################################################################
    referenceZa = np.asarray([40])
    vol_crown = 4.0 * np.pi / 3 * rcr_crown * rcr_crown * hcr_crown
    area_crown = np.pi * rcr_crown * rcr_crown
    density = lai_crown / (std_crown * vol_crown)
    tgthx = np.tan(np.deg2rad(referenceZa))
    cthetx = np.cos(np.deg2rad(referenceZa))
    hc_crown = 2 * hcr_crown
    [upArea, upVol] = crosscutting_ellipsoid(np.asarray([0]), referenceZa, hc_crown, hcr_crown, rcr_crown)
    pl = upVol / upArea / cthetx
    gapx = np.exp(-density * G_crown * pl)
    LAIeff = cthetx * std_crown * upArea * (1 - gapx)
    iLAI = LAIeff / number_crown
    iPo = (Pcs + Pch) / number_crown

    layer1 = layertemp
    layer2 = layertemp + number_crown
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_crown, agl_crown, 1-ec, 0, sza, vza, vsa)
    vbi[layer1:layer2,:] = vb[layer1:layer2,:]
    vfi[layer1:layer2,:] = vf[layer1:layer2,:]
    Ems[layer1:layer2,:] = Emci

    fHc[layer1:layer2,:] = (iLAI * m *(1-rinf))*ec
    Po[layer1:layer2,:] = iPo
    mi[layer1:layer2,:] = m
    LAIi[layer1:layer2,:] = iLAI
    rinfi[layer1:layer2,:] = rinf
    layertemp = layertemp + number_crown

    #################################################################
    #### Trunk
    #################################################################
    # iLAI = lai_crown / number_crown
    # iPo = Pc / number_crown
    # layer1 = layertemp
    # layer2 = layertemp + number_crown
    # [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
    #     scattering_coefficient(lai_crown, agl_crown, rc, tc, sza, vza, vsa)
    # vbi[layer1:layer2] = vb
    # vfi[layer1:layer2] = vf
    # Ems[layer1:layer2] = 10.0
    # fHc[layer1:layer2] = (iLAI * m *(1-rinf))*(1-rc-tc)
    # Po[layer1:layer2] = iPo
    # mi[layer1:layer2] = m
    # LAIi[layer1:layer2] = iLAI
    # rinfi[layer1:layer2] = rinf
    layertemp = layertemp + number_trunk

    #################################################################
    #### unveg
    #################################################################
    iLAI = lai_unveg / number_unveg
    iPo = (Pvs+Pvh) / number_unveg
    layer1 = layertemp
    layer2 = layertemp + number_unveg
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_unveg, agl_unveg, 1-ev, 0, sza, vza, vsa)
    vbi[layer1:layer2, :] = vb[layer1:layer2, :]
    vfi[layer1:layer2, :] = vf[layer1:layer2, :]
    Ems[layer1:layer2, :] = Emci

    fHc[layer1:layer2, :] = (iLAI * m * (1 - rinf)) * ec
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = iLAI
    rinfi[layer1:layer2, :] = rinf
    layertemp = layertemp + number_unveg

    ###################################################################
    ###  SOIL
    ###################################################################

    layer1 = layertemp
    layer2 = layertemp + number_soil
    fHs = (1 - rinf * rinf) * es / (1 - rinf * (1-es)) * es


    Ems[layer1:layer2, :] = Emsi
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = 1.0
    rinfi[layer1:layer2, :] = rinf

    layertemp = layertemp + number_soil

    ######################################################################
    cont = 1
    rs = 1-es
    fbottom = (rs - rinf) / (1 - rinf * rs)
    F1 = np.zeros([number_layer,number_angle])
    F2 = np.zeros([number_layer,number_angle])
    F1top = 0
    F1[number_layer - 1,:] = 0
    count = 0
    number_grid = number_layer - 1
    while (cont):
        F1topn = -rinf * F2[0,:]
        F1[0,:] = F1topn
        for j in range(number_grid):
            F1[j + 1,:] = F1[j,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        F2[number_grid,:] = fbottom * F1[number_grid,:] + Ems[number_grid,:] * fHs
        for j in range(number_grid - 1, -1, -1):
            F2[j,:] = F2[j + 1,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        count = count + 1
        if count > 10: break
        cont = np.max(np.abs(F1topn - F1top)) > 0.01

    Emini = (F1 + rinfi * F2) / (1 - rinfi * rinfi)
    Eplui = (F2 + rinfi * F1) / (1 - rinfi * rinfi)


    diffuseRad0 = np.sum(LAIi[1:number_grid] *(Emini[1:number_grid] * vbi[1:number_grid] + Eplui[1:number_grid] * vfi[1:number_grid]) * Po[1:number_grid])
    diffuseRad1 = diffuseRad0 + (LAIi[0] * (Emini[0] * vbi[0] + Eplui[0] * vfi[0]) * Po[0])
    diffuseRad2 = diffuseRad1 + Emini[number_grid] * rs * Po[number_grid]

    return diffuseRad2/10


def multiple_scattering_hom_sail_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                           lai_crown, lai_trunk, lai_unveg,
                                           ec, et, ev, es,
                                           vza0, sza0, vaa0,
                                           saa0=0, agl_crown = 54, agl_unveg = 54):
    '''通过4SAIL的方法计算多次散射项，通过JACOB方法进行了方程求解，是数值方法'''
    eps = 0.001
    number_crown = 10
    number_trunk = 10
    number_unveg = 10
    number_soil = 1
    number_layer = number_crown + number_trunk + number_unveg + number_soil

    vsa0 = np.abs(vaa0 - saa0)
    vza = np.tile(vza0,[number_layer,1])
    sza = np.tile(sza0,[number_layer,1])
    vsa = np.tile(vsa0,[number_layer,1])
    number_angle = np.size(vza0)

    Ems = np.zeros([number_layer,number_angle])  # thermal emission for each layer
    Eplui = np.zeros([number_layer,number_angle] ) # upward emission for each layer
    Emini = np.zeros([number_layer,number_angle])  # downward emission for each layer
    mi = np.zeros([number_layer,number_angle])
    LAIi = np.zeros([number_layer,number_angle])  # LAI for each layer, like lavd
    rinfi = np.zeros([number_layer,number_angle])
    fHc = np.zeros([number_layer,number_angle])   # probability for vegetation in each layer
    vbi = np.zeros([number_layer,number_angle])   # probability for backward hemisphere space
    vfi = np.zeros([number_layer,number_angle])    # probability for forward hemisphere space
    Po = np.zeros([number_layer,number_angle])  # visible proportions for each layer

    Rss = 10
    Rsh = 10
    Rcs = 10
    Rch = 10
    Rvs = 10
    Rvh = 10
    Rts = 10
    Rth = 10

    ks = Pss/(Pss+Psh)
    Emsi = Rss * ks + (1-ks) * Rsh
    kc = Pcs/(Pcs+Pch)
    Emci = Rcs * kc + (1-kc) * Rch
    kv = Pvs/(Pvs+Pvh)
    Emvi = Rvs * kv + (1-kv) * Rvh

    kt = np.zeros(np.size(kv))
    ind = Pth + Pts >0
    if np.sum(ind) >0 : kt[ind] = Pts[ind] / (Pts[ind] + Pth[ind])
    Emti = Rts * kt + (1-kt) * Rth

    Emsi = np.tile(Emsi,[number_soil,1])
    Emci = np.tile(Emci,[number_crown,1])
    Emvi = np.tile(Emvi,[number_unveg,1])
    Emti = np.tile(Emti, [number_trunk,1])

    layertemp = 0
    #################################################################
    #### Crown
    #################################################################
    iLAI = lai_crown / number_crown
    iPo = (Pcs+Pch) / number_crown
    layer1 = layertemp
    layer2 = layertemp + number_crown
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_crown, agl_crown, 1-ec, 0, sza, vza, vsa)
    vbi[layer1:layer2,:] = vb[layer1:layer2,:]
    vfi[layer1:layer2,:] = vf[layer1:layer2,:]
    Ems[layer1:layer2,:] = Emci

    fHc[layer1:layer2,:] = (iLAI * m *(1-rinf))*ec
    Po[layer1:layer2,:] = iPo
    mi[layer1:layer2,:] = m
    LAIi[layer1:layer2,:] = iLAI
    rinfi[layer1:layer2,:] = rinf
    layertemp = layertemp + number_crown

    #################################################################
    #### Trunk
    #################################################################
    # iLAI = lai_crown / number_crown
    # iPo = Pc / number_crown
    # layer1 = layertemp
    # layer2 = layertemp + number_crown
    # [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
    #     scattering_coefficient(lai_crown, agl_crown, rc, tc, sza, vza, vsa)
    # vbi[layer1:layer2] = vb
    # vfi[layer1:layer2] = vf
    # Ems[layer1:layer2] = 10.0
    # fHc[layer1:layer2] = (iLAI * m *(1-rinf))*(1-rc-tc)
    # Po[layer1:layer2] = iPo
    # mi[layer1:layer2] = m
    # LAIi[layer1:layer2] = iLAI
    # rinfi[layer1:layer2] = rinf
    layertemp = layertemp + number_trunk

    #################################################################
    #### unveg
    #################################################################
    iLAI = lai_unveg / number_unveg
    iPo = (Pvs+Pvh) / number_unveg
    layer1 = layertemp
    layer2 = layertemp + number_unveg
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_unveg, agl_unveg, 1-ev, 0, sza, vza, vsa)
    vbi[layer1:layer2, :] = vb[layer1:layer2, :]
    vfi[layer1:layer2, :] = vf[layer1:layer2, :]
    Ems[layer1:layer2, :] = Emci

    fHc[layer1:layer2, :] = (iLAI * m * (1 - rinf)) * ec
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = iLAI
    rinfi[layer1:layer2, :] = rinf
    layertemp = layertemp + number_unveg

    ###################################################################
    ###  SOIL
    ###################################################################

    layer1 = layertemp
    layer2 = layertemp + number_soil
    fHs = (1 - rinf * rinf) * es / (1 - rinf * (1-es)) * es


    Ems[layer1:layer2, :] = Emsi
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = 1.0
    rinfi[layer1:layer2, :] = rinf

    layertemp = layertemp + number_soil

    ######################################################################
    cont = 1
    rs = 1-es
    fbottom = (rs - rinf) / (1 - rinf * rs)
    F1 = np.zeros([number_layer,number_angle])
    F2 = np.zeros([number_layer,number_angle])
    F1top = 0
    F1[number_layer - 1,:] = 0
    count = 0
    number_grid = number_layer - 1
    while (cont):
        F1topn = -rinf * F2[0,:]
        F1[0,:] = F1topn
        for j in range(number_grid):
            F1[j + 1,:] = F1[j,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        F2[number_grid,:] = fbottom * F1[number_grid,:] + Ems[number_grid,:] * fHs
        for j in range(number_grid - 1, -1, -1):
            F2[j,:] = F2[j + 1,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        count = count + 1
        if count > 10: break
        cont = np.max(np.abs(F1topn - F1top)) > 0.01

    Emini = (F1 + rinfi * F2) / (1 - rinfi * rinfi)
    Eplui = (F2 + rinfi * F1) / (1 - rinfi * rinfi)


    diffuseRad0 = np.sum(LAIi[1:number_grid] *(Emini[1:number_grid] * vbi[1:number_grid] + Eplui[1:number_grid] * vfi[1:number_grid]) * Po[1:number_grid])
    diffuseRad1 = diffuseRad0 + (LAIi[0] * (Emini[0] * vbi[0] + Eplui[0] * vfi[0]) * Po[0])
    diffuseRad2 = diffuseRad1 + Emini[number_grid] * rs * Po[number_grid]

    return diffuseRad2/10


def multiple_scattering_crown_sail_rad_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                                 lai_crown, lai_trunk, lai_unveg, std_crown, hcr_crown, rcr_crown,
                                                 ec, et, ev, es, Tss, Tsh, Tcs, Tch, Tvs, Tvh, Tts, Tth,
                                                 vza0, sza0, vaa0,
                                                 saa0=0, G_crown = 0.5, agl_crown = 54, agl_unveg = 54):
    '''通过4SAIL的方法计算多次散射项，通过JACOB方法进行了方程求解，是数值方法'''
    eps = 0.001
    number_crown = 10
    number_trunk = 10
    number_unveg = 10
    number_soil = 1
    number_layer = number_crown + number_trunk + number_unveg + number_soil

    vsa0 = np.abs(vaa0 - saa0)
    vza = np.tile(vza0,[number_layer,1])
    sza = np.tile(sza0,[number_layer,1])
    vsa = np.tile(vsa0,[number_layer,1])
    number_angle = np.size(vza0)

    Ems = np.zeros([number_layer,number_angle])  # thermal emission for each layer
    Eplui = np.zeros([number_layer,number_angle] ) # upward emission for each layer
    Emini = np.zeros([number_layer,number_angle])  # downward emission for each layer
    mi = np.zeros([number_layer,number_angle])
    LAIi = np.zeros([number_layer,number_angle])  # LAI for each layer, like lavd
    rinfi = np.zeros([number_layer,number_angle])
    fHc = np.zeros([number_layer,number_angle])   # probability for vegetation in each layer
    vbi = np.zeros([number_layer,number_angle])   # probability for backward hemisphere space
    vfi = np.zeros([number_layer,number_angle])    # probability for forward hemisphere space
    Po = np.zeros([number_layer,number_angle])  # visible proportions for each layer
    Rss = planck(Tss)
    Rsh = planck(Tsh)
    Rcs = planck(Tcs)
    Rch = planck(Tch)
    Rvs = planck(Tvs)
    Rvh = planck(Tvh)
    Rts = planck(Tts)
    Rth = planck(Tth)


    ks = Pss/(Pss+Psh)
    Emsi = Rss * ks + (1-ks) * Rsh
    kc = Pcs/(Pcs+Pch)
    Emci = Rcs * kc + (1-kc) * Rch
    kv = Pvs/(Pvs+Pvh)
    Emvi = Rvs * kv + (1-kv) * Rvh

    kt = np.zeros(np.size(kv))
    ind = Pth + Pts >0
    if np.sum(ind) >0 : kt[ind] = Pts[ind] / (Pts[ind] + Pth[ind])
    Emti = Rts * kt + (1-kt) * Rth

    Emsi = np.tile(Emsi,[number_soil,1])
    Emci = np.tile(Emci,[number_crown,1])
    Emvi = np.tile(Emvi,[number_unveg,1])
    Emti = np.tile(Emti, [number_trunk,1])

    layertemp = 0
    #################################################################
    #### Crown
    #################################################################
    referenceZa = np.asarray([40])
    vol_crown = 4.0 * np.pi / 3 * rcr_crown * rcr_crown * hcr_crown
    area_crown = np.pi * rcr_crown * rcr_crown
    density = lai_crown / (std_crown * vol_crown)
    tgthx = np.tan(np.deg2rad(referenceZa))
    cthetx = np.cos(np.deg2rad(referenceZa))
    hc_crown = 2 * hcr_crown
    [upArea, upVol] = crosscutting_ellipsoid(np.asarray([0]), referenceZa, hc_crown, hcr_crown, rcr_crown)
    pl = upVol / upArea / cthetx
    gapx = np.exp(-density * G_crown * pl)
    LAIeff = cthetx * std_crown * upArea * (1 - gapx)
    iLAI = LAIeff / number_crown
    iPo = (Pcs + Pch) / number_crown

    layer1 = layertemp
    layer2 = layertemp + number_crown
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_crown, agl_crown, 1-ec, 0, sza, vza, vsa)
    vbi[layer1:layer2,:] = vb[layer1:layer2,:]
    vfi[layer1:layer2,:] = vf[layer1:layer2,:]
    Ems[layer1:layer2,:] = Emci

    fHc[layer1:layer2,:] = (iLAI * m *(1-rinf))*ec
    Po[layer1:layer2,:] = iPo
    mi[layer1:layer2,:] = m
    LAIi[layer1:layer2,:] = iLAI
    rinfi[layer1:layer2,:] = rinf
    layertemp = layertemp + number_crown

    #################################################################
    #### Trunk
    #################################################################
    # iLAI = lai_crown / number_crown
    # iPo = Pc / number_crown
    # layer1 = layertemp
    # layer2 = layertemp + number_crown
    # [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
    #     scattering_coefficient(lai_crown, agl_crown, rc, tc, sza, vza, vsa)
    # vbi[layer1:layer2] = vb
    # vfi[layer1:layer2] = vf
    # Ems[layer1:layer2] = 10.0
    # fHc[layer1:layer2] = (iLAI * m *(1-rinf))*(1-rc-tc)
    # Po[layer1:layer2] = iPo
    # mi[layer1:layer2] = m
    # LAIi[layer1:layer2] = iLAI
    # rinfi[layer1:layer2] = rinf
    layertemp = layertemp + number_trunk

    #################################################################
    #### unveg
    #################################################################
    iLAI = lai_unveg / number_unveg
    iPo = (Pvs+Pvh) / number_unveg
    layer1 = layertemp
    layer2 = layertemp + number_unveg
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_unveg, agl_unveg, 1-ev, 0, sza, vza, vsa)
    vbi[layer1:layer2, :] = vb[layer1:layer2, :]
    vfi[layer1:layer2, :] = vf[layer1:layer2, :]
    Ems[layer1:layer2, :] = Emvi

    fHc[layer1:layer2, :] = (iLAI * m * (1 - rinf)) * ec
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = iLAI
    rinfi[layer1:layer2, :] = rinf
    layertemp = layertemp + number_unveg

    ###################################################################
    ###  SOIL
    ###################################################################

    layer1 = layertemp
    layer2 = layertemp + number_soil
    fHs = (1 - rinf * rinf) * es / (1 - rinf * (1-es)) * es


    Ems[layer1:layer2, :] = Emsi
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = 1.0
    rinfi[layer1:layer2, :] = rinf

    layertemp = layertemp + number_soil

    ######################################################################
    cont = 1
    rs = 1-es
    fbottom = (rs - rinf) / (1 - rinf * rs)
    F1 = np.zeros([number_layer,number_angle])
    F2 = np.zeros([number_layer,number_angle])
    F1top = 0
    F1[number_layer - 1,:] = 0
    count = 0
    number_grid = number_layer - 1
    while (cont):
        F1topn = -rinf * F2[0,:]
        F1[0,:] = F1topn
        for j in range(number_grid):
            F1[j + 1,:] = F1[j,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        F2[number_grid,:] = fbottom * F1[number_grid,:] + Ems[number_grid,:] * fHs
        for j in range(number_grid - 1, -1, -1):
            F2[j,:] = F2[j + 1,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        count = count + 1
        if count > 10: break
        cont = np.max(np.abs(F1topn - F1top)) > 0.01

    Emini = (F1 + rinfi * F2) / (1 - rinfi * rinfi)
    Eplui = (F2 + rinfi * F1) / (1 - rinfi * rinfi)


    diffuseRad0 = np.sum(LAIi[1:number_grid] *(Emini[1:number_grid] * vbi[1:number_grid] + Eplui[1:number_grid] * vfi[1:number_grid]) * Po[1:number_grid])
    diffuseRad1 = diffuseRad0 + (LAIi[0] * (Emini[0] * vbi[0] + Eplui[0] * vfi[0]) * Po[0])
    diffuseRad2 = diffuseRad1 + Emini[number_grid] * rs * Po[number_grid]

    return diffuseRad2


def multiple_scattering_hom_sail_rad_endmember(Pss, Psh, Pcs, Pch, Pvs, Pvh, Pts, Pth,
                                               lai_crown, lai_trunk, lai_unveg,
                                               ec, et, ev, es, Rss, Rsh, Rcs, Rch, Rvs, Rvh, Rts, Rth,
                                               vza0, sza0, vaa0,
                                               saa0=0, agl_crown = 54, agl_unveg = 54):
    '''通过4SAIL的方法计算多次散射项，通过JACOB方法进行了方程求解，是数值方法'''
    eps = 0.001
    number_crown = 10
    number_trunk = 10
    number_unveg = 10
    number_soil = 1
    number_layer = number_crown + number_trunk + number_unveg + number_soil

    vsa0 = np.abs(vaa0 - saa0)
    vza = np.tile(vza0,[number_layer,1])
    sza = np.tile(sza0,[number_layer,1])
    vsa = np.tile(vsa0,[number_layer,1])
    number_angle = np.size(vza0)

    Ems = np.zeros([number_layer,number_angle])  # thermal emission for each layer
    Eplui = np.zeros([number_layer,number_angle] ) # upward emission for each layer
    Emini = np.zeros([number_layer,number_angle])  # downward emission for each layer
    mi = np.zeros([number_layer,number_angle])
    LAIi = np.zeros([number_layer,number_angle])  # LAI for each layer, like lavd
    rinfi = np.zeros([number_layer,number_angle])
    fHc = np.zeros([number_layer,number_angle])   # probability for vegetation in each layer
    vbi = np.zeros([number_layer,number_angle])   # probability for backward hemisphere space
    vfi = np.zeros([number_layer,number_angle])    # probability for forward hemisphere space
    Po = np.zeros([number_layer,number_angle])  # visible proportions for each layer



    ks = Pss/(Pss+Psh)
    Emsi = Rss * ks + (1-ks) * Rsh
    kc = Pcs/(Pcs+Pch)
    Emci = Rcs * kc + (1-kc) * Rch
    kv = Pvs/(Pvs+Pvh)
    Emvi = Rvs * kv + (1-kv) * Rvh

    kt = np.zeros(np.size(kv))
    ind = Pth + Pts >0
    if np.sum(ind) >0 : kt[ind] = Pts[ind] / (Pts[ind] + Pth[ind])
    Emti = Rts * kt + (1-kt) * Rth

    Emsi = np.tile(Emsi,[number_soil,1])
    Emci = np.tile(Emci,[number_crown,1])
    Emvi = np.tile(Emvi,[number_unveg,1])
    Emti = np.tile(Emti, [number_trunk,1])

    layertemp = 0
    #################################################################
    #### Crown
    #################################################################
    iLAI = lai_crown / number_crown
    iPo = (Pcs+Pch) / number_crown
    layer1 = layertemp
    layer2 = layertemp + number_crown
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_crown, agl_crown, 1-ec, 0, sza, vza, vsa)
    vbi[layer1:layer2,:] = vb[layer1:layer2,:]
    vfi[layer1:layer2,:] = vf[layer1:layer2,:]
    Ems[layer1:layer2,:] = Emci

    fHc[layer1:layer2,:] = (iLAI * m *(1-rinf))*ec
    Po[layer1:layer2,:] = iPo
    mi[layer1:layer2,:] = m
    LAIi[layer1:layer2,:] = iLAI
    rinfi[layer1:layer2,:] = rinf
    layertemp = layertemp + number_crown

    #################################################################
    #### Trunk
    #################################################################
    # iLAI = lai_crown / number_crown
    # iPo = Pc / number_crown
    # layer1 = layertemp
    # layer2 = layertemp + number_crown
    # [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
    #     scattering_coefficient(lai_crown, agl_crown, rc, tc, sza, vza, vsa)
    # vbi[layer1:layer2] = vb
    # vfi[layer1:layer2] = vf
    # Ems[layer1:layer2] = 10.0
    # fHc[layer1:layer2] = (iLAI * m *(1-rinf))*(1-rc-tc)
    # Po[layer1:layer2] = iPo
    # mi[layer1:layer2] = m
    # LAIi[layer1:layer2] = iLAI
    # rinfi[layer1:layer2] = rinf
    layertemp = layertemp + number_trunk

    #################################################################
    #### unveg
    #################################################################
    iLAI = lai_unveg / number_unveg
    iPo = (Pvs+Pvh) / number_unveg
    layer1 = layertemp
    layer2 = layertemp + number_unveg
    [dso, pdd, tdd, tsd, psd, tdo, pdo, gammasdf, gammasdb, vf, vb, m, rinf] = \
        scattering_coefficient(lai_unveg, agl_unveg, 1-ev, 0, sza, vza, vsa)
    vbi[layer1:layer2, :] = vb[layer1:layer2, :]
    vfi[layer1:layer2, :] = vf[layer1:layer2, :]
    Ems[layer1:layer2, :] = Emvi

    fHc[layer1:layer2, :] = (iLAI * m * (1 - rinf)) * ec
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = iLAI
    rinfi[layer1:layer2, :] = rinf
    layertemp = layertemp + number_unveg

    ###################################################################
    ###  SOIL
    ###################################################################

    layer1 = layertemp
    layer2 = layertemp + number_soil
    fHs = (1 - rinf * rinf) * es / (1 - rinf * (1-es)) * es


    Ems[layer1:layer2, :] = Emsi
    Po[layer1:layer2, :] = iPo
    mi[layer1:layer2, :] = m
    LAIi[layer1:layer2, :] = 1.0
    rinfi[layer1:layer2, :] = rinf

    layertemp = layertemp + number_soil

    ######################################################################
    cont = 1
    rs = 1-es
    fbottom = (rs - rinf) / (1 - rinf * rs)
    F1 = np.zeros([number_layer,number_angle])
    F2 = np.zeros([number_layer,number_angle])
    F1top = 0
    F1[number_layer - 1,:] = 0
    count = 0
    number_grid = number_layer - 1
    while (cont):
        F1topn = -rinf * F2[0,:]
        F1[0,:] = F1topn
        for j in range(number_grid):
            F1[j + 1,:] = F1[j,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        F2[number_grid,:] = fbottom * F1[number_grid,:] + Ems[number_grid,:] * fHs
        for j in range(number_grid - 1, -1, -1):
            F2[j,:] = F2[j + 1,:] * (1 - mi[j,:] * LAIi[j,:]) + Ems[j,:] * fHc[j,:]
        count = count + 1
        if count > 10: break
        cont = np.max(np.abs(F1topn - F1top)) > 0.01

    Emini = (F1 + rinfi * F2) / (1 - rinfi * rinfi)
    Eplui = (F2 + rinfi * F1) / (1 - rinfi * rinfi)


    diffuseRad0 = np.sum(LAIi[1:number_grid] *(Emini[1:number_grid] * vbi[1:number_grid] + Eplui[1:number_grid] * vfi[1:number_grid]) * Po[1:number_grid])
    diffuseRad1 = diffuseRad0 + (LAIi[0] * (Emini[0] * vbi[0] + Eplui[0] * vfi[0]) * Po[0])
    diffuseRad2 = diffuseRad1 + Emini[number_grid] * rs * Po[number_grid]

    return diffuseRad2
