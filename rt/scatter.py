import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rt.gap import *
from rt.hotspot import *
from rt.proportion import *

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
    vza_index_bottom = np.asarray(vza / 10, dtype=int)
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
    vza_index = np.asarray(vza / 10, dtype=np.int_)
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
    # scheme = quadpy.sn.stroud_1967_7_c(3)
    # voxels = scheme.points
    # vo = scheme.weights * vol_crown
    voxels,number_voxel = ellipsoid_grid_crown_FRT()
    vo = voxels[3,:]/np.sum(voxels[3,:])*vol_crown

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
    number_voxel_width = np.int_(row_width/(row_height/number_voxel_height))
    number_voxel_blank = np.int_(row_blank/(row_height/number_voxel_height))
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

