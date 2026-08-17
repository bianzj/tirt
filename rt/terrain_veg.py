import numpy as np

from rt.crown import *
from base.util_data import *
###----------------------------------
### 这里考虑了山地本身的遮挡和子像元内植被的复合效应
### 1. determine the mountainous surface
### 2. determine the forest strcture
###----------------------------------



class Terrain_Veg:

    ### h, r, density
    shapes = np.asarray([[10,10,0.003],[10,5,0.001],[10,15,0.001]])
    ### lai, hspot, stand, radi_horizontal, radi_vertical
    forestshape = np.asarray([1.0,0.001,3,3,0.2])

    crown = Crown()

    n_shape,n_dim = np.shape(shapes)

    Es = 0.955
    Em = 0.975

    Tss = 320
    Tsh = 305
    Tvs = 303
    Tvh = 300



    vza_ = np.asarray([0,55])
    vaa_ = np.asarray([90,90])
    sza = 20
    saa = 250
    n_part = 10
    n_angle = np.size(vza_)
    canopy_effective_emissivity = []
    component_effective_emissivity = []
    canopy_brightness_temperature = []


    ####----------------------------------------------
    ### PLANE AND MOUNTAIN surface first
    ###-----------------------------------------------
    def calculate_component_direct_emissivity(self, kangle,ifP = 0):
        shapes = self.shapes
        Esoil = self.Es
        n_shape = self.n_shape
        n_dem = self.n_dim
        n_part = self.n_part
        vza = self.vza_[kangle]
        vaa = self.vaa_[kangle]
        sza = self.sza
        saa = self.saa
        raa = np.abs(vaa - saa)

        if raa > 180: raa = 360- raa



        rd = np.pi / 180.0
        ui = np.cos(sza * rd)
        uv = np.cos(vza * rd)
        si = np.sin(sza * rd)
        sv = np.sin(vza * rd)
        up = np.cos(raa * rd)
        tantv = np.tan(vza * rd)
        tants = np.tan(sza * rd)

        ####----------------------------------------------
        ### PLANE surface
        ### 确定 山地的结构的情况
        ###-----------------------------------------------
        height1r = 0
        projv = 0
        projs = 0
        poccupied = 0
        projvs = 0
        for kshape2 in range(n_shape):
            shape2 = shapes[kshape2]
            radius2 = shape2[1]
            height2 = shape2[0]
            density2 = shape2[2]
            poccupied = poccupied + density2 * np.pi * (radius2*radius2)
            dh = height2 - height1r
            alpha2 = np.arctan(radius2/height2)
            L2_v = dh * tantv
            if L2_v < radius2: L2_v = radius2*1.0
            theta2_v = np.arctan(L2_v*np.tan(alpha2)/radius2)
            gamma2_v = np.arcsin(radius2/L2_v)
            if (dh > 0) + (L2_v > radius2):
                projv_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v - np.pi/2.0)*radius2*radius2
                projv = projv + projv_mount
            L2_s = dh * tants
            if L2_s < radius2: L2_s = radius2*1.0
            theta2_s = np.arctan(L2_s*np.tan(alpha2)/radius2)
            gamma2_s = np.arcsin(radius2/L2_s)
            if (dh > 0) + (L2_s > radius2):
                projs_mount = density2 * (1.0/np.tan(gamma2_s) + gamma2_s - np.pi/2.0)*radius2*radius2
                projs = projs + projs_mount

        projv = projv/(1-poccupied)
        projs = projs/(1-poccupied)
        Overlapping = np.sqrt(tantv * tantv + tants * tants - 2 * tantv * tants * up) / (tantv + tants)
        projvs = projv + projs * Overlapping

        # partmax = np.max([projv, projs])
        # partmin = np.min([projv, projs])
        # # if (raa == 0) | (raa==180):
        # #     Overlapping = 0
        # projvs = partmax + partmin * Overlapping

        ### 地形导致的平坦下垫面的光照和阴影情况
        pPlaneV = np.exp(-projv)*(1-poccupied)
        pPlaneS = np.exp(-projs)*(1-poccupied)
        pPlaneV_sunlit = np.exp(-projvs)*(1-poccupied)
        pPlaneV_shaded = pPlaneV - pPlaneV_sunlit

        ### 森林的属性
        forestshape = self.forestshape
        lai = forestshape[0]
        std = forestshape[1]
        hspot = forestshape[4]
        hcr = forestshape[3]
        rcr = forestshape[2]
        Pcomtemp = proportion_bidirectional_crown_one(lai, std, hspot, hcr, rcr, np.asarray([vza]),
                                                      np.asarray([sza]), np.asarray([raa]))
        ### 只有光照山坡的光照植被和光照土壤才是光照
        ### 在阴影山坡的所有植被和土壤全是阴影
        pPlaneV_veg_sunlit = pPlaneV_sunlit*Pcomtemp[2][0]
        pPlaneV_veg_shaded = pPlaneV_sunlit*Pcomtemp[3][0] + pPlaneV_shaded*Pcomtemp[2][0] + pPlaneV_shaded*Pcomtemp[3][0]
        pPlaneV_soil_sunlit = pPlaneV_sunlit*Pcomtemp[0][0]
        pPlaneV_soil_shaded = pPlaneV_sunlit*Pcomtemp[1][0] + pPlaneV_shaded*Pcomtemp[0][0] + pPlaneV_shaded*Pcomtemp[1][0]
        ### 测试集合
        alla1 = pPlaneV_veg_sunlit+pPlaneV_veg_shaded+pPlaneV_soil_sunlit+pPlaneV_soil_shaded
        ####----------------------------------------------
        ### Mountain surface
        ###-----------------------------------------------

        pMountV = 0
        pMountV_sunlit = 0
        pMComV = 0

        Pcom = 0
        pMount_veg = 0
        pMount_soil = 0
        pMount_veg_sunlit = 0
        pMount_soil_sunlit = 0
        for kshape1 in range(n_shape):
            shape1 = shapes[kshape1]

            radius1 = shape1[1]
            height1 = shape1[0]
            density1 = shape1[2]
            height1r = height1
            dheight = height1 / n_part
            projv = 0
            projs = 0

            ### 及圆锥的角度情况
            alpha1 = np.arctan(radius1 / height1)
            ### 最高点在观测方向倾斜投影的距离
            L1_v = height1 * tantv
            if L1_v < radius1: L1_v = radius1 * 1.0
            theta1_v = np.arctan(L1_v * np.tan(alpha1) / radius1)
            gamma1_v = np.arcsin(radius1 / L1_v)
            ### 最高点在太阳方向倾斜投影的距离，及圆锥的角度情况
            L1_s = height1 * tants
            if L1_s < radius1: L1_s = radius1 * 1.0
            theta1_s = np.arctan(L1_s * np.tan(alpha1) / radius1)
            gamma1_s = np.arcsin(radius1 / L1_s)

            weight = 0
            ### 这里计算的投影只是圆球方向出来的部分
            ### 首先不管哪个角度都会有最少的等效圆的投影面积；
            ### 在大于圆锥角度后，会有侵占平坦地表的圆锥面积，该面积才是导致圆锥被挡住的原因，因此只计算该部分
            for kh in range(n_part):
                height_temp = dheight * (kh + 0.5)
                dh = shapes[:, 0] - height_temp
                alpha2 = np.arctan(shapes[:,1] / shapes[:,0])
                r1 = height_temp / height1 * radius1
                w = (r1*np.pi*dheight+r1*2+dheight)  ### 高度并不是等效的，这里是不同高度的权重

                L2_v = dh * tantv
                ind =  L2_v < shapes[:,1]
                L2_v[ind] = shapes[ind,1]
                theta2_v = np.arctan(L2_v * np.tan(alpha2) / shapes[:,1])
                gamma2_v = np.arcsin(shapes[:,1] / L2_v)
                projv_mount = (shapes[:,2] * (1.0/np.tan(gamma2_v) + gamma2_v - np.pi / 2.0) * shapes[:,1] * shapes[:,1])
                ind = (dh < 0) + (L2_v <= shapes[:,1])
                projv_mount[ind] = 0
                projv = projv + np.sum(projv_mount)*(w)

                L2_s = dh * tants
                ind =  (L2_s < shapes[:,1])
                L2_s[ind] = shapes[ind,1]
                theta2_s = np.arctan(L2_s * np.tan(alpha2) / shapes[:,1])
                gamma2_s = np.arcsin(shapes[:,1] / L2_s)
                projs_mount = ( shapes[:,2] * (1.0/np.tan(gamma2_s) + gamma2_s - np.pi / 2.0) * shapes[:,1] * shapes[:,1])
                ind = (dh < 0) + (L2_s <= (shapes[:,1]))
                projs_mount[ind] = 0
                projs = projs + np.sum(projs_mount)*w

                weight = weight + (w)

            projv = projv / n_part/weight/(1-poccupied)  ### 因为没有考虑本身圆形情况，延伸出来的部分
            projs = projs / n_part/weight/(1-poccupied)
            Overlapping = np.sqrt(tantv * tantv + tants * tants - 2 * tantv * tants * up) / (tantv + tants)
            projvs = projv + projs * Overlapping

            # partmax = np.max([projv, projs])
            # partmin = np.min([projv, projs])
            # # if (raa == 0) | (raa==180):
            # #     Overlapping = 0
            # projvs = partmax + partmin * Overlapping

            ### 计算山地的平均光照和可视比例
            gapv_mount = np.exp(-projv)
            gapvs_mount = np.exp(-projvs)

            slope = alpha1 * 180 / np.pi
            pMountV_temp = 0
            pMountV_sunlit_temp = 0

            if L1_v <= radius1:
                ### 如果没有投影出来，那可视就都是本圆
                pMountV_temp =  density1 * gapv_mount * (( np.pi) * radius1 * radius1)
                if sza > slope:
                    ### 如果此时光照角度大于圆锥倾斜，则光照部分类似于椭球的方式计算
                    cosphi = uv * ui + sv * si * up
                    pMountV_sunlit_temp = density1 * gapvs_mount *(( np.pi) * radius1 * radius1) * (1 + cosphi) * 0.5
                else:
                    ### 如果此时光照角度小于圆锥倾斜，则所有可视，所有光照
                    gapvs_mount = 1.0
                    pMountV_sunlit_temp = density1 * gapv_mount * ((np.pi) * radius1 * radius1)
            else:
                ### 投影出本圆，那可视就都是本圆+投影出来的部分，是Li圆锥的计算方法
                pMountV_temp = density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1)
                if sza > slope:
                    ### 如果此时光照角度大于圆锥倾斜，则光照是类椭球计算方式
                    cosphi = uv * ui + sv * si * up
                    pMountV_sunlit_temp = density1 * gapvs_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1) * (1 + cosphi) * 0.5
                else:
                    ### 如果此时光照角度大于圆锥倾斜，则光照是完全的光照
                    pMountV_sunlit_temp = density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1)
            pMountV = pMountV + pMountV_temp
            pMountV_sunlit = pMountV_sunlit + pMountV_sunlit_temp

            ### 这里不进行方位角的积分，选择了采用方位角的中心进行计算
            forestshape = self.forestshape
            lai = forestshape[0]
            std = forestshape[1]
            hspot = forestshape[4]
            hcr = forestshape[3]
            rcr = forestshape[2]


            ### 计算每个离散360方向上的角度贡献
            n = 360
            lza = np.zeros(n) + slope
            laa = np.arange(0, n, 1)
            vzatemp = np.zeros(n) + vza
            vaatemp = np.zeros(n) + vaa
            vlatemp = np.abs(vaa - laa)
            vsatemp = np.abs(vaatemp - saa)
            szatemp = np.zeros(n) + sza
            stdtemp = np.zeros(n) + std
            hspotemp = np.zeros(n) + hspot
            hcrtemp = np.zeros(n) + hcr

            uii = np.cos(lza * rd)
            uvv = np.cos(vzatemp * rd)
            sii = np.sin(lza * rd)
            svv = np.sin(vzatemp * rd)
            upp = np.cos(vlatemp * rd)
            # tantv = np.tan(vzatemp * rd)
            tantl = np.tan(vlatemp * rd)
            cosang = uvv * uii + svv * sii * upp
            ang = np.arccos(cosang) / rd

            ### 当大于0的时候，才是能看到的情况，否则是看不到的，也是没有贡献的
            ind = cosang >= 0
            weight_all = np.sum(cosang[ind])
            weight = cosang[ind]/weight_all
            vzatemp = vzatemp[ind]
            vsatemp = vsatemp[ind]
            szatemp = szatemp[ind]
            hspotemp = hspotemp[ind]
            # stdtemp = stdtemp[ind]
            # hcrtemp = hcrtemp[ind]

            stdtemp = stdtemp[ind]*uii[ind]
            hcrtemp = hcrtemp[ind]*uii[ind]

            ### 在这些有效角度下的贡献，其贡献的权重是其方向上投影有关系的

            Pcomtemp = proportion_bidirectional_crown_one(lai, stdtemp, hspotemp, hcrtemp, rcr, np.asarray(vzatemp), np.asarray(szatemp), np.asarray(vsatemp))

            Pcom = np.sum(Pcomtemp * weight,axis=1)
            pVeg = Pcom[2]+Pcom[3]
            pSoil = Pcom[0] + Pcom[1]
            pVeg_sunlit = Pcom[2]
            pSoil_sunlit = Pcom[0]
            pMount_veg = pMount_veg + pMountV_temp * pVeg                                ### 植被的贡献
            pMount_soil = pMount_soil + pMountV_temp * pSoil                             ### 土壤的贡献
            pMount_soil_sunlit = pMount_soil_sunlit + pMountV_sunlit_temp * pSoil_sunlit ### 阳坡 * 光照土壤
            pMount_veg_sunlit = pMount_veg_sunlit + pMountV_sunlit_temp * pVeg_sunlit    ### 阳坡 * 光照植被




        ### 进行了坡地与平坦地表间的归一化操作，temp是坡地的变化比例
        pMountV_sunlit_fraction = pMountV_sunlit/pMountV
        pMountVnew = 1 - pPlaneV
        temp = pMountVnew/pMountV
        pMountV_sunlit = temp * pMountV_sunlit
        pMountV_shaded = pMountVnew - pMountV_sunlit

        pMountV_veg_sunlit = pMount_veg_sunlit * temp
        pMountV_veg_shaded = pMount_veg * temp - pMountV_veg_sunlit
        pMountV_soil_sunlit = pMount_soil_sunlit * temp
        pMountV_soil_shaded = pMount_soil * temp - pMountV_soil_sunlit

        ### 测试集，加成1和加成2
        alla2 = pMountV_veg_sunlit+pMountV_veg_shaded+pMountV_soil_shaded+pMountV_soil_sunlit
        alla = pPlaneV_veg_sunlit + pMountV_veg_sunlit+ pPlaneV_veg_shaded + pMountV_veg_shaded+ pPlaneV_soil_sunlit + pMountV_soil_sunlit+ pPlaneV_soil_shaded + pMountV_soil_shaded
        # print(pPlaneV_sunlit,pPlaneV_shaded,pMountV_sunlit,pMountV_shaded)
        # print(Pcom,np.sum(Pcom))
        if ifP ==1:
            return pPlaneV, pMountV
        elif ifP == 2:
            return pPlaneV_soil_sunlit + pMountV_soil_sunlit, pPlaneV_soil_shaded + pMountV_soil_shaded,\
                   pPlaneV_veg_sunlit + pMountV_veg_sunlit, pPlaneV_veg_shaded + pMountV_veg_shaded

        elif ifP ==3:
            return 0, 0
        else:
            return pPlaneV+ pMountV


    def calculate_component_scatter_emissivity(self, kangle,ifP = 0):
        shapes = self.shapes
        Es = self.Es
        n_shape = self.n_shape
        n_dem = self.n_dim
        n_part = self.n_part
        vza = self.vza_[kangle]
        vaa = self.vaa_[kangle]
        sza = self.sza
        saa = self.saa
        raa = vaa - saa

        if raa > 180:raa = 360-raa


        rd = np.pi / 180.0
        ui = np.cos(sza * rd)
        uv = np.cos(vza * rd)
        si = np.sin(sza * rd)
        sv = np.sin(vza * rd)
        up = np.cos(raa * rd)
        tantv = np.tan(vza * rd)
        tants = np.tan(sza * rd)

        ####----------------------------------------------
        ### PLANE surface
        ###-----------------------------------------------
        height1r = 0
        projv = 0
        projs = 0
        poccupied = 0
        projvs = 0
        for kshape2 in range(n_shape):
            shape2 = shapes[kshape2]
            radius2 = shape2[1]
            height2 = shape2[0]
            density2 = shape2[2]
            poccupied = poccupied + density2 * np.pi * (radius2*radius2)
            dh = height2 - height1r
            alpha2 = np.arctan(radius2/height2)

            L2_v = dh * tantv
            if L2_v < radius2: L2_v = radius2*1.0
            theta2_v = np.arctan(L2_v*np.tan(alpha2)/radius2)
            gamma2_v = np.arcsin(radius2/L2_v)
            if dh <= 0 or L2_v <= radius2: continue
            projv_mount = density2 * (np.cos(gamma2_v) + gamma2_v + np.pi/2.0)*radius2*radius2
            projv = projv + projv_mount


            L2_s = dh * tants
            if L2_s < radius2: L2_s = radius2*1.0
            theta2_s = np.arctan(L2_s*np.tan(alpha2)/radius2)
            gamma2_s = np.arcsin(radius2/L2_s)
            if dh <= 0 or L2_s <= radius2: continue
            projs_mount = density2 * (1.0/np.tan(gamma2_s) + gamma2_s + np.pi/2.0)*radius2*radius2
            projs = projs + projs_mount

        Overlapping = np.sqrt(tantv * tantv + tants * tants - 2 * tantv * tants * up) / (tantv + tants)
        projvs = projv + projs * Overlapping
        pPlaneV = np.exp(-projv)*(1-poccupied)
        pPlaneS = np.exp(-projs)*(1-poccupied)
        pPlaneV_sunlit = np.exp(-projvs)*(1-poccupied)
        pPlaneV_shaded = pPlaneV - pPlaneV_sunlit

        ### 半球投影的归一化
        if2 = 1.0
        ### 上半球的方向
        n_hza0 = 90
        n_haa0 = 60
        hza0 = np.linspace(0, 89, n_hza0) # hemisphere space
        haa0 = np.linspace(0, 360, n_haa0)
        hza,haa = np.meshgrid(hza0,haa0)
        hza = np.reshape(hza,-1)  ### to linear
        haa = np.reshape(haa,-1)  ### to linear

        ### step of hza in rad unit rather zhan degree unit
        dangle = np.pi/2.0/n_hza0 # 角度的步长
        tantemp = np.tan(hza * rd)   ### 半球的正切
        # fweiplus = np.sin(hza*rd)*dangle *np.cos(hza*rd) # 角度的归一化方式, 一般是0.5
        fweiplus = np.sin(hza * rd) * dangle  ##### 一般是1.0
        fweiplusSum = np.sum(if2*fweiplus)  ### 权重累计，如果没有meshgrid，应该是0.5，现在总和是0.5*n_haa0

        alpha2 = np.arctan(shapes[:, 1] / shapes[:, 0])

        L2_v = shapes[:,0] * tantv
        ind = L2_v < shapes[:, 1]
        L2_v[ind] = shapes[ind, 1]
        theta2_v = np.arctan(L2_v * np.tan(alpha2) / shapes[:, 1])
        gamma2_v = np.arcsin(shapes[:, 1] / L2_v)
        projv_mount = (shapes[:, 2] * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * shapes[:, 1] * shapes[:, 1])
        tempSv = np.sum(projv_mount)
        f0 = np.exp(-tempSv)
        i0v = 1-f0


        eu = 0.0
        ed = 0.0
        au = 0.5
        ad = 0.5

        for kshape1 in range(n_shape):
            shape1 = shapes[kshape1]

            radius1 = shape1[1]
            height1 = shape1[0]
            density1 = shape1[2]
            height1r = height1
            dheight = height1 / n_part
            projv = 0
            projs = 0

            alpha1 = np.arctan(radius1 / height1)
            L1_v = height1 * tantv
            if L1_v < radius1: L1_v = radius1 * 1.0
            theta1_v = np.arctan(L1_v * np.tan(alpha1) / radius1)
            gamma1_v = np.arcsin(radius1 / L1_v)


            for kh in range(n_part):
                heighttemp = dheight * (kh + 0.5)
                tempdS = ( dheight * density1)* tantv *(heighttemp/height1)*(2*radius1)

                projStemp = 0
                proj_roof = 0
                tempSv = 0
                for kshape2 in range(n_shape):
                    shape2 = shapes[kshape2]
                    radius2 = shape2[1]
                    height2 = shape2[0]
                    density2 = shape2[2]

                    if (height2 - heighttemp > 0):

                        dh = height2 - height1r
                        alpha2 = np.arctan(radius2 / height2)
                        L2_v = dh * tantemp
                        ind = L2_v < radius2
                        L2_v[ind] = radius2
                        theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                        gamma2_v = np.arcsin(radius2 / L2_v)
                        projtemp_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2
                        ind = (dh<=0)+(L2_v<=radius2)
                        projtemp_mount[ind] = 0
                        projStemp = projStemp + projtemp_mount

                        L2_v = dh * tantv
                        if L2_v < radius2: L2_v = radius2 * 1.0
                        theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                        gamma2_v = np.arcsin(radius2 / L2_v)
                        if dh <= 0 or L2_v <= radius2: continue
                        projv_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2
                        tempSv = tempSv + projv_mount

                ftempv = np.exp(-(tempSv))
                ftemp = np.exp(-(projStemp))
                eu = eu + np.sum(if2 * ftemp * fweiplus / fweiplusSum) * tempdS * ftempv

                projStemp = 0
                tempSv = 0
                for kshape2 in range(n_shape):
                    shape2 = shapes[kshape2]
                    radius2 = shape2[1]
                    height2 = shape2[0]
                    density2 = shape2[2]

                    if (height2 - heighttemp < 0):

                        dh = height2
                        alpha2 = np.arctan(radius2 / height2)
                        L2_v = dh * tantemp
                        ind = L2_v < radius2
                        L2_v[ind] = radius2
                        theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                        gamma2_v = np.arcsin(radius2 / L2_v)
                        # if dh <= 0 or L2_v <= radius2: continue
                        projtemp_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2

                        ind = (dh <= 0) + (L2_v <= radius2)
                        projtemp_mount[ind] = 0

                        projStemp = projStemp + projtemp_mount
                    if (height2 - heighttemp > 0):

                        dh = heighttemp
                        alpha2 = np.arctan(radius2 / heighttemp)
                        L2_v = dh * tantemp
                        ind = L2_v < radius2
                        L2_v[ind] = radius2
                        theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                        gamma2_v = np.arcsin(radius2 / L2_v)
                        # if dh <= 0 or L2_v <= radius2: continue
                        projtemp_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2
                        ind = (dh <= 0) + (L2_v <= radius2)
                        projtemp_mount[ind] = 0
                        projStemp = projStemp + projtemp_mount


                ftemp = np.exp(-(projStemp))
                ed = ed + np.sum(if2 * ftemp * fweiplus/ fweiplusSum) * tempdS * ftempv

        eu = eu / i0v*0.5
        ed = ed / i0v*0.5
        p = 1 - eu - ed
        if i0v == 0:
            p = 0
            eu = 0
            ed = 0

        ### 到达街道，然后街道的比例，然后街道离开地面到达天空的概率
        heighttemp = 0
        ws = 0
        C = 1.05
        projStemp = 0
        poccupied = 0
        tempSv = 0
        for kshape2 in range(n_shape):
            shape2 = shapes[kshape2]
            radius2 = shape2[1]
            height2 = shape2[0]
            density2 = shape2[2]
            poccupied = poccupied + density2 * np.pi * (radius2*radius2)
            if (height2 - heighttemp > 0):

                dh = height2 - height1r
                alpha2 = np.arctan(radius2 / height2)
                L2_v = dh * tantemp
                ind = L2_v < radius2
                L2_v[ind] = radius2
                theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                gamma2_v = np.arcsin(radius2 / L2_v)
                if dh <= 0 or L2_v <= radius2: continue
                projtemp_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2
                projStemp = projStemp + projtemp_mount

                L2_v = dh * tantv
                if L2_v < radius2: L2_v = radius2 * 1.0
                theta2_v = np.arctan(L2_v * np.tan(alpha2) / radius2)
                gamma2_v = np.arcsin(radius2 / L2_v)
                if dh <= 0 or L2_v <= radius2: continue
                projv_mount = density2 * (1.0/np.tan(gamma2_v) + gamma2_v + np.pi / 2.0) * radius2 * radius2
                tempSv = tempSv + projv_mount

        ftemp = np.exp(-(projStemp))
        ftempv = np.exp(-(tempSv)) ### 到达该体元的概率
        tempds = 1-poccupied
        ws = (1-np.sum(if2*(ftemp) * fweiplus/ fweiplusSum)) * ftempv * tempds



        eww = Es * (1 - Es) * p * i0v
        ews = Es * (1-Es) * ws
        esw = Es *(1-Es) * ed * i0v

        vsraa = (180.0 - np.abs(raa)) / 180.0  # sunlit part becasue of raa
        vhraa = 1 - vsraa
        ## 场景中有多少光照组分

        area = 0
        fts = 1.0
        for kshape1 in range(n_shape):
            shape1 = shapes[kshape1]
            radius1 = shape1[1]
            height1 = shape1[0]
            density1 = shape1[2]
            height1r = height1
            L1_s = height1 * tants
            if L1_s < radius1: L2_v = radius1 * 1.0
            gamma1_s = np.arcsin(radius1 / L1_s)
            projs_mount = density1 * (1.0/np.tan(gamma1_s) + gamma1_s + np.pi / 2.0) * radius1 * radius1
            area = area + projs_mount
        if area > 0:
            fts = (1 - np.exp(-area)) / (area)


        emw = ews + eww
        emws = emw * fts*vsraa
        emwh = emw - emws
        ems = esw
        emss = ems * pPlaneS
        emsh = ems - emss

        if ifP == 1:
            return ems,emw
        elif ifP == 2:
            return emss, emsh, emws, emwh
        elif ifP == 3:
            return 0, 0, 0, 0
        else:
            return esw + ews + eww

    def set_strcutural_input(self,shapes):
        self.shapes = shapes
        self.n_shape, self.n_dim = np.shape(shapes)

    def set_structural_input(self, shapes):
        return self.set_strcutural_input(shapes)

    def set_thermal_input(self,Tss,Tsh,Tms,Tmh):
        self.Ts_sunlit = Tss
        self.Ts_shaded = Tsh
        self.Tm_sunlit = Tms
        self.Tm_shaded = Tmh


    def set_angular_input(self,vza_,vaa_,sza,saa):
        self.vaa_ = vaa_
        self.vza_ = vza_
        self.sza = sza
        self.saa = saa
        self.n_angle = np.size(vza_)

    def set_spectral_input(self,es,em):
        self.Es = es
        self.Em = em

    def calculate_effective_component_emissivity(self,ifP = 0):
        emissivity_ = []
        for kangle in range(self.n_angle):
            direct = self.calculate_component_direct_emissivity(kangle,ifP)
            # scatter = self.calculate_component_scatter_emissivity(kangle,ifP)
            scatter = 0
            direct = np.asarray(direct)
            scatter = np.asarray(scatter)
            # scatter[:] = 0
            # direct[:] = 0
            emissivity_.append(direct+scatter)

        self.canopy_effective_emissivity = np.asarray(emissivity_ )
        return np.asarray(emissivity_ )

