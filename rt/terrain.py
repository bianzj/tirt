import numpy as np
from base.util_plot import *
from base.util_data import *
from base.util_plot import *
from rt.physicsF import *

###----------------------------------
### 这里考虑了山地本身，相互遮挡
###----------------------------------


class Terrain:

    ### h, r, density
    shapes = np.asarray([[10,10,0.001],[10,10,0.001]])
    n_shape,n_dim = np.shape(shapes)



    Em = 0.975
    Es = 0.955
    Ev = 0.975

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
        n_shape = self.n_shape
        n_dem = self.n_dim
        n_part = self.n_part
        vza = self.vza_[kangle]
        vaa = self.vaa_[kangle]
        sza = self.sza
        saa = self.saa
        raa = vaa - saa

        if raa > 180: raa = 360- raa

        Es = self.Es
        Em = self.Em

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

        pPlaneV = np.exp(-projv)*(1-poccupied)
        pPlaneS = np.exp(-projs)*(1-poccupied)
        pPlaneV_sunlit = np.exp(-projvs)*(1-poccupied)
        pPlaneV_shaded = pPlaneV - pPlaneV_sunlit


        ####----------------------------------------------
        ### Mountain surface
        ###-----------------------------------------------

        pMountV = 0
        pMountV_sunlit = 0
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

            L1_s = height1 * tants
            if L1_s < radius1: L1_s = radius1 * 1.0
            theta1_s = np.arctan(L1_s * np.tan(alpha1) / radius1)
            gamma1_s = np.arcsin(radius1 / L1_s)


            weight = 0
            for kh in range(n_part):
                height_temp = dheight * (kh + 0.5)
                dh = shapes[:, 0] - height_temp
                alpha2 = np.arctan(shapes[:,1] / shapes[:,0])
                r1 = height_temp / height1 * radius1
                w = (r1*np.pi*dheight+r1*2+dheight)


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

            projv = projv / n_part/weight/(1-poccupied)
            projs = projs / n_part/weight/(1-poccupied)
            Overlapping = np.sqrt(tantv * tantv + tants * tants - 2 * tantv * tants * up) / (tantv + tants)
            projvs = projv + projs * Overlapping

            # partmax = np.max([projv, projs])
            # partmin = np.min([projv, projs])
            # # if (raa == 0) | (raa==180):
            # #     Overlapping = 0
            # projvs = partmax + partmin * Overlapping

            gapv_mount = np.exp(-projv)
            gapvs_mount = np.exp(-projvs)

            slope = alpha1 * 180 / np.pi
            if L1_v <= radius1:
                pMountV = pMountV + density1 * gapv_mount * (( np.pi) * radius1 * radius1)
                if sza > slope:
                    cosphi = uv * ui + sv * si * up
                    pMountV_sunlit = pMountV_sunlit + density1 *(( np.pi) * radius1 * radius1) * (1 + cosphi) * 0.5
                else:
                    pMountV_sunlit = pMountV_sunlit + density1 * gapvs_mount * ((np.pi) * radius1 * radius1)
            else:
                pMountV = pMountV + density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1)
                if sza > slope:
                    cosphi = uv * ui + sv * si * up
                    pMountV_sunlit = pMountV_sunlit + density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1) * (1 + cosphi) * 0.5
                else:
                    pMountV_sunlit = pMountV_sunlit + density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1)





        pMountV_sunlit_fraction = pMountV_sunlit/pMountV
        pMountV = 1 - pPlaneV
        pMountV_sunlit = pMountV * pMountV_sunlit_fraction
        pMountV_shaded = pMountV - pMountV_sunlit

        # print(pPlaneV_sunlit,pPlaneV_shaded,pMountV_sunlit,pMountV_shaded)
        if ifP ==1:
            return pPlaneV, pMountV
        elif ifP == 2:
            return pPlaneV_sunlit*Es,pPlaneV_shaded*Es,pMountV_sunlit*Em,pMountV_shaded*Em
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
            scatter = self.calculate_component_scatter_emissivity(kangle,ifP)
            # scatter = 0
            direct = np.asarray(direct)
            scatter = np.asarray(scatter)
            # scatter[:] = 0
            # direct[:] = 0
            emissivity_.append(direct+scatter)

        self.canopy_effective_emissivity = np.asarray(emissivity_ )
        return np.asarray(emissivity_ )

