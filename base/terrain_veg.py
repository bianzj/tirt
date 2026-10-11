import numpy as np

from .crown import *
from .physicsF import projection_overlap
from .terrain import Terrain, _self_lit_fraction, _unit
###----------------------------------
### 这里考虑了山地本身的遮挡和子像元内植被的复合效应
### 1. determine the mountainous surface
### 2. determine the forest strcture
###----------------------------------



class Terrain_Veg:

    ### h, r, density
    shapes = np.asarray([[10,10,0.003],[10,5,0.001],[10,15,0.001]])
    ### lai, stand, radi_horizontal, radi_vertical, hspot, trunk height (crown base)
    forestshape = np.asarray([1.0,0.001,3,3,0.2,3])

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
        Overlapping = projection_overlap(tantv, tants, up)
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
                                                      np.asarray([sza]), np.asarray([raa]), hcb=forestshape[5])
        ### 只有光照山坡的光照植被和光照土壤才是光照
        ### 在阴影山坡的所有植被和土壤全是阴影
        self._leaf_view_ground = Pcomtemp[2][0] + Pcomtemp[3][0]
        pPlaneV_veg_sunlit = pPlaneV_sunlit*Pcomtemp[2][0]
        pPlaneV_veg_shaded = pPlaneV_sunlit*Pcomtemp[3][0] + pPlaneV_shaded*Pcomtemp[2][0] + pPlaneV_shaded*Pcomtemp[3][0]
        pPlaneV_soil_sunlit = pPlaneV_sunlit*Pcomtemp[0][0]
        pPlaneV_soil_shaded = pPlaneV_sunlit*Pcomtemp[1][0] + pPlaneV_shaded*Pcomtemp[0][0] + pPlaneV_shaded*Pcomtemp[1][0]
        ### 森林与背景之间的散射（同平地树冠模型），在地形阴影中的部分全部归为阴影叶片
        scat = np.asarray(emissivity_scatter_analytical(lai, self.Em, Esoil, np.asarray([vza]), np.asarray([sza]))[0], dtype=float).reshape(4, -1)[:, 0]
        pPlaneV_scat_sunlit = pPlaneV_sunlit * scat[2]
        pPlaneV_scat_shaded = pPlaneV_sunlit * scat[3] + pPlaneV_shaded * (scat[2] + scat[3])
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
        pMount_scat = 0
        pMount_scat_sunlit = 0
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
                r1 = (1 - height_temp / height1) * radius1
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
            Overlapping = projection_overlap(tantv, tants, up)
            projvs = projv + projs * Overlapping

            # partmax = np.max([projv, projs])
            # partmin = np.min([projv, projs])
            # # if (raa == 0) | (raa==180):
            # #     Overlapping = 0
            # projvs = partmax + partmin * Overlapping

            ### 计算山地的平均光照和可视比例
            # exp(-gap_shift) is common to all shapes and cancels in the normalisation by temp; it avoids 0/0 as vza -> 90
            if kshape1 == 0 or projv < gap_shift:
                rescale = 1.0 if kshape1 == 0 else np.exp(projv - gap_shift)
                pMountV, pMountV_sunlit = pMountV * rescale, pMountV_sunlit * rescale
                pMount_veg, pMount_soil = pMount_veg * rescale, pMount_soil * rescale
                pMount_veg_sunlit, pMount_soil_sunlit = pMount_veg_sunlit * rescale, pMount_soil_sunlit * rescale
                pMount_scat, pMount_scat_sunlit = pMount_scat * rescale, pMount_scat_sunlit * rescale
                gap_shift = projv
            gapv_mount = np.exp(-(projv - gap_shift))
            gapvs_mount = np.exp(-(projvs - gap_shift))

            slope = alpha1 * 180 / np.pi
            pMountV_temp = 0
            pMountV_sunlit_temp = 0

            if L1_v <= radius1:
                ### 如果没有投影出来，那可视就都是本圆
                pMountV_temp =  density1 * gapv_mount * (( np.pi) * radius1 * radius1)
            else:
                ### 投影出本圆，那可视就都是本圆+投影出来的部分，是Li圆锥的计算方法
                pMountV_temp = density1 * gapv_mount * ((1.0/np.tan(gamma1_v) + gamma1_v + np.pi / 2.0) * radius1*radius1)
            # other cones' shadows only; self-shading is resolved per facet orientation below
            mutual = np.exp(-projs * Overlapping)
            beta1 = np.pi / 2 - alpha1
            v_dir = _unit(np.asarray(vza * rd), np.asarray(vaa * rd))
            s_dir = _unit(np.asarray(sza * rd), np.asarray(saa * rd))
            pMountV_sunlit_temp = pMountV_temp * mutual * _self_lit_fraction(beta1, v_dir, s_dir)
            pMountV = pMountV + pMountV_temp
            pMountV_sunlit = pMountV_sunlit + pMountV_sunlit_temp

            forestshape = self.forestshape
            lai = forestshape[0]
            std = forestshape[1]
            hspot = forestshape[4]
            hcr = forestshape[3]
            rcr = forestshape[2]

            ### 计算每个离散360方向上的角度贡献
            # facet normals: zenith = side slope beta1 = atan(h/r), azimuth laa
            laa = np.arange(0, 360, 1) * rd
            normals = _unit(np.full(laa.size, beta1), laa)
            cosang = normals @ v_dir
            coss = normals @ s_dir
            ind = cosang > 1e-9
            weight = cosang[ind] / np.sum(cosang[ind])
            normals, cosang, coss = normals[ind], cosang[ind], coss[ind]
            # facet-local geometry (cf. slope_veg/slope1); relative azimuth between the in-plane projections of v and s
            v_t = v_dir - cosang[:, None] * normals
            s_t = s_dir - coss[:, None] * normals
            cos_raa = np.sum(v_t * s_t, 1) / np.maximum(np.linalg.norm(v_t, axis=1) * np.linalg.norm(s_t, axis=1), 1e-12)
            vza_loc = np.minimum(np.arccos(np.clip(cosang, -1, 1)) / rd, 89.999)
            sza_loc = np.minimum(np.arccos(np.clip(coss, -1, 1)) / rd, 89.0)
            raa_loc = np.arccos(np.clip(cos_raa, -1, 1)) / rd
            # Boolean crowns standing on a facet: densities per facet area are lai*cos(beta), std*cos(beta); the path
            # factor is 1/cos(local zenith) (gap MC on tilted planes: 0.736 vs 0.738; old global-angle form 0.746/0.776 vs 0.603)
            cos_beta = np.cos(beta1)
            Pcomtemp = np.asarray(proportion_bidirectional_crown_one(lai * cos_beta, std * cos_beta, hspot, hcr, rcr, vza_loc, sza_loc, raa_loc, hcb=forestshape[5]))
            # facets turned away from the sun are entirely shaded
            lit = coss > 0
            Pcomtemp[1] = Pcomtemp[1] + Pcomtemp[0] * ~lit
            Pcomtemp[3] = Pcomtemp[3] + Pcomtemp[2] * ~lit
            Pcomtemp[0] = Pcomtemp[0] * lit
            Pcomtemp[2] = Pcomtemp[2] * lit

            Pcom = np.sum(Pcomtemp * weight,axis=1)
            pVeg = Pcom[2]+Pcom[3]
            pSoil = Pcom[0] + Pcom[1]
            pVeg_sunlit = Pcom[2]
            pSoil_sunlit = Pcom[0]
            pMount_veg = pMount_veg + pMountV_temp * pVeg                                ### 植被的贡献
            pMount_soil = pMount_soil + pMountV_temp * pSoil                             ### 土壤的贡献
            pMount_soil_sunlit = pMount_soil_sunlit + pMountV_temp * mutual * pSoil_sunlit
            pMount_veg_sunlit = pMount_veg_sunlit + pMountV_temp * mutual * pVeg_sunlit
            # forest-background scattering of the facet canopy in facet-local angles; unlit facets are all shaded
            scat = np.asarray(emissivity_scatter_analytical(lai * cos_beta, self.Em, Esoil, vza_loc, sza_loc)[0], dtype=float).reshape(4, -1)
            pMount_scat = pMount_scat + pMountV_temp * np.sum((scat[2] + scat[3]) * weight)
            pMount_scat_sunlit = pMount_scat_sunlit + pMountV_temp * mutual * np.sum(scat[2] * lit * weight)




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
        pMountV_scat_sunlit = pMount_scat_sunlit * temp
        pMountV_scat_shaded = pMount_scat * temp - pMountV_scat_sunlit
        # composition of the facet surface seen by the sensor, used as reflector in the scatter term
        self._leaf_view_facet = pMount_veg / (pMount_veg + pMount_soil)
        ### 测试集，加成1和加成2
        alla2 = pMountV_veg_sunlit+pMountV_veg_shaded+pMountV_soil_shaded+pMountV_soil_sunlit
        alla = pPlaneV_veg_sunlit + pMountV_veg_sunlit+ pPlaneV_veg_shaded + pMountV_veg_shaded+ pPlaneV_soil_sunlit + pMountV_soil_sunlit+ pPlaneV_soil_shaded + pMountV_soil_shaded
        # print(pPlaneV_sunlit,pPlaneV_shaded,pMountV_sunlit,pMountV_shaded)
        # print(Pcom,np.sum(Pcom))
        if ifP ==1:
            return pPlaneV, pMountVnew
        elif ifP == 2:
            return (pPlaneV_soil_sunlit + pMountV_soil_sunlit) * self.Es, (pPlaneV_soil_shaded + pMountV_soil_shaded) * self.Es,\
                   (pPlaneV_veg_sunlit + pMountV_veg_sunlit) * self.Em + pPlaneV_scat_sunlit + pMountV_scat_sunlit,\
                   (pPlaneV_veg_shaded + pMountV_veg_shaded) * self.Em + pPlaneV_scat_shaded + pMountV_scat_shaded

        elif ifP ==3:
            return 0, 0
        else:
            return pPlaneV+ pMountV


    def calculate_component_scatter_emissivity(self, kangle,ifP = 0):
        """Terrain-scale single scattering between canopy-covered facets and canopy-covered ground.

        Geometry (i0v, eu/ed/p, ws, path sunlit shares) is shared with the bare Terrain model. Every surface is a
        soil+crown mixture: it emits Es*G (soil) and Em*(1-G) (leaves), G = hemispherical canopy gap, and reflects with
        the composition the sensor sees on it. The crown kernel used here has no within-canopy scattering, so nothing
        is counted twice.
        """
        t = Terrain.scatter_terms(self, kangle)
        Es, Em = self.Es, self.Em
        lai, std, rcr, hcr, hspot, hcb = self.forestshape[:6]
        beta = np.mean(np.arctan(self.shapes[:, 0] / self.shapes[:, 1]))
        sza = self.sza

        # 4-point Gauss-Legendre in mu = cos(theta) on (0,1) for hemispherical (cosine-weighted) averages
        x, w = np.polynomial.legendre.leggauss(4)
        mu, w = 0.5 * (x + 1), 0.5 * w
        rel = np.arange(0, 360, 60.0)
        MU, REL = np.meshgrid(mu, np.minimum(rel, 360 - rel), indexing="ij")
        WT = (np.repeat(w * mu, rel.size) / rel.size).reshape(MU.shape) * 2
        zen = np.degrees(np.arccos(MU)).reshape(-1)

        def hemispherical(lai_x, std_x):
            P = np.asarray(proportion_bidirectional_crown_one(lai_x, std_x, hspot, hcr, rcr, zen, np.full(zen.size, min(sza, 89.0)), REL.reshape(-1), hcb=hcb))
            wt = WT.reshape(-1)
            gap = np.sum(wt * (P[0] + P[1]))
            soil_lit = np.sum(wt * P[0]) / max(np.sum(wt * (P[0] + P[1])), 1e-12)
            leaf_lit = np.sum(wt * P[2]) / max(np.sum(wt * (P[2] + P[3])), 1e-12)
            return gap, soil_lit, leaf_lit

        G_g, sig_s, sig_l = hemispherical(lai, std)
        G_f = hemispherical(lai * np.cos(beta), std * np.cos(beta))[0]

        rho_f = self._leaf_view_facet * (1 - Em) + (1 - self._leaf_view_facet) * (1 - Es)
        rho_g = self._leaf_view_ground * (1 - Em) + (1 - self._leaf_view_ground) * (1 - Es)
        x_ff = t["p"] * t["i0v"] * rho_f
        x_fg = t["ws"] * rho_g
        x_gf = t["ed"] * t["i0v"] * rho_f
        lit_f = x_ff * t["s_ww"] + x_fg * t["s_ws"]
        lit_g = x_gf * t["ground_sunlit"]

        soil = Es * (G_f * (x_ff + x_fg) + G_g * x_gf)
        soil_sunlit = Es * sig_s * (G_f * lit_f + G_g * lit_g)
        leaf = Em * ((1 - G_f) * (x_ff + x_fg) + (1 - G_g) * x_gf)
        leaf_sunlit = Em * sig_l * ((1 - G_f) * lit_f + (1 - G_g) * lit_g)

        if ifP == 1:
            return soil, leaf
        elif ifP == 2:
            return soil_sunlit, soil - soil_sunlit, leaf_sunlit, leaf - leaf_sunlit
        elif ifP == 3:
            return 0, 0, 0, 0
        else:
            return soil + leaf

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
            scatter = self.calculate_component_scatter_emissivity(kangle,ifP)
            direct = np.asarray(direct)
            scatter = np.asarray(scatter)
            # scatter[:] = 0
            # direct[:] = 0
            emissivity_.append(direct+scatter)

        self.canopy_effective_emissivity = np.asarray(emissivity_ )
        return np.asarray(emissivity_ )
