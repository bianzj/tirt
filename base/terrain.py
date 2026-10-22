import numpy as np
from .physicsF import *

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
        raa = np.abs((vaa - saa + 180.0) % 360.0 - 180.0)

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
        """Single-bounce terrain/ground interreflection with each emitter's own emissivity."""
        from .terrain_exchange import TerrainExchange
        t = TerrainExchange.scatter_terms(self, kangle)
        Es, Em = self.Es, self.Em
        mount_via_mount = Em * (1 - Em) * t["p"] * t["i0v"]
        mount_via_ground = Em * (1 - Es) * t["ws"]
        ground_via_mount = Es * (1 - Em) * t["ed"] * t["i0v"]
        emw = mount_via_mount + mount_via_ground
        emws = mount_via_mount * t["s_ww"] + mount_via_ground * t["s_ws"]
        emwh = emw - emws
        ems = ground_via_mount
        emss = ems * t["ground_sunlit"]
        emsh = ems - emss

        if ifP == 1:
            return ems,emw
        elif ifP == 2:
            return emss, emsh, emws, emwh
        elif ifP == 3:
            return 0, 0, 0, 0
        else:
            return ems + emw

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
            # scatter = 0
            direct = np.asarray(direct)
            scatter = np.asarray(scatter)
            # scatter[:] = 0
            # direct[:] = 0
            emissivity_.append(direct+scatter)

        self.canopy_effective_emissivity = np.asarray(emissivity_ )
        return np.asarray(emissivity_ )
