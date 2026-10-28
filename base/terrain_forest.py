"""Crown-covered terrain geometry and exchange restored from GitHub 86c3829."""
import numpy as np

from .emissivity import canopy_scatter
from .terrain_exchange import TerrainExchange as Terrain, _self_lit_fraction, _unit, projection_overlap
from .terrain_veg import Terrain_Veg
from .thermal_crown import component_fractions, geometry_mode
from .vegetation import vegetation_proportions

_HEMISPHERICAL_CACHE = {}


class TerrainForest(Terrain_Veg):
    forestshape = np.asarray([1., .001, 3., 3., .2, 3.])
    thermal_crown_geometry = 'legacy'
    thermal_hotspot_inner = None
    thermal_hotspot_outer = None
    thermal_forest_transport = 'legacy'
    # Optional canopy closures (canopy_closure.install): fractions(kind, lai, hspot, shape, vza, vaa, sza, saa,
    # row_azimuth, base) -> 4 components and scatter(kind, lai, el, es, vza, sza, density, vertical, horizontal, ...).
    canopy_fraction_model = None
    canopy_scatter_model = None
    _collect_optical_geometry = False
    # Coarse quadrature: 3 height bands and 10 facet azimuths change the 100 m hill scene by at most
    # 0.08 K (12 azimuths: 0.03 K) while running about 4x faster.
    n_part = 3
    facet_azimuths = 10

    def set_radiative_transport(self,mode='legacy'):
        if mode not in ('legacy','reciprocal_exchange'):
            raise ValueError('Unknown thermal forest transport: '+str(mode))
        if mode=='reciprocal_exchange' and self.thermal_crown_geometry!='clumping_dual':
            raise ValueError('Reciprocal forest transport requires clumping_dual.')
        self.thermal_forest_transport=mode

    def set_structural_input(self, shapes):
        array = np.array(shapes, dtype=float, copy=True)
        if (array.ndim != 2 or array.shape[1] != 3 or len(array) == 0
                or not np.isfinite(array).all() or np.any(array[:, :2] <= 0)
                or np.any(array[:, 2] < 0)):
            raise ValueError('Terrain requires finite positive heights/radii and nonnegative densities.')
        if np.sum(array[:, 2]*np.pi*array[:, 1]**2) >= 1:
            raise ValueError('Cone base fractions must sum to less than one.')
        self.shapes = array
        self.n_shape, self.n_dim = array.shape

    def set_hotspot_geometry(self, mode='legacy', inner=None, outer=None):
        geometry_mode(dict(thermal_crown_geometry=mode))
        if self.thermal_forest_transport=='reciprocal_exchange' and mode!='clumping_dual':
            raise ValueError('Reciprocal forest transport requires clumping_dual.')
        self.thermal_crown_geometry = mode
        self.thermal_hotspot_inner = inner
        self.thermal_hotspot_outer = outer

    def _forest_proportions(self, lai, density, hspot, height, width, view, solar, relative, hcb=0., azimuths=None):
        kind = self.vegetation_model
        if self.canopy_fraction_model is not None and kind in ('hom', 'row'):
            view_azimuth, solar_azimuth, row_azimuth = self._canopy_azimuths(relative, azimuths)
            view, solar, view_azimuth, solar_azimuth, row_azimuth = (a.ravel() for a in np.broadcast_arrays(
                np.atleast_1d(view).astype(float), np.atleast_1d(solar).astype(float),
                np.atleast_1d(view_azimuth).astype(float), np.atleast_1d(solar_azimuth).astype(float),
                np.atleast_1d(row_azimuth).astype(float)))
            return np.asarray([self.canopy_fraction_model(kind, float(lai), float(hspot), self.row_shape[:3], v, a, s, sa, r)
                               for v, a, s, sa, r in zip(view, view_azimuth, solar, solar_azimuth, row_azimuth)]).T
        if kind in ('hom', 'row'):
            row_width, row_blank, row_height, row_azimuth = self.row_shape
            if kind == 'row' and azimuths is not None:
                view_azimuth, solar_azimuth, row_azimuth = azimuths
            else:
                view_azimuth, solar_azimuth = relative, 0.
            return vegetation_proportions(kind, lai, hspot, view, view_azimuth, solar, solar_azimuth,
                                          row_width=row_width, row_blank=row_blank, row_height=row_height,
                                          row_azimuth=row_azimuth, row_geometry=self.row_geometry)
        if self.canopy_fraction_model is not None and kind == 'crown':
            view, solar, relative = np.broadcast_arrays(np.atleast_1d(view), np.atleast_1d(solar), np.atleast_1d(relative))
            return np.asarray([self.canopy_fraction_model('crown', lai, hspot, (density, height, width), v, r, s, 0., 0., hcb)
                               for v, s, r in zip(view.ravel(), solar.ravel(), relative.ravel())]).T
        if self.thermal_crown_geometry == 'legacy':
            return vegetation_proportions('crown', lai, hspot, view, relative, solar, 0.,
                                          density, height, width)
        outer=self.thermal_hotspot_outer
        if self.thermal_forest_transport=='reciprocal_exchange' and outer is None:
            from .forest_exchange import outer_width
            outer=outer_width(width,height,hcb)
        settings = dict(lai=lai, stand_density=density, hspot=hspot, height=height, width=width,
                        thermal_crown_geometry=self.thermal_crown_geometry,
                        thermal_hotspot_inner=self.thermal_hotspot_inner,
                        thermal_hotspot_outer=outer, crown_base=hcb)
        return component_fractions(settings, solar, view, relative)

    def _canopy_azimuths(self, relative, azimuths):
        """(view azimuth, solar azimuth, row azimuth) in the frame of the evaluated surface."""
        if self.vegetation_model == 'row' and azimuths is not None:
            return azimuths
        return relative, 0., self.row_shape[3]

    def _canopy_scatter(self, lai, vza, sza, std, hcr, rcr, relative, azimuths=None):
        if self.canopy_scatter_model is None:
            return canopy_scatter(self.vegetation_model, lai, self.Em, self.Es, vza, sza, std, hcr, rcr)
        view_azimuth, solar_azimuth, row_azimuth = self._canopy_azimuths(relative, azimuths)
        if np.ndim(row_azimuth) and np.size(row_azimuth) > 1:
            return np.column_stack([self.canopy_scatter_model(
                self.vegetation_model, lai, self.Em, self.Es, v, s, std, hcr, rcr, a, sa, (*self.row_shape[:3], float(r)))
                for v, s, a, sa, r in zip(*np.broadcast_arrays(vza, sza, view_azimuth, solar_azimuth, row_azimuth))])
        return self.canopy_scatter_model(self.vegetation_model, lai, self.Em, self.Es, vza, sza, std, hcr, rcr,
                                         view_azimuth, solar_azimuth, (*self.row_shape[:3], float(np.ravel(row_azimuth)[0])))

    def calculate_effective_component_emissivity(self, ifP=0):
        if ifP not in (0, 1, 2, 3):
            raise ValueError('Component mode must be 0, 1, 2 or 3.')
        values = []
        for angle in range(self.n_angle):
            direct = np.asarray(self.calculate_component_direct_emissivity(angle, 2), dtype=float)
            scatter = np.asarray(self.calculate_component_scatter_emissivity(angle, 2), dtype=float)
            budget = max(0., 1.-direct.sum())
            if scatter.sum() > budget:
                scatter *= budget / scatter.sum()
            components = direct + scatter
            if ifP == 0:
                values.append(components.sum())
            elif ifP == 1:
                values.append([components[:2].sum(), components[2:].sum()])
            elif ifP == 3:
                values.append(np.zeros(4))
            else:
                values.append(components)
        self.canopy_effective_emissivity = np.asarray(values)
        return self.canopy_effective_emissivity

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
        ### ȷ�� ɽ�صĽṹ�����
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

        ### ���ε��µ�ƽ̹�µ���Ĺ��պ���Ӱ���
        pPlaneV = np.exp(-projv)*(1-poccupied)
        pPlaneS = np.exp(-projs)*(1-poccupied)
        pPlaneV_sunlit = np.exp(-projvs)*(1-poccupied)
        pPlaneV_shaded = pPlaneV - pPlaneV_sunlit

        ### ɭ�ֵ�����
        forestshape = self.forestshape
        lai = forestshape[0]
        std = forestshape[1]
        hspot = forestshape[4]
        hcr = forestshape[3]
        rcr = forestshape[2]
        Pcomtemp = self._forest_proportions(lai, std, hspot, hcr, rcr, np.asarray([vza]),
                                                      np.asarray([sza]), np.asarray([raa]), hcb=forestshape[5],
                                                      azimuths=(np.asarray([vaa]), np.asarray([saa]), self.row_shape[3]))
        ### ֻ�й���ɽ�µĹ���ֲ���͹����������ǹ���
        ### ����Ӱɽ�µ�����ֲ��������ȫ����Ӱ
        self._leaf_view_ground = Pcomtemp[2][0] + Pcomtemp[3][0]
        self._reciprocal_ground_states=pPlaneV*np.array([1-self._leaf_view_ground,self._leaf_view_ground])
        pPlaneV_veg_sunlit = pPlaneV_sunlit*Pcomtemp[2][0]
        pPlaneV_veg_shaded = pPlaneV_sunlit*Pcomtemp[3][0] + pPlaneV_shaded*Pcomtemp[2][0] + pPlaneV_shaded*Pcomtemp[3][0]
        pPlaneV_soil_sunlit = pPlaneV_sunlit*Pcomtemp[0][0]
        pPlaneV_soil_shaded = pPlaneV_sunlit*Pcomtemp[1][0] + pPlaneV_shaded*Pcomtemp[0][0] + pPlaneV_shaded*Pcomtemp[1][0]
        if self._collect_optical_geometry:
            from .vnir import spherical_coefficients
            _,ko,sob,sof=spherical_coefficients(float(sza),float(vza),float(raa))
            optical_flat=np.array([sob/ko*pPlaneV_veg_sunlit,sof/ko*pPlaneV_veg_sunlit,
                                   pPlaneV_soil_sunlit])
            optical_slope=np.zeros(3)
            # visible canopy parts: [sunlit view weight x incidence, view weight, on facet, vza, sza, raa,
            # view azimuth, solar azimuth, row azimuth] in the frame of the part
            optical_parts=[[float(pPlaneV_sunlit),float(pPlaneV),0.,float(vza),float(sza),float(raa),
                            float(vaa),float(saa),float(self.row_shape[3])]]
        ### ɭ���뱳��֮���ɢ�䣨ͬƽ������ģ�ͣ����ڵ�����Ӱ�еĲ���ȫ����Ϊ��ӰҶƬ
        if self.thermal_forest_transport=='reciprocal_exchange':
            scat=np.zeros(4)
        else:
            scat = self._canopy_scatter(lai, np.asarray([vza]), np.asarray([sza]), std, hcr, rcr, np.asarray([raa]),
                                        (np.asarray([vaa]), np.asarray([saa]), self.row_shape[3]))[:, 0]
        pPlaneV_scat_sunlit = pPlaneV_sunlit * scat[2]
        pPlaneV_scat_shaded = pPlaneV_sunlit * scat[3] + pPlaneV_shaded * (scat[2] + scat[3])
        pPlaneV_scat_soil_sunlit = pPlaneV_sunlit * scat[0]
        pPlaneV_scat_soil_shaded = pPlaneV_sunlit * scat[1] + pPlaneV_shaded * (scat[0] + scat[1])
        ### ���Լ���
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
        pMount_scat_soil = 0
        pMount_scat_soil_sunlit = 0
        for kshape1 in range(n_shape):
            shape1 = shapes[kshape1]

            radius1 = shape1[1]
            height1 = shape1[0]
            density1 = shape1[2]
            height1r = height1
            dheight = height1 / n_part
            projv = 0
            projs = 0

            ### ��Բ׶�ĽǶ����
            alpha1 = np.arctan(radius1 / height1)
            ### ��ߵ��ڹ۲ⷽ����бͶӰ�ľ���
            L1_v = height1 * tantv
            if L1_v < radius1: L1_v = radius1 * 1.0
            theta1_v = np.arctan(L1_v * np.tan(alpha1) / radius1)
            gamma1_v = np.arcsin(radius1 / L1_v)
            ### ��ߵ���̫��������бͶӰ�ľ��룬��Բ׶�ĽǶ����
            L1_s = height1 * tants
            if L1_s < radius1: L1_s = radius1 * 1.0
            theta1_s = np.arctan(L1_s * np.tan(alpha1) / radius1)
            gamma1_s = np.arcsin(radius1 / L1_s)

            weight = 0
            ### ��������ͶӰֻ��Բ��������Ĳ���
            ### ���Ȳ����ĸ��Ƕȶ��������ٵĵ�ЧԲ��ͶӰ�����
            ### �ڴ���Բ׶�ǶȺ󣬻�����ռƽ̹�ر���Բ׶�������������ǵ���Բ׶����ס��ԭ�����ֻ����ò���
            for kh in range(n_part):
                height_temp = dheight * (kh + 0.5)
                dh = shapes[:, 0] - height_temp
                alpha2 = np.arctan(shapes[:,1] / shapes[:,0])
                r1 = (1 - height_temp / height1) * radius1
                w = (r1*np.pi*dheight+r1*2+dheight)  ### �߶Ȳ����ǵ�Ч�ģ������ǲ�ͬ�߶ȵ�Ȩ��

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

            projv = projv / n_part/weight/(1-poccupied)  ### ��Ϊû�п��Ǳ���Բ���������������Ĳ���
            projs = projs / n_part/weight/(1-poccupied)
            Overlapping = projection_overlap(tantv, tants, up)
            projvs = projv + projs * Overlapping

            # partmax = np.max([projv, projs])
            # partmin = np.min([projv, projs])
            # # if (raa == 0) | (raa==180):
            # #     Overlapping = 0
            # projvs = partmax + partmin * Overlapping

            ### ����ɽ�ص�ƽ�����պͿ��ӱ���
            # exp(-gap_shift) is common to all shapes and cancels in the normalisation by temp; it avoids 0/0 as vza -> 90
            if kshape1 == 0 or projv < gap_shift:
                rescale = 1.0 if kshape1 == 0 else np.exp(projv - gap_shift)
                pMountV, pMountV_sunlit = pMountV * rescale, pMountV_sunlit * rescale
                pMount_veg, pMount_soil = pMount_veg * rescale, pMount_soil * rescale
                pMount_veg_sunlit, pMount_soil_sunlit = pMount_veg_sunlit * rescale, pMount_soil_sunlit * rescale
                pMount_scat, pMount_scat_sunlit = pMount_scat * rescale, pMount_scat_sunlit * rescale
                if self._collect_optical_geometry:
                    optical_slope*=rescale
                    for part in optical_parts:
                        if part[2]:
                            part[0]*=rescale; part[1]*=rescale
                gap_shift = projv
            gapv_mount = np.exp(-(projv - gap_shift))
            gapvs_mount = np.exp(-(projvs - gap_shift))

            slope = alpha1 * 180 / np.pi
            pMountV_temp = 0
            pMountV_sunlit_temp = 0

            if L1_v <= radius1:
                ### ���û��ͶӰ�������ǿ��ӾͶ��Ǳ�Բ
                pMountV_temp =  density1 * gapv_mount * (( np.pi) * radius1 * radius1)
            else:
                ### ͶӰ����Բ���ǿ��ӾͶ��Ǳ�Բ+ͶӰ�����Ĳ��֣���LiԲ׶�ļ��㷽��
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

            # facet normals: zenith = side slope beta1 = atan(h/r), azimuth laa
            laa = np.arange(self.facet_azimuths) * (2.0 * np.pi / self.facet_azimuths)
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
            local_azimuths = None
            if self.vegetation_model == 'row':
                # in-plane basis of each facet: t1 = d(normal)/d(beta), t2 = vertical x normal direction
                phi = laa[ind]
                t1 = np.stack([np.cos(beta1) * np.cos(phi), np.cos(beta1) * np.sin(phi), -np.sin(beta1) * np.ones_like(phi)], -1)
                t2 = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], -1)
                row_dir = _unit(np.asarray(np.pi / 2), np.asarray(self.row_shape[3] * rd))
                r_t = row_dir - (normals @ row_dir)[:, None] * normals
                local = lambda w: np.degrees(np.arctan2(np.sum(w * t2, 1), np.sum(w * t1, 1))) % 360.0
                local_azimuths = (local(v_t), local(s_t), local(r_t))
            # Boolean crowns standing on a facet: densities per facet area are lai*cos(beta), std*cos(beta); the path
            # factor is 1/cos(local zenith) (gap MC on tilted planes: 0.736 vs 0.738; old global-angle form 0.746/0.776 vs 0.603)
            cos_beta = np.cos(beta1)
            Pcomtemp = np.asarray(self._forest_proportions(lai * cos_beta, std * cos_beta, hspot, hcr, rcr, vza_loc, sza_loc, raa_loc, hcb=forestshape[5],
                                                           azimuths=local_azimuths))
            # facets turned away from the sun are entirely shaded
            lit = coss > 0
            Pcomtemp[1] = Pcomtemp[1] + Pcomtemp[0] * ~lit
            Pcomtemp[3] = Pcomtemp[3] + Pcomtemp[2] * ~lit
            Pcomtemp[0] = Pcomtemp[0] * lit
            Pcomtemp[2] = Pcomtemp[2] * lit
            if self._collect_optical_geometry:
                incident=np.maximum(coss,0.)/ui
                phases=np.array([spherical_coefficients(float(ss),float(vv),float(aa))[1:]
                    for ss,vv,aa in zip(sza_loc,vza_loc,raa_loc)])
                optical_slope+=pMountV_temp*mutual*np.array([
                    np.sum(weight*incident*Pcomtemp[2]*phases[:,1]/phases[:,0]),
                    np.sum(weight*incident*Pcomtemp[2]*phases[:,2]/phases[:,0]),
                    np.sum(weight*incident*Pcomtemp[0])])
                view_az,solar_az,row_az=self._canopy_azimuths(raa_loc,local_azimuths)
                view_az,solar_az,row_az=np.broadcast_arrays(view_az,solar_az,row_az)
                optical_parts.extend([float(pMountV_temp*mutual*w*l*i),float(pMountV_temp*w),1.,
                                      float(vv),float(ss),float(aa),float(va),float(sa),float(ra)]
                                     for w,l,i,vv,ss,aa,va,sa,ra in zip(weight,lit,incident,vza_loc,sza_loc,raa_loc,
                                                                         view_az,solar_az,row_az))
            Pcom = np.sum(Pcomtemp * weight,axis=1)
            pVeg = Pcom[2]+Pcom[3]
            pSoil = Pcom[0] + Pcom[1]
            pVeg_sunlit = Pcom[2]
            pSoil_sunlit = Pcom[0]
            pMount_veg = pMount_veg + pMountV_temp * pVeg                                ### ֲ���Ĺ���
            pMount_soil = pMount_soil + pMountV_temp * pSoil                             ### �����Ĺ���
            pMount_soil_sunlit = pMount_soil_sunlit + pMountV_temp * mutual * pSoil_sunlit
            pMount_veg_sunlit = pMount_veg_sunlit + pMountV_temp * mutual * pVeg_sunlit
            # forest-background scattering of the facet canopy in facet-local angles; unlit facets are all shaded
            if self.thermal_forest_transport=='reciprocal_exchange':
                scat=np.zeros((4,len(vza_loc)))
            else:
                scat = self._canopy_scatter(lai * cos_beta, vza_loc, sza_loc, std * cos_beta, hcr, rcr, raa_loc,
                                            local_azimuths)
            pMount_scat = pMount_scat + pMountV_temp * np.sum((scat[2] + scat[3]) * weight)
            pMount_scat_sunlit = pMount_scat_sunlit + pMountV_temp * mutual * np.sum(scat[2] * lit * weight)
            pMount_scat_soil = pMount_scat_soil + pMountV_temp * np.sum((scat[0] + scat[1]) * weight)
            pMount_scat_soil_sunlit = pMount_scat_soil_sunlit + pMountV_temp * mutual * np.sum(scat[0] * lit * weight)




        ### �������µ���ƽ̹�ر���Ĺ�һ��������temp���µصı仯����
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
        pMountV_scat_soil_sunlit = pMount_scat_soil_sunlit * temp
        pMountV_scat_soil_shaded = pMount_scat_soil * temp - pMountV_scat_soil_sunlit
        # composition of the facet surface seen by the sensor, used as reflector in the scatter term
        self._leaf_view_facet = pMount_veg / (pMount_veg + pMount_soil)
        self._reciprocal_first_hit_states=np.r_[pMountVnew*np.array([
            1-self._leaf_view_facet,self._leaf_view_facet]),self._reciprocal_ground_states]
        if self._collect_optical_geometry:
            self._optical_first_coefficients=optical_flat+temp*optical_slope
            parts=np.asarray(optical_parts)
            parts[:,:2]*=np.where(parts[:,2:3]>0,temp,1.)
            self._optical_canopy_parts=parts
        ### ���Լ����ӳ�1�ͼӳ�2
        alla2 = pMountV_veg_sunlit+pMountV_veg_shaded+pMountV_soil_shaded+pMountV_soil_sunlit
        alla = pPlaneV_veg_sunlit + pMountV_veg_sunlit+ pPlaneV_veg_shaded + pMountV_veg_shaded+ pPlaneV_soil_sunlit + pMountV_soil_sunlit+ pPlaneV_soil_shaded + pMountV_soil_shaded
        # print(pPlaneV_sunlit,pPlaneV_shaded,pMountV_sunlit,pMountV_shaded)
        # print(Pcom,np.sum(Pcom))
        if ifP ==1:
            return pPlaneV, pMountVnew
        elif ifP == 2:
            return (pPlaneV_soil_sunlit + pMountV_soil_sunlit) * self.Es + pPlaneV_scat_soil_sunlit + pMountV_scat_soil_sunlit,\
                   (pPlaneV_soil_shaded + pMountV_soil_shaded) * self.Es + pPlaneV_scat_soil_shaded + pMountV_scat_soil_shaded,\
                   (pPlaneV_veg_sunlit + pMountV_veg_sunlit) * self.Em + pPlaneV_scat_sunlit + pMountV_scat_sunlit,\
                   (pPlaneV_veg_shaded + pMountV_veg_shaded) * self.Em + pPlaneV_scat_shaded + pMountV_scat_shaded

        elif ifP ==3:
            return 0, 0
        else:
            return pPlaneV+ pMountV

    def calculate_component_scatter_emissivity(self, kangle,ifP = 0):
        """Legacy terrain reflection or the opt-in reciprocal joint closure."""
        if self.thermal_forest_transport=='reciprocal_exchange':
            from .thermal_forest_exchange import scatter_components
            parts=scatter_components(self)
            if ifP==2:
                return tuple(parts)
            if ifP==1:
                return parts[:2].sum(),parts[2:].sum()
            if ifP==3:
                return (0.,0.,0.,0.)
            return float(parts.sum())
        # Retain the legacy additive terrain reflection for existing inputs.
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
            # Independent of the view direction: computed once per sun position and canopy.
            key = (self.vegetation_model, self.thermal_crown_geometry, str(self.thermal_hotspot_inner),
                   str(self.thermal_hotspot_outer), self.thermal_forest_transport, tuple(self.row_shape),
                   self.row_geometry, self.canopy_fraction_model, float(lai_x), float(std_x), float(hspot),
                   float(hcr), float(rcr), float(hcb), float(sza), float(self.saa))
            if key in _HEMISPHERICAL_CACHE:
                return _HEMISPHERICAL_CACHE[key]
            azimuths = (self.saa + REL.reshape(-1), self.saa, self.row_shape[3])
            P = np.asarray(self._forest_proportions(lai_x, std_x, hspot, hcr, rcr, zen, np.full(zen.size, min(sza, 89.0)), REL.reshape(-1), hcb=hcb,
                                                    azimuths=azimuths))
            wt = WT.reshape(-1)
            gap = np.sum(wt * (P[0] + P[1]))
            soil_lit = np.sum(wt * P[0]) / max(np.sum(wt * (P[0] + P[1])), 1e-12)
            leaf_lit = np.sum(wt * P[2]) / max(np.sum(wt * (P[2] + P[3])), 1e-12)
            if len(_HEMISPHERICAL_CACHE) > 4096:
                _HEMISPHERICAL_CACHE.clear()
            _HEMISPHERICAL_CACHE[key] = (gap, soil_lit, leaf_lit)
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
