

#####################################################
#### Homgeneous vegetation scene
#####################################################

import numpy as np
from .hotspot import *
from .gap import *
from .scatter import *
from .proportion import *
from .radiance import *
from .emissivity import *
from .vegetation import vegetation_proportions
from .thermal_crown import geometry_mode, thermal_emission

class Crown():

    ### structual variables
    leaf_area_index = 3.0
    clumping_index = 1.0
    stand_density = 0.05
    hcr = 1.0
    rcr = 1.0
    leaf_average_inclination_angle = 53
    hspot = 0.2
    G = 0.5
    thermal_crown_geometry = "legacy"
    thermal_hotspot_inner = None
    thermal_hotspot_outer = None
    wavelength = 10.5

    ### optical and thermal variables
    number_band = 1
    emissivity_leaf = np.asarray([0.985])
    emissivity_soil = np.asarray([0.955])

    temperature_leaf_sunlit = np.asarray([303])
    temperature_leaf_shaded = np.asarray([300])
    temperature_soil_sunlit = np.asarray([303])
    temperature_soil_shaded = np.asarray([300])

    ### viewing and solar gemoetry
    number_angle = 1
    viewing_zenith_angle = np.asarray([10])
    solar_zenith_angle = np.asarray([25])
    relative_azimuth_angle = np.asarray([0])



    def set_structure(self,lai,hspot,stand_density,hcr,rcr):
        self.leaf_area_index = lai
        self.hspot = hspot
        self.stand_density = stand_density
        self.rcr = rcr
        self.hcr = hcr

    def set_hotspot_geometry(self, mode="legacy", inner=None, outer=None):
        self.thermal_crown_geometry = geometry_mode({"thermal_crown_geometry": mode})
        self.thermal_hotspot_inner = inner
        self.thermal_hotspot_outer = outer

    def set_thermal(self,temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded):
        self.temperature_soil_shaded = np.asarray(temperature_soil_shaded)
        self.temperature_soil_sunlit = np.asarray(temperature_soil_sunlit)
        self.temperature_leaf_sunlit = np.asarray(temperature_leaf_sunlit)
        self.temperature_leaf_shaded = np.asarray(temperature_leaf_shaded)

    def set_optical(self,wavelength,emissivity_soil,emissivity_leaf):
        self.wavelength = wavelength
        self.emissivity_soil = emissivity_soil
        self.emissivity_leaf = emissivity_leaf

    def set_angle(self,view_zenith_angle,solar_zenith_angle,relative_azimuth_angle):
        if type(view_zenith_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(view_zenith_angle)])
        if type(solar_zenith_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
        if type(relative_azimuth_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(relative_azimuth_angle)])

        if (self.number_angle > 1) and type(view_zenith_angle) != numpy.ndarray:
            view_zenith_angle = np.repeat(view_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(solar_zenith_angle) != numpy.ndarray:
            solar_zenith_angle = np.repeat(solar_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(relative_azimuth_angle) != numpy.ndarray:
            relative_azimuth_angle = np.repeat(relative_azimuth_angle, self.number_angle)

        self.viewing_zenith_angle = np.asarray(view_zenith_angle)
        self.solar_zenith_angle = np.asarray(solar_zenith_angle)
        self.relative_azimuth_angle = np.asarray(relative_azimuth_angle)


    def run(self,ifradiance = 0):

        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tls = self.temperature_leaf_sunlit
        Tlh = self.temperature_leaf_shaded
        el = self.emissivity_leaf
        es = self.emissivity_soil
        wl = self.wavelength

        vza = self.viewing_zenith_angle
        sza = self.solar_zenith_angle
        vsa = self.relative_azimuth_angle
        lai = self.leaf_area_index
        hspot = self.hspot
        std = self.stand_density
        hcr = self.hcr
        rcr = self.rcr

        if self.thermal_crown_geometry != "legacy":
            settings = dict(lai=lai, stand_density=std, height=hcr, width=rcr,
                            hspot=hspot, soil_emissivity=float(np.asarray(es).item()),
                            leaf_emissivity=float(np.asarray(el).item()),
                            thermal_crown_geometry=self.thermal_crown_geometry,
                            thermal_hotspot_inner=self.thermal_hotspot_inner,
                            thermal_hotspot_outer=self.thermal_hotspot_outer,
                            crown_base=getattr(self, "crown_base", None))
            self.radiance, temperature, self.Pcom = thermal_emission(
                settings, [Tss, Tsh, Tls, Tlh], wl, sza, vza, vsa)
            return self.radiance if ifradiance == 1 else temperature

        ### Pleaf 和 Psoil 是植被和土壤组分的可视比例
        ### Ksoil 和 Kleaf 是植被和土壤组分的可视比例的光照比例 K = Psunlitvisible / Pvisible

        Ecom = np.asarray([es, es, el, el])
        Tcom = np.asarray([Tss, Tsh, Tls, Tlh])
        Rcom = planck(wl,Tcom)

        Pcom = vegetation_proportions(
            "crown", lai, hspot, vza, vsa, sza, 0.0,
            density=std, crown_height=hcr, crown_radius=rcr)
        Ecom_direct, Ecom_direct_sum = emissivity_direct(Pcom, Ecom)
        Rcom_direct = radiance_direct(Ecom_direct, Rcom)

        Ecom_scatter = canopy_scatter("crown", lai, el, es, vza, sza, std, hcr, rcr)
        # Limit scatter to the remaining emissivity budget, as in thermal_crown.
        budget = np.maximum(0.0, 1.0 - np.sum(Ecom_direct, axis=0))
        total = np.sum(Ecom_scatter, axis=0)
        Ecom_scatter = Ecom_scatter * np.minimum(1.0, np.divide(budget, total, out=np.ones_like(total), where=total > 0))
        Rcom_scatter = radiance_scatter(Ecom_scatter, Rcom)

        self.radiance = Rcom_direct + Rcom_scatter

        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(self.wavelength,self.radiance)
