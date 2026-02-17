

#####################################################
#### crown-shape scene
#####################################################

import numpy as np
from rt.hotspot import *
from rt.gap import *
from rt.scatter import *
from rt.proportion import *
from rt.radiance import *
from rt.emissivity import *


class Crown_Voxel():
    ### structual variables
    leaf_area_index = 3.0
    clumping_index = 1.0
    stand_density = 0.05
    hcr = 1.0
    rcr = 1.0
    offz = 0.0
    leaf_average_inclination_angle = 53
    hspot = 0.2
    G = 0.5
    Pcom = 0
    Ecom_scatter = 0

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
    view_zenith_angle = np.asarray([10])
    solar_zenith_angle = np.asarray([25])
    view_azimuth_angle = np.asarray([0])
    solar_azimuth_angle = np.asarray([0])

    def set_structure(self, lai, hspot, stand_density, hcr, rcr,offz = 1.0):
        self.leaf_area_index = lai
        self.hspot = hspot
        self.stand_density = stand_density
        self.rcr = rcr
        self.hcr = hcr
        self.hc = 2* hcr + offz
        self.offz = offz

    def set_thermal(self,temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded,):
        self.temperature_soil_shaded = np.asarray(temperature_soil_shaded)
        self.temperature_soil_sunlit = np.asarray(temperature_soil_sunlit)
        self.temperature_leaf_sunlit = np.asarray(temperature_leaf_sunlit)
        self.temperature_leaf_shaded = np.asarray(temperature_leaf_shaded)

    def set_optical(self,wavelength,emissivity_soil,emissivity_leaf):
        self.wavelength = wavelength
        self.emissivity_soil = emissivity_soil
        self.emissivity_leaf = emissivity_leaf

    def set_angle(self,view_zenith_angle,solar_zenith_angle,view_azimuth_angle,solar_azimuth_angle=0):

        self.number_angle = np.max([self.number_angle, np.size(view_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(view_azimuth_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_azimuth_angle)])

        if (self.number_angle > 1) and type(view_zenith_angle) != numpy.ndarray:
            view_zenith_angle = np.resize(view_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(solar_zenith_angle) != numpy.ndarray:
            solar_zenith_angle = np.resize(solar_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(view_azimuth_angle) != numpy.ndarray:
            view_azimuth_angle = np.resize(view_azimuth_angle, self.number_angle)
        if (self.number_angle > 1) and type(solar_azimuth_angle) != numpy.ndarray:
            solar_azimuth_angle = np.resize(solar_azimuth_angle, self.number_angle)
        self.view_zenith_angle = np.asarray(view_zenith_angle)
        self.solar_zenith_angle = np.asarray(solar_zenith_angle)
        self.view_azimuth_angle = np.asarray(view_azimuth_angle)
        self.solar_azimuth_angle = np.asarray(solar_azimuth_angle)


    def run(self,ifradiance = 0):

        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tls = self.temperature_leaf_sunlit
        Tlh = self.temperature_leaf_shaded
        el = self.emissivity_leaf
        es = self.emissivity_soil

        vza = self.view_zenith_angle
        sza = self.solar_zenith_angle
        vaa = self.view_azimuth_angle
        saa = self.solar_azimuth_angle
        lai = self.leaf_area_index
        hspot = self.hspot
        std = self.stand_density
        offz = self.offz

        hcr = self.hcr
        rcr = self.rcr
        vsa = np.abs(vaa - saa)
        ### Pleaf 和 Psoil 是植被和土壤组分的可视比例
        ### Ksoil 和 Kleaf 是植被和土壤组分的可视比例的光照比例 K = Psunlitvisible / Pvisible

        Ecom = np.asarray([es, es, el, el])
        Tcom = np.asarray([Tss, Tsh, Tls, Tlh])
        Rcom = planck(Tcom)

        hc = hcr * 2 + offz
        Pcom = proportion_bidirectional_crown_voxel_one(lai, std, hspot,hc, hcr, rcr, vza, sza, vaa,saa)
        Ecom_direct, Ecom_direct_sum = emissivity_direct(Pcom, Ecom)
        Rcom_direct = radiance_direct(Ecom_direct, Rcom)

        Rcom_scatter = 0
        # Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vza, sza)
        Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_crown_voxel_one(lai,std,hc,hcr,rcr,el,es,vza,sza,vaa,saa)
        Rcom_scatter = radiance_scatter(Ecom_scatter, Rcom)

        self.radiance = Rcom_direct + Rcom_scatter
        self.Ecom_scatter = Ecom_scatter_sum
        self.Pcom = Pcom

        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(self.radiance, self.wavelength)

    def runInv(self,ifradiance = 0):


        el = self.emissivity_leaf
        es = self.emissivity_soil

        vza = self.view_zenith_angle
        sza = self.solar_zenith_angle
        vaa = self.view_azimuth_angle
        saa = self.solar_azimuth_angle
        lai = self.leaf_area_index
        hspot = self.hspot
        std = self.stand_density
        offz = self.offz

        hcr = self.hcr
        rcr = self.rcr
        vsa = np.abs(vaa - saa)
        ### Pleaf 和 Psoil 是植被和土壤组分的可视比例
        ### Ksoil 和 Kleaf 是植被和土壤组分的可视比例的光照比例 K = Psunlitvisible / Pvisible

        Ecom = np.asarray([es, es, el, el])


        hc = hcr * 2 + offz
        Pcom = proportion_bidirectional_crown_voxel_one(lai, std, hspot,hc, hcr, rcr, vza, sza, vaa,saa)
        Ecom_direct, Ecom_direct_sum = emissivity_direct(Pcom, Ecom)

        # Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vza, sza)
        Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_crown_voxel_one(lai,std,hc,hcr,rcr,el,es,vza,sza,vaa,saa)

        return Ecom_direct + Ecom_scatter


    def fetch_Pcom(self):
        return self.Pcom


    def fetch_Scatter(self):
        return self.Ecom_scatter