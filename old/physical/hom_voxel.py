
#####################################################
#### Homgeneous vegetation scene
#####################################################

import numpy as np
from old.utils import *
from old.hotspot import *
from old.gap import *
from old.scatter import *
from old.proportion import *
from old.emissivity import *
from old.radiance import *

class Hom_Voxel():

    ### structual variables
    leaf_area_index = 1.0
    clumping_index = 1.0
    leaf_average_inclination_angle = 53
    number_voxel = 10
    number_vegetation = 10
    number_background = 1
    hotspot = 0.2
    G = 0.5

    ### optical and thermal variables
    number_band = 1
    emissivity_leaf = np.asarray([0.985])
    emissivity_soil = np.asarray([0.955])

    temperature_sunlit_leaf = np.asarray([303])
    temperature_shaded_leaf = np.asarray([300])
    temperature_sunlit_soil = np.asarray([303])
    temperature_shaded_soil = np.asarray([300])

    ### viewing and solar gemoetry
    number_angle = 1
    view_zenith_angle = np.asarray([10])
    solar_zenith_angle = np.asarray([25])
    azimuth_angle = np.asarray([0])


    def set_structure(self,leaf_area_index,hspot,number_voxel = 50):
        self.leaf_area_index = leaf_area_index
        self.hspot = hspot
        self.number_voxel = number_voxel
        self.number_vegetation = number_voxel

    def set_thermal(self,temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded,):
        self.temperature_soil_shaded = np.asarray(temperature_soil_shaded)
        self.temperature_soil_sunlit = np.asarray(temperature_soil_sunlit)
        self.temperature_leaf_sunlit = np.asarray(temperature_leaf_sunlit)
        self.temperature_leaf_shaded = np.asarray(temperature_leaf_shaded)

    def set_optical(self,wavelength,emissivity_soil,emissivity_leaf):
        self.wavelength = wavelength
        self.emissivity_soil = emissivity_soil
        self.emissivity_leaf = emissivity_leaf


    def set_angle(self,view_zenith_angle,solar_zenith_angle,relative_azimuth_angle):

        self.number_angle = np.max([self.number_angle, np.size(view_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
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

        vza = self.viewing_zenith_angle
        sza = self.solar_zenith_angle
        vsa = self.relative_azimuth_angle
        lai = self.leaf_area_index
        hspot = self.hspot

        ### Pleaf 和 Psoil 是植被和土壤组分的可视比例
        ### Ksoil 和 Kleaf 是植被和土壤组分的可视比例的光照比例 K = Psunlitvisible / Pvisible

        Ecom = np.asarray([es,es,el,el])
        Tcom = np.asarray([Tss,Tsh,Tls,Tlh])
        Rcom = planck(Tcom)

        Pcom = proportion_bidirectional_hom_voxel_one(lai,hspot,vza,sza,vsa)
        Ecom_direct,Ecom_direct_sum = emissivity_direct(Pcom,Ecom)
        Rcom_direct = radiance_direct(Ecom_direct,Rcom)

        # Ecom_scatter,Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vza, sza)
        Ecom_scatter,Ecom_scatter_sum = emissivity_scatter_hom_voxel_one(lai,el,es,vza,sza,vsa)
        Rcom_scatter = radiance_scatter(Ecom_scatter,Rcom)
        self.radiance = Rcom_direct + Rcom_scatter

        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(self.radiance, self.wavelength)