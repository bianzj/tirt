

#####################################################
#### Homgeneous vegetation scene
#####################################################

import numpy as np
from rt.hotspot import *
from rt.gap import *
from rt.scatter import *
from rt.proportion import *
from rt.radiance import *
from rt.emissivity import *


class Slope_Veg():

    ### structual variables
    rd = np.pi/180.0
    leaf_area_index = 3.0
    clumping_index = 1.0
    stand_density = 0.05
    hcr = 1.0
    rcr = 1.0
    leaf_average_inclination_angle = 53
    pza = 0
    paa = 0
    hspot = 0.2
    G = 0.5
    wl = 10.5

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
    viewing_azimuth_angle = np.asarray([0])
    solar_azimth_angle = np.asarray([0])

    viewing_zenith_angle_norm = np.asarray([10])
    solar_zenith_angle_norm = np.asarray([25])
    viewing_azimuth_angle_norm = np.asarray([0])
    solar_azimuth_angle_norm = np.asarray([0])


    def set_slope(self,pza=0,paa=0):
        self.pza = pza
        self.paa = paa


    def set_structure(self,lai,hspot,stand_density,hcr,rcr):
        self.leaf_area_index = lai
        self.hspot = hspot
        self.stand_density = stand_density
        self.rcr = rcr
        self.hcr = hcr

    def set_thermal(self,temperature_soil_sunlit,temperature_soil_shaded,temperature_leaf_sunlit,temperature_leaf_shaded,):
        self.temperature_soil_shaded = np.asarray(temperature_soil_shaded)
        self.temperature_soil_sunlit = np.asarray(temperature_soil_sunlit)
        self.temperature_leaf_sunlit = np.asarray(temperature_leaf_sunlit)
        self.temperature_leaf_shaded = np.asarray(temperature_leaf_shaded)

    def set_optical(self,wavelength,emissivity_soil,emissivity_leaf):
        self.wavelength = wavelength
        self.emissivity_soil = emissivity_soil
        self.emissivity_leaf = emissivity_leaf

    def set_angle(self,view_zenith_angle,solar_zenith_angle,view_azimuth_angle,solar_azimuth_angle):
        self.number_angle = 1
        if type(view_zenith_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(view_zenith_angle)])
        if type(solar_zenith_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
        if type(view_azimuth_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(view_azimuth_angle)])
        if type(solar_azimuth_angle) == numpy.ndarray:
            self.number_angle = np.max([self.number_angle, np.size(solar_azimuth_angle)])

        if (self.number_angle > 1) and type(view_zenith_angle) != numpy.ndarray:
            view_zenith_angle = np.repeat(view_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(solar_zenith_angle) != numpy.ndarray:
            solar_zenith_angle = np.repeat(solar_zenith_angle, self.number_angle)
        if (self.number_angle > 1) and type(view_azimuth_angle) != numpy.ndarray:
            view_azimuth_angle = np.repeat(view_azimuth_angle, self.number_angle)
        if (self.number_angle > 1) and type(solar_azimuth_angle) != numpy.ndarray:
            solar_azimuth_angle = np.repeat(solar_azimuth_angle, self.number_angle)

        self.viewing_zenith_angle = np.asarray(view_zenith_angle)
        self.solar_zenith_angle = np.asarray(solar_zenith_angle)
        self.viewing_azimuth_angle = np.asarray(view_azimuth_angle)
        self.solar_azimuth_angle = np.asarray(solar_azimuth_angle)

    def angle_transfer(self):
        self.viewing_zenith_angle_norm,self.viewing_azimuth_angle_norm = \
            slope1(self.viewing_zenith_angle,self.viewing_azimuth_angle,self.pza,self.paa)
        self.solar_zenith_angle_norm,self.solar_azimuth_angle_norm = \
            slope1(self.solar_zenith_angle,self.solar_azimuth_angle,self.pza,self.paa)

        # hcr = self.hcr
        # rcr = self.rcr
        # rd = self.rd
        # ###----------------
        # ### tree crown
        # ###----------------
        # vza = self.viewing_zenith_angle
        # sza = self.solar_zenith_angle
        # theta = np.arctan(hcr / rcr * np.tan(np.deg2rad(vza)))
        # vzanew = theta / rd
        # theta = np.arctan(hcr / rcr * np.tan(np.deg2rad(sza)))
        # szanew = theta / rd
        # self.viewing_zenith_angle_norm,self.viewing_azimuth_angle_norm = \
        #     slope1(vzanew,self.viewing_azimuth_angle,self.pza,self.paa)
        # self.solar_zenith_angle_norm,self.solar_azimuth_angle_norm = \
        #     slope1(szanew,self.solar_azimuth_angle,self.pza,self.paa)

    def run(self,ifradiance = 0):

        self.angle_transfer()
        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tls = self.temperature_leaf_sunlit
        Tlh = self.temperature_leaf_shaded
        el = self.emissivity_leaf
        es = self.emissivity_soil
        wl = self.wavelength

        ###----------------
        ### slope structure
        ###----------------
        vza = self.viewing_zenith_angle_norm
        sza = self.solar_zenith_angle_norm
        saa = self.solar_azimuth_angle_norm
        vaa = self.viewing_azimuth_angle_norm



        sza0 = self.solar_zenith_angle
        saa0 = self.solar_azimuth_angle
        vsa = np.abs(vaa - saa)

        pza = self.pza
        paa = self.paa
        lai = self.leaf_area_index
        hspot = self.hspot
        std = self.stand_density
        hcr = self.hcr
        rcr = self.rcr
        rd = self.rd

        ###----------------
        ### tree crown
        ###----------------
        theta = np.arctan(hcr / rcr * np.tan(np.deg2rad(vza)))
        vzanew = theta / rd
        theta = np.arctan(hcr / rcr * np.tan(np.deg2rad(sza)))
        szanew = theta / rd

        rd = np.pi/180.0
        stdnew = std * np.cos(pza * rd)
        hcrnew = hcr * np.cos(pza * rd)
        hspotnew = hspot / np.cos(pza * rd)

        stdnew = std
        hspotnew = hspot
        hcrnew = hcr
        vzanew = vza
        szanew = sza

        ### Pleaf 和 Psoil 是植被和土壤组分的可视比例
        ### Ksoil 和 Kleaf 是植被和土壤组分的可视比例的光照比例 K = Psunlitvisible / Pvisible

        Ecom = np.asarray([es, es, el, el])
        Tcom = np.asarray([Tss, Tsh, Tls, Tlh])
        Rcom = planck(wl,Tcom)

        Pcom = proportion_bidirectional_crown_one(lai, stdnew, hspotnew, hcrnew,rcr, vzanew, szanew, vsa)
        Ecom_direct, Ecom_direct_sum = emissivity_direct(Pcom, Ecom)

        ############---------------------------
        ### to check if the  slope is all in the shadow
        ### vector i and vector s under vector p
        rd = np.pi/180.0
        uii = np.cos(sza0 * rd)
        uvv = np.cos(pza * rd)
        sii = np.sin(sza0 * rd)
        svv = np.sin(pza * rd)
        upp = np.cos(np.abs(paa - saa) * rd)
        cosang = uvv * uii + svv * sii * upp
        ang = np.arccos(cosang) / rd
        ###------------------------------------
        if np.min(cosang)<0:
            Ecom_direct_new = np.vstack([np.zeros(self.number_angle),Ecom_direct[0]+Ecom_direct[1],np.zeros(self.number_angle),Ecom_direct[2]+Ecom_direct[3]])
        else:
            Ecom_direct_new = Ecom_direct
        Rcom_direct = radiance_direct(Ecom_direct_new, Rcom)
        Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vzanew, szanew)

        if np.min(cosang) < 0:
            Ecom_scatter_new = np.vstack([np.zeros(self.number_angle), Ecom_scatter[0] + Ecom_scatter[1], np.zeros(self.number_angle),
                 Ecom_scatter[2] + Ecom_scatter[3]])
        else:
            Ecom_scatter_new = Ecom_scatter
        Rcom_scatter = radiance_scatter(Ecom_scatter_new, Rcom)
        self.radiance = Rcom_direct + Rcom_scatter

        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(wl,self.radiance)



