import numpy as np
import numpy
from base.hotspot import *
from base.gap import *
from base.scatter import *
from base.proportion import *
from base.emissivity import *
from base.radiance import *

class Row_Voxel():

    leaf_area_index = 1.0
    clumping_index = 1.0
    leaf_average_inclination_angle = 53
    number_voxel_width = 50
    number_voxel_blank = 50
    number_voxel_height = 100
    number_voxel_spacing = 100
    hotspot = 0.2
    G = 0.5
    Pcom = 0
    Ecom_scatter = 0

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
    view_azimuth_angle = np.asarray([0])
    solar_azimuth_angle = np.asarray([0])

    def set_structure(self, lai, hspot, row_width, row_blank, row_height, number_voxel_spacing = 150):
        self.leaf_area_index = lai
        self.hspot = hspot
        self.row_width = row_width
        self.row_blank = row_blank
        self.row_height = row_height
        row_spacing = row_width + row_blank
        scale = row_spacing / number_voxel_spacing
        self.number_voxel_height = np.int_(row_height/scale)
        self.number_voxel_width = np.int_(row_width /scale)
        self.number_voxel_blank = number_voxel_spacing - self.number_voxel_width

    def set_thermal(self,temperature_soil_sunlit,temperature_soil_shaded,
                    temperature_leaf_sunlit,temperature_leaf_shaded,):
        self.temperature_soil_shaded = np.asarray(temperature_soil_shaded)
        self.temperature_soil_sunlit = np.asarray(temperature_soil_sunlit)
        self.temperature_leaf_sunlit = np.asarray(temperature_leaf_sunlit)
        self.temperature_leaf_shaded = np.asarray(temperature_leaf_shaded)

    def set_optical(self,wavelength,emissivity_soil,emissivity_leaf):
        self.wavelength = wavelength
        self.emissivity_soil = emissivity_soil
        self.emissivity_leaf = emissivity_leaf


    def set_angle(self,view_zenith_angle,solar_zenith_angle,
                  view_azimuth_angle,solar_azimuth_angle,row_azimuth_angle):

        self.number_angle = np.max([self.number_angle,np.size(view_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
        self.number_angle = np.max([self.number_angle,np.size(view_azimuth_angle)])
        self.number_angle = np.max([self.number_angle,np.size(solar_azimuth_angle)])
        self.number_angle = np.max([self.number_angle,np.size(row_azimuth_angle)])

        if (self.number_angle >= 1) and type(view_zenith_angle) != numpy.ndarray:
            view_zenith_angle = np.resize(view_zenith_angle,self.number_angle)
        if (self.number_angle >= 1) and type(solar_zenith_angle) != numpy.ndarray:
            solar_zenith_angle = np.resize(solar_zenith_angle, self.number_angle)
        if (self.number_angle >= 1) and type(view_azimuth_angle) != numpy.ndarray:
            view_azimuth_angle = np.resize(view_azimuth_angle, self.number_angle)
        if (self.number_angle >= 1) and type(solar_azimuth_angle) != numpy.ndarray:
            solar_azimuth_angle = np.resize(solar_azimuth_angle, self.number_angle)
        if (self.number_angle >= 1) and type(row_azimuth_angle) != numpy.ndarray:
            row_azimuth_angle = np.resize(row_azimuth_angle, self.number_angle)


        self.view_zenith_angle = np.asarray(view_zenith_angle)
        self.solar_zenith_angle = np.asarray(solar_zenith_angle)
        self.view_azimuth_angle = np.asarray(view_azimuth_angle)
        self.solar_azimuth_angle = np.asarray(solar_azimuth_angle)
        self.row_azimuth_angle = np.asarray(row_azimuth_angle)


    def run(self, ifradiance = 0):

        vza = self.view_zenith_angle
        sza = self.solar_zenith_angle
        raa = self.row_azimuth_angle
        saa = self.solar_azimuth_angle
        vaa = self.view_azimuth_angle
        vsa = np.abs(saa - vaa)
        sra = np.abs(saa - raa) % 180
        vra = np.abs(vaa - raa) % 180
        vsa[vsa>180] = 360 - vsa[vsa>180]


        lai = self.leaf_area_index
        G = self.G
        hspot = self.hspot
        row_width = self.row_width
        row_blank = self.row_blank
        row_height = self.row_height
        nvw = self.number_voxel_width
        nvh = self.number_voxel_height
        nvb = self.number_voxel_blank

        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tls = self.temperature_leaf_sunlit
        Tlh = self.temperature_leaf_shaded
        el = self.emissivity_leaf
        es = self.emissivity_soil


        Ecom = np.asarray([es,es,el,el])
        Tcom = np.asarray([Tss,Tsh,Tls,Tlh])
        Rcom = planck(self.wavelength,Tcom)

        Pcom = proportion_bidirectional_row_voxel_one(
            lai,hspot,row_width,row_blank,row_height,vza,vaa,sza,saa,raa,nvw,nvh,nvb)
        Ecom_direct,Ecom_direct_sum = emissivity_direct(Pcom,Ecom)
        Rcom_direct = radiance_direct(Ecom_direct,Rcom)

        # Ecom_scatter,Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vza, sza)
        Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_row_voxel_one(
            lai,row_width,row_blank,row_height, el, es, vza, sza,vaa,saa,raa)
        Rcom_scatter = radiance_scatter(Ecom_scatter,Rcom)

        self.radiance = Rcom_direct + Rcom_scatter
        self.Ecom_scatter = Ecom_scatter_sum
        self.Pcom = Pcom


        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(self.wavelength, self.radiance)

    def runInv(self, ifradiance = 0):

        vza = self.view_zenith_angle
        sza = self.solar_zenith_angle
        raa = self.row_azimuth_angle
        saa = self.solar_azimuth_angle
        vaa = self.view_azimuth_angle
        vsa = np.abs(saa - vaa)
        sra = np.abs(saa - raa) % 180
        vra = np.abs(vaa - raa) % 180
        vsa[vsa>180] = 360 - vsa[vsa>180]


        lai = self.leaf_area_index
        G = self.G
        hspot = self.hspot
        row_width = self.row_width
        row_blank = self.row_blank
        row_height = self.row_height
        nvw = self.number_voxel_width
        nvh = self.number_voxel_height
        nvb = self.number_voxel_blank

        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tls = self.temperature_leaf_sunlit
        Tlh = self.temperature_leaf_shaded
        el = self.emissivity_leaf
        es = self.emissivity_soil


        Ecom = np.asarray([es,es,el,el])
        Tcom = np.asarray([Tss,Tsh,Tls,Tlh])
        Rcom = planck(self.wavelength,Tcom)

        Pcom = proportion_bidirectional_row_voxel_one(
            lai,hspot,row_width,row_blank,row_height,vza,vaa,sza,saa,raa,nvw,nvh,nvb)
        Ecom_direct,Ecom_direct_sum = emissivity_direct(Pcom,Ecom)


        # Ecom_scatter,Ecom_scatter_sum = emissivity_scatter_analytical(lai, el, es, vza, sza)
        Ecom_scatter, Ecom_scatter_sum = emissivity_scatter_row_voxel_one(
            lai,row_width,row_blank,row_height, el, es, vza, sza,vaa,saa,raa)

        return Ecom_direct + Ecom_scatter



    def fetch_Pcom(self):
        return self.Pcom


    def fetch_Scatter(self):
        return self.Ecom_scatter
