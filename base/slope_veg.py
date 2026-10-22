

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
from .vegetation import CanopyGeometry


class Slope_Veg(CanopyGeometry):
    """Slope-local canopy geometry with the analytical scattering approximation."""

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
    wavelength = 10.5
    forestshape = np.asarray([3.0, 0.05, 1.0, 1.0, 0.2])

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
        self.forestshape = np.asarray([lai, stand_density, rcr, hcr, hspot])

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

    def run(self, ifradiance=0):
        vza, vaa = self.viewing_zenith_angle, self.viewing_azimuth_angle
        sza, saa = self.solar_zenith_angle, self.solar_azimuth_angle
        pza, paa = self.pza, self.paa
        local_view_cos = (np.cos(np.deg2rad(pza)) * np.cos(np.deg2rad(vza)) +
                          np.sin(np.deg2rad(pza)) * np.sin(np.deg2rad(vza)) *
                          np.cos(np.deg2rad(vaa - paa)))
        local_sun_cos = (np.cos(np.deg2rad(pza)) * np.cos(np.deg2rad(sza)) +
                         np.sin(np.deg2rad(pza)) * np.sin(np.deg2rad(sza)) *
                         np.cos(np.deg2rad(saa - paa)))
        local_vza = np.clip(np.rad2deg(np.arccos(np.clip(local_view_cos, 0, 1))), 0, 89.999)
        local_sza = np.clip(np.rad2deg(np.arccos(np.clip(local_sun_cos, 0, 1))), 0, 89.999)
        proportions = self.local_canopy_proportions(vza, vaa, sza, saa, pza, paa)
        emissivities = np.asarray([self.emissivity_soil, self.emissivity_soil,
                                   self.emissivity_leaf, self.emissivity_leaf])
        temperatures = np.asarray([self.temperature_soil_sunlit, self.temperature_soil_shaded,
                                   self.temperature_leaf_sunlit, self.temperature_leaf_shaded])
        component_radiance = planck(self.wavelength, temperatures)
        direct, _ = emissivity_direct(proportions, emissivities)
        lai = self.leaf_area_index
        empty = lai == 0.0 or (self.vegetation_model == "crown" and self.stand_density == 0.0)
        if empty:
            scatter = np.zeros_like(proportions)
        else:
            _, density, crown_radius, crown_height = self.forestshape[:4]
            scatter = canopy_scatter(self.vegetation_model, lai, self.emissivity_leaf, self.emissivity_soil,
                                     local_vza, local_sza, density, crown_height, crown_radius)
            budget = np.maximum(0.0, 1.0 - np.sum(direct, axis=0))
            total = np.sum(scatter, axis=0)
            scatter = scatter * np.minimum(1.0, np.divide(budget, total, out=np.ones_like(total), where=total > 0))
        shadow = np.asarray(local_sun_cos).reshape(-1) <= 0.0
        scatter[1, shadow] += scatter[0, shadow]
        scatter[3, shadow] += scatter[2, shadow]
        scatter[0, shadow] = scatter[2, shadow] = 0.0
        self.radiance = (radiance_direct(direct, component_radiance) +
                         radiance_scatter(scatter, component_radiance))
        if not np.all(np.isfinite(self.radiance)):
            raise FloatingPointError(f"Non-finite slope+{self.vegetation_model} radiance")
        return self.radiance if ifradiance == 1 else inv_planck(self.wavelength, self.radiance)
