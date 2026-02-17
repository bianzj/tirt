import numpy as np
from old.proportion import *
from old.emissivity import *
from old.radiance import *

class Endmember_crown_hom():
    Esky = 0
    Esun = 0

    # lai_scene = 1.34
    # ifcrown = 0
    # lai_crown = 0.51
    # std_crown= 0.068
    # hcr_crown = 6.5
    # rcr_crown = 1.61
    # G_crown = 0.5
    # hspot_crown = 0.2
    # CI_crown = 1.0
    # alg_crown = 53
    #
    # lai_trunk = 0.5
    # dbh_trunk = 0.39
    # hc_trunk = 12.5
    #
    # hspot_veg = 0.1
    # lai_veg = 0.48
    # agl_veg = 53
    # CI_veg = 1.0
    #
    # wavelength = 10.5
    # number_angle = 1
    # solar_zenith_angle = np.asarray([0])
    # view_zenith_angle = np.asarray([0])
    # solar_azimuth_angle = np.asarray([0])
    # view_azimuth_angle = np.asarray([0])
    #
    # emissivity_crown = 0.985
    # emissivity_trunk = 0.930
    # emissivity_veg = 0.985
    # emissivity_soil = 0.950
    #
    # temperature_soil_sunlit = 320
    # temperature_soil_shaded = 305
    # temperature_crown_sunlit = 303
    # temperature_crown_shaded = 300
    # temperature_trunk_sunlit = 318
    # temperature_trunk_shaded = 305
    # temperature_veg_sunlit = 303
    # temperature_veg_shaded = 300


    def set_structure(self,variables):
        self.lai_crown = variables['lai_crown']
        self.lai_veg = variables['lai_veg']
        self.std_crown = variables['std_crown']
        self.hcr_crown = variables['hcr_crown']
        self.rcr_crown = variables['rcr_crown']
        self.hspot_crown = variables['hspot_crown']
        self.hspot_veg = variables['hspot_veg']
        self.dbh_trunk = variables['dbt_trunk']
        self.hc_trunk = variables['hc_trunk']
        self.CI_crown = variables['CI_crown']
        self.CI_veg = variables['CI_veg']
        self.alg_crown = variables['alg_crown']
        self.alg_veg = variables['alg_veg']
        self.lai_trunk = self.hc_trunk * self.dbh_trunk

    def set_thermal(self,variables):
        self.temperature_soil_sunlit = variables['temperature_soil_sunlit']
        self.temperature_soil_shaded = variables['temperature_soil_shaded']
        self.temperature_crown_sunlit = variables['temperature_crown_sunlit']
        self.temperature_crown_shaded = variables['temperature_crown_shaded']
        self.temperature_veg_shaded = variables['temperature_veg_sunlit']
        self.temperature_veg_sunlit = variables['temperature_veg_shaded']
        self.temperature_trunk_sunlit = variables['temperature_trunk_sunlit']
        self.temperature_trunk_shaded = variables['temperature_trunk_shaded']

    def set_optical(self,variables):
        self.wavelength =  variables['wavelength']
        self.emissivity_soil =  variables['emissivity_soil']
        self.emissivity_crown =  variables['emissivity_crown']
        self.emissivity_trunk =  variables['emissivity_trunk']
        self.emissivity_veg =  variables['emissivity_veg']


    def set_angle(self,variables):
        view_zenith_angle = variables['view_zenith_angle']
        view_azimuth_angle = variables['view_azimuth_angle']
        solar_zenith_angle = variables['solar_zenith_angle']
        solar_azimuth_angle = variables['solar_azimuth_angle']

        self.number_angle = 1
        self.number_angle = np.max([self.number_angle, np.size(view_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_zenith_angle)])
        self.number_angle = np.max([self.number_angle, np.size(view_azimuth_angle)])
        self.number_angle = np.max([self.number_angle, np.size(solar_azimuth_angle)])

        if (self.number_angle > 1) :
            view_zenith_angle = np.resize(view_zenith_angle, self.number_angle)
        if (self.number_angle > 1) :
            solar_zenith_angle = np.resize(solar_zenith_angle, self.number_angle)
        if (self.number_angle > 1) :
            view_azimuth_angle = np.resize(view_azimuth_angle, self.number_angle)
        if (self.number_angle > 1):
            solar_azimuth_angle = np.resize(solar_azimuth_angle, self.number_angle)

        self.view_zenith_angle = np.asarray(view_zenith_angle)
        self.solar_zenith_angle = np.asarray(solar_zenith_angle)
        self.view_azimuth_angle = np.asarray(view_azimuth_angle)
        self.solar_azimuth_angle = np.asarray(solar_azimuth_angle)

    def run(self,ifradiance = 0):


        std_crown = self.std_crown
        hcr_crown = self.hcr_crown
        rcr_crown = self.rcr_crown



        lai_crown = self.lai_crown
        hspot_crown = self.hspot_crown
        hc_trunk = self.hc_trunk
        dbh_trunk = self.dbh_trunk
        std_trunk = std_crown
        lai_unveg = self.lai_veg
        hspot_unveg = self.hspot_veg
        vza = self.view_zenith_angle
        sza = self.solar_zenith_angle
        vaa = self.view_azimuth_angle
        saa = self.solar_azimuth_angle
        ec = self.emissivity_crown
        es = self.emissivity_soil
        ev = self.emissivity_veg
        et = self.emissivity_trunk
        Tss = self.temperature_soil_sunlit
        Tsh = self.temperature_soil_shaded
        Tcs = self.temperature_crown_sunlit
        Tch = self.temperature_crown_shaded
        Tvs = self.temperature_veg_sunlit
        Tvh = self.temperature_veg_shaded
        Tts = self.temperature_trunk_sunlit
        Tth = self.temperature_trunk_shaded

        hc_crown = hcr_crown * 2 + hc_trunk
        lai_trunk = hc_trunk * dbh_trunk

        Ecom = np.asarray([es, es, ec, ec,ev,ev,et,et])
        Tcom = np.asarray([Tss, Tsh, Tcs, Tch,Tvs,Tvh,Tts,Tth])
        Rcom = planck(Tcom)


        Rcom_scatter = 0
        Pcom = proportion_bidirectional_crown_voxel_one_endmember(lai_crown, std_crown, hspot_crown, hc_crown, hcr_crown, rcr_crown,
                                                               hc_trunk, dbh_trunk, lai_unveg, hspot_unveg, vza, sza, vaa, saa)
        Ecom_direct, Ecom_direct_sum = emissivity_direct(Pcom, Ecom)
        Rcom_direct = radiance_direct(Ecom_direct, Rcom)

        # Pcomsum = np.sum(Pcom,axis = 0)
        # plt.plot(Pcomsum)
        # plt.show()


        ###################################
        ##### 光谱不变多次散射
        ###################################
        Ecom_scatter,Ecom_scatter_sum = \
            emissivity_scatter_crown_endmember(lai_crown, std_crown, hc_crown, hcr_crown, rcr_crown, lai_unveg, hc_trunk, dbh_trunk,
                                               ec, et, ev, es,
                                               vza, sza, vaa, saa)
        Rcom_scatter = radiance_scatter(Ecom_scatter,Rcom)

        ###################################
        ##### SAIL多次散射
        ###################################
        # Rcom_scatter = radiance_scatter_crown_endmember_sail(Pcom,
        #                                           lai_crown,lai_trunk,lai_unveg,std_crown,hcr_crown,rcr_crown,
        #                                           ec,et,ev,es,Tss,Tsh,Tcs,Tch,Tvs,Tvh,Tts,Tth,
        #                                           vza,sza,vaa)
        # self.radiance = radiance_hom_endmember(Ess, Esh, Ecs, Ech, Evs, Evh, Ets, Eth,
        #                                        Tss, Tsh, Tcs, Tch, Tvs, Tvh, Tts, Tth)

        self.radiance = Rcom_direct + Rcom_scatter

        if ifradiance == 1:
            return self.radiance
        else:
            return inv_planck(self.radiance, self.wavelength)