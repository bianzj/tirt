import numpy as np


def sample():
     variables ={
         'lai_crown' : 0.5,
         'lai_veg' : 0.1,

         'hcr_crown':3.0,
         'rcr_crown':1.0,
         'hspot_crown':0.1,
         'std_crown':0.075,
         'G_crown':0.5,
         'CI_crown':1.0,
         'alg_crown':54,

         'hc_trunk':10.11,
         'dbt_trunk':0.3,

         'hspot_veg':0.1,
         'G_veg':0.5,
         'CI_veg':1.0,
         'alg_veg':54,

         'wavelength':10.5,
         'emissivity_soil':0.950,
         'emissivity_crown':0.975,
         'emissivity_veg':0.975,
         'emissivity_trunk':0.930,


         'view_zenith_angle':np.hstack([np.linspace(60,1,60),np.linspace(0,60,61)]),
         'view_azimuth_angle':np.hstack([np.repeat(0,61),np.repeat(180,60)]),
         'solar_zenith_angle':np.asarray([45]),
         'solar_azimuth_angle':np.asarray([0]),


         'temperature_soil_sunlit':320,
         'temperature_soil_shaded':305,
         'temperature_crown_sunlit':303,
         'temperature_crown_shaded':300,
         'temperature_trunk_sunlit':320,
         'temperature_trunk_shaded':305,
         'temperature_veg_sunlit':303,
         'temperature_veg_shaded':300,
     }

     return variables


def sample_layer():
    variables = {
        'lai_crown': 1.1,
        'lai_veg': 0.35,

        'hcr_crown': 6.0,
        'rcr_crown': 1.4,
        'hspot_crown': 0.1,
        'std_crown': 0.09,
        'G_crown': 0.5,
        'CI_crown': 1.0,
        'alg_crown': 54,

        'hc_trunk': 10.11,
        'dbt_trunk': 0.4,

        'hspot_veg': 0.1,
        'G_veg': 0.5,
        'CI_veg': 1.0,
        'alg_veg': 54,

        'wavelength': 10.5,
        'emissivity_soil': 0.950,
        'emissivity_crown': 0.985,
        'emissivity_veg': 0.985,
        'emissivity_trunk': 0.930,

        'view_zenith_angle': np.hstack([np.linspace(60, 1, 60), np.linspace(0, 60, 61)]),
        'view_azimuth_angle': np.hstack([np.repeat(0, 61), np.repeat(180, 60)]),
        'solar_zenith_angle': np.asarray([45]),
        'solar_azimuth_angle': np.asarray([0]),

        'temperature_soil_sunlit': 320,
        'temperature_soil_shaded': 305,
        'temperature_crown_sunlit': np.asarray([303, 303, 303, 303, 303]),
        'temperature_crown_shaded': np.asarray([300, 300, 300, 300, 300]),
        'temperature_trunk_sunlit': 318,
        'temperature_trunk_shaded': 305,
        'temperature_veg_sunlit': 303,
        'temperature_veg_shaded': 300,
    }

    return variables