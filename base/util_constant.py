import numpy as np

'''
util_constant
常量
以后常量尽量的用大写，多个名词间用下划线明确
'''

rd = np.pi/180.0
ndvi_max_avhrr = 0.95
ndvi_min_avhrr = 0.05
ndvi_max_modis = 0.95
ndvi_min_modis = 0.05

'''modis'''
band31_modis = 10.765
band32_modis = 12.021

'''viirs'''
band15_viirs = 10.710
band16_viirs = 11.834


'''ahi'''
Band1_AHI =    0.47
Band2_AHI  =   0.51
Band3_AHI   =  0.64
Band4_AHI    = 0.86
Band5_AHI    = 1.6
Band6_AHI   =  2.3
Band7_AHI   =  3.9
Band8_AHI   =  6.2
Band9_AHI   =  6.9
Band10_AHI  =  7.3
Band11_AHI  =  8.6
Band12_AHI  =  9.6
Band13_AHI  = 10.4
Band14_AHI  = 11.2
Band15_AHI  = 12.4
Band16_AHI  = 13.3

'''slstr'''

ns_nadir_g = 1500
nl_nadir_g = 1200
ns_obliq_g = 900
nl_obliq_g = 1200
ndvi_min_g = 0.05
ndvi_max_g = 0.99
wl8_g = 10.8
wl9_g = 12.0
# to resize for lat and lon
adjust_xmin_n_g = 0
adjust_xmax_n_g = 95 + 1
adjust_xmin_o_g = 36
adjust_xmax_o_g = 94 + 1
# to vary from oblique to nadir
adjust_xmin_no_g = 548
adjust_xmax_no_g = 548 + 900
adjust_ymin_no_g = 2
adjust_ymax_no_g = 1200

# aodvalue = 0.05 ### low
aodvalue = 0.12 ### normal
# aodvalue = 0.15 ### high

'''
util_constant
'''
rd = np.pi/180.0
planck_c1 = 11910.439340652
planck_c2 = 14388.291040407
temperature_threshold = 100
temperature_zero = 273.15


