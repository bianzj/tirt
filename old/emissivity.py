import numpy as np
import scipy.integrate as sci
from gap import *
from hotspot import *
from utils import *
from scatter import *
import matplotlib.pyplot as plt
import quadpy


def emissivity_direct(Pcom,Ecom):
    number_component = np.size(Ecom)
    number_angle = np.size(Pcom[0])
    # Pcomnew = np.zeros([number_component,number_angle])
    Pcomnew = []
    Pcomnew_sum = 0
    for k in range(number_component):
        Pcomtemp = Pcom[k]
        Pcomnew.append(Pcomtemp * Ecom[k])
        if np.size(Pcomtemp) > number_angle:
            Pcomnew_sum_temp = np.sum(Pcomtemp * Ecom[k],axis=0)
        else:
            Pcomnew_sum_temp = Pcomtemp * Ecom[k]
        Pcomnew_sum = Pcomnew_sum + Pcomnew_sum_temp

    return Pcomnew,np.asarray(Pcomnew_sum)



def emissivity_scatter_analytical(lai, emis_l, emis_s, vza, sza):
    '''
    组分有效发射率，组分观测比例*组分发射率 + 多次散射项
    :param lai:  叶面积指数
    :param Psoil:  土壤观测比例
    :param Pleaf:  叶片观测比例
    :param Psoil_sunlit:  光照土壤观测比例
    :param Pleaf_sunlit:  光照叶片观测比例
    :param emis_l: 叶片发射率
    :param emis_s: 土壤发射率
    :param vza:    观测天顶角
    :param sza:    太阳天顶角
    :return:
    '''

    refl_s = 1-emis_s
    refl_l = 1-emis_l
    Mleaf_sunlit, Mleaf_shaded = multiple_scattering_analytical_sunlit(
        lai, vza, sza, refl_s, refl_l)
    number_angle = np.size(vza)
    eemis_ss = np.zeros(number_angle)
    eemis_sh = np.zeros(number_angle)
    eemis_ls =  Mleaf_sunlit
    eemis_lh =  Mleaf_shaded
    eemis_sum = eemis_lh + eemis_ls + eemis_sh + eemis_ss
    eemis = np.asarray([eemis_ss,eemis_sh,eemis_ls,eemis_lh])

    return eemis,eemis_sum



def emissivity_scatter_hom_endmember(lai_crown, lai_unveg, hc_trunk, dbh_trunk, std_trunk,
                    ec, et, ev, es,vza0, sza0, vaa0, saa0=0):
    '''
    组分有效发射率，组分观测比例*组分发射率 + 多次散射项
    :param lai:  叶面积指数
    :param Psoil:  土壤观测比例
    :param Pleaf:  叶片观测比例
    :param Psoil_sunlit:  光照土壤观测比例
    :param Pleaf_sunlit:  光照叶片观测比例
    :param emis_l: 叶片发射率
    :param emis_s: 土壤发射率
    :param vza:    观测天顶角
    :param sza:    太阳天顶角
    :return:
    '''

    emiss,emish,emics,emich,emivs,emivh,emits,emith,memis = \
        multiple_scattering_hom_spectral_invariance_endmember(
        lai_crown, lai_unveg, hc_trunk, dbh_trunk, std_trunk,
        ec, et, ev, es,
        vza0, sza0, vaa0, saa0)

    eemis_ = np.asarray([ emiss,emish,emics,emich,emivs,emivh,emits,emith])

    # plt.plot(eemis_ss,'r-')
    # plt.plot(eemis_sh,'b-')
    # plt.plot(eemis_cs,'ro-')
    # plt.plot(eemis_ch,'bo-')
    # plt.plot(eemis_vs,'r^-')
    # plt.plot(eemis_vh,'b^-')
    # plt.plot(eemis_ts,'rs-')
    # plt.plot(eemis_th,'bs-')
    # plt.plot(eemis,'k-')
    # plt.show()
    return eemis_, memis


def emissivity_scatter_crown_endmember(lai_crown,std_crown, hc,hcr,rcr,lai_unveg,hc_trunk,dbh_trunk,
                    ec,et,ev,es,vza0,sza0,vaa0,saa0=0):
    '''
    组分有效发射率，组分观测比例*组分发射率 + 多次散射项
    :param lai:  叶面积指数
    :param Psoil:  土壤观测比例
    :param Pleaf:  叶片观测比例
    :param Psoil_sunlit:  光照土壤观测比例
    :param Pleaf_sunlit:  光照叶片观测比例
    :param emis_l: 叶片发射率
    :param emis_s: 土壤发射率
    :param vza:    观测天顶角
    :param sza:    太阳天顶角
    :return:
    '''


    emiss,emish,emics,emich,emivs,emivh,emits,emith,memis = \
        multiple_scattering_crown_spectral_invariance_endmember(
        lai_crown,std_crown,hc,hcr,rcr, lai_unveg, hc_trunk, dbh_trunk,
        ec, et, ev, es,
        vza0, sza0, vaa0, saa0)

    emis_ = np.asarray([emiss,emish,emics,emich,emivs,emivh,emits,emith])


    return emis_,memis


def emissivity_scatter_hom_voxel_one(lai_crown, ec, es, vza0, sza0, vaa0, saa0=0, G_crown=0.5):
    number_angle = np.size(vza0)
    eemis_ss = np.zeros(number_angle)
    eemis_sh = np.zeros(number_angle)
    emics, emich, memis = multiple_scattering_hom_spectral_invariance(lai_crown, ec, es, vza0, sza0, vaa0, saa0, G_crown)
    eemis = np.asarray([eemis_ss,eemis_sh,emics,emich])

    return eemis,memis

def emissivity_scatter_row_voxel_one(lai_crown, row_width,row_blank,row_height,
                                                  ec, es,vza0, sza0, vaa0, saa0,raa0):
    number_angle = np.size(vza0)
    eemis_ss = np.zeros(number_angle)
    eemis_sh = np.zeros(number_angle)
    emics, emich, memis = multiple_scattering_row_spectral_invariance(lai_crown, row_width,row_blank,row_height,
                                                  ec, es,vza0, sza0, vaa0, saa0,raa0)
    eemis = np.asarray([eemis_ss,eemis_sh,emics,emich])

    return eemis,memis


def emissivity_scatter_crown_voxel_one(lai_crown, std_crown, hc, hcr, rcr,
                                       ec,  es,vza0, sza0, vaa0, saa0=0, G_crown=0.5):
    number_angle = np.size(vza0)
    eemis_ss = np.zeros(number_angle)
    eemis_sh = np.zeros(number_angle)
    emics, emich, memis = multiple_scattering_crown_spectral_invariance(
        lai_crown, std_crown, hc, hcr, rcr, ec,  es,vza0, sza0, vaa0, saa0, G_crown)
    eemis = np.asarray([eemis_ss,eemis_sh,emics,emich])

    return eemis,memis