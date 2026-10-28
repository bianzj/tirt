import numpy as np
import scipy.integrate as sci
from .gap import *
from .hotspot import *
from .scatter import *
import matplotlib.pyplot as plt


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



def emissivity_scatter_analytical(lai, emis_l, emis_s, vza, sza, view_gap=None, hemi_gap=None, sunlit_volume=None):
    '''
    组分有效发射率，组分观测比例*组分发射率 + 多次散射项
    :param lai:  叶面积指数
    :param emis_l: 叶片发射率
    :param emis_s: 土壤发射率
    :param vza:    观测天顶角
    :param sza:    太阳天顶角
    :param view_gap: 实际冠层的视向间隙率；缺省为均质冠层
    :param hemi_gap: 实际冠层的半球间隙率；缺省为均质冠层
    :param sunlit_volume: 实际冠层的光照叶片比例；缺省为均质冠层
    :return:
    '''

    refl_s = 1-emis_s
    refl_l = 1-emis_l
    Mleaf_sunlit, Mleaf_shaded = multiple_scattering_analytical_sunlit(
        lai, vza, sza, refl_s, refl_l, bv=view_gap, M=hemi_gap, Vsunlit=sunlit_volume)
    number_angle = np.size(vza)
    eemis_ss = np.zeros(number_angle)
    eemis_sh = np.zeros(number_angle)
    eemis_ls =  Mleaf_sunlit
    eemis_lh =  Mleaf_shaded
    eemis_sum = eemis_lh + eemis_ls + eemis_sh + eemis_ss
    eemis = np.asarray([eemis_ss,eemis_sh,eemis_ls,eemis_lh])

    return eemis,eemis_sum


_MU_NODES, _MU_WEIGHTS = np.polynomial.legendre.leggauss(12)
_MU_NODES, _MU_WEIGHTS = 0.5 * (_MU_NODES + 1.0), 0.5 * _MU_WEIGHTS
_MU_ZENITH = np.degrees(np.arccos(_MU_NODES))
_HOM_LAI_TABLE = np.linspace(0.0, 30.0, 3001)
_HOM_HEMI_TABLE = np.asarray([2.0 * np.sum(_MU_WEIGHTS * _MU_NODES * gap_probability_hom_analytical(x, _MU_ZENITH))
                              for x in _HOM_LAI_TABLE])


def _equivalent_hemispherical_gap(hemi_numeric):
    """Map a cos-weighted hemispherical gap onto the fitted homogeneous closure exp(-0.825 LAI_eq)."""
    lai_eq = np.interp(-np.asarray(hemi_numeric, float), -_HOM_HEMI_TABLE, _HOM_LAI_TABLE)
    return gap_probability_hom_hemisphere_analytical(lai_eq)


def canopy_scatter(kind, lai, emis_l, emis_s, vza, sza, density=None, crown_height=None, crown_radius=None):
    """Leaf-source multiple scattering (4, N) using the actual canopy's view and hemispherical gaps.

    Homogeneous and row canopies keep the established closure. Crowns use the
    crown gap function, so the term vanishes as crown density goes to zero.
    """
    vza = np.atleast_1d(np.asarray(vza, float))
    sza = np.atleast_1d(np.asarray(sza, float))
    if lai <= 0.0 or (kind == "crown" and (density is None or density <= 0.0)):
        return np.zeros((4, max(vza.size, sza.size)))
    if kind != "crown":
        return np.asarray(emissivity_scatter_analytical(lai, emis_l, emis_s, vza, sza)[0], float).reshape(4, -1)
    view_gap = gap_probability_crown_analytical(lai, density, crown_height, crown_radius, vza)
    hemi_numeric = 2.0 * np.sum(_MU_WEIGHTS * _MU_NODES *
                                gap_probability_crown_analytical(lai, density, crown_height, crown_radius, _MU_ZENITH))
    hemi_gap = _equivalent_hemispherical_gap(hemi_numeric)
    return np.asarray(emissivity_scatter_analytical(lai, emis_l, emis_s, vza, sza, view_gap, hemi_gap)[0],
                      float).reshape(4, -1)




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