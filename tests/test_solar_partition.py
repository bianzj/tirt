import math
import pytest
from base.radiation import _shortwave_diffuse_fraction
from base.tirteb_radiation import prepare_shortwave


def test_default_partition_preserved_when_direct_beam_is_admissible():
    assert _shortwave_diffuse_fraction({'shortwave_diffuse_fraction':.2},800.,30.) == .2


def test_grazing_sun_does_not_create_super_solar_direct_normal_irradiance():
    settings={'shortwave_diffuse_fraction':.2,'solar_day_of_year':214}
    fraction=_shortwave_diffuse_fraction(settings,131.1,88.93829568731898)
    normal=131.1*(1-fraction)/math.cos(math.radians(88.93829568731898))
    maximum=1361.*(1+.033*math.cos(2*math.pi*214/365.))
    assert normal <= maximum*(1+1e-12)
    assert fraction > .8


def test_below_horizon_positive_irradiance_is_diffuse():
    assert _shortwave_diffuse_fraction({'shortwave_diffuse_fraction':.2},10.,95.) == 1.


def test_small_positive_diffuse_irradiance_is_not_discarded():
    settings=dict(surface_model='plane',vegetation_model='bare',lai=0.,
                  leaf_reflectance=.015,leaf_transmittance=0.,soil_reflectance=.1,
                  shortwave_diffuse_fraction=1.)
    result=prepare_shortwave(settings,dict(shortwave=5.),85.,0.)
    assert result['diffuseVrad_soil_shaded'] == pytest.approx(4.5)


def test_direct_beam_is_capped_by_meinel_clear_sky_dni_and_total_is_preserved():
    from base.radiation import _solar_direct_normal_limit
    settings={'shortwave_diffuse_fraction':.2,'solar_day_of_year':214}
    for sza, total in ((80., 300.), (85., 150.), (30., 900.)):
        fraction=_shortwave_diffuse_fraction(settings,total,sza)
        cosine=math.cos(math.radians(sza))
        direct, diffuse = total*(1-fraction), total*fraction
        assert direct+diffuse == pytest.approx(total)
        assert direct/cosine <= _solar_direct_normal_limit(settings,sza)*(1+1e-12)
    assert _shortwave_diffuse_fraction(settings,900.,30.) == .2
    air_mass=1/(math.cos(math.radians(60.))+.50572*(96.07995-60.)**-1.6364)
    e0=1+.033*math.cos(2*math.pi*214/365.)
    assert _solar_direct_normal_limit(settings,60.) == pytest.approx(1367*e0*.7**(air_mass**.678))
    assert _solar_direct_normal_limit(settings,95.) == 0.


def test_clearness_index_above_one_is_flagged_and_summarised():
    from base.radiation import clearness_index, clearness_index_warning
    settings=dict(surface_model='plane',vegetation_model='bare',lai=0.,leaf_reflectance=.015,
                  leaf_transmittance=0.,soil_reflectance=.1,shortwave_diffuse_fraction=.2,solar_day_of_year=214)
    e0=1+.033*math.cos(2*math.pi*214/365.)
    assert clearness_index(settings,600.,60.) == pytest.approx(600/(1361*e0*.5))
    bright=prepare_shortwave(settings,dict(shortwave=300.),80.,270.)
    normal=prepare_shortwave(settings,dict(shortwave=100.),80.,270.)
    assert bright['clearness_index_exceeds_one'] and not normal['clearness_index_exceeds_one']
    states=[dict(time=214.76,**bright),dict(time=214.5,**normal)]
    message=clearness_index_warning(states)
    assert '1/2' in message and '214.7600' in message
    assert clearness_index_warning([dict(time=214.5,**normal)]) is None
