import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from base.aerodynamic import _stability_functions
from base.aerodynamic import aerodynamic_resistance
from base.io import read_defined
from base.modelio import read_input
from base.pipeline import run_node
from base.eb import saturation_vapor_pressure, SIGMASB
import tirteb_run


@pytest.fixture(scope='module')
def settings(tmp_path_factory):
    root = Path(__file__).resolve().parents[1]
    config = read_input(root / 'input_tirteb.csv')
    config.update(_input_base=str(root), surface_model='plane', vegetation_model='bare',
                  lai='0', soil_temperature_c='25', soil_heat_method='direct')
    empty = tmp_path_factory.mktemp('optics') / 'coefficients.csv'
    empty.write_text('',encoding='utf-8')
    config['optical_coefficients_file'] = str(empty)
    return dict(read_defined(None, config))


def test_stable_similarity_is_not_neutral():
    psim, psih, _ = _stability_functions(100.)
    assert psim(10.) == pytest.approx(-0.5)
    assert psih(10.) == pytest.approx(-0.5)


def test_small_negative_obukhov_length_has_unstable_correction():
    psim, psih, _ = _stability_functions(-2.)
    assert psim(1.) > 0
    assert psih(1.) > 0


def test_wall_flux_is_per_horizontal_ground(settings):
    st = dict(N_wallsunlit=1., N_wallshaded=1., soil_sunlit_fraction=1.,
              soil_shaded_fraction=0., lai_sunlit=0., lai_shaded=0.)
    urban = dict(settings, urban_density=.0015, urban_building_length=10.,
                 urban_building_width=10., urban_building_height=20.)
    assert tirteb_run._urban_surface_totals(urban, st)['Rn'] == pytest.approx(1.2)


@pytest.mark.parametrize('shortwave,longwave,initial,sign', [(800.,400.,25.,-1), (0.,250.,15.,1)])
def test_temperature_and_stability_converge_together(settings, shortwave, longwave, initial, sign):
    meteo = dict(time=214.5, air_temperature_c=25., vapor_pressure=15., pressure=1013.,
                 wind_speed=3., shortwave=shortwave, longwave=longwave)
    s = dict(settings, soil_temperature_c=initial)
    state = run_node(s, meteo, 0, 2024, 214, 39.9,116.4, max_iter=150,utc_offset_hours=8.)
    assert state['L'] * sign > 0
    assert abs(state['L']) < 1e8
    assert state['closed']
    assert state['stability_closed']
    assert state['stability_residual_m_inv'] < 1e-5
    assert state['dt_seconds'] == s['dtime']


def test_neutral_state_is_a_solution_instead_of_a_forced_default(settings):
    meteo = dict(time=214.5, air_temperature_c=26.85,
                 vapor_pressure=saturation_vapor_pressure(26.85), pressure=1013.,
                 wind_speed=2., shortwave=0., longwave=SIGMASB*300.**4)
    st = run_node(dict(settings, soil_temperature_c=26.85), meteo,0,2024,214,39.9,116.4,
                  max_iter=100,utc_offset_hours=8.)
    assert st['closed'] and st['stability_closed']
    assert abs(st['scene_sensible_heat_W_m2']) < 1e-8
    assert st['stability_residual_m_inv'] < 1e-8


def test_stable_warm_start_can_cross_to_unstable(settings):
    meteo = dict(time=214.5, air_temperature_c=25., vapor_pressure=15.,pressure=1013.,
                 wind_speed=2.,shortwave=800.,longwave=400.)
    st = run_node(settings, meteo,0,2024,214,39.9,116.4,initial_state={'L':0.5},
                  max_iter=150,utc_offset_hours=8.)
    assert st['closed'] and st['L'] < 0.


def test_very_stable_bulk_transfer_does_not_return_to_neutral(settings):
    meteo = dict(wind_speed=1.,air_temperature_c=25.)
    states = []
    for length in (float('inf'),15.,0.01,0.001):
        st = {'L':length}
        aerodynamic_resistance(dict(settings, max_aerodynamic_resistance=0.),meteo,st)
        states.append(st)
    assert states[2]['raa_soil'] >= states[1]['raa_soil']
    assert states[3]['raa_soil'] >= states[2]['raa_soil']
    assert states[3]['raa_soil'] > states[0]['raa_soil']


def test_each_resistance_segment_is_capped_as_in_scope(settings):
    meteo = dict(wind_speed=0.05, air_temperature_c=25.)
    canopy = dict(settings, vegetation_model='hom', lai=2.)
    for model, keys, segments in ((dict(settings), ('raa_soil',), 1), (canopy, ('raa_leaf', 'raa_soil'), 2)):
        free, capped = {'L': 0.01}, {'L': 0.01}
        aerodynamic_resistance(dict(model, max_aerodynamic_resistance=0.), meteo, free)
        aerodynamic_resistance(model, meteo, capped)
        assert settings['max_aerodynamic_resistance'] == 400.
        for key in keys:
            assert free[key] > segments*400.
            assert capped[key] <= segments*400. + 1e-9
    windy = {'L': float('inf')}
    calm = {'L': float('inf')}
    aerodynamic_resistance(canopy, dict(wind_speed=3., air_temperature_c=25.), windy)
    aerodynamic_resistance(dict(canopy, max_aerodynamic_resistance=0.), dict(wind_speed=3., air_temperature_c=25.), calm)
    assert windy['raa_leaf'] == pytest.approx(calm['raa_leaf'])


def test_street_conduction_commits_the_material_used_by_its_flux(settings):
    from base.evapo import _depth_differential,update_thermal_history
    from base.heat_conduction import conduction_profile
    st=dict(Tss=310.,Tsh=300.,Troof_sunlit=300.,Troof_shaded=300.,
            Twall_sunlit=300.,Twall_shaded=300.,dt_seconds=1800.,
            depth_profile_street_sunlit=[300.,300.,300.,300.])
    expected,_,_=conduction_profile(settings,st,310.,'depth_profile_street_sunlit',urban=False)
    flux=_depth_differential(settings,st,310.,'depth_profile_street_sunlit',urban=False)
    update_thermal_history(settings,st)
    assert st['depth_profile_street_sunlit'] == pytest.approx(expected)
    from base.heat_conduction import column_layers, thermal_properties
    layers = column_layers()
    widths = [layers[0]/2]+[(a+b)/2 for a, b in zip(layers, layers[1:])]+[layers[-1]/2]
    _, capacity = thermal_properties(settings)
    stored = capacity*sum(w*(t-300.) for w, t in zip(widths, expected))
    assert flux*1800. == pytest.approx(stored)


def test_air_density_uses_pressure_and_virtual_temperature():
    from base.eb import air_density
    rho = air_density(dict(pressure=837., vapor_pressure=15., air_temperature_c=30.))
    tv = 303.15/(1-0.378*15./837.)
    assert rho == pytest.approx(83700./(287.05*tv))
    assert 0.94 < rho < 0.97


def test_fluxes_and_stability_use_step_air_density(settings):
    from base.eb import CP
    meteo = dict(time=214.5, air_temperature_c=30., vapor_pressure=15., pressure=837.,
                 wind_speed=3., shortwave=800., longwave=400.)
    st = run_node(settings, meteo, 0, 2024, 214, 39.9, 116.4, max_iter=150, utc_offset_hours=8.)
    rho = st['air_density_kg_m3']
    assert rho == pytest.approx(83700./(287.05*303.15/(1-0.378*15./837.)))
    assert st['H_soilsunlit'] == pytest.approx(rho*CP/st['raa_soil']*(st['Tss']-303.15))
    target = -0.4*9.81*st['scene_sensible_heat_W_m2']/(rho*CP*303.15*st['ustar']**3)
    assert 1/st['L'] == pytest.approx(target, rel=1e-5)


def test_every_step_reports_convergence_and_stability_limit(settings):
    meteo = dict(time=214.5, air_temperature_c=25., vapor_pressure=15., pressure=1013.,
                 wind_speed=1., shortwave=0., longwave=300.)
    st = run_node(dict(settings, soil_temperature_c=15.), meteo, 0, 2024, 214, 39.9, 116.4,
                  max_iter=150, utc_offset_hours=8.)
    for key in ('closed', 'energy_balance_residual_W_m2', 'fallback_solver_used', 'stability_limited'):
        assert key in st
    assert isinstance(st['fallback_solver_used'], bool)
    assert st['stability_limited'] == (st['stability_reference_height']/st['L'] > 1.)
    limited, free = {'L': 0.5}, {'L': 1000.}
    aerodynamic_resistance(settings, dict(wind_speed=1., air_temperature_c=25.), limited)
    aerodynamic_resistance(settings, dict(wind_speed=1., air_temperature_c=25.), free)
    assert limited['stability_limited'] and not free['stability_limited']


def test_single_slope_stability_uses_flux_per_horizontal_area(settings):
    from base.scene_fluxes import horizontal_area_factor, scene_totals
    slope = dict(settings, surface_model='terrain', terrain_type=1, terrain_slope=40., terrain_aspect=180.)
    assert horizontal_area_factor(slope) == pytest.approx(1/math.cos(math.radians(40.)))
    assert horizontal_area_factor(dict(slope, terrain_type=2)) == 1.
    assert horizontal_area_factor(dict(settings, surface_model='plane')) == 1.
    meteo = dict(time=214.5, air_temperature_c=25., vapor_pressure=15., pressure=1013.,
                 wind_speed=3., shortwave=800., longwave=400.)
    st = run_node(slope, meteo, 0, 2024, 214, 39.9, 116.4, max_iter=150, utc_offset_hours=8.)
    slope_h = scene_totals(slope, st)['H']
    assert st['scene_sensible_heat_W_m2'] == pytest.approx(slope_h/math.cos(math.radians(40.)))
    assert st['stability_closed']


def test_sparse_canopy_exchange_converges_continuously_to_bare_soil(settings):
    meteo = dict(time=214.5, air_temperature_c=25., vapor_pressure=15., pressure=1013.,
                 wind_speed=3., shortwave=800., longwave=400.)
    bare = run_node(settings, meteo, 0, 2024, 214, 39.9, 116.4, max_iter=150, utc_offset_hours=8.)
    gaps = {key: [] for key in ('raa_soil', 'ustar', 'Tss')}
    for lai in (0.01, 1e-4, 1e-6):
        st = run_node(dict(settings, vegetation_model='hom', lai=lai), meteo, 0, 2024, 214, 39.9, 116.4,
                      max_iter=150, utc_offset_hours=8.)
        assert st['closed']
        for key in gaps:
            gaps[key].append(abs(st[key]-bare[key]))
    for key, values in gaps.items():
        assert values[0] >= values[1] >= values[2], (key, values)
        assert values[-1] <= 1e-3*max(abs(bare[key]), 1.), (key, values)
    sparse, dense = {'L': float('inf')}, {'L': float('inf')}
    aerodynamic_resistance(dict(settings, lai=1e-4), dict(wind_speed=3., air_temperature_c=25.), sparse)
    aerodynamic_resistance(dict(settings, lai=2.), dict(wind_speed=3., air_temperature_c=25.), dense)
    assert sparse['aerodynamic_roughness_length'] == pytest.approx(0.01, rel=1e-3)
    assert sparse['aerodynamic_displacement'] < 1e-6
    assert dense['aerodynamic_canopy_weight'] == 1.


def test_min_wind_speed_is_the_only_wind_floor(settings, tmp_path):
    from base.io import read_meteo
    meteo_file = tmp_path / 'meteo.txt'
    meteo_file.write_text('1\n214.5 25 15 1013 0.08 0 400\n', encoding='utf-8')
    assert read_meteo(meteo_file)[0]['wind_speed'] == pytest.approx(0.08)
    assert settings['min_wind_speed'] == pytest.approx(0.05)

    def ustar(wind, floor):
        state = {'L': float('inf')}
        aerodynamic_resistance(dict(settings, min_wind_speed=floor), dict(wind_speed=wind, air_temperature_c=25.), state)
        return state['ustar']

    assert ustar(0.08, 0.05) < ustar(0.15, 0.05) < ustar(0.3, 0.05)
    assert ustar(0.08, 1.0) == pytest.approx(ustar(1.0, 0.05))
