"""Exercise actual temporal calculations through either source or frozen HTTP APIs."""
import json
import math

METEO = '2\n214.5 25 15 1013 2 700 350\n214.5208 26 15 1013 2 750 355\n'


def check_dynamic_api(request):
    reports = []
    for surface, vegetation in (('plane', 'hom'), ('plane', 'row'), ('plane', 'crown'),
                                ('slope', 'crown'), ('terrain', 'crown'), ('urban', 'crown')):
        overrides = dict(surface_model=surface, vegetation_model=vegetation, lai=2,
                         canopy_height=3, crown_rad_a=2, crown_rad_b=2, crown_height=3,
                         crown_density=.02, slope_angle=20, slope_aspect=110,
                         terrain_height=100, terrain_radius=100, terrain_density=.00001,
                         start_time=214.5, end_time=214.6, observation_time_mode='all',
                         spectral_source='direct', wavelengths='0.67;0.86;10.5',
                         geometry_pair_mode='zip', vza='15;30', vaa='77;160',
                         emissivity_leaf=.975, emissivity_soil=.955, output_intermediate='yes')
        data = json.loads(request('/api/run_time', dict(overrides=overrides, meteo_text=METEO))[1])
        assert data['ok'] and data['count'] == 2, data
        assert data['observation_count'] == 4 and data['intermediate_count'] == 2
        assert {(r['vza'], r['vaa']) for r in data['observations']} == {(15, 77), (30, 160)}
        for row in data['observations']:
            for key in ('0p67', '0p86', '10p5'):
                assert row[key] is not None and math.isfinite(row[key]), (surface, vegetation, row)
            assert 150 < row['10p5'] < 400
        for row in data['series']:
            assert all(row[k] is not None and math.isfinite(row[k]) for k in ('Tv_C', 'LE', 'H', 'G'))
        reports.append(dict(surface=surface, vegetation=vegetation, time_steps=data['count'],
                            observations=data['observation_count'],
                            closed_steps=sum(r['closed'] == 'True' for r in data['series']),
                            max_residual_W_m2=max(r['energy_balance_residual_W_m2'] for r in data['series'])))
        print('Dynamic API passed:', surface, vegetation, flush=True)
    # Uploads must override defaults, including view pairs and wavelength files.
    data = json.loads(request('/api/run_time', dict(overrides=dict(overrides, surface_model='plane',
                              vegetation_model='hom'), meteo_text=METEO, files={
        'directions': '22 123\n', 'wavelengths': '0.67\n0.86\n10.5\n',
        'vegetation_spectrum': '400 .1 .1\n860 .4 .4\n2400 .1 .1\n',
        'soil_spectrum': '400 .08\n860 .12\n2400 .2\n'}))[1])
    assert data['observation_count'] == 2, data
    assert {(r['vza'], r['vaa']) for r in data['observations']} == {(22, 123)}
    assert any(abs(a['0p86'] - b['0p86']) > 1e-5
               for a, b in zip(data['observations'], json.loads(request('/api/run_time', dict(
                   overrides=dict(overrides, surface_model='plane', vegetation_model='hom'),
                   meteo_text=METEO, files={'directions': '22 123\n'}))[1])['observations']))
    return dict(cases=reports, uploads_passed=True, manual_direction_pairs_passed=True)
