"""Actual meteorology regression for the evening physiological transition."""
from pathlib import Path
from base.io import read_defined,read_meteo
from base.modelio import read_input
from base.pipeline import run_node


def test_real_half_hour_series_closes_through_evening_transition():
    root=Path(__file__).resolve().parents[1]
    config=read_input(root/'input_tirteb.csv')
    config.update(_input_base=str(root),surface_model='plane',vegetation_model='hom')
    settings=dict(read_defined(None,config))
    records=[r for r in read_meteo(root/config['meteo']) if 214<=r['time']<=215]
    settings['soil_temperature_c']=sum(r['air_temperature_c'] for r in records)/len(records)
    previous=None
    for i,meteo in enumerate(records):
        st=run_node(settings,meteo,i,2024,172,39.9,116.4,initial_state=previous,max_iter=50,utc_offset_hours=8.)
        assert st['closed'], (meteo['time'],st['energy_balance_residual_W_m2'],st['stability_residual_m_inv'])
        previous=st
