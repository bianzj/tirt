"""Urban radiation must obey blackbody equilibrium and whole-ground budgets."""
import math

import numpy as np
import pytest

from base.eb import SIGMASB
from base.evapo import _urban_component_incoming, _urban_longwave_inputs, urban_evapotranspiration
from base.scene_fluxes import urban_weights
from base import tirteb_radiation as tr


def configuration(**changes):
    config = dict(surface_model="urban", vegetation_model="hom", lai=2.0, hspot=0.1,
                  leaf_emissivity=.96, soil_emissivity=.91,
                  leaf_reflectance=.12, leaf_transmittance=.08, soil_reflectance=.2,
                  urban_emissivity_street=.91, urban_emissivity_wall=.87, urban_emissivity_roof=.9,
                  urban_absorption_street=.8, urban_absorption_wall=.7, urban_absorption_roof=.75,
                  urban_building_length=12., urban_building_width=6., urban_building_height=15.,
                  urban_density=.002, urban_building_azimuth=0.,
                  soil_heat_method="direct", urban_heat_method="direct")
    config.update(changes)
    return config


def thermal_state(**changes):
    state = dict(sza=40., saa=55., sunlit_surface=True,
                 Tss=300., Tsh=300., Tls=300., Tlh=300.,
                 Troof_sunlit=300., Troof_shaded=300., Twall_sunlit=300., Twall_shaded=300.,
                 raa_leaf=50., raa_roof=50., raa_wall=50., raa_street=50.)
    state.update(changes)
    return state


def sunlight(direct=700., diffuse=100.):
    return dict(directVrad_incident=direct, diffuseVrad_incident=diffuse,
                diffuseVrad_soil=diffuse, longwave_down=SIGMASB * 300.**4)


@pytest.mark.parametrize("lai", [0., 2.])
def test_isothermal_blackbody_environment_has_zero_component_longwave(lai):
    config = configuration(lai=lai)
    incoming = _urban_longwave_inputs(config, thermal_state(), SIGMASB * 300.**4)
    for component, emissivity, faces in (("roof", .9, 1.), ("wall", .87, 1.),
                                         ("street", .91, 1.), ("leaf", .96, 2.)):
        if component == "leaf" and lai == 0:
            continue
        for suffix in ("sunlit", "shaded"):
            assert incoming[f"{component}_{suffix}"] == pytest.approx(emissivity * faces * SIGMASB * 300.**4, abs=1e-8)


@pytest.mark.parametrize("changes", [{"urban_density": 0.}, {"urban_building_height": 0.}])
def test_no_walls_reduces_to_plane_canopy_longwave(changes):
    config = configuration(**changes)
    state = thermal_state(Tss=315., Tsh=295., Tls=310., Tlh=292.)
    actual = _urban_longwave_inputs(config, state, 340.)
    expected = tr.longwave(config, 340., (315., 295., 310., 292.), 40., 55., True, .91)
    assert [actual[k] for k in ("street_sunlit", "street_shaded", "leaf_sunlit", "leaf_shaded")] == pytest.approx(expected)
    assert np.isfinite(list(actual.values())).all()


def test_radiation_rotates_with_rectangular_building_bearing():
    a, b = thermal_state(saa=30.), thermal_state(saa=95.)
    first = _urban_component_incoming(configuration(), a, sunlight())
    rotated = _urban_component_incoming(configuration(urban_building_azimuth=65.), b, sunlight())
    for component in ("roof", "wall", "street"):
        assert first[component] == pytest.approx(rotated[component], abs=1e-8)


def test_overhead_sun_does_not_create_a_sunlit_vertical_wall():
    state = thermal_state(sza=0.)
    _urban_component_incoming(configuration(), state, sunlight())
    assert state["wall_sunlit_fraction"] == 0.
    assert state["SW_direct_wall_sunlit"] == 0.


@pytest.mark.parametrize("lai", [0., 2.])
def test_city_shortwave_absorption_and_escape_equal_input(lai):
    config, state = configuration(lai=lai), thermal_state()
    _urban_component_incoming(config, state, sunlight())
    soil_sun, soil_shade, leaf_sun, leaf_shade = tr.flux_weights(config, 40., 55.)
    state.update(soil_sunlit_fraction=soil_sun, lai_sunlit=leaf_sun, lai_shaded=leaf_shade)
    absorbed = 0.
    for name, weight in urban_weights(config, state).items():
        component = "leaf" if name.startswith("leaf") else name[:-6]
        suffix = "sunlit" if name.endswith("sunlit") else "shaded"
        absorbed += weight * (state.get(f"SW_direct_{component}_{suffix}", 0.)
                              + state.get(f"SW_diffuse_{component}_{suffix}", 0.)) * (lai if component == "leaf" else 1.)
    assert absorbed <= 800. + 1e-7
    assert absorbed + state.get("urban_SW_escaped", 0.) == pytest.approx(800., abs=1e-7)


def test_urban_flux_uses_same_canopy_and_street_longwave():
    config, state = configuration(), thermal_state()
    meteo = dict(pressure=1013., vapor_pressure=15., air_temperature_c=26.85)
    urban_evapotranspiration(config, meteo, state, sunlight(0., 0.))
    for component in ("roof", "wall", "street", "leaf"):
        for suffix in ("sunlit", "shaded"):
            assert state[f"N_{component}{suffix}"] == pytest.approx(0., abs=1e-7)


def test_zero_area_wall_slots_remain_finite_and_in_equilibrium():
    state = thermal_state()
    config = configuration(urban_density=0.)
    meteo = dict(pressure=1013., vapor_pressure=15., air_temperature_c=26.85)
    urban_evapotranspiration(config, meteo, state, sunlight(0., 0.))
    assert state["N_wallsunlit"] == pytest.approx(0., abs=1e-7)
    assert state["N_wallshaded"] == pytest.approx(0., abs=1e-7)


@pytest.mark.parametrize("sza", [40., 100.])
def test_nonisothermal_component_longwave_closes_at_city_boundary(sza):
    config = configuration()
    state = thermal_state(sza=sza, Tss=320., Tsh=295., Tls=312., Tlh=289.,
                          Twall_sunlit=330., Twall_shaded=296., Troof_sunlit=325., Troof_shaded=290.)
    incoming = _urban_longwave_inputs(config, state, 330.)
    state.update(incoming)
    weights = tr.flux_weights(config, sza, 55., sza < 90.)
    state.update(soil_sunlit_fraction=weights[0], lai_sunlit=weights[2], lai_shaded=weights[3])
    temperatures = dict(roof=(325., 290.), wall=(330., 296.), street=(320., 295.), leaf=(312., 289.))
    emissivity = dict(roof=.9, wall=.87, street=.91, leaf=.96)
    net = 0.
    for name, area in urban_weights(config, state).items():
        component = "leaf" if name.startswith("leaf") else name[:-6]
        index = 0 if name.endswith("sunlit") else 1
        suffix = "sunlit" if index == 0 else "shaded"
        outgoing = emissivity[component] * SIGMASB * temperatures[component][index]**4 * (2. if component == "leaf" else 1.)
        net += area * (incoming[f"{component}_{suffix}"] - outgoing) * (2. if component == "leaf" else 1.)
    assert net == pytest.approx(330. - incoming["urban_LW_escaped"], abs=1e-7)


@pytest.mark.parametrize("reflectance", [.2, 1.])
def test_spectral_shortwave_conserves_with_canopy_and_reflections(reflectance):
    config = configuration(leaf_reflectance_spectrum=[reflectance] * 2001,
                           leaf_transmittance_spectrum=[0.] * 2001,
                           urban_reflectance_roof_spectrum=[reflectance] * 2001,
                           urban_reflectance_wall_spectrum=[reflectance] * 2001,
                           urban_reflectance_street_spectrum=[reflectance] * 2001,
                           wavelengths=np.arange(400., 2401.))
    state = thermal_state()
    radiation = dict(sunlight(), directVrad_incident_spectrum=[700. / 2001.] * 2001,
                     diffuseVrad_incident_spectrum=[100. / 2001.] * 2001)
    _urban_component_incoming(config, state, radiation)
    weights = tr.flux_weights(config, 40., 55.)
    state.update(soil_sunlit_fraction=weights[0], lai_sunlit=weights[2], lai_shaded=weights[3])
    absorbed = 0.
    for name, weight in urban_weights(config, state).items():
        component = "leaf" if name.startswith("leaf") else name[:-6]
        suffix = "sunlit" if name.endswith("sunlit") else "shaded"
        absorbed += weight * (state[f"SW_direct_{component}_{suffix}"] + state[f"SW_diffuse_{component}_{suffix}"]) * (2. if component == "leaf" else 1.)
    assert absorbed + state["urban_SW_escaped"] == pytest.approx(800., abs=1e-7)
    if reflectance == 1.:
        assert absorbed == pytest.approx(0., abs=1e-7)
    else:
        assert state["directPrad_leaf"] > state["diffusePrad_leaf_shaded"]
