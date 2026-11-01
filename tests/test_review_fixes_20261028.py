"""Regression tests for the 2026-10-28 code-review fixes."""
import math
import warnings

import numpy as np
import pytest

import run
from base import eb
from base.inputio import read_input
from base.leafopt import fluspect
from base.path import crosscutting_ellipsoid
from base.pipeline import _ambient_vapour_pressure
from base.terrain_rt import compound_terrain_proportions


@pytest.mark.parametrize("za", [20.0, 50.0, 70.0])
def test_ellipsoid_branch_switch_is_at_silhouette_tangent_height(za):
    """z0 = b / sqrt(a^2/(b tan)^2 + 1); the old code used a instead of a^2."""
    hc, hcr, rcr = 10.0, 2.0, 3.0
    zk = np.linspace(hc - 2 * hcr + 1e-3, hc - 1e-3, 8001)
    area, _ = crosscutting_ellipsoid(zk, np.full(zk.size, za), hc, hcr, rcr)
    switch = zk[int(np.argmax(np.abs(np.diff(area))))] - (hc - hcr)
    z0 = hcr / math.sqrt(rcr ** 2 / (hcr * math.tan(math.radians(za))) ** 2 + 1.0)
    assert switch == pytest.approx(z0, abs=2 * (zk[1] - zk[0]))


def _lit_vegetation(sza, saa, vza, vaa, sunlit):
    return (0.6, 0.4, 0.0, 0.0)


def test_compound_terrain_mutual_shadow_changes_sunlit_share():
    shares = []
    for sza in (10.0, 40.0, 70.0):
        result = compound_terrain_proportions({}, sza, 0.0, 0.0, 0.0, _lit_vegetation)
        assert result.sum() == pytest.approx(1.0) and np.all(result >= 0.0)
        shares.append(result[0])
    # Before the fix mountain_sun + mountain_shade == mountain, so the shadowing was lost.
    assert shares[0] > shares[1] > shares[2]


def test_zip_pairs_solar_angles_instead_of_expanding_scenarios():
    config = {"sza": "20;30;40", "vza": "0;30;60", "vaa": "0;0;0",
              "geometry_pair_mode": "zip", "geometry_mode": "0"}
    assert len(run._scenario_configs(config)) == 1
    geometries = run._geometry(config, run.DEFAULT_INPUT)
    assert [(g[0], g[2]) for g in geometries] == [(0.0, 20.0), (30.0, 30.0), (60.0, 40.0)]
    cross = dict(config, geometry_pair_mode="cross")
    assert len(run._scenario_configs(cross)) == 3


def test_stomatal_humidity_follows_hourly_vapour_pressure():
    assert _ambient_vapour_pressure({"ea": 15.0}, {"vapor_pressure": 22.5}) == 22.5
    assert _ambient_vapour_pressure({"ea": 15.0}, {"vapor_pressure": float("nan")}) == 15.0


def test_urban_temperature_method_is_honoured():
    settings = {"temperature_method": "direct", "urban_temperature_method": "bogus"}
    with pytest.raises(ValueError, match="temperature_method"):
        eb.rebalance_urban(settings, {"pressure": 1013.0, "vapor_pressure": 10.0,
                                      "air_temperature_c": 20.0}, {})


def test_robust_denominator_includes_leaf_emission_of_both_faces():
    rho, ra, rs, p, t = 1.2, 50.0, 200.0, 1013.0, 300.0
    one = eb._jacobian(rho, ra, rs, p, t, 0.98, 1.0, 0.0)
    six = eb._jacobian(rho, ra, rs, p, t, 0.98, 6.0, 0.0)
    assert six - one == pytest.approx(5.0 * 4.0 * 0.98 * eb.SIGMASB * t ** 3)


def _coefficients(row_count=3, kca=0.0):
    # wavelength, n, Kab, Kca, Ks, Kw, Kdm, ..., Kant (index 11)
    return [[400.0 + i, 1.45, 0.0, kca, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0.0] for i in range(row_count)]


def test_non_absorbing_leaf_does_not_divide_by_zero():
    leaf = {"Cab": 0.0, "Cw": 0.0, "Cdm": 0.0, "Cs": 0.0, "N": 1.5}
    reflectance, transmittance = fluspect(_coefficients(), leaf)
    assert all(math.isfinite(r) and math.isfinite(t) for r, t in zip(reflectance, transmittance))
    assert reflectance[0] + transmittance[0] == pytest.approx(1.0, abs=1e-6)


def test_carotenoids_are_optional_and_absorb():
    base = {"Cab": 0.0, "Cw": 0.0, "Cdm": 0.001, "Cs": 0.0, "N": 1.5}
    coefficients = _coefficients(kca=0.05)
    r0, t0 = fluspect(coefficients, base)
    r1, t1 = fluspect(coefficients, dict(base, Cca=10.0))
    assert r1[0] + t1[0] < r0[0] + t0[0]
    assert fluspect(coefficients, [0.0, 0.0, 0.001, 0.0, 1.5]) == (r0, t0)


def test_comma_separated_value_warns(tmp_path):
    path = tmp_path / "input.csv"
    path.write_text("section,parameter,value,description\ngeometry,vza,0,30,60\n", encoding="utf-8")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        config = read_input(path)
    assert config["vza"] == "0"
    assert any("split by ','" in str(item.message) for item in caught)
