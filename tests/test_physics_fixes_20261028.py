"""Regression tests for the second round of 2026-10-28 physics fixes."""
import math
import warnings

import numpy as np
import pytest

from base.physicsF import inv_planck, planck
from base.soilopt import _soilwat
from base.urban import validate_urban_shapes
from base.vegetation import vegetation_proportions


def test_negative_view_zenith_mirrors_instead_of_nadir():
    kwargs = dict(kind="hom", lai=2.0, hspot=0.1)
    mirrored = vegetation_proportions(vza=-40.0, vaa=0.0, sza=30.0, saa=0.0, **kwargs)
    expected = vegetation_proportions(vza=40.0, vaa=180.0, sza=30.0, saa=0.0, **kwargs)
    nadir = vegetation_proportions(vza=0.0, vaa=0.0, sza=30.0, saa=0.0, **kwargs)
    np.testing.assert_allclose(mirrored, expected)
    assert not np.allclose(mirrored, nadir)


def test_urban_shapes_reject_zero_height_and_full_cover():
    validate_urban_shapes([[10.0, 10.0, 20.0, 0.001, 0.0, 90.0]])
    with pytest.raises(ValueError, match="positive"):
        validate_urban_shapes([[10.0, 10.0, 0.0, 0.001, 0.0, 90.0]])
    with pytest.raises(ValueError, match="cover"):
        validate_urban_shapes([[10.0, 10.0, 20.0, 0.01, 0.0, 90.0]])


def test_inv_planck_of_nonpositive_radiance_is_nan():
    assert math.isnan(float(inv_planck(10.5, 0.0)))
    assert float(inv_planck(10.5, planck(10.5, 300.0))) == pytest.approx(300.0)


def test_planck_warns_on_celsius_guess():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert float(planck(10.5, 26.85)) == pytest.approx(float(planck(10.5, 300.0)))
    assert any("kelvin" in str(item.message) for item in caught)


def test_soilwat_dry_limit_and_darkening():
    args = dict(rw=0.02, p=0.5, rbac=0.3, kw=50.0, film=0.015)
    assert _soilwat(0.3, mu=0.0, **args) == 0.3
    # Stronger water absorption (two-way film path) darkens the wet soil.
    wet = _soilwat(0.3, mu=1.0, **args)
    assert wet < 0.3
    assert _soilwat(0.3, mu=1.0, **dict(args, kw=100.0)) < wet


def test_night_slope_is_not_sunlit():
    from base.solar import slope_direction
    _, _, facing = slope_direction(95.0, 180.0, 30.0, 180.0)
    assert facing  # incidence > 0 even though the sun is below the horizon
    # pipeline combines it with the global elevation
    assert not (facing and 95.0 < 90.0)
