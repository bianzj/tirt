"""Degenerate limits of the shortwave closure: vanishing canopy and single-slope illumination."""
import math
from pathlib import Path

import pytest

import run
from base import canopy_closure as cc
from base.inputio import read_input
from base.solar import slope_direction

ROOT = Path(__file__).resolve().parents[1]
BAND = dict(leaf_reflectance=.45, leaf_transmittance=.40, soil_reflectance=.25,
            roof_reflectance=.20, wall_reflectance=.25, street_reflectance=.18)
HOT, FORWARD = (40., 130., 40., 130.), (40., 310., 40., 130.)
CASES = (("hom", cc.HOMOGENEOUS), ("row", (.5, .5, 1.)), ("crown", (.02, 3., 3.)))


def _converges(values, limit):
    errors = [abs(v - limit) for v in values]
    assert all(b < a for a, b in zip(errors, errors[1:])), errors
    assert errors[-1] < 1.e-7, errors


@pytest.mark.parametrize("kind,shape", CASES)
@pytest.mark.parametrize("angles", (HOT, FORWARD))
def test_vanishing_lai_converges_to_bare_soil(kind, shape, angles):
    limit = cc.reflectance(kind, 0., .15, shape, .45, .40, .25, *angles)
    assert limit == .25
    _converges([cc.reflectance(kind, lai, .15, shape, .45, .40, .25, *angles)
                for lai in (1.e-4, 1.e-6, 1.e-8)], limit)


@pytest.mark.parametrize("angles", (HOT, FORWARD))
def test_vanishing_crown_density_converges_to_bare_soil(angles):
    limit = cc.reflectance("crown", 1.5, .15, (0., 3., 3.), .45, .40, .25, *angles)
    _converges([cc.reflectance("crown", 1.5, .15, (d, 3., 3.), .45, .40, .25, *angles)
                for d in (1.e-3, 1.e-6, 1.e-9)], limit)


def _brf(surface, veg, vza, vaa, sza, saa, **over):
    config = dict(read_input(ROOT / "input.csv"))
    config.update({k: str(v) for k, v in over.items()})
    config.update(surface_model=surface, vegetation_model=veg, lai="2", hspot="0.1")
    return run._optical_band(config, surface, veg, .86, (vza, vaa, sza, saa), BAND)


@pytest.mark.parametrize("veg", ("bare", "hom", "row", "crown"))
def test_zero_slope_is_the_plane(veg):
    # A nonzero aspect rotates the local azimuth frame; the plane hotspot moments use a
    # 96-azimuth grid, so rotation invariance holds to quadrature (~1e-10), not roundoff.
    for aspect, tolerance in ((0., 1.e-12), (73., 1.e-9)):
        for fraction in (0., .2, 1.):
            slope = _brf("slope", veg, 50., 20., 35., 200., slope_angle=0., slope_aspect=aspect,
                         row_azimuth=17., shortwave_diffuse_fraction=fraction)
            plane = _brf("plane", veg, 50., 20., 35., 200., row_azimuth=17.,
                         shortwave_diffuse_fraction=fraction)
            assert slope == pytest.approx(plane, abs=tolerance)


@pytest.mark.parametrize("veg", ("bare",
                                 *(pytest.param(v, marks=pytest.mark.xfail(strict=True, reason=(
                                     "known limitation, TiRT theory doc 13: the first-order term does not "
                                     "vanish with the grazing incident cosine, so BRF is not proportional to "
                                     "cos(i) near the slope horizon")))
                                   for v in ("hom", "row", "crown"))))
def test_sun_crossing_the_slope_horizon_is_continuous(veg):
    # Slope 40 facing north, sun due south: the local incidence vanishes at sza = 50.
    def value(sza, fraction):
        return _brf("slope", veg, 30., 0., sza, 180., slope_angle=40., slope_aspect=0.,
                    shortwave_diffuse_fraction=fraction)
    for fraction in (0., .2):
        values = [value(s, fraction) for s in (49., 49.9, 49.99, 50.01, 50.1)]
        assert abs(values[2] - values[3]) < 2.e-4
        assert values[3] == values[4]
        assert values[0] > values[1] > values[2] > values[3]
    assert value(50.01, 0.) == 0.


@pytest.mark.parametrize("veg", ("hom", "row"))
def test_horizontal_irradiance_normalization(veg):
    slope, aspect, sza, saa = 25., 160., 35., 200.
    local_view, local_vaa, _ = slope_direction(50., 20., slope, aspect)
    local_sun, local_saa, _ = slope_direction(sza, saa, slope, aspect)
    row = slope_direction(90., 17., slope, aspect)[1]
    incidence = math.cos(math.radians(local_sun)) / math.cos(math.radians(sza))
    for fraction, weight in ((0., incidence), (1., (1. + math.cos(math.radians(slope))) / 2.)):
        tilted = _brf("slope", veg, 50., 20., sza, saa, slope_angle=slope, slope_aspect=aspect,
                      row_azimuth=17., shortwave_diffuse_fraction=fraction)
        local = _brf("plane", veg, local_view, local_vaa, local_sun, local_saa, row_azimuth=row,
                     shortwave_diffuse_fraction=fraction)
        assert tilted == pytest.approx(weight * local, abs=1.e-12)


@pytest.mark.parametrize("veg", ("hom", "row", "crown"))
def test_slope_rotation_equivalence(veg):
    def value(turn):
        return _brf("slope", veg, 50., 20. + turn, 35., 200. + turn, slope_angle=20., slope_aspect=10. + turn,
                    row_azimuth=10. + turn, shortwave_diffuse_fraction=.2)
    base = value(0.)
    for turn in (37., 360.):
        assert value(turn) == pytest.approx(base, abs=1.e-10)


def test_shaded_bare_soil_is_not_counted_as_sunlit():
    from base.radiation import _vegetation_proportions
    settings = dict(lai=0., surface_type=0, distribution=0, vegetation_model="bare")
    assert list(_vegetation_proportions(settings, 30., 0., 0., 0., True)) == [1., 0., 0., 0.]
    assert list(_vegetation_proportions(settings, 30., 0., 0., 0., False)) == [0., 1., 0., 0.]


def test_bare_compound_terrain_shortwave_sees_the_terrain():
    def value(vza, vaa, fraction):
        return _brf("terrain", "bare", vza, vaa, 30., 0.,
                    terrain_height=20., terrain_radius=10., terrain_density=0.001, terrain_type=2,
                    lai=0., shortwave_diffuse_fraction=fraction)
    # The terrain fraction of the view changes with direction; a flat soil would not depend on it.
    values = {value(vza, vaa, 0.) for vza, vaa in ((0., 0.), (30., 180.), (60., 180.))}
    assert len(values) == 3
    assert max(values) < 0.25 + 1.e-9
    # A fully diffuse sky sees only the Lambertian soil albedo.
    assert value(30., 180., 1.) == pytest.approx(0.25, abs=1.e-12)


def test_shortwave_terrain_reads_the_same_cone_classes_as_the_thermal_band():
    from base.terrain_rt import _shape_values
    assert _shape_values({}) == ((20., 10., 0.001),)
    assert _shape_values(dict(terrain_height="10;20", terrain_radius="5;10", terrain_density="0.002;0.001")) == (
        (10., 5., 0.002), (20., 10., 0.001))
    assert _shape_values(dict(terrain_height="10;20", terrain_radius=10., terrain_density=0.001)) == (
        (10., 10., 0.001), (20., 10., 0.001))
