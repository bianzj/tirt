"""Regression: a global row direction must be rotated into each terrain facet's frame."""
import numpy as np

from base.terrain_rt import compound_terrain_proportions


def _settings(row_azimuth=0.0):
    return dict(
        surface_model="terrain", terrain_type=2,
        terrain_height="20", terrain_radius="10", terrain_density="0.001",
        vegetation_model="row", row_azimuth=row_azimuth, lai=1.5, hspot=0.15,
        row_width=0.8, row_blank=0.8, row_height=1.5,
    )


def test_row_azimuth_is_rotated_per_facet():
    seen = []

    def proportions(sza, saa, vza, vaa, sunlit, row_azimuth=None):
        seen.append(row_azimuth)
        return np.array([0.25, 0.25, 0.25, 0.25])

    compound_terrain_proportions(_settings(0.0), 40.0, 180.0, 60.0, 90.0, proportions)
    values = [value for value in seen if value is not None]
    assert values, "compound terrain did not pass a local row azimuth"
    assert len({round(float(value), 6) for value in values}) > 1, "row azimuth was the same on every facet"
    assert all(0.0 <= float(value) < 360.0 for value in values)


def test_legacy_callback_without_keyword_still_works():
    calls = []

    def proportions(sza, saa, vza, vaa, sunlit):
        calls.append(1)
        return np.array([0.25, 0.25, 0.25, 0.25])

    result = compound_terrain_proportions(_settings(45.0), 40.0, 180.0, 60.0, 90.0, proportions)
    assert calls
    assert np.isclose(result.sum(), 1.0)
