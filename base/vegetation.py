"""Shared canopy geometry for the coupled surface models.

The analytical four-component kernels are the existing TiRT hom/row/crown
models. Geometry is independent of material emissivity and temperature.
"""

import numpy as np

from .physicsF import slope1
from .proportion import (
    proportion_bidirectional_crown_one,
    proportion_bidirectional_hom_one,
    proportion_bidirectional_row_one,
)


ROW_GEOMETRIES = ("legacy", "conditioned_corridor", "parametric")


def vegetation_proportions(kind, lai, hspot, vza, vaa, sza, saa,
                          density=0.02, crown_height=3.0, crown_radius=3.0,
                          row_width=0.5, row_blank=0.5, row_height=1.0,
                          row_azimuth=0.0, row_geometry="legacy"):
    """Return soil sun/shade and leaf sun/shade probabilities, shape (4, N)."""
    if kind not in {"bare", "hom", "row", "crown"}:
        raise ValueError(f"Unsupported vegetation model: {kind}")
    if row_geometry not in ROW_GEOMETRIES:
        raise ValueError(f"Unknown thermal row geometry: {row_geometry}")
    vza, vaa, sza, saa, row_azimuth = [a.reshape(-1) for a in np.broadcast_arrays(
        np.asarray(vza, float), np.asarray(vaa, float), np.asarray(sza, float),
        np.asarray(saa, float), np.asarray(row_azimuth, float))]
    if lai < 0.0 or hspot < 0.0:
        raise ValueError("LAI and hotspot must be nonnegative")
    if (kind == "bare" or lai == 0.0 or
            (kind == "crown" and np.all(np.asarray(density) == 0.0))):
        result = np.zeros((4, vza.size))
        result[0] = 1.0
        return result
    # The horizon and zero-hotspot values are evaluated as one-sided limits.
    vza = np.clip(vza, 0.0, 89.999)
    sza = np.clip(sza, 0.0, 89.999)
    relative = np.abs((vaa - saa + 180.0) % 360.0 - 180.0)
    hotspot = max(hspot, 1.0e-12)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if kind == "hom" or (kind == "row" and row_blank == 0.0):
            components = proportion_bidirectional_hom_one(lai, hotspot, vza, sza, relative)
        elif kind == "row" and row_geometry == "parametric":
            if row_width <= 0.0 or row_blank < 0.0 or row_height <= 0.0:
                raise ValueError("Row width/height must be positive; row gap must be nonnegative")
            from .row_parametric import row_fractions
            components = np.asarray([row_fractions(lai, hotspot, row_width, row_blank, row_height, v, a, s, sa, az)
                                     for v, a, s, sa, az in zip(vza, vaa, sza, saa, row_azimuth)]).T
        elif kind == "row":
            if row_width <= 0.0 or row_blank < 0.0 or row_height <= 0.0:
                raise ValueError("Row width/height must be positive; row gap must be nonnegative")
            components = proportion_bidirectional_row_one(
                lai, hotspot, row_width, row_blank, row_height,
                vza, vaa % 360.0, sza, saa % 360.0, row_azimuth % 360.0)
        else:
            components = proportion_bidirectional_crown_one(
                lai, density, hotspot, crown_height, crown_radius, vza, sza, relative)
    result = np.asarray(components, float).reshape(4, -1)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(f"Non-finite {kind} canopy proportions")
    # The empirical hotspot approximation must respect its directional
    # marginal probabilities; enforce probability bounds, not a scene fallback.
    soil = np.clip(result[0] + result[1], 0.0, 1.0)
    leaf = 1.0 - soil
    soil_sun = np.clip(result[0], 0.0, soil)
    leaf_sun = np.clip(result[2], 0.0, leaf)
    result = np.asarray([soil_sun, soil - soil_sun, leaf_sun, leaf - leaf_sun])
    if kind == "row" and row_geometry == "conditioned_corridor" and row_blank > 0.0:
        from .thermal_row import conditioned_proportions
        result = conditioned_proportions(lai, hspot, row_width, row_blank, row_height,
                                         vza, vaa, sza, saa, row_azimuth, result)
    return result


class CanopyGeometry:
    vegetation_model = "crown"
    row_shape = (0.5, 0.5, 1.0, 0.0)
    row_geometry = "parametric"

    def set_vegetation_input(self, kind, row_width=0.5, row_blank=0.5,
                             row_height=1.0, row_azimuth=0.0, row_geometry="parametric"):
        if kind not in {"bare", "hom", "row", "crown"}:
            raise ValueError(f"Unsupported vegetation model: {kind}")
        if row_geometry not in ROW_GEOMETRIES:
            raise ValueError(f"Unknown thermal row geometry: {row_geometry}")
        self.vegetation_model = kind
        self.row_shape = (row_width, row_blank, row_height, row_azimuth)
        self.row_geometry = row_geometry

    def canopy_proportions(self, vza, vaa, sza, saa, density=None,
                           crown_height=None, row_azimuth=None, row_geometry=None):
        lai, std, rcr, hcr, hspot = self.forestshape[:5]
        width, blank, height, azimuth = self.row_shape
        return vegetation_proportions(
            self.vegetation_model, lai, hspot, vza, vaa, sza, saa,
            std if density is None else density,
            hcr if crown_height is None else crown_height, rcr,
            width, blank, height, azimuth if row_azimuth is None else row_azimuth,
            self.row_geometry if row_geometry is None else row_geometry)

    def local_canopy_proportions(self, vza, vaa, sza, saa, slope, aspect,
                                 density=None, crown_height=None):
        """Project global directions and row bearing into each slope plane."""
        vza, vaa, sza, saa, slope, aspect = np.broadcast_arrays(
            np.atleast_1d(vza), np.atleast_1d(vaa), np.atleast_1d(sza),
            np.atleast_1d(saa), np.atleast_1d(slope), np.atleast_1d(aspect))
        with np.errstate(divide="ignore", invalid="ignore"):
            local_vza, local_vaa = slope1(vza, vaa, slope, aspect)
            local_sza, local_saa = slope1(sza, saa, slope, aspect)
            row = np.full(vza.shape, self.row_shape[3])
            _, local_row = slope1(np.full(vza.shape, 90.0), row, slope, aspect)
        # Avoid the legacy slope1's artificial 0.81-degree nadir offset.
        cosine_view = (np.cos(np.deg2rad(slope)) * np.cos(np.deg2rad(vza)) +
                       np.sin(np.deg2rad(slope)) * np.sin(np.deg2rad(vza)) *
                       np.cos(np.deg2rad(vaa - aspect)))
        cosine_sun = (np.cos(np.deg2rad(slope)) * np.cos(np.deg2rad(sza)) +
                      np.sin(np.deg2rad(slope)) * np.sin(np.deg2rad(sza)) *
                      np.cos(np.deg2rad(saa - aspect)))
        local_vza = np.rad2deg(np.arccos(np.clip(cosine_view, 0.0, 1.0)))
        local_sza = np.rad2deg(np.arccos(np.clip(cosine_sun, 0.0, 1.0)))
        proportions = self.canopy_proportions(
            local_vza, local_vaa, local_sza, local_saa, density,
            crown_height, local_row)
        shadow = cosine_sun.reshape(-1) <= 0.0
        proportions[1, shadow] += proportions[0, shadow]
        proportions[3, shadow] += proportions[2, shadow]
        proportions[0, shadow] = 0.0
        proportions[2, shadow] = 0.0
        return proportions
