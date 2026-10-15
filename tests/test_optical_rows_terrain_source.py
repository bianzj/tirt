"""Corridor limits, source-depth transport and optical/thermal isolation."""
import math
import unittest
from unittest.mock import patch
import numpy as np

from base.optical_row import clear_corridor, conditioned_first
from base.optical_crown import transport, geometry, crown_reflectance
from base.optical_crown_profile import beam_down_escape
from base.optical_terrain import exchange_geometry, terrain_crown_reflectance
from base.radiation import directional_reflectance
from base.vnir import joint_transmission, spherical_coefficients

ROW = dict(surface_model='plane',vegetation_model='row',distribution=2,surface_type=0,
           lai=1.5,hspot=.15,row_width=.5,row_blank=.5,row_height=1.,row_azimuth=0.)
FOREST = dict(surface_model='terrain',vegetation_model='crown',surface_type=1,
              lai=1.5,hspot=.025,height=3.,width=3.,stand_density=.02,
              terrain_type=2,terrain_height=20.,terrain_radius=30.,terrain_density=.00015,
              optical_crown_geometry='clumping_dual',optical_crown_source='diffusion_profile',
              optical_crown_quadrature=4)


class OpticalRowsTerrainSourceTests(unittest.TestCase):
    def test_corridor_tracks_signed_shadows_and_rotation(self):
        self.assertAlmostEqual(clear_corridor(.5,.5,1.,30.,0.,75.,180.,0.), .5)
        self.assertEqual(clear_corridor(.5,.5,1.,30.,0.,75.,90.,0.), 0.)
        a = clear_corridor(.5,.5,1.,15.,12.,20.,25.,0.)
        b = clear_corridor(.5,.5,1.,15.,85.,20.,98.,73.)
        self.assertAlmostEqual(a,b,places=14)
        self.assertLess(clear_corridor(.5,.5,1.,15.,12.,20.,205.,0.),a)

    def test_along_row_first_scattering_is_conditional_slab(self):
        for view, azimuth in ((0.,0.),(30.,0.),(75.,0.),(75.,180.)):
            ks,ko,sob,sof = spherical_coefficients(30.,view,azimuth)
            gs,gv = .5+.5*math.exp(-ks*3),.5+.5*math.exp(-ko*3)
            joint,path = joint_transmission(3.,ks,ko,.15,30.,view,azimuth)
            actual = conditioned_first(1.5,.15,.45,.42,.3,30.,view,azimuth,gs,gv,.5)
            expected = .15+.5*((sob*.45+sof*.42)*path+.3*joint)
            self.assertAlmostEqual(actual,expected,places=12)
            # Jointly sunlit empty soil survives even at grazing angles.
            soil_only = conditioned_first(1.5,.15,0.,0.,1.,30.,view,azimuth,gs,gv,.5)
            self.assertGreaterEqual(soil_only,.5)
            self.assertLessEqual(soil_only,min(gs,gv)+1e-12)

    def test_closed_corridor_and_row90_retain_diffuse_closure(self):
        for settings,view,azimuth in ((ROW,75.,90.),(dict(ROW,row_azimuth=90.),75.,0.)):
            actual = directional_reflectance(settings,30.,0.,view,azimuth,.45,.42,.3)
            old = directional_reflectance(dict(settings,optical_row_method='effective_slab'),
                                          30.,0.,view,azimuth,.45,.42,.3)
            self.assertEqual(actual,old)
        # The first correction is linear in materials; subtracting its
        # low-albedo limit leaves exactly the same higher-order term.
        values = []
        for method in ('conditioned_corridor','effective_slab'):
            settings = dict(ROW,optical_row_method=method)
            def r(a):
                return directional_reflectance(settings,30.,0.,75.,0.,.45*a,.42*a,.3*a)
            values.append(r(1.)-(2*r(.02)/.02-r(.04)/.04))
        self.assertAlmostEqual(*values,places=12)

    def test_row_rotation_continuity_and_illumination_linearity(self):
        beam = directional_reflectance(ROW,30.,0.,75.,0.,.45,.42,.3)
        rotated = directional_reflectance(dict(ROW,row_azimuth=73.),30.,73.,75.,73.,.45,.42,.3)
        self.assertAlmostEqual(beam,rotated,places=12)
        sky = directional_reflectance(ROW,30.,0.,75.,0.,.45,.42,.3,diffuse_fraction=1.)
        mix = directional_reflectance(ROW,30.,0.,75.,0.,.45,.42,.3,diffuse_fraction=.4)
        self.assertAlmostEqual(mix,.6*beam+.4*sky,places=12)
        edge = math.degrees(math.atan(.5))
        left = directional_reflectance(ROW,30.,0.,edge-1e-6,90.,.45,.42,.3)
        right = directional_reflectance(ROW,30.,0.,edge+1e-6,90.,.45,.42,.3)
        self.assertLess(abs(left-right),1e-6)
        self.assertEqual(directional_reflectance(ROW,30.,0.,75.,0.,0.,0.,0.),0.)
        self.assertEqual(directional_reflectance(ROW,30.,0.,75.,0.,.45,.42,.3,sunlit=False),0.)

    def test_beam_primary_escape_and_energy_without_spatial_grid(self):
        shape = (1.5,.02,3.,3.,.025)
        self.assertLess(beam_down_escape(shape,30.,12),geometry(shape)[0])
        self.assertLess(beam_down_escape(shape,75.,12),beam_down_escape(shape,30.,12))
        with patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('spatial grid')):
            for lai in (1e-7,1.5,6.):
                for rho,tau,soil in ((.06,.03,.18),(.45,.42,.3),(.5,.5,1.)):
                    _,budget = transport((lai,.02,3.,3.,.025),30.,rho,tau,soil,4,
                                         'clumping_dual','beam_diffusion_profile')
                    self.assertLess(abs(budget['balance_residual']),1e-10)
                    if rho+tau == 1:
                        self.assertAlmostEqual(budget['reflected'],1.,places=10)

    def test_terrain_exchange_area_reciprocity(self):
        for height,radius,density in ((20.,30.,.00015),(5.,30.,.00015),(40.,20.,.0001)):
            matrix = exchange_geometry(((height,radius,density),))
            areas = np.array([density*math.pi*radius*math.hypot(height,radius),
                              1-density*math.pi*radius**2])
            exchange = areas[:,None]*matrix
            np.testing.assert_allclose(exchange,exchange.T,atol=1e-12)
            self.assertTrue((matrix >= 0).all())
            self.assertTrue((matrix.sum(axis=1) <= 1).all())

    def test_new_terrain_source_reduces_to_plane_when_no_cones(self):
        settings = dict(FOREST,terrain_density=0.)
        actual = terrain_crown_reflectance(settings,30.,0.,25.,0.,.45,.42,.3)
        expected = crown_reflectance(dict(settings,surface_model='plane'),30.,25.,0.,.45,.42,.3)
        self.assertEqual(actual,expected)


if __name__ == '__main__':
    unittest.main()
