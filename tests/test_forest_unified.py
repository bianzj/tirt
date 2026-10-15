"""Shared geometry/closure with distinct optical and thermal sources."""
import json
import unittest
from unittest.mock import patch
import numpy as np
import run
from base import forest_exchange as common
from base import thermal_forest_exchange as thermal
from base import optical_terrain_reciprocal as optical
from base.radiation import directional_reflectance


SETTINGS=dict(surface_model='terrain',vegetation_model='crown',surface_type=1,terrain_type=2,
    lai=1.5,stand_density=.02,height=2.,width=2.,hspot=.0375,crown_bottom=3.,
    terrain_height=100.,terrain_radius=100.,terrain_density=1.35e-5,
    optical_crown_geometry='clumping_dual',forest_transport='reciprocal_exchange')


class UnifiedForestTests(unittest.TestCase):
    def test_thermal_and_optical_import_identical_geometry_and_solver(self):
        self.assertIs(thermal.exchange_operator,common.exchange_operator)
        self.assertIs(thermal.higher_orders,common.higher_orders)
        self.assertIs(optical.exchange_operator,common.exchange_operator)
        self.assertIs(optical.higher_orders,common.higher_orders)
        self.assertIs(thermal.outer_width,optical.outer_width)

    def test_common_switch_selects_both_entry_points(self):
        with patch('base.optical_terrain_reciprocal.reflectance',return_value=.123) as function:
            value=directional_reflectance(dict(SETTINGS,optical_terrain_method='finite_crown_ports'),
                30.,0.,30.,90.,.45,.42,.3)
            self.assertEqual(value,.123);function.assert_called_once()
        config=dict(lai=1.5,hspot=.0375,crown_density=.02,crown_rad_a=2.,crown_rad_b=2.,
            crown_height=3.,terrain_height=100.,terrain_radius=100.,terrain_density=1.35e-5,
            thermal_crown_geometry='clumping_dual',forest_transport='reciprocal_exchange',
            soil_sunlit_k=318.15,soil_shaded_k=303.15,leaf_sunlit_k=306.15,leaf_shaded_k=303.15)
        spectrum=dict(soil_e=.95,leaf_e=.985,terrain_e=.95)
        common_result=run._surface_result(config,'terrain','crown',10.5,(30.,90.,30.,0.),spectrum)
        config.pop('forest_transport');config['thermal_forest_transport']='reciprocal_exchange'
        explicit_result=run._surface_result(config,'terrain','crown',10.5,(30.,90.,30.,0.),spectrum)
        np.testing.assert_array_equal(common_result,explicit_result)

    def test_optical_reflection_and_transmission_stay_distinct(self):
        reflected=optical.reflectance(SETTINGS,30.,0.,30.,0.,.6,0.,.3)
        transmitted=optical.reflectance(SETTINGS,30.,0.,30.,0.,0.,.6,.3)
        self.assertGreater(abs(reflected-transmitted),1e-4)
        direct=optical.reflectance(SETTINGS,30.,0.,30.,90.,.45,.42,.3)
        changed=optical.reflectance(dict(SETTINGS,soil_emissivity=.1,leaf_emissivity=.1),
            30.,0.,30.,90.,.45,.42,.3)
        self.assertEqual(direct,changed)

    def test_discrete_state_flux_balance_and_material_extremes(self):
        for rho,tau,soil in ((0.,0.,0.),(.06,.03,.18),(.45,.42,.3),(.5,.5,1.)):
            for f in (0.,.4,1.):
                value,budget=optical.reflectance(SETTINGS,30.,0.,50.,90.,rho,tau,soil,f,return_budget=True)
                self.assertGreaterEqual(value,0.)
                self.assertLess(abs(budget['residual']),1e-10)
                if rho+tau==1 and soil==1:
                    self.assertAlmostEqual(budget['absorbed'],0.,places=11)
                if rho+tau+soil==0:
                    self.assertEqual(value,0.)
        value,budget=optical.reflectance(SETTINGS,30.,0.,30.,90.,.45,.42,.3,
            sunlit=False,return_budget=True)
        self.assertEqual(value,0.);self.assertEqual(budget['incident'],0.)

    def test_rotation_diffuse_linearity_and_no_crown_cubature(self):
        with patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('volume integration')):
            direct=optical.reflectance(SETTINGS,30.,0.,45.,90.,.45,.42,.3)
            rotated=optical.reflectance(SETTINGS,30.,63.,45.,153.,.45,.42,.3)
            sky=optical.reflectance(SETTINGS,30.,0.,45.,90.,.45,.42,.3,1.)
            mixed=optical.reflectance(SETTINGS,30.,0.,45.,90.,.45,.42,.3,.4)
        self.assertAlmostEqual(direct,rotated,places=12)
        self.assertAlmostEqual(mixed,.6*direct+.4*sky,places=12)

    def test_material_and_geometry_validation(self):
        with self.assertRaises(ValueError):
            optical.reflectance(SETTINGS,30.,0.,30.,90.,.8,.4,.3)
        with self.assertRaises(ValueError):
            common.outer_width(2.,2.,-1.)
        with self.assertRaises(ValueError):
            common.higher_orders(np.ones((4,4)),np.ones(4),np.ones(4))


if __name__=='__main__':
    unittest.main()
