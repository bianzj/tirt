"""Physical invariants of the opt-in four-state thermal forest closure."""
import unittest
from unittest.mock import patch
import numpy as np
import run
from base.terrain_forest import TerrainForest
from base.thermal_forest_exchange import exchange_operator,canopy_block
from base.physicsF import planck


def model(emissivity=(.95,.985),bottom=3.,outer=None,lai=1.5,density=.02):
    m=TerrainForest()
    m.set_hotspot_geometry('clumping_dual',outer=outer)
    m.set_radiative_transport('reciprocal_exchange')
    m.forestshape=np.array([lai,density,2.,2.,.0375,bottom])
    m.set_structural_input([[100.,100.,1.35e-5]])
    m.set_spectral_input(*emissivity)
    m.set_angular_input(np.array([30.]),np.array([90.]),30.,0.)
    return m


class ThermalForestExchangeTests(unittest.TestCase):
    def test_probability_budget_and_area_reciprocity(self):
        for h,r,lai,density in ((100.,100.,1.5,.02),(20.,30.,1.5,.02),(30.,80.,.1,.001)):
            f,a=exchange_operator(((h,r,.4/(np.pi*r*r)),),lai,density,2.,2.)
            self.assertGreaterEqual(f.min(),0.)
            self.assertLessEqual(f.sum(axis=1).max(),1+1e-12)
            np.testing.assert_allclose(a[:,None]*f,(a[:,None]*f).T,atol=1e-12)

    def test_blackbody_and_no_emission_limits(self):
        black=model((1.,1.))
        parts=black.calculate_effective_component_emissivity(2)[0]
        self.assertAlmostEqual(parts.sum(),1.,places=12)
        self.assertAlmostEqual(float(parts@planck(10.5,np.full(4,300.))),float(planck(10.5,300.)),places=12)
        mirror=model((0.,0.))
        np.testing.assert_allclose(mirror.calculate_effective_component_emissivity(2),0.,atol=1e-12)

    def test_joint_closure_matches_scattering_orders(self):
        m=model()
        effective=m.calculate_effective_component_emissivity(2)[0]
        b=m.thermal_exchange_budget
        current=b['transport']@b['source']
        series=current.copy()
        for _ in range(30):
            current=b['transport']@current
            series+=current
        np.testing.assert_allclose(b['first_hit']@series,b['scattered'],atol=1e-13)
        direct=np.asarray(m.calculate_component_direct_emissivity(0,2))
        np.testing.assert_allclose(effective,direct+b['scattered'],atol=1e-13)
        self.assertLessEqual(effective.sum(),1.)

    def test_isothermal_directional_kirchhoff_identity(self):
        m=model();m.calculate_effective_component_emissivity(2)
        b=m.thermal_exchange_budget
        eps=np.array([m.Es,m.Em,m.Es,m.Em])
        expected=b['first_hit']@np.linalg.solve(np.eye(4)-b['transport'],eps)
        self.assertAlmostEqual(float(m.canopy_effective_emissivity.sum()),float(expected),places=12)

    def test_height_scale_and_explicit_override(self):
        args=(1.5,.02,.0375,2.,2.,np.array([30.]),np.array([30.]),np.array([90.]))
        m=model()
        auto=m._forest_proportions(*args,hcb=3.)
        explicit=model(outer=.4)._forest_proportions(*args,hcb=3.)
        np.testing.assert_array_equal(auto,explicit)
        changed=m._forest_proportions(*args,hcb=0.)
        self.assertGreater(np.max(abs(auto-changed)),1e-4)
        fixed=model(outer=.7)
        np.testing.assert_array_equal(fixed._forest_proportions(*args,hcb=0.),fixed._forest_proportions(*args,hcb=3.))

    def test_no_old_canopy_scatter_or_crown_volume_integration(self):
        m=model()
        with patch('base.terrain_forest.emissivity_scatter_analytical',side_effect=AssertionError('double local closure')),\
             patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('volume integration')):
            self.assertTrue(np.isfinite(m.calculate_effective_component_emissivity(2)).all())

    def test_absent_foliage_keeps_macro_ground_reflection(self):
        for options in (dict(lai=0.),dict(density=0.)):
            m=model(**options)
            parts=m.calculate_effective_component_emissivity(2)[0]
            self.assertGreater(parts[:2].sum(),.95)
            np.testing.assert_allclose(parts[2:],0.,atol=1e-12)
        local,escape,entry=canopy_block(0.,1.)
        np.testing.assert_array_equal(local,np.zeros((2,2)))

    def test_invalid_geometry_mode_and_multiple_cones_rejected(self):
        m=TerrainForest()
        with self.assertRaises(ValueError):
            m.set_radiative_transport('reciprocal_exchange')
        m=model()
        with self.assertRaises(ValueError):
            m.set_hotspot_geometry('volume')
        with self.assertRaises(ValueError):
            exchange_operator(((20.,30.,.0001),(10.,20.,.0001)),1.,.02,2.,2.)


if __name__=='__main__':
    unittest.main()
