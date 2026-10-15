"""Physical limits and explicit finite-terrain source-depth isolation."""
import math
import unittest
import numpy as np
from base import optical_crown as crown,optical_crown_intrinsic as intrinsic
from base.optical_terrain_ports import reflectance


class IntrinsicCrownTests(unittest.TestCase):
    shape=(1.5,.02,3.,3.,.025)
    scales=(.025,1.)
    settings=dict(surface_model='terrain',vegetation_model='crown',lai=1.5,
        stand_density=.02,height=3.,width=3.,hspot=.025,optical_crown_geometry='clumping_dual',
        optical_crown_source='diffusion_profile',terrain_radius=30.,terrain_height=20.,
        terrain_density=.00015,terrain_type=2,crown_bottom=3.,optical_crown_quadrature=12)

    def test_intrinsic_depth_preserves_lad_under_facet_area_conversion(self):
        cs=math.cos(math.atan2(20.,30.))
        facet=(self.shape[0]*cs,self.shape[1]*cs,*self.shape[2:])
        self.assertAlmostEqual(intrinsic.intrinsic_depth(self.shape),3.9788735772973833,12)
        self.assertAlmostEqual(intrinsic.intrinsic_depth(facet),intrinsic.intrinsic_depth(self.shape),12)

    def test_go_first_scattering_is_unchanged(self):
        args=(self.shape,30.,.45,.42,.3,12,'clumping_dual')
        old=crown.primary(*args,'beam_diffusion_profile',self.scales)
        new=intrinsic.primary(*args,self.scales)
        np.testing.assert_array_equal(old[0],new[0])
        np.testing.assert_array_equal(old[1],new[1])
        self.assertEqual(old[3],new[3])

    def test_local_flux_integrates_directional_radiance_independently(self):
        # Higher-order BRF is azimuthally symmetric. Its independent dense
        # quadrature is compared with the operator's exited source energy.
        for rho,tau,soil in ((.45,.42,0.),(.45,.42,.3),(.5,.5,1.),(.06,.03,.18)):
            higher,budget=intrinsic.transport(self.shape,30.,rho,tau,soil,12,'clumping_dual',self.scales)
            first=crown.primary(self.shape,30.,rho,tau,soil,12,'clumping_dual','beam_diffusion_profile',self.scales)[1].sum()
            mu,w,_=crown.quadrature(64)
            from base.gap import gap_probability_crown_analytical as gap
            from base.optical_crown_profile import escape_ratios
            depth,att,_,_=intrinsic.source_profile(self.shape,rho,tau,12)
            flux=0.
            for m,weight in zip(mu,w):
                view=float(np.degrees(np.arccos(m)))
                gv=float(gap(*self.shape[:4],view))
                brf=higher[0]*(1-gv)/self.shape[0]*escape_ratios(view,depth,att)[0]+higher[2]*gv
                flux+=2*m*weight*brf
            self.assertAlmostEqual(first+flux,budget['reflected'],delta=3e-5)
            self.assertAlmostEqual(budget['balance_residual'],0.,12)
            if rho+tau==1 and soil==1:
                self.assertAlmostEqual(budget['reflected'],1.,12)

    def test_sparse_forest_invalid_escape_is_rejected_explicitly(self):
        with self.assertRaisesRegex(ValueError,'invalid escaped flux'):
            intrinsic.source_profile((1.5,.005,3.,3.,.025),.06,.03,12)

    def test_thin_crown_and_black_material_limits(self):
        shape=(.01,.02,3.,3.,.025)
        for rho,tau,soil in ((0.,0.,0.),(.06,.03,.18),(.45,.42,.3)):
            _,budget=intrinsic.transport(shape,30.,rho,tau,soil,12,'clumping_dual',self.scales)
            self.assertGreaterEqual(budget['reflected'],0.)
            self.assertLessEqual(budget['reflected'],1.)
            self.assertAlmostEqual(budget['balance_residual'],0.,12)
            if rho+tau+soil==0.:
                self.assertEqual(budget['reflected'],0.)

    def test_explicit_selection_has_no_cache_or_shared_profile_side_effect(self):
        settings=dict(self.settings)
        baseline=reflectance(settings,30.,0.,0.,0.,.45,.42,.3)
        old_profile=crown.transport(self.shape,30.,.45,.42,.3,12,'clumping_dual','beam_diffusion_profile',self.scales)
        candidate=reflectance(dict(settings,optical_terrain_source_depth='intrinsic_crown'),30.,0.,0.,0.,.45,.42,.3)
        self.assertGreater(candidate,baseline)
        self.assertEqual(baseline,reflectance(settings,30.,0.,0.,0.,.45,.42,.3))
        new_profile=crown.transport(self.shape,30.,.45,.42,.3,12,'clumping_dual','beam_diffusion_profile',self.scales)
        np.testing.assert_array_equal(old_profile[0],new_profile[0])
        self.assertEqual(old_profile[1],new_profile[1])

    def test_zero_terrain_limit_is_continuous_for_new_mode(self):
        settings=dict(self.settings,optical_terrain_source_depth='intrinsic_crown',terrain_density=0.)
        plane=reflectance(settings,30.,0.,0.,0.,.45,.42,.3)
        tiny=reflectance(dict(settings,terrain_density=1e-10),30.,0.,0.,0.,.45,.42,.3)
        self.assertAlmostEqual(plane,tiny,delta=1e-8)

    def test_bare_terrain_does_not_use_intrinsic_source(self):
        settings=dict(self.settings,lai=0.)
        old=reflectance(settings,30.,0.,30.,0.,.45,.42,.3)
        new=reflectance(dict(settings,optical_terrain_source_depth='intrinsic_crown'),30.,0.,30.,0.,.45,.42,.3)
        self.assertEqual(old,new)

    def test_zero_terrain_and_zero_lai_direct_diffuse_dark_limits(self):
        settings=dict(self.settings,lai=0.,terrain_density=0.,optical_terrain_source_depth='intrinsic_crown')
        for diffuse,sunlit in ((0.,True),(.4,True),(.4,False),(0.,False)):
            value=reflectance(settings,30.,0.,30.,0.,.45,.42,.3,diffuse,sunlit)
            self.assertAlmostEqual(value,.3*((1-diffuse)*sunlit+diffuse),15)


if __name__=='__main__':
    unittest.main()
