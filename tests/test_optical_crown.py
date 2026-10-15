"""Physical contracts of the optical-only ellipsoidal crown transport."""
import unittest
from unittest.mock import patch

import numpy as np

from base.optical_crown import crown_reflectance, fractions, transport
from base.radiation import directional_reflectance
from base.brightness import vegetation_brightness


def settings(**changes):
    result = dict(surface_model='plane',vegetation_model='crown',surface_type=1,
                  distribution=0,lai=3.,stand_density=1/36,height=1.5,width=2.,hspot=.2/3,
                  leaf_emissivity=.98,soil_emissivity=.95,leaf_reflectance=.02,
                  leaf_transmittance=0.,soil_reflectance=.05,hapke_enable=0)
    result.update(changes)
    return result


class OpticalCrownTests(unittest.TestCase):
    def test_bare_and_black_limits(self):
        self.assertAlmostEqual(crown_reflectance(settings(lai=0),30,60,90,.45,.42,.3),.3)
        self.assertEqual(crown_reflectance(settings(),30,60,90,0,0,0),0)
        self.assertAlmostEqual(crown_reflectance(settings(lai=1e-7),30,0,0,.45,.42,.3),.3,places=6)

    def test_conservative_and_absorbing_energy_balance(self):
        for lai in (.1,1.,3.,6.):
            shape=(lai,1/36,1.5,2.,.2/3)
            for rho,tau,rs in ((.06,.03,.18),(.45,.42,.3),(.5,.5,1.),(0,0,0)):
                _, audit=transport(shape,30.,rho,tau,rs)
                self.assertLess(abs(audit['balance_residual']),1e-10)
                self.assertGreaterEqual(audit['absorbed'],-1e-10)
                if rho+tau==1 and rs==1:
                    self.assertAlmostEqual(audit['reflected'],1.,places=10)

    def test_joint_probabilities_and_exact_hotspot(self):
        shape=(3.,1/36,1.5,2.,.2/3)
        for v in (0.,30.,60.,89.5):
            for az in (0.,90.,180.):
                f=fractions(*shape,30.,v,az)
                self.assertGreaterEqual(min(f),0.)
                self.assertAlmostEqual(sum(f),1.)
        f=fractions(*shape,30.,30.,0.)
        self.assertAlmostEqual(f[1],0.,places=10)
        self.assertAlmostEqual(f[3],0.,places=10)
        np.testing.assert_allclose(f,fractions(*shape,30.+1e-12,30.,1e-12),atol=1e-12)
        np.testing.assert_allclose(f,fractions(*shape,29.999999999999996,30.,8.537736562175269e-7),atol=1e-12)

    def test_hotspot_retained_in_both_bands(self):
        for rho,tau,rs in ((.06,.03,.18),(.45,.42,.3)):
            vals=[crown_reflectance(settings(),30,v,0,rho,tau,rs) for v in (25,30,35)]
            self.assertGreater(vals[1],max(vals[0],vals[2]))

    def test_isotropic_sky_and_illumination_superposition(self):
        s=settings()
        direct=crown_reflectance(s,30,30,0,.45,.42,.3)
        sky=crown_reflectance(s,30,30,0,.45,.42,.3,1)
        self.assertEqual(sky,crown_reflectance(s,70,30,170,.45,.42,.3,1))
        self.assertAlmostEqual(crown_reflectance(s,30,30,0,.45,.42,.3,.3),.7*direct+.3*sky)
        self.assertAlmostEqual(crown_reflectance(s,30,30,0,.45,.42,.3,.3,False),.3*sky)
        self.assertEqual(crown_reflectance(s,30,30,0,.45,.42,.3,0,False),0.)

    def test_transmission_changes_directional_scattering(self):
        # Same single-scattering albedo, different reflection/transmission.
        a=crown_reflectance(settings(),30,30,0,.8,0,.3)
        b=crown_reflectance(settings(),30,30,0,0,.8,.3)
        self.assertGreater(a,b)

    def test_statistical_forest_has_no_azimuth_handedness(self):
        a=crown_reflectance(settings(),30,60,90,.45,.42,.3)
        b=crown_reflectance(settings(),30,60,270,.45,.42,.3)
        self.assertEqual(a,b)

    def test_optical_dispatch_and_legacy_option(self):
        s=settings()
        expected=crown_reflectance(s,30,30,0,.45,.42,.3)
        self.assertEqual(directional_reflectance(s,30,0,30,0,.45,.42,.3),expected)
        old=directional_reflectance(settings(optical_crown_method='effective_slab'),30,0,30,0,.45,.42,.3)
        self.assertNotAlmostEqual(old,expected,places=5)

    def test_thermal_path_never_calls_optical_transport(self):
        state=dict(Tss=315.,Tsh=295.,Tls=305.,Tlh=298.,sza=30.,saa=0.,sunlit_surface=True)
        with patch('base.optical_crown.crown_reflectance',side_effect=AssertionError('Thermal entered optical kernel')):
            for wave in (3.7,8.,10.5,12.):
                a=vegetation_brightness(settings(),state,wave,30.,0.)
                b=vegetation_brightness(settings(optical_crown_method='effective_slab'),state,wave,30.,0.)
                np.testing.assert_array_equal(a,b)

    def test_invalid_materials(self):
        for rho,tau,rs in ((.6,.6,.3),(-.1,.2,.3),(.4,.4,1.1),(np.nan,.2,.3)):
            with self.assertRaises(ValueError):
                crown_reflectance(settings(),30,0,0,rho,tau,rs)

    def test_clumping_uses_original_parameterized_geometry(self):
        from base.radiation import component_proportions
        s=settings(lai=2.9999961956235524)
        shape=tuple(s[k] for k in ('lai','stand_density','height','width','hspot'))
        for view in (0.,25.,30.,35.,60.):
            for az in (0.,90.,180.):
                np.testing.assert_allclose(fractions(*shape,30.,view,az,'clumping'),
                    component_proportions(s,30.,0.,view,az),atol=1e-14)

    def test_clumping_uniform_has_no_crown_spatial_integral(self):
        # A new geometry avoids cache hits concealing an unwanted volume call.
        s=settings(lai=2.713, optical_crown_geometry='clumping',optical_crown_quadrature=4)
        with patch('base.optical_crown.proportion_bidirectional_crown_voxel_one',side_effect=AssertionError('volume fractions')), \
             patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('spatial grid')), \
             patch('base.optical_crown.sunlit_down_escape',side_effect=AssertionError('spatial source')):
            direct=crown_reflectance(s,30.,30.,0.,.45,.42,.3)
            sky=crown_reflectance(s,30.,30.,0.,.45,.42,.3,1.)
            self.assertGreater(direct,0.)
            self.assertGreater(sky,0.)
            self.assertEqual(sky,crown_reflectance(s,60.,30.,90.,.45,.42,.3,1.))
            self.assertAlmostEqual(crown_reflectance(s,30.,30.,0.,.45,.42,.3,.3),.7*direct+.3*sky)

    def test_clumping_energy_and_hotspot(self):
        for source in ('uniform','sunlit_volume'):
            for lai in (.1,1.,3.,6.):
                shape=(lai,1/36,1.5,2.,.2/3)
                for rho,tau,rs in ((.06,.03,.18),(.45,.42,.3),(.5,.5,1.)):
                    _,audit=transport(shape,30.,rho,tau,rs,12,'clumping',source)
                    self.assertLess(abs(audit['balance_residual']),1e-10)
                    if rho+tau==1 and rs==1:
                        self.assertAlmostEqual(audit['reflected'],1.,places=10)
            s=settings(optical_crown_geometry='clumping',optical_crown_source=source)
            for rho,tau,rs in ((.06,.03,.18),(.45,.42,.3)):
                values=[crown_reflectance(s,30.,v,0.,rho,tau,rs) for v in (25.,30.,35.)]
                self.assertGreater(values[1],max(values[0],values[2]))

    def test_clumping_options_do_not_affect_thermal(self):
        state=dict(Tss=315.,Tsh=295.,Tls=305.,Tlh=298.,sza=30.,saa=0.,sunlit_surface=True)
        for source in ('uniform','sunlit_volume'):
            for wave in (3.7,8.,10.5,12.):
                np.testing.assert_array_equal(vegetation_brightness(settings(),state,wave,30.,0.),
                    vegetation_brightness(settings(optical_crown_geometry='clumping',optical_crown_source=source),state,wave,30.,0.))

    def test_invalid_geometry_options(self):
        for options in (dict(optical_crown_geometry='typo'),dict(optical_crown_source='typo')):
            with self.assertRaises(ValueError):
                crown_reflectance(settings(**options),30.,0.,0.,.45,.42,.3)

    def test_dual_hotspot_probability_limits_and_marginals(self):
        from base.gap import gap_probability_crown_analytical
        for lai in (.001,1.,3.,10.):
            shape=(lai,1/36,1.5,2.,.2/3)
            for sun in (0.,30.,75.,89.5):
                for view in (0.,30.,75.,89.5):
                    for az in (0.,90.,180.):
                        p=fractions(*shape,sun,view,az,'clumping_dual')
                        gv=float(gap_probability_crown_analytical(*shape[:4],view))
                        gs=float(gap_probability_crown_analytical(*shape[:4],sun))
                        self.assertGreaterEqual(min(p),-1e-15)
                        self.assertAlmostEqual(sum(p),1.,places=14)
                        self.assertAlmostEqual(p[0]+p[1],gv,places=14)
                        self.assertAlmostEqual(p[2]+p[3],1-gv,places=14)
                        self.assertGreaterEqual(p[0]+1e-15,max(0.,gv+gs-1))
                        self.assertLessEqual(p[0],min(gv,gs)+1e-15)
                        reverse=fractions(*shape,view,sun,az,'clumping_dual')
                        self.assertAlmostEqual(p[0],reverse[0],places=14)
                p=fractions(*shape,sun,sun,0.,'clumping_dual')
                self.assertAlmostEqual(p[1],0.,places=14)
                self.assertAlmostEqual(p[3],0.,places=14)

    def test_dual_hotspot_inner_and_outer_have_independent_effects(self):
        shape=(2.9999961956,1/36,1.5,2.,.2/3)
        base=np.array(fractions(*shape,30.,25.,0.,'clumping_dual',(.2/3,2/1.5)))
        for scales in ((.2/6,2/1.5),(.2/3,1/1.5)):
            changed=np.array(fractions(*shape,30.,25.,0.,'clumping_dual',scales))
            self.assertGreater(abs(base-changed).max(),1e-5)
            self.assertAlmostEqual(changed[0]+changed[1],base[0]+base[1])
            self.assertAlmostEqual(changed[2]+changed[3],base[2]+base[3])
        for scales in ((0.,0.),(0.,1.),(.1,0.)):
            p=fractions(*shape,30.,25.,0.,'clumping_dual',scales)
            self.assertTrue(np.isfinite(p).all())
            self.assertGreaterEqual(min(p),0.)

    def test_dual_hotspot_no_spatial_integral_or_old_hotspot(self):
        s=settings(lai=2.819,optical_crown_geometry='clumping_dual',optical_crown_quadrature=4)
        with patch('base.optical_crown.proportion_bidirectional_crown_voxel_one',side_effect=AssertionError('spatial fractions')), \
             patch('base.optical_crown.proportion_bidirectional_crown_one',side_effect=AssertionError('old hotspot')), \
             patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('spatial grid')), \
             patch('base.optical_crown.sunlit_down_escape',side_effect=AssertionError('spatial source')):
            beam=crown_reflectance(s,30.,25.,0.,.45,.42,.3)
            sky=crown_reflectance(s,30.,25.,0.,.45,.42,.3,1.)
            self.assertEqual(sky,crown_reflectance(s,70.,25.,90.,.45,.42,.3,1.))
            self.assertAlmostEqual(crown_reflectance(s,30.,25.,0.,.45,.42,.3,.4),.6*beam+.4*sky)
            self.assertGreater(beam,0.)

    def test_dual_hotspot_energy_balance(self):
        for lai in (.1,1.,3.,6.):
            for sun in (0.,30.,75.):
                for rho,tau,rs in ((.06,.03,.18),(.45,.42,.3),(.5,.5,1.)):
                    _,budget=transport((lai,1/36,1.5,2.,.2/3),sun,rho,tau,rs,12,'clumping_dual','uniform')
                    self.assertLess(abs(budget['balance_residual']),1e-10)
                    if rho+tau==1 and rs==1:self.assertAlmostEqual(budget['reflected'],1.,places=10)

    def test_dual_hotspot_dispatch_parameters_and_thermal_isolation(self):
        s=settings(optical_crown_geometry='clumping_dual')
        a=crown_reflectance(s,30.,25.,0.,.45,.42,.3)
        self.assertEqual(a,directional_reflectance(s,30.,0.,25.,0.,.45,.42,.3))
        for option in ('optical_hotspot_inner','optical_hotspot_outer'):
            b=crown_reflectance(dict(s,**{option:.02}),30.,25.,0.,.45,.42,.3)
            self.assertNotAlmostEqual(a,b,places=6)
        state=dict(Tss=315.,Tsh=295.,Tls=305.,Tlh=298.,sza=30.,saa=0.,sunlit_surface=True)
        with patch('base.optical_crown.dual_hotspot_fractions',side_effect=AssertionError('Thermal entered dual hotspot')):
            for wave in (3.7,8.,10.5,12.):
                np.testing.assert_array_equal(vegetation_brightness(settings(),state,wave,30.,0.),
                                              vegetation_brightness(s,state,wave,30.,0.))

    def test_dual_hotspot_invalid_options(self):
        for options in (dict(optical_hotspot_inner=-1),dict(optical_hotspot_outer=np.nan),
                        dict(optical_hotspot_outer=float('inf')),dict(height=0),
                        dict(optical_crown_source='sunlit_volume')):
            with self.assertRaises(ValueError):
                crown_reflectance(settings(optical_crown_geometry='clumping_dual',**options),30.,0.,0.,.45,.42,.3)

    def test_dual_hotspot_overlap_and_depth_closed_forms(self):
        from base.optical_hotspot import attenuation_mean,disk_overlap
        self.assertEqual(attenuation_mean(0.),1.)
        self.assertAlmostEqual(disk_overlap(np.pi,np.pi,0.),np.pi)
        self.assertEqual(disk_overlap(np.pi,np.pi,2.),0.)
        self.assertAlmostEqual(disk_overlap(np.pi,np.pi,1.),2*np.pi/3-np.sqrt(3)/2)
        self.assertAlmostEqual(disk_overlap(4*np.pi,np.pi,.5),np.pi)

    def test_independent_hemisphere_integral(self):
        # A different integration grid checks the directional observation
        # against the independently balanced flux system, including grazing rays.
        mu,w=np.polynomial.legendre.leggauss(20)
        integral=0.
        for m,weight in zip((mu+1)/2,w/2):
            for az in (np.arange(48)+.5)*360/48:
                integral+=2*m*weight/48*crown_reflectance(settings(),30,float(np.rad2deg(np.arccos(m))),az,.5,.5,1.)
        self.assertLess(abs(integral-1),1e-3)


if __name__=='__main__':
    unittest.main()
