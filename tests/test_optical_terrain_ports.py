"""Physical controls for the optical-only whole-crown port experiment."""
import math
import unittest
from unittest.mock import patch
import numpy as np
from base import optical_crown as crown
from base import optical_terrain_ports as model

SETTINGS = dict(surface_model='terrain',vegetation_model='crown',surface_type=1,
    lai=1.5,hspot=.025,height=3.,width=3.,stand_density=.02,crown_bottom=3.,
    terrain_type=2,terrain_height=20.,terrain_radius=30.,terrain_density=.00015,
    optical_crown_geometry='clumping_dual',optical_crown_source='diffusion_profile',
    optical_crown_quadrature=4)


class FiniteCrownPortsTests(unittest.TestCase):
    def test_finite_cone_intersections(self):
        self.assertAlmostEqual(model.cone_hit(np.array([[0.,0.,25.]]),np.array([0.,0.,-1.]),20.,30.)[0],5.)
        # x=16-.6t and cone side x=1.5*(20-(10-.8t)); t=5/9.
        self.assertAlmostEqual(model.cone_hit(np.array([[16.,0.,10.]]),np.array([-.6,0.,-.8]),20.,30.)[0],5/9)
        # Above the apex, a globally upward ray can point below its slope's
        # tangent while seeing no finite terrain: the audit counterexample.
        origin = np.array([[312.25329,344.61922,14.97801]])
        direction = np.array([.479380141,-.721245229,.5])
        normal = np.array([-.307049028,.461966668,.832050294])
        radial = normal[:2]/np.linalg.norm(normal[:2])
        centre = origin[0,:2]-18.349641746*radial
        origin[0,:2] -= centre
        self.assertLess(direction@normal,0.)
        self.assertTrue(np.isinf(model.cone_hit(origin,direction,20.,30.)[0]))

    def test_conditioned_ports_partition_and_black_ground_escape(self):
        cs = math.cos(math.atan2(20,30))
        shape = (1.5*cs,.02*cs,3.,3.,.025)
        geo = (20.,30.,6.,3.)
        for attenuation in (0.,.4,2.):
            ports = model.ports(shape,*geo,attenuation,12)
            self.assertAlmostEqual(sum(ports),1.,places=12)
            self.assertTrue(all(0<=value<=1 for value in ports))
            self.assertGreater(ports[1],0.)
            self.assertGreater(ports[2],0.)
        for rho,tau,soil in ((0.,0.,0.),(.45,.42,0.),(.45,.42,.3),(.5,.5,1.)):
            _,_,exits,absorbed,_ = model.state(shape,30.,rho,tau,soil,4,'clumping_dual',(.025,1.),geo)
            self.assertAlmostEqual(sum(exits)+absorbed,1.,places=10)
            if soil==1:
                self.assertAlmostEqual(sum(exits),1.,places=10)
            if rho:
                self.assertGreater(exits[1],0.)
                self.assertGreater(exits[2],0.)

    def test_closed_ports_recover_original_local_transport(self):
        shape = (1.25,.016666666666,3.,3.,.025)
        args = (shape,30.,.45,.42,.3,4,'clumping_dual',(.025,1.))
        expected,budget = crown.transport(*args[:7],'beam_diffusion_profile',args[7])
        model.state.cache_clear()
        try:
            with patch.object(model,'ports',return_value=(1.,0.,0.)):
                actual,_,out,absorbed,_ = model.state(*args,(20.,30.,6.,3.))
                np.testing.assert_allclose(actual,expected,atol=1e-12)
                self.assertAlmostEqual(out[0],budget['reflected'],places=12)
                self.assertEqual(out[1],0.)
                self.assertEqual(out[2],0.)
        finally:
            model.state.cache_clear()

    def test_global_flux_conservation_with_external_canopy_exchange(self):
        with patch('base.optical_crown.ellipsoid_grid_crown_FRT',side_effect=AssertionError('volume integration')):
            for rho,tau,soil in ((.06,.03,.18),(.45,.42,0.),(.45,.42,.3),(.5,.5,1.)):
                for f in (0.,.4,1.):
                    value,budget = model.reflectance(SETTINGS,30.,0.,75.,90.,rho,tau,soil,f,return_budget=True)
                    self.assertTrue(math.isfinite(value) and value>=0)
                    self.assertLess(abs(budget['residual']),1e-10)
                    if rho+tau==1:
                        self.assertAlmostEqual(budget['escaped'],budget['incident'],places=10)
                        self.assertAlmostEqual(budget['absorbed'],0.,places=10)
                    if f==0:
                        self.assertAlmostEqual(budget['incident'],1.,places=10)

    def test_rotation_diffuse_linearity_and_zero_light(self):
        args = (SETTINGS,30.,0.,75.,90.,.45,.42,.3)
        beam = model.reflectance(*args)
        rotated = model.reflectance(SETTINGS,30.,73.,75.,163.,.45,.42,.3)
        self.assertAlmostEqual(beam,rotated,places=12)
        sky = model.reflectance(*args,1.)
        mixed = model.reflectance(*args,.4)
        self.assertAlmostEqual(mixed,.6*beam+.4*sky,places=12)
        self.assertEqual(model.reflectance(*args,sunlit=False),0.)
        self.assertEqual(model.reflectance(SETTINGS,30.,0.,75.,90.,0.,0.,0.),0.)

    def test_no_cones_and_invalid_material(self):
        settings = dict(SETTINGS,terrain_density=0.)
        expected = crown.crown_reflectance(settings,30.,25.,90.,.45,.42,.3)
        actual = model.reflectance(settings,30.,0.,25.,90.,.45,.42,.3)
        self.assertEqual(actual,expected)
        for rho,tau,soil in ((.6,.5,.3),(.4,.3,1.1),(-.1,.3,.3),(float('nan'),.3,.3)):
            with self.assertRaises(ValueError):
                model.reflectance(SETTINGS,30.,0.,25.,90.,rho,tau,soil)


if __name__=='__main__':
    unittest.main()
