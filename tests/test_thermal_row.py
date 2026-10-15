"""Empty corridor, probability, thermal and optical isolation controls."""
import math
import unittest
import numpy as np
from base.thermal_row import conditioned_proportions
from base.vegetation import vegetation_proportions
from base.row import Row
from base.physicsF import planck


def probabilities(view,azimuth,solar=30.,sun_azimuth=0.,row=0.,blank=.5):
    angles = (np.atleast_1d(view),np.atleast_1d(azimuth),np.atleast_1d(solar),np.atleast_1d(sun_azimuth))
    old = vegetation_proportions('row',1.5,.15,*angles,row_width=.5,row_blank=blank,row_height=1.,row_azimuth=row)
    return old,conditioned_proportions(1.5,.15,.5,blank,1.,*angles,row,old)


class ThermalRowTests(unittest.TestCase):
    def test_along_row_restores_empty_soil_at_grazing(self):
        for view,azimuth in ((0.,0.),(30.,0.),(75.,0.),(75.,180.)):
            old,new = probabilities(view,azimuth)
            gs = .5+.5*math.exp(-.5*3/math.cos(math.radians(30)))
            gv = .5+.5*math.exp(-.5*3/math.cos(math.radians(view)))
            self.assertGreaterEqual(new[0,0],.5)
            self.assertLessEqual(new[0,0],min(gs,gv)+1e-12)
            self.assertAlmostEqual(new[:2,0].sum(),gv,places=12)
            self.assertAlmostEqual(new[2:,0].sum(),1-gv,places=12)
            self.assertTrue((new>=0).all())
            self.assertAlmostEqual(new[:,0].sum(),1.,places=12)

    def test_coincident_rays_and_closed_corridors(self):
        _,p = probabilities(30.,0.)
        self.assertLess(abs(p[1,0])+abs(p[3,0]),1e-12)
        for view,azimuth,row in ((75.,90.,0.),(0.,0.,90.),(75.,0.,90.)):
            old,new = probabilities(view,azimuth,row=row)
            np.testing.assert_array_equal(new,old)
        old,new = probabilities(75.,0.,blank=0.)
        np.testing.assert_array_equal(new,old)

    def test_rotation_batch_and_isothermal_invariance(self):
        old,p = probabilities([0.,30.,75.],[0.,0.,180.])
        _,rotated = probabilities([0.,30.,75.],[73.,73.,253.],sun_azimuth=73.,row=73.)
        np.testing.assert_allclose(p,rotated,atol=1e-12)
        model = Row()
        model.set_structure(1.5,.15,.5,.5,1.)
        model.set_angle(np.array([0.,30.,75.]),np.array([30.]),np.array([0.,0.,180.]),np.array([0.]),np.array([0.]))
        model.set_thermal(300.,300.,300.,300.)
        for emissivity in (1.,.95):
            model.set_optical(10.5,emissivity,emissivity)
            model.set_geometry('legacy'); before = model.run(ifradiance=1).copy()
            model.set_geometry('conditioned_corridor'); after = model.run(ifradiance=1)
            np.testing.assert_allclose(before,after,atol=1e-12)
            if emissivity==1:
                np.testing.assert_allclose(after,planck(10.5,300.),atol=1e-12)

    def test_geometry_is_material_independent_and_legacy_replay(self):
        model = Row()
        model.set_structure(1.5,.15,.5,.5,1.)
        model.set_angle(np.array([75.]),np.array([30.]),np.array([0.]),np.array([0.]),np.array([0.]))
        model.set_thermal(318.15,303.15,306.15,303.15)
        model.set_optical(10.5,.955,.985)
        model.set_geometry('legacy'); old = model.run(ifradiance=1).copy()
        model.set_geometry('conditioned_corridor'); new = model.run(ifradiance=1)
        self.assertGreater(float(new[0]),float(old[0]))
        model.set_geometry('legacy')
        np.testing.assert_array_equal(model.run(ifradiance=1),old)
        with self.assertRaises(ValueError):
            model.set_geometry('unknown')


if __name__=='__main__':
    unittest.main()
