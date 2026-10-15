"""Physical boundary and unit checks for the explicit VNIR interface."""
import unittest

import numpy as np

from base.optical import reflectance
from optical_run import simulate
from base.radiation import directional_reflectance
from base.vnir import joint_transmission, diffuse_field


def bare():
    return dict(surface_model="plane", vegetation_model="bare", surface_type=0,
                distribution=0, lai=0, wavelengths=[670,860],
                leaf_reflectance_spectrum=[0.06,0.45],
                leaf_transmittance_spectrum=[0.03,0.42],
                soil_reflectance_spectrum=[0.18,0.30])


class OpticalTests(unittest.TestCase):
    def canopy(self, **changes):
        settings = dict(bare(),lai=3,hspot=0.1,vegetation_model='hom')
        settings.update(changes)
        return settings

    def test_hotspot_in_both_bands(self):
        for wave in (670,860):
            s=self.canopy()
            peak=reflectance(s,wave,30,0,30,0)
            shoulders=[reflectance(s,wave,v,0,30,0) for v in (25,35)]
            self.assertGreater(peak,max(shoulders)+0.02)

    def test_hotspot_width_responds_to_parameter(self):
        near=[reflectance(self.canopy(hspot=h),670,28,0,30,0) for h in (0.01,0.1,0.5)]
        self.assertTrue(near[0]<near[1]<near[2])
        peaks=[reflectance(self.canopy(hspot=h),670,30,0,30,0) for h in (0.01,0.1,0.5)]
        np.testing.assert_allclose(peaks,peaks[0],atol=1e-12)

    def test_joint_path_limits_and_reciprocity(self):
        k=0.5/np.cos(np.deg2rad(30))
        bottom,path=joint_transmission(3,k,k,0.1,30,30,0)
        self.assertAlmostEqual(bottom,np.exp(-k*3))
        self.assertAlmostEqual(path,(1-np.exp(-k*3))/k)
        bottom,path=joint_transmission(3,k,k,0,30,30,0)
        self.assertAlmostEqual(bottom,np.exp(-2*k*3))
        self.assertAlmostEqual(path,(1-np.exp(-2*k*3))/(2*k))
        a=joint_transmission(3,0.6,0.8,0.1,30,50,20)
        b=joint_transmission(3,0.8,0.6,0.1,50,30,20)
        np.testing.assert_allclose(a,b,atol=1e-12)

    def test_diffuse_sky_has_no_solar_hotspot(self):
        s=self.canopy(shortwave_diffuse_fraction=1)
        values=[reflectance(s,860,30,az,30,0) for az in (0,90,180)]
        np.testing.assert_allclose(values,values[0],atol=1e-12)
        self.assertAlmostEqual(values[0],reflectance(s,860,30,0,60,90),places=10)

    def test_illumination_superposition_and_shade(self):
        s=self.canopy()
        direct=directional_reflectance(s,30,0,30,0,.45,.42,.3,diffuse_fraction=0)
        sky=directional_reflectance(s,30,0,30,0,.45,.42,.3,diffuse_fraction=1)
        mixed=directional_reflectance(s,30,0,30,0,.45,.42,.3,diffuse_fraction=.3)
        self.assertAlmostEqual(mixed,.7*direct+.3*sky)
        shaded=directional_reflectance(s,30,0,30,0,.45,.42,.3,diffuse_fraction=.3,sunlit=False)
        self.assertAlmostEqual(shaded,.3*sky)

    def test_black_canopy_and_bare_limits(self):
        s=self.canopy(leaf_reflectance_spectrum=[0,0],leaf_transmittance_spectrum=[0,0],soil_reflectance_spectrum=[0,0])
        self.assertEqual(reflectance(s,670,30,0,30,0),0)
        s=self.canopy(lai=1e-9)
        self.assertAlmostEqual(reflectance(s,670,30,0,30,0),.18,places=8)

    def test_conservative_flux_and_absorbing_soil(self):
        # Exact hemispheric flux balance, independent of directional BRF clipping.
        for lai in (.1,3,20):
            for sky in (0,.4,1):
                *_,top=diffuse_field(lai,.6,.5,.5,1,1-sky,sky)
                self.assertAlmostEqual(top,1,places=10)
                *_,ground,top=diffuse_field(lai,.6,.5,.5,0,1-sky,sky)
                absorbed_soil=ground+(1-sky)*np.exp(-.6*lai)
                self.assertAlmostEqual(top+absorbed_soil,1,places=10)

    def test_diffuse_field_layer_convergence(self):
        a=diffuse_field(6,.6,.45,.42,.3,1,0,64)
        b=diffuse_field(6,.6,.45,.42,.3,1,0,256)
        np.testing.assert_allclose(a[-2:],b[-2:],atol=1e-10)
        def outgoing(field):
            z,d,u,_,_=field
            return np.sum((.4*d+.3*u)*np.exp(-.7*z))*6/len(z)
        self.assertLess(abs(outgoing(a)-outgoing(b)),2e-4)

    def test_direct_single_scattering_reciprocity(self):
        # Negligible multiple scattering isolates the reciprocal phase/path.
        s=self.canopy(leaf_reflectance_spectrum=[1e-8]*2,
                      leaf_transmittance_spectrum=[0]*2,
                      soil_reflectance_spectrum=[0]*2)
        a=reflectance(s,670,50,20,30,0)
        b=reflectance(s,670,30,0,50,20)
        self.assertAlmostEqual(a/b,1,places=6)

    def test_lambertian_bare_control(self):
        for wave,rho in ((670,0.18),(860,0.30)):
            for vza in (0,30,60):
                for vaa in (0,90,180,270):
                    self.assertAlmostEqual(reflectance(bare(),wave,vza,vaa,30,0),rho)

    def test_black_ground(self):
        settings=bare()
        settings["soil_reflectance_spectrum"]=[0,0]
        self.assertEqual(reflectance(settings,670,45,180,30,0),0)

    def test_nm_and_thermal_separation(self):
        for wave in (0.67,10.5,10500,float("nan")):
            with self.assertRaises(ValueError):
                reflectance(bare(),wave,0,0,30,0)

    def test_conservation_and_missing_spectrum(self):
        settings=bare()
        settings["leaf_transmittance_spectrum"]=[0.96,0.7]
        with self.assertRaises(ValueError):
            reflectance(settings,670,0,0,30,0)
        settings=bare()
        settings["soil_reflectance_spectrum"]=[0.18]
        with self.assertRaises(ValueError):
            reflectance(settings,670,0,0,30,0)

    def test_local_geometry_scope(self):
        for surface in ("slope","terrain"):
            settings=bare()
            settings["surface_model"]=surface
            with self.assertRaises(ValueError):
                reflectance(settings,670,0,0,30,0)
        with self.assertRaises(ValueError):
            reflectance(bare(),670,90,0,30,0)

    def test_interpolation_does_not_extrapolate(self):
        self.assertAlmostEqual(reflectance(bare(),765,0,0,30,0),0.24)
        with self.assertRaises(ValueError):
            reflectance(bare(),500,0,0,30,0)

    def test_two_planes_signed_angles(self):
        rows=list(simulate(dict(settings=bare(),sza=30,saa=355,max_vza=60,vza_step=5)))
        self.assertEqual(len(rows),100)
        self.assertEqual(len({(r['plane'],r['signed_vza'],r['wavelength_nm']) for r in rows}),100)
        positive=next(r for r in rows if r['plane']=='parallel' and r['signed_vza']==30)
        self.assertEqual(positive['vaa'],355)
        negative=next(r for r in rows if r['plane']=='perpendicular' and r['signed_vza']==-30)
        self.assertEqual(negative['vaa'],265)


if __name__=='__main__':
    unittest.main()
