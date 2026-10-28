"""Physical controls for heat capacity, implicit conduction and interval history."""
import math
import unittest

from base.evapo import (
    _depth_differential, _envelope_differential, _heat_flux,
    _surface_differential, _thermal_inertia, update_thermal_history,
)
from base.heat_conduction import COLUMN_DEPTH, column_layers


def column_storage(profile, old, capacity):
    """Energy stored by the finite-volume column relative to the old profile [J m-2]."""
    layers = column_layers()
    widths = [layers[0] / 2] + [(a + b) / 2 for a, b in zip(layers, layers[1:])] + [layers[-1] / 2]
    return capacity * sum(w * (new - before) for w, new, before in zip(widths, profile, old))


def state_at(surface=330.0, dt=10.0):
    state = dict(Tss=surface, Tsh=surface, Troof_sunlit=surface,
                 Troof_shaded=surface, Twall_sunlit=surface,
                 Twall_shaded=surface, dt_seconds=dt,
                 thermal_history_dt_seconds=[])
    for name in ("sunlit", "shaded", "roof_sunlit", "roof_shaded",
                 "wall_sunlit", "wall_shaded", "street_sunlit", "street_shaded"):
        state["depth_profile_" + name] = [300.0] * 4
        state["surface_history_" + name] = [300.0]
    return state


class HeatConductionTests(unittest.TestCase):
    def setUp(self):
        self.settings = dict(soil_cs=1000.0, soil_rhos=1.0, soil_lambda=1.0,
                             soil_depth_step=0.1, urban_cs=1000.0,
                             urban_rhos=1800.0, urban_lambda=1.0,
                             urban_depth_step=0.1, urban_wall_layers=3,
                             urban_wall_thickness=0.3,
                             urban_internal_temperature_c=26.85)

    def test_urban_inertia_uses_volumetric_capacity_without_density(self):
        self.assertAlmostEqual(_thermal_inertia(self.settings, urban=True), math.sqrt(1000))
        self.settings["urban_rhos"] = 3600
        self.assertAlmostEqual(_thermal_inertia(self.settings, urban=True), math.sqrt(1000))

    def test_soil_density_changes_depth_diffusion_and_surface_inertia(self):
        light, heavy = state_at(), state_at()
        update_thermal_history(self.settings, light)
        update_thermal_history(dict(self.settings, soil_rhos=2.0), heavy)
        self.assertGreater(light["depth_profile_sunlit"][1], heavy["depth_profile_sunlit"][1])
        self.assertGreater(heavy["depth_profile_sunlit"][1], 300.0)
        self.assertAlmostEqual(_thermal_inertia(dict(self.settings, soil_rhos=2.0)) /
                               _thermal_inertia(self.settings), math.sqrt(2))

    def test_flux_uses_same_new_boundary_profile_that_update_commits(self):
        state = state_at()
        original = state["depth_profile_sunlit"][:]
        flux = _depth_differential(self.settings, state, 330.0, "depth_profile_sunlit")
        self.assertEqual(state["depth_profile_sunlit"], original)
        update_thermal_history(self.settings, state)
        profile = state["depth_profile_sunlit"]
        top = column_layers()[0]
        self.assertEqual(profile[0], 330.0)
        self.assertAlmostEqual(flux, (330-profile[1])/top + 1000*top/2*(330-300)/10)

    def test_long_step_solves_actual_diffusion_without_fourier_clipping(self):
        state = state_at(dt=1.0e8)
        update_thermal_history(self.settings, state)
        profile = state["depth_profile_sunlit"]
        self.assertTrue(all(300 <= value <= 330 for value in profile))
        # A zero-flux bottom lets a long step equilibrate the whole column.
        for value in profile:
            self.assertAlmostEqual(value, 330, delta=1e-2)

    def test_larger_time_step_advances_subsurface_further(self):
        short, long = state_at(dt=1), state_at(dt=100)
        update_thermal_history(self.settings, short)
        update_thermal_history(self.settings, long)
        self.assertGreater(long["depth_profile_sunlit"][1], short["depth_profile_sunlit"][1])

    def test_constant_profile_remains_constant(self):
        state = state_at(surface=300, dt=1e8)
        update_thermal_history(self.settings, state)
        for key, profile in state.items():
            if key.startswith("depth_profile_"):
                for value in profile:
                    self.assertAlmostEqual(value, 300)
        self.assertAlmostEqual(_depth_differential(self.settings, state, 300, "depth_profile_sunlit"), 0)

    def test_wall_flux_and_committed_profile_share_implicit_indoor_boundary(self):
        state = state_at()
        flux = _envelope_differential(self.settings, state, 330, "depth_profile_wall_sunlit")
        self.assertAlmostEqual(flux, 187.5)
        self.assertAlmostEqual(state["G_wall_inner_sunlit"], 37.5)
        update_thermal_history(self.settings, state)
        self.assertAlmostEqual(state["depth_profile_wall_sunlit"][1], 311.25)
        self.assertAlmostEqual(state["depth_profile_wall_sunlit"][-1], 300)

    def test_profile_storage_equals_net_conductive_energy(self):
        state = state_at()
        surface_flux = _depth_differential(self.settings, state, 330, "depth_profile_sunlit")
        update_thermal_history(self.settings, state)
        profile = state["depth_profile_sunlit"]
        self.assertAlmostEqual(column_storage(profile, [300.0] * len(profile), 1000) / 10, surface_flux)

    def test_steady_wall_recovers_conductivity_over_full_thickness(self):
        state = state_at(dt=1e8)
        state["depth_profile_wall_sunlit"] = [330, 320, 310, 300]
        flux = _envelope_differential(self.settings, state, 330, "depth_profile_wall_sunlit")
        self.assertAlmostEqual(flux, 100)
        update_thermal_history(self.settings, state)
        for actual, expected in zip(state["depth_profile_wall_sunlit"], [330, 320, 310, 300]):
            self.assertAlmostEqual(actual, expected)

    def test_single_layer_wall_connects_current_surface_to_indoor_boundary(self):
        state = state_at()
        settings = dict(self.settings, urban_wall_layers=1)
        flux = _envelope_differential(settings, state, 330, "depth_profile_wall_sunlit")
        self.assertAlmostEqual(flux, 100)
        update_thermal_history(settings, state)
        self.assertEqual(state["depth_profile_wall_sunlit"], [330, 300])

    def test_street_uses_same_soil_capacity_for_flux_and_profile(self):
        settings = dict(self.settings, soil_rhos=100)
        state = state_at()
        flux = _heat_flux(settings, state, 330, 100, "depth", "surface_history_street_sunlit",
                          "depth_profile_street_sunlit", 0.25, urban=False)
        soil_flux = _depth_differential(settings, state, 330, "depth_profile_sunlit")
        urban_flux = _depth_differential(settings, state, 330, "depth_profile_sunlit", urban=True)
        update_thermal_history(settings, state)
        self.assertAlmostEqual(flux, soil_flux)
        self.assertNotAlmostEqual(flux, urban_flux)
        self.assertEqual(state["depth_profile_street_sunlit"], state["depth_profile_sunlit"])
        profile = state["depth_profile_street_sunlit"]
        self.assertAlmostEqual(column_storage(profile, [300.0] * len(profile), 100 * 1000) / 10, flux)

    def test_column_is_two_metres_deep_and_refines_toward_the_surface(self):
        layers = column_layers()
        self.assertAlmostEqual(sum(layers), COLUMN_DEPTH)
        self.assertAlmostEqual(COLUMN_DEPTH, 2.0)
        self.assertLessEqual(layers[0], 0.01)
        self.assertTrue(all(b >= a for a, b in zip(layers, layers[1:])))
        self.assertLessEqual(len(layers) + 1, 25)

    def test_sinusoidal_surface_has_zero_mean_flux_and_analytic_amplitude(self):
        settings = dict(soil_cs=1180.0, soil_rhos=1100.0, soil_lambda=1.55)
        inertia = math.sqrt(1.55 * 1180.0 * 1100.0)
        omega, amplitude, mean, dt = 2 * math.pi / 86400, 10.0, 298.15, 1800.0
        state = state_at(surface=mean, dt=dt)
        for name in ("sunlit", "shaded", "street_sunlit", "street_shaded"):
            state["depth_profile_" + name] = [mean]
        fluxes = []
        for k in range(1, 8 * 48 + 1):
            surface = mean + amplitude * math.sin(omega * k * dt)
            flux = _depth_differential(settings, state, surface, "depth_profile_sunlit")
            state["Tss"] = state["Tsh"] = surface
            update_thermal_history(settings, state)
            if k > 6 * 48:
                fluxes.append(flux)
        self.assertLess(abs(sum(fluxes) / len(fluxes)), 0.5)
        exact = inertia * amplitude * math.sqrt(omega)
        self.assertAlmostEqual((max(fluxes) - min(fluxes)) / 2 / exact, 1.0, delta=0.02)

    def test_bottom_is_zero_flux_not_the_initial_temperature(self):
        state = state_at(surface=330, dt=3600)
        for _ in range(2000):
            update_thermal_history(self.settings, state)
        self.assertAlmostEqual(state["depth_profile_sunlit"][-1], 330, places=3)

    def test_surface_convolution_uses_previous_interval_lengths(self):
        state = dict(dt_seconds=1.0, thermal_history_dt_seconds=[4.0],
                     surface_history_sunlit=[302.0, 300.0])
        settings = dict(soil_cs=1, soil_rhos=1, soil_lambda=1)
        flux = _surface_differential(settings, state, 304, "surface_history_sunlit")
        expected = 2/math.sqrt(math.pi) * (2 + 0.5*(math.sqrt(5)-1))
        self.assertAlmostEqual(flux, expected)

    def test_surface_convolution_truncates_at_7_hours_in_elapsed_time(self):
        state = dict(dt_seconds=47*3600, thermal_history_dt_seconds=[3600, 3600],
                     surface_history_sunlit=[300, 290, 280])
        settings = dict(soil_cs=1, soil_rhos=1, soil_lambda=1)
        flux = _surface_differential(settings, state, 304, "surface_history_sunlit")
        expected = 2/math.sqrt(math.pi)*(4/(47*3600))*math.sqrt(7*3600)
        self.assertAlmostEqual(flux, expected)

    def test_history_retains_7_hours_and_interval_boundary_sample(self):
        state = state_at(surface=300, dt=1800)
        for _ in range(100):
            update_thermal_history(self.settings, state)
        self.assertEqual(len(state["surface_history_sunlit"]), 15)
        self.assertEqual(len(state["thermal_history_dt_seconds"]), 14)
        self.assertEqual(sum(state["thermal_history_dt_seconds"]), 7*3600)

    def test_variable_interval_history_preserves_the_true_timing(self):
        state = state_at(surface=300, dt=3600)
        for dt in (3600, 7200, 900):
            state["dt_seconds"] = dt
            update_thermal_history(self.settings, state)
        self.assertEqual(state["thermal_history_dt_seconds"], [900, 7200, 3600])
        self.assertEqual(len(state["surface_history_sunlit"]), 4)


if __name__ == "__main__":
    unittest.main()
