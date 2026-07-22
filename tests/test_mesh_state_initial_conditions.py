import unittest

import numpy as np

from sidm_bh.halos import NFWProfile
from sidm_bh.initial_conditions import (
    hydrostatic_pressure_code,
    hydrostatic_state_from_profile,
    isothermal_state_from_profile,
)
from sidm_bh.halos import SingularIsothermalSphere
from sidm_bh.mesh import SphericalGrid
from sidm_bh.state import FluidState
from sidm_bh.units import SimulationScales


class SphericalGridTest(unittest.TestCase):
    def test_volume_sum_matches_domain_volume(self):
        grid = SphericalGrid.from_log_spacing(1.0e-3, 10.0, 64)

        self.assertAlmostEqual(
            grid.cell_volumes_code.sum(),
            (10.0**3 - 1.0e-3**3) / 3.0,
        )
        self.assertEqual(grid.num_cells, 64)


class FluidStateTest(unittest.TestCase):
    def test_conservative_variables(self):
        state = FluidState(
            density=np.array([2.0, 3.0]),
            radial_velocity=np.array([0.5, -1.0]),
            velocity_dispersion=np.array([4.0, 2.0]),
        )

        conservative = state.conservative
        self.assertEqual(conservative.shape, (3, 2))
        np.testing.assert_allclose(conservative[0], [2.0, 3.0])
        np.testing.assert_allclose(conservative[1], [1.0, -3.0])
        np.testing.assert_allclose(conservative[2], state.density * state.specific_energy)


class HydrostaticInitialConditionTest(unittest.TestCase):
    def test_natural_sis_scales_give_constant_isothermal_state(self):
        profile = SingularIsothermalSphere(sound_speed_km_s=4.2)
        scales = SimulationScales.for_singular_isothermal_sphere(4.2)
        grid = SphericalGrid.from_log_spacing(1.0e-3, 1.0e3, 64)

        state = isothermal_state_from_profile(profile, grid, scales, 4.2)

        self.assertAlmostEqual(scales.velocity_scale_km_s, np.sqrt(2.0) * 4.2)
        np.testing.assert_allclose(
            state.velocity_dispersion,
            np.full(grid.num_cells, 1.0 / np.sqrt(2.0)),
        )
        np.testing.assert_allclose(
            state.density,
            1.0 / grid.centers_code**2,
        )

    def test_isothermal_power_law_can_be_recovered_with_outer_pressure(self):
        grid = SphericalGrid.from_log_spacing(1.0, 100.0, 1024)
        radius = grid.centers_code
        velocity_squared = 3.0
        density = 2.0 * velocity_squared / radius**2
        enclosed_mass = 2.0 * velocity_squared * radius
        outer_radius = grid.interfaces_code[-1]
        outer_pressure = 2.0 * velocity_squared**2 / outer_radius**2

        pressure = hydrostatic_pressure_code(
            grid,
            density,
            enclosed_mass,
            outer_pressure_code=outer_pressure,
        )

        np.testing.assert_allclose(pressure / density, velocity_squared, rtol=5.0e-3)

    def test_nfw_hydrostatic_state_has_positive_velocity_dispersion(self):
        profile = NFWProfile(scale_density_msun_pc3=3.7, scale_radius_pc=30.0)
        scales = SimulationScales(radius_scale_pc=30.0, density_scale_msun_pc3=3.7)
        grid = SphericalGrid.from_log_spacing(0.01, 100.0, 128)

        state = hydrostatic_state_from_profile(profile, grid, scales)

        self.assertEqual(len(state.density), grid.num_cells)
        self.assertTrue(np.all(state.velocity_dispersion > 0.0))


if __name__ == "__main__":
    unittest.main()
