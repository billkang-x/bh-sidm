import math
import unittest

import numpy as np

from sidm_bh.conduction import cell_conductivity_code
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sidm import CONDUCTIVITY_A
from sidm_bh.state import FluidState
from sidm_bh.timescales import (
    inward_flux_median_radius_code,
    local_timescale_profiles_code,
)


class TimescaleProfilesTest(unittest.TestCase):
    def setUp(self):
        self.grid = SphericalGrid.from_log_spacing(0.1, 10.0, 8)
        self.state = FluidState(
            density=np.full(8, 2.0),
            radial_velocity=np.zeros(8),
            velocity_dispersion=np.full(8, 3.0),
        )

    def test_constant_temperature_timescales_match_definitions(self):
        sigma = 0.4
        black_hole_mass = 5.0
        result = local_timescale_profiles_code(
            self.state,
            self.grid,
            sigma,
            black_hole_mass_code=black_hole_mass,
        )
        conductivity = cell_conductivity_code(self.state, sigma)
        expected_collision = 1.0 / (CONDUCTIVITY_A * 2.0 * sigma * 3.0)

        np.testing.assert_allclose(result.collision_code, expected_collision)
        np.testing.assert_allclose(
            result.conduction_radius_code,
            1.5 * 2.0 * self.grid.centers_code**2 / conductivity,
        )
        np.testing.assert_allclose(
            result.conduction_gradient_code,
            result.conduction_radius_code,
        )
        np.testing.assert_allclose(
            result.thermal_length_code,
            self.grid.centers_code,
        )
        self.assertTrue(np.all(np.isinf(result.inflow_code)))

    def test_inward_flux_median_uses_log_radius_weights(self):
        velocity = np.zeros(8)
        velocity[2] = -1.0
        velocity[5] = -100.0
        state = FluidState(
            density=np.ones(8),
            radial_velocity=velocity,
            velocity_dispersion=np.ones(8),
        )
        median = inward_flux_median_radius_code(state, self.grid)
        self.assertEqual(median, self.grid.centers_code[5])

    def test_no_inward_flux_returns_nan(self):
        radius = inward_flux_median_radius_code(self.state, self.grid)
        self.assertTrue(math.isnan(radius))


if __name__ == "__main__":
    unittest.main()
