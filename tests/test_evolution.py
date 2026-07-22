import unittest

import numpy as np

from sidm_bh.accretion import inner_boundary_sidm_accretion_rate_code
from sidm_bh.evolution import evolve_sidm, fluid_mass_code
from sidm_bh.mesh import SphericalGrid
from sidm_bh.solver import stable_timestep_code
from sidm_bh.state import FluidState


class EvolutionDriverTest(unittest.TestCase):
    def setUp(self):
        self.grid = SphericalGrid(np.array([1.0, 1.5, 2.0]))

    def test_inner_accretion_updates_black_hole_and_closes_mass_budget(self):
        state = FluidState(
            density=np.array([2.0, 1.0]),
            radial_velocity=np.array([-0.01, 0.0]),
            velocity_dispersion=np.ones(2),
        )
        black_hole_mass = 0.2
        end_time = 0.5 * stable_timestep_code(
            state,
            self.grid,
            black_hole_mass_code=black_hole_mass,
        )
        expected_accreted = (
            inner_boundary_sidm_accretion_rate_code(state, self.grid) * end_time
        )

        result = evolve_sidm(
            state,
            self.grid,
            end_time,
            initial_black_hole_mass_code=black_hole_mass,
        )

        self.assertEqual(result.history.num_steps, 1)
        self.assertAlmostEqual(result.history.times_code[-1], end_time)
        self.assertAlmostEqual(
            result.final_black_hole_mass_code,
            black_hole_mass + expected_accreted,
        )
        self.assertAlmostEqual(
            result.history.max_absolute_mass_budget_residual_code,
            0.0,
            places=13,
        )

    def test_signed_outer_flux_is_included_in_budget(self):
        state = FluidState(
            density=np.array([1.0, 2.0]),
            radial_velocity=np.array([0.0, 0.02]),
            velocity_dispersion=np.ones(2),
        )
        initial_mass = fluid_mass_code(state, self.grid)
        end_time = 1.0e-3

        result = evolve_sidm(state, self.grid, end_time)

        self.assertLess(fluid_mass_code(result.final_state, self.grid), initial_mass)
        self.assertGreater(result.history.cumulative_outer_flux_code[-1], 0.0)
        self.assertAlmostEqual(
            result.history.max_absolute_mass_budget_residual_code,
            0.0,
            places=13,
        )

    def test_inner_domain_inflow_is_not_counted_as_black_hole_loss(self):
        state = FluidState(
            density=np.array([1.0, 1.0]),
            radial_velocity=np.array([0.01, 0.0]),
            velocity_dispersion=np.ones(2),
        )
        end_time = 1.0e-3

        result = evolve_sidm(
            state,
            self.grid,
            end_time,
            initial_black_hole_mass_code=0.5,
        )

        self.assertEqual(result.final_black_hole_mass_code, 0.5)
        self.assertGreater(result.history.cumulative_inner_inflow_code[-1], 0.0)
        self.assertAlmostEqual(
            result.history.max_absolute_mass_budget_residual_code,
            0.0,
            places=13,
        )

    def test_optional_conduction_changes_only_thermal_substep(self):
        grid = SphericalGrid.from_log_spacing(0.5, 4.0, 8)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.linspace(0.5, 1.5, grid.num_cells),
        )
        end_time = 1.0e-3

        adiabatic = evolve_sidm(state, grid, end_time)
        conducting = evolve_sidm(
            state,
            grid,
            end_time,
            sigma_over_m_code=2.0,
        )

        np.testing.assert_allclose(
            conducting.final_state.density,
            adiabatic.final_state.density,
        )
        np.testing.assert_allclose(
            conducting.final_state.radial_velocity,
            adiabatic.final_state.radial_velocity,
        )
        self.assertFalse(
            np.allclose(
                conducting.final_state.velocity_dispersion,
                adiabatic.final_state.velocity_dispersion,
            )
        )

    def test_max_steps_failure_is_explicit(self):
        state = FluidState(
            density=np.ones(2),
            radial_velocity=np.zeros(2),
            velocity_dispersion=np.ones(2),
        )

        with self.assertRaises(RuntimeError):
            evolve_sidm(state, self.grid, end_time_code=100.0, max_steps=1)


if __name__ == "__main__":
    unittest.main()
