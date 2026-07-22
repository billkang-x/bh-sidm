import unittest

import numpy as np

from sidm_bh.conduction import (
    cell_conductivity_code,
    conduction_system_code,
    implicit_conduction_step,
    interface_conductivity_code,
    solve_tridiagonal,
)
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sidm import conductivity_code
from sidm_bh.state import FluidState


class TridiagonalSolverTest(unittest.TestCase):
    def test_matches_dense_linear_solve(self):
        lower = np.array([-1.0, -2.0, -1.0])
        diagonal = np.array([4.0, 5.0, 6.0, 4.0])
        upper = np.array([-0.5, -1.0, -1.5])
        rhs = np.array([1.0, 2.0, 3.0, 4.0])
        matrix = np.diag(diagonal) + np.diag(lower, -1) + np.diag(upper, 1)

        expected = np.linalg.solve(matrix, rhs)

        np.testing.assert_allclose(
            solve_tridiagonal(lower, diagonal, upper, rhs),
            expected,
        )


class ImplicitConductionTest(unittest.TestCase):
    def setUp(self):
        self.grid = SphericalGrid.from_log_spacing(0.1, 10.0, 32)

    def test_constant_temperature_is_unchanged(self):
        state = FluidState(
            density=np.geomspace(10.0, 0.1, self.grid.num_cells),
            radial_velocity=np.linspace(-0.2, 0.2, self.grid.num_cells),
            velocity_dispersion=np.full(self.grid.num_cells, 1.7),
        )

        updated = implicit_conduction_step(
            state,
            self.grid,
            dt_code=3.0,
            sigma_over_m_code=2.0,
        )

        np.testing.assert_allclose(updated.density, state.density)
        np.testing.assert_allclose(updated.radial_velocity, state.radial_velocity)
        np.testing.assert_allclose(
            updated.velocity_dispersion,
            state.velocity_dispersion,
        )

    def test_vectorized_cell_conductivity_matches_scalar_closure(self):
        state = FluidState(
            density=np.array([0.5, 1.0, 2.0]),
            radial_velocity=np.zeros(3),
            velocity_dispersion=np.array([0.7, 1.1, 1.8]),
        )
        sigma = 2.3

        expected = np.array(
            [
                conductivity_code(rho, velocity, sigma)
                for rho, velocity in zip(
                    state.density,
                    state.velocity_dispersion,
                    strict=True,
                )
            ]
        )

        np.testing.assert_allclose(
            cell_conductivity_code(state, sigma),
            expected,
        )

    def test_cell_conductivity_accepts_local_cross_section(self):
        state = FluidState(
            density=np.array([0.8, 1.2, 2.0]),
            radial_velocity=np.zeros(3),
            velocity_dispersion=np.array([0.7, 1.1, 1.8]),
        )
        sigma = np.array([0.4, 1.0, 2.3])
        expected = np.array(
            [
                conductivity_code(rho, velocity, local_sigma)
                for rho, velocity, local_sigma in zip(
                    state.density,
                    state.velocity_dispersion,
                    sigma,
                    strict=True,
                )
            ]
        )
        np.testing.assert_allclose(cell_conductivity_code(state, sigma), expected)

    def test_no_flux_boundaries_conserve_thermal_energy(self):
        state = FluidState(
            density=np.geomspace(5.0, 0.2, self.grid.num_cells),
            radial_velocity=np.zeros(self.grid.num_cells),
            velocity_dispersion=np.linspace(0.5, 2.0, self.grid.num_cells),
        )
        initial_thermal_energy = np.sum(
            1.5
            * state.density
            * state.velocity_dispersion**2
            * self.grid.cell_volumes_code
        )

        updated = implicit_conduction_step(
            state,
            self.grid,
            dt_code=10.0,
            sigma_over_m_code=3.0,
        )

        final_thermal_energy = np.sum(
            1.5
            * updated.density
            * updated.velocity_dispersion**2
            * self.grid.cell_volumes_code
        )
        self.assertAlmostEqual(final_thermal_energy, initial_thermal_energy)

    def test_conduction_reduces_temperature_range(self):
        state = FluidState(
            density=np.ones(self.grid.num_cells),
            radial_velocity=np.zeros(self.grid.num_cells),
            velocity_dispersion=np.sqrt(
                np.where(
                    np.arange(self.grid.num_cells) == self.grid.num_cells // 2,
                    9.0,
                    1.0,
                )
            ),
        )
        initial_temperature_range = np.ptp(state.velocity_dispersion**2)

        updated = implicit_conduction_step(
            state,
            self.grid,
            dt_code=0.1,
            sigma_over_m_code=1.0,
        )

        self.assertLess(
            np.ptp(updated.velocity_dispersion**2),
            initial_temperature_range,
        )

    def test_assembled_system_matches_appendix_b_diagonal(self):
        grid = SphericalGrid(np.array([1.0, 2.0, 3.0, 4.0]))
        state = FluidState(
            density=np.array([1.0, 2.0, 3.0]),
            radial_velocity=np.zeros(3),
            velocity_dispersion=np.array([1.0, 2.0, 3.0]),
        )
        interface_kappa = interface_conductivity_code(np.array([2.0, 4.0, 6.0]))

        lower, diagonal, upper, _ = conduction_system_code(
            state,
            grid,
            dt_code=0.2,
            conductivity_at_interfaces_code=interface_kappa,
        )

        expected_diagonal = state.density.copy()
        expected_diagonal[1:] -= lower
        expected_diagonal[:-1] -= upper
        np.testing.assert_allclose(diagonal, expected_diagonal)


if __name__ == "__main__":
    unittest.main()
