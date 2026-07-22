import unittest

import numpy as np

from sidm_bh.fluxes import convective_flux_code, roe_flux_code
from sidm_bh.mesh import SphericalGrid
from sidm_bh.reconstruction import mc_reconstruct_primitive, mc_slopes, minmod
from sidm_bh.solver import advance_hyperbolic_step
from sidm_bh.state import FluidState


class MCLimiterTest(unittest.TestCase):
    def test_minmod_requires_common_sign(self):
        first = np.array([1.0, -3.0, 2.0, 0.0])
        second = np.array([2.0, -2.0, -1.0, 4.0])
        third = np.array([0.5, -4.0, 3.0, 5.0])

        np.testing.assert_allclose(
            minmod(first, second, third),
            [0.5, -2.0, 0.0, 0.0],
        )

    def test_nonuniform_mc_is_exact_for_linear_data(self):
        grid = SphericalGrid(np.array([1.0, 1.3, 2.0, 3.7, 6.0, 10.0]))
        radius = grid.centers_code
        density = 2.0 + 0.1 * radius
        velocity = -0.4 + 0.03 * radius
        pressure = 1.0 + 0.2 * radius
        state = FluidState.from_pressure(density, velocity, pressure)

        left, right = mc_reconstruct_primitive(state, grid)
        interfaces = grid.interfaces_code[2:-2]

        np.testing.assert_allclose(left.density[2:-2], 2.0 + 0.1 * interfaces)
        np.testing.assert_allclose(right.density[2:-2], 2.0 + 0.1 * interfaces)
        np.testing.assert_allclose(
            left.radial_velocity[2:-2],
            -0.4 + 0.03 * interfaces,
        )
        np.testing.assert_allclose(
            right.radial_velocity[2:-2],
            -0.4 + 0.03 * interfaces,
        )
        np.testing.assert_allclose(left.pressure[2:-2], 1.0 + 0.2 * interfaces)
        np.testing.assert_allclose(right.pressure[2:-2], 1.0 + 0.2 * interfaces)

    def test_mc_reconstruction_creates_no_new_extrema(self):
        grid = SphericalGrid.from_log_spacing(1.0, 10.0, 6)
        values = np.array([[1.0, 1.0, 2.0, 4.0, 4.0, 3.0]])
        slopes = mc_slopes(values, grid.centers_code)
        reconstructed_left = values[:, :-1] + slopes[:, :-1] * (
            grid.interfaces_code[1:-1] - grid.centers_code[:-1]
        )
        reconstructed_right = values[:, 1:] + slopes[:, 1:] * (
            grid.interfaces_code[1:-1] - grid.centers_code[1:]
        )
        local_min = np.minimum(values[:, :-1], values[:, 1:])
        local_max = np.maximum(values[:, :-1], values[:, 1:])

        self.assertTrue(np.all(reconstructed_left >= local_min))
        self.assertTrue(np.all(reconstructed_left <= local_max))
        self.assertTrue(np.all(reconstructed_right >= local_min))
        self.assertTrue(np.all(reconstructed_right <= local_max))


class RoeFluxTest(unittest.TestCase):
    def test_identical_states_reduce_to_physical_flux(self):
        state = FluidState(
            density=np.array([1.0, 2.0]),
            radial_velocity=np.array([-0.2, 0.3]),
            velocity_dispersion=np.array([0.7, 1.1]),
        )

        np.testing.assert_allclose(
            roe_flux_code(state, state),
            convective_flux_code(state),
        )

    def test_stationary_contact_is_preserved(self):
        pressure = 3.0
        left = FluidState.from_pressure(
            np.array([1.0]),
            np.array([0.0]),
            np.array([pressure]),
        )
        right = FluidState.from_pressure(
            np.array([4.0]),
            np.array([0.0]),
            np.array([pressure]),
        )

        np.testing.assert_allclose(
            roe_flux_code(left, right)[:, 0],
            [0.0, pressure, 0.0],
            atol=1.0e-14,
        )

    def test_supersonic_flow_uses_left_state_flux(self):
        left = FluidState(
            density=np.array([1.0]),
            radial_velocity=np.array([5.0]),
            velocity_dispersion=np.array([0.5]),
        )
        right = FluidState(
            density=np.array([0.7]),
            radial_velocity=np.array([4.5]),
            velocity_dispersion=np.array([0.6]),
        )

        np.testing.assert_allclose(
            roe_flux_code(left, right),
            convective_flux_code(left),
        )

    def test_mc_roe_step_conserves_mass_through_boundary_fluxes(self):
        grid = SphericalGrid.from_log_spacing(1.0, 4.0, 8)
        radius = grid.centers_code
        state = FluidState(
            density=1.0 + 0.1 * radius,
            radial_velocity=-0.02 + 0.01 * radius,
            velocity_dispersion=0.8 + 0.02 * radius,
        )
        dt = 1.0e-4
        initial_mass = np.sum(state.density * grid.cell_volumes_code)

        updated = advance_hyperbolic_step(
            state,
            grid,
            dt,
            dark_matter_enclosed_mass_code=np.zeros(grid.num_cells),
            reconstruction="mc",
            riemann_solver="roe",
            positivity_fallback=False,
        )

        final_mass = np.sum(updated.density * grid.cell_volumes_code)
        inner_flux = state.density[0] * state.radial_velocity[0]
        outer_flux = state.density[-1] * state.radial_velocity[-1]
        expected_change = -dt * (
            grid.interface_areas_code[-1] * outer_flux
            - grid.interface_areas_code[0] * inner_flux
        )
        self.assertAlmostEqual(final_mass - initial_mass, expected_change)


if __name__ == "__main__":
    unittest.main()
