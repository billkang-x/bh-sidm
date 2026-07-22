import unittest

import numpy as np

from sidm_bh.accretion import inner_boundary_sidm_accretion_rate_code
from sidm_bh.fluxes import convective_flux_code, rusanov_flux_code, sound_speed_code
from sidm_bh.mesh import SphericalGrid
from sidm_bh.solver import (
    apply_gravity_kick,
    advance_hyperbolic_step,
    cfl_timestep_code,
    enclosed_dark_matter_mass_code,
    gravity_timestep_code,
    stable_timestep_code,
    zero_gradient_interface_states,
)
from sidm_bh.state import FluidState


class RusanovFluxTest(unittest.TestCase):
    def test_identical_states_reduce_to_physical_flux(self):
        state = FluidState(
            density=np.array([1.0, 2.0]),
            radial_velocity=np.array([-0.3, 0.4]),
            velocity_dispersion=np.array([0.7, 0.8]),
        )

        np.testing.assert_allclose(
            rusanov_flux_code(state, state),
            convective_flux_code(state),
        )

    def test_zero_gradient_states_align_with_interfaces(self):
        state = FluidState(
            density=np.array([1.0, 2.0, 3.0]),
            radial_velocity=np.array([-1.0, 0.0, 1.0]),
            velocity_dispersion=np.ones(3),
        )

        left, right = zero_gradient_interface_states(state)

        np.testing.assert_allclose(left.density, [1.0, 1.0, 2.0, 3.0])
        np.testing.assert_allclose(right.density, [1.0, 2.0, 3.0, 3.0])


class HyperbolicSolverTest(unittest.TestCase):
    def test_exact_gravity_kick_preserves_density_and_velocity_dispersion(self):
        state = FluidState(
            density=np.array([1.0, 2.0]),
            radial_velocity=np.array([-0.3, 0.4]),
            velocity_dispersion=np.array([0.7, 0.8]),
        )
        acceleration = np.array([-2.0, 3.0])

        kicked = apply_gravity_kick(state, acceleration, dt_code=0.25)

        np.testing.assert_allclose(kicked.density, state.density)
        np.testing.assert_allclose(
            kicked.radial_velocity,
            state.radial_velocity + 0.25 * acceleration,
        )
        np.testing.assert_allclose(
            kicked.velocity_dispersion,
            state.velocity_dispersion,
        )

    def test_cfl_timestep_uses_local_width_and_signal_speed(self):
        grid = SphericalGrid(np.array([1.0, 2.0, 4.0]))
        state = FluidState(
            density=np.ones(2),
            radial_velocity=np.array([0.2, -0.1]),
            velocity_dispersion=np.array([0.5, 2.0]),
        )
        signal_speed = np.abs(state.radial_velocity) + sound_speed_code(state)

        expected = 0.3 * np.min(grid.widths_code / signal_speed)

        self.assertAlmostEqual(cfl_timestep_code(state, grid, 0.3), expected)

    def test_uniform_density_enclosed_mass_matches_analytic_integral(self):
        grid = SphericalGrid.from_log_spacing(0.2, 10.0, 32)
        state = FluidState(
            density=np.full(grid.num_cells, 2.5),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )

        enclosed = enclosed_dark_matter_mass_code(
            state,
            grid,
            inner_enclosed_mass_code=0.7,
        )
        expected = 0.7 + 2.5 * (
            grid.centers_code**3 - grid.interfaces_code[0] ** 3
        ) / 3.0

        np.testing.assert_allclose(enclosed, expected)

    def test_stable_timestep_limits_strong_point_mass_acceleration(self):
        grid = SphericalGrid(np.array([0.1, 0.2, 0.4]))
        state = FluidState(
            density=np.ones(2),
            radial_velocity=np.zeros(2),
            velocity_dispersion=np.full(2, 0.1),
        )
        total_mass = np.full(2, 10.0)

        gravity_dt = gravity_timestep_code(state, grid, total_mass, 0.2)
        stable_dt = stable_timestep_code(
            state,
            grid,
            0.2,
            black_hole_mass_code=10.0,
            dark_matter_enclosed_mass_code=np.zeros(2),
        )

        self.assertLess(gravity_dt, cfl_timestep_code(state, grid, 0.2))
        self.assertAlmostEqual(stable_dt, gravity_dt)
        updated = advance_hyperbolic_step(
            state,
            grid,
            stable_dt,
            black_hole_mass_code=10.0,
            dark_matter_enclosed_mass_code=np.zeros(2),
        )
        self.assertTrue(np.all(updated.velocity_dispersion > 0.0))

    def test_gravity_limit_tightens_when_gas_cools_during_inflow(self):
        grid = SphericalGrid(np.array([0.1, 0.2]))
        hot = FluidState(
            density=np.ones(1),
            radial_velocity=np.array([-1.0]),
            velocity_dispersion=np.array([0.2]),
        )
        cold = FluidState(
            density=np.ones(1),
            radial_velocity=np.array([-1.0]),
            velocity_dispersion=np.array([0.02]),
        )
        total_mass = np.array([5.0])

        hot_dt = gravity_timestep_code(hot, grid, total_mass, 0.2)
        cold_dt = gravity_timestep_code(cold, grid, total_mass, 0.2)

        self.assertAlmostEqual(cold_dt / hot_dt, 0.1)

    def test_finite_volume_mass_change_equals_boundary_flux(self):
        grid = SphericalGrid(np.array([1.0, 1.5, 2.2, 3.0]))
        state = FluidState(
            density=np.array([1.0, 1.2, 1.4]),
            radial_velocity=np.array([-0.02, 0.01, 0.03]),
            velocity_dispersion=np.array([0.8, 0.7, 0.6]),
        )
        dt = 1.0e-3
        initial_mass = np.sum(state.density * grid.cell_volumes_code)

        updated = advance_hyperbolic_step(
            state,
            grid,
            dt,
            dark_matter_enclosed_mass_code=np.zeros(grid.num_cells),
        )

        final_mass = np.sum(updated.density * grid.cell_volumes_code)
        inner_flux = state.density[0] * state.radial_velocity[0]
        outer_flux = state.density[-1] * state.radial_velocity[-1]
        expected_change = -dt * (
            grid.interface_areas_code[-1] * outer_flux
            - grid.interface_areas_code[0] * inner_flux
        )
        self.assertAlmostEqual(final_mass - initial_mass, expected_change)

    def test_inner_accretion_balances_fluid_loss_when_outer_flux_is_zero(self):
        grid = SphericalGrid(np.array([1.0, 1.5, 2.0]))
        state = FluidState(
            density=np.array([2.0, 1.0]),
            radial_velocity=np.array([-0.01, 0.0]),
            velocity_dispersion=np.array([1.0, 1.0]),
        )
        dt = 1.0e-4
        initial_mass = np.sum(state.density * grid.cell_volumes_code)
        accreted_mass = inner_boundary_sidm_accretion_rate_code(state, grid) * dt

        updated = advance_hyperbolic_step(
            state,
            grid,
            dt,
            dark_matter_enclosed_mass_code=np.zeros(grid.num_cells),
        )

        final_mass = np.sum(updated.density * grid.cell_volumes_code)
        self.assertAlmostEqual(initial_mass - final_mass, accreted_mass)


if __name__ == "__main__":
    unittest.main()
