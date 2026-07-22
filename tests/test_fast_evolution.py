import importlib.util
import unittest

import numpy as np

from sidm_bh.evolution import evolve_sidm
from sidm_bh.halos import NFWProfile
from sidm_bh.baryons import HernquistBaryons
from sidm_bh.initial_conditions import hydrostatic_state_from_profile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.state import FluidState
from sidm_bh.units import SimulationScales
from sidm_bh.stage3 import static_baryon_equilibrium_state


@unittest.skipUnless(importlib.util.find_spec("numba"), "numba is not installed")
class FastEvolutionTest(unittest.TestCase):
    def test_fast_mc_roe_matches_reference_implementation(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        scales = SimulationScales(30.0, 3.7)
        profile = NFWProfile(3.7, 30.0)
        grid = SphericalGrid.from_log_spacing(0.005 / 30.0, 10.0, 16)
        state = hydrostatic_state_from_profile(profile, grid, scales)
        end_time = scales.time_to_code(1.0e-4)
        black_hole_mass = scales.mass_to_code(100.0)
        sigma = scales.sigma_over_m_to_code(50.0)
        common = dict(
            initial_black_hole_mass_code=black_hole_mass,
            sigma_over_m_code=sigma,
            cfl_number=0.15,
            entropy_fix=0.1,
            source_integration="euler",
            max_steps=10000,
        )

        reference = evolve_sidm(
            state,
            grid,
            end_time,
            reconstruction="mc",
            riemann_solver="roe",
            positivity_fallback=False,
            **common,
        )
        accelerated = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, end_time]),
            **common,
        )

        self.assertEqual(accelerated.num_steps, reference.history.num_steps)
        self.assertAlmostEqual(
            accelerated.final_black_hole_mass_code,
            reference.final_black_hole_mass_code,
        )
        np.testing.assert_allclose(
            accelerated.final_state.conservative,
            reference.final_state.conservative,
            rtol=2.0e-13,
            atol=1.0e-14,
        )

    def test_fast_static_baryon_run_matches_reference(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        scales = SimulationScales(30.0, 3.7)
        profile = NFWProfile(3.7, 30.0)
        grid = SphericalGrid.from_log_spacing(0.005 / 30.0, 10.0, 16)
        baryons = HernquistBaryons(1.0e4, 0.3)
        state, baryon_mass = static_baryon_equilibrium_state(
            profile,
            grid,
            scales,
            100.0,
            baryons,
        )
        end_time = scales.time_to_code(1.0e-4)
        common = dict(
            initial_black_hole_mass_code=scales.mass_to_code(100.0),
            sigma_over_m_code=scales.sigma_over_m_to_code(50.0),
            cfl_number=0.15,
            entropy_fix=0.1,
            baryon_enclosed_mass_code=baryon_mass,
            source_integration="euler",
            max_steps=10000,
        )

        reference = evolve_sidm(
            state,
            grid,
            end_time,
            reconstruction="mc",
            riemann_solver="roe",
            positivity_fallback=False,
            **common,
        )
        accelerated = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, end_time]),
            **common,
        )

        self.assertAlmostEqual(
            accelerated.final_black_hole_mass_code,
            reference.final_black_hole_mass_code,
        )
        np.testing.assert_allclose(
            accelerated.final_state.conservative,
            reference.final_state.conservative,
            rtol=3.0e-13,
            atol=2.0e-14,
        )

    def test_fast_assembling_baryon_run_matches_reference(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        scales = SimulationScales(30.0, 3.7)
        profile = NFWProfile(3.7, 30.0)
        grid = SphericalGrid.from_log_spacing(0.005 / 30.0, 10.0, 16)
        state, _ = static_baryon_equilibrium_state(
            profile, grid, scales, 100.0, baryons=None
        )
        _, full_baryon_mass = static_baryon_equilibrium_state(
            profile, grid, scales, 100.0, HernquistBaryons(1.0e4, 0.3)
        )
        end_time = scales.time_to_code(2.0e-4)
        common = dict(
            initial_black_hole_mass_code=scales.mass_to_code(100.0),
            sigma_over_m_code=scales.sigma_over_m_to_code(50.0),
            cfl_number=0.15,
            entropy_fix=0.1,
            baryon_enclosed_mass_code=full_baryon_mass,
            baryon_assembly_time_code=scales.time_to_code(1.0e-4),
            source_integration="euler",
            max_steps=10000,
        )
        reference = evolve_sidm(
            state,
            grid,
            end_time,
            reconstruction="mc",
            riemann_solver="roe",
            positivity_fallback=False,
            **common,
        )
        accelerated = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, end_time]),
            **common,
        )

        self.assertEqual(accelerated.num_steps, reference.history.num_steps)
        self.assertAlmostEqual(
            accelerated.final_black_hole_mass_code,
            reference.final_black_hole_mass_code,
        )
        np.testing.assert_allclose(
            accelerated.final_state.conservative,
            reference.final_state.conservative,
            rtol=3.0e-13,
            atol=2.0e-14,
        )

    def test_evolving_baryon_reservoir_separates_growth_channels(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        initial_black_hole_mass = 0.2
        result = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=initial_black_hole_mass,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=0.05,
            baryon_radiative_efficiency=0.1,
        )

        self.assertAlmostEqual(
            result.final_black_hole_mass_code,
            initial_black_hole_mass
            + result.dark_matter_accreted_masses_code[-1]
            + result.baryon_accreted_masses_code[-1],
        )
        self.assertAlmostEqual(
            result.baryon_accreted_masses_code[-1],
            0.9 * result.baryon_gas_consumed_masses_code[-1],
        )
        self.assertAlmostEqual(
            result.baryon_remaining_masses_code[-1]
            + result.baryon_gas_consumed_masses_code[-1],
            1.0,
        )
        self.assertLess(result.max_mass_budget_residual_code, 1.0e-12)

    def test_dark_bondi_reservoir_separates_supply_from_capture(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=-0.1 * np.ones(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        initial_black_hole_mass = 0.2
        result = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=initial_black_hole_mass,
            sigma_over_m_code=0.0,
            dark_capture_model="bondi_reservoir",
            dark_bondi_lambda=0.0,
        )

        self.assertGreater(result.dark_matter_supplied_masses_code[-1], 0.0)
        self.assertEqual(result.dark_matter_accreted_masses_code[-1], 0.0)
        self.assertAlmostEqual(
            result.inner_dark_matter_reservoir_masses_code[-1],
            result.dark_matter_supplied_masses_code[-1],
        )
        self.assertAlmostEqual(
            result.final_black_hole_mass_code,
            initial_black_hole_mass,
        )
        self.assertLess(result.max_mass_budget_residual_code, 1.0e-12)

    def test_dark_bondi_capture_is_limited_by_available_reservoir(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=-0.1 * np.ones(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        initial_black_hole_mass = 0.2
        initial_reservoir = 0.01
        result = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=initial_black_hole_mass,
            sigma_over_m_code=0.0,
            dark_capture_model="bondi_reservoir",
            dark_bondi_lambda=1.0e9,
            initial_dark_matter_reservoir_code=initial_reservoir,
        )

        captured = result.dark_matter_accreted_masses_code[-1]
        supplied = result.dark_matter_supplied_masses_code[-1]
        reservoir = result.inner_dark_matter_reservoir_masses_code[-1]
        self.assertGreater(captured, 0.0)
        self.assertGreaterEqual(reservoir, 0.0)
        self.assertLessEqual(captured, initial_reservoir + supplied + 1.0e-14)
        self.assertAlmostEqual(captured + reservoir, initial_reservoir + supplied)
        self.assertAlmostEqual(
            result.final_black_hole_mass_code,
            initial_black_hole_mass + captured,
        )
        self.assertLess(result.max_mass_budget_residual_code, 1.0e-12)

    def test_influence_gated_capture_matches_resolved_boundary_flux(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=-0.1 * np.ones(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
        )
        direct = evolve_mc_roe_fast(**common)
        gated = evolve_mc_roe_fast(
            **common,
            dark_capture_model="influence_gated",
            dark_flux_capture_mass_threshold_code=0.1,
        )

        np.testing.assert_allclose(
            gated.black_hole_masses_code,
            direct.black_hole_masses_code,
            rtol=2.0e-13,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            gated.dark_matter_accreted_masses_code,
            direct.dark_matter_accreted_masses_code,
            rtol=2.0e-13,
            atol=1.0e-14,
        )
        self.assertLess(gated.max_mass_budget_residual_code, 1.0e-12)

    def test_feedback_expands_evolving_baryon_profile(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        result = evolve_mc_roe_fast(
            state,
            grid,
            np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=0.05,
            feedback_ratio_per_consumed_mass_code=100.0,
            feedback_expansion_exponent=0.5,
        )

        self.assertGreater(result.feedback_to_binding_ratios[-1], 0.0)
        self.assertGreater(result.baryon_scale_radii_code[-1], 0.1)

    def test_bondi_limit_caps_eddington_baryon_inflow(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=1.0,
        )
        eddington_only = evolve_mc_roe_fast(**common)
        bondi_limited = evolve_mc_roe_fast(
            **common,
            baryon_bondi_inflow_coefficient_code=1.0e-3,
        )

        self.assertLess(
            bondi_limited.baryon_accreted_masses_code[-1],
            eddington_only.baryon_accreted_masses_code[-1],
        )
        self.assertLess(
            bondi_limited.baryon_accretion_rates_code[-1],
            eddington_only.baryon_accretion_rates_code[-1],
        )

    def test_expansion_and_heating_reduce_evolving_bondi_inflow(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-2]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=1.0,
            baryon_bondi_inflow_coefficient_code=0.1,
            feedback_ratio_per_consumed_mass_code=1.0e4,
            feedback_expansion_exponent=0.5,
        )
        fixed_ambient = evolve_mc_roe_fast(**common)
        evolving_ambient = evolve_mc_roe_fast(
            **common,
            evolve_bondi_ambient=True,
            baryon_initial_sound_speed_code=0.1,
            feedback_thermal_energy_per_consumed_mass_code=1.0e5,
        )

        self.assertLess(
            evolving_ambient.baryon_accreted_masses_code[-1],
            fixed_ambient.baryon_accreted_masses_code[-1],
        )
        self.assertLess(
            evolving_ambient.baryon_accretion_rates_code[-1],
            fixed_ambient.baryon_accretion_rates_code[-1],
        )

    def test_finite_cooling_restores_bondi_inflow(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-2]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=1.0,
            baryon_bondi_inflow_coefficient_code=0.1,
            evolve_bondi_ambient=True,
            baryon_initial_sound_speed_code=0.1,
            feedback_thermal_energy_per_consumed_mass_code=1.0e5,
        )
        no_cooling = evolve_mc_roe_fast(**common)
        infinite_cooling = evolve_mc_roe_fast(
            **common,
            baryon_cooling_time_code=np.inf,
        )
        fast_cooling = evolve_mc_roe_fast(
            **common,
            baryon_cooling_time_code=1.0e-5,
        )

        self.assertAlmostEqual(
            no_cooling.baryon_thermal_energies_code[-1],
            1.0e5 * no_cooling.baryon_gas_consumed_masses_code[-1],
        )
        np.testing.assert_allclose(
            infinite_cooling.baryon_thermal_energies_code,
            no_cooling.baryon_thermal_energies_code,
        )
        np.testing.assert_allclose(
            infinite_cooling.baryon_accreted_masses_code,
            no_cooling.baryon_accreted_masses_code,
        )
        self.assertLess(
            fast_cooling.baryon_thermal_energies_code[-1],
            no_cooling.baryon_thermal_energies_code[-1],
        )
        self.assertGreater(
            fast_cooling.baryon_accreted_masses_code[-1],
            no_cooling.baryon_accreted_masses_code[-1],
        )

    def test_cloudy_cooling_restores_bondi_inflow(self):
        from sidm_bh.constants import M_SUN_CGS, MYR_CGS, PC_CGS
        from sidm_bh.cooling import load_cloudy_cooling_table
        from sidm_bh.fast_evolution import evolve_mc_roe_fast

        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.ones(grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.ones(grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-2]),
            initial_black_hole_mass_code=0.2,
            sigma_over_m_code=0.0,
            evolving_baryon_total_mass_code=1.0,
            evolving_baryon_scale_radius_code=0.1,
            baryon_eddington_inflow_coefficient_code=1.0,
            baryon_bondi_inflow_coefficient_code=0.1,
            evolve_bondi_ambient=True,
            baryon_initial_sound_speed_code=0.1,
            feedback_thermal_energy_per_consumed_mass_code=1.0e5,
        )
        no_cooling = evolve_mc_roe_fast(**common)
        cloudy = evolve_mc_roe_fast(
            **common,
            baryon_cloudy_cooling_table=load_cloudy_cooling_table(),
            baryon_initial_gas_density_cgs=(
                300.0 * M_SUN_CGS / PC_CGS**3
            ),
            baryon_velocity_unit_cgs=1.0e7,
            baryon_time_unit_s=MYR_CGS,
        )

        self.assertLess(
            cloudy.baryon_thermal_energies_code[-1],
            no_cooling.baryon_thermal_energies_code[-1],
        )
        self.assertGreater(
            cloudy.baryon_accreted_masses_code[-1],
            no_cooling.baryon_accreted_masses_code[-1],
        )

    def test_evolving_hernquist_without_accretion_matches_fixed_profile(self):
        from sidm_bh.fast_evolution import evolve_mc_roe_fast
        from sidm_bh.sources import enclosed_baryon_mass_code

        scales = SimulationScales(30.0, 3.7)
        profile = NFWProfile(3.7, 30.0)
        grid = SphericalGrid.from_log_spacing(0.005 / 30.0, 10.0, 16)
        state, _ = static_baryon_equilibrium_state(
            profile, grid, scales, 100.0, baryons=None
        )
        baryons = HernquistBaryons(5.0e4, 0.3)
        fixed_mass = enclosed_baryon_mass_code(baryons, grid, scales)
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, scales.time_to_code(1.0e-4)]),
            initial_black_hole_mass_code=scales.mass_to_code(100.0),
            sigma_over_m_code=scales.sigma_over_m_to_code(50.0),
        )

        fixed = evolve_mc_roe_fast(
            **common,
            baryon_enclosed_mass_code=fixed_mass,
        )
        evolving = evolve_mc_roe_fast(
            **common,
            evolving_baryon_total_mass_code=scales.mass_to_code(5.0e4),
            evolving_baryon_scale_radius_code=scales.radius_to_code(0.3),
        )

        np.testing.assert_allclose(
            evolving.final_state.conservative,
            fixed.final_state.conservative,
            rtol=3.0e-13,
            atol=2.0e-14,
        )
        self.assertAlmostEqual(
            evolving.final_black_hole_mass_code,
            fixed.final_black_hole_mass_code,
        )


if __name__ == "__main__":
    unittest.main()
