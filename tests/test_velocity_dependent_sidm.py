"""Tests for the velocity-dependent SIDM transport closure."""

from __future__ import annotations

import unittest

import numpy as np

from sidm_bh.fast_evolution import evolve_mc_roe_fast
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sidm import (
    effective_cross_section_ratio_tables,
    maxwellian_viscosity_cross_section_ratio,
    rutherford_momentum_transfer_ratio,
    rutherford_viscosity_cross_section_ratio,
)
from sidm_bh.state import FluidState


class VelocityDependentCrossSectionTests(unittest.TestCase):
    def test_rutherford_low_velocity_limit_and_monotonicity(self):
        ratios = rutherford_momentum_transfer_ratio(
            np.array([0.0, 1.0e-4, 0.1, 1.0, 10.0])
        )
        self.assertAlmostEqual(ratios[0], 1.0)
        self.assertTrue(np.all(np.diff(ratios) < 0.0))

    def test_viscosity_cross_section_limit_and_monotonicity(self):
        ratios = rutherford_viscosity_cross_section_ratio(
            np.array([0.0, 1.0e-4, 0.1, 1.0, 10.0])
        )
        self.assertAlmostEqual(ratios[0], 1.0)
        self.assertTrue(np.all(np.diff(ratios) < 0.0))

    def test_maxwellian_averages_low_velocity_limit_and_monotonicity(self):
        dispersion_ratios = np.logspace(-5, 2, 40)
        for velocity_power in (3, 5):
            effective = maxwellian_viscosity_cross_section_ratio(
                dispersion_ratios,
                velocity_power=velocity_power,
            )
            self.assertAlmostEqual(effective[0], 1.0, places=7)
            self.assertTrue(np.all(np.diff(effective) < 0.0))
            self.assertLess(effective[-1], 1.0e-5)

    def test_lookup_table_resolves_direct_quadrature(self):
        log_ratio, k3_table, k5_table = effective_cross_section_ratio_tables()
        probes = np.logspace(-5.5, 3.5, 100)
        for velocity_power, table in ((3, k3_table), (5, k5_table)):
            interpolated = np.interp(np.log10(probes), log_ratio, table)
            direct = maxwellian_viscosity_cross_section_ratio(
                probes,
                velocity_power=velocity_power,
            )
            np.testing.assert_allclose(interpolated, direct, rtol=2.0e-4)

    def test_large_transition_speed_recovers_constant_solver(self):
        grid = SphericalGrid.from_log_spacing(0.05, 5.0, 16)
        state = FluidState(
            density=np.linspace(1.0, 0.2, grid.num_cells),
            radial_velocity=np.zeros(grid.num_cells),
            velocity_dispersion=np.linspace(0.4, 0.8, grid.num_cells),
        )
        common = dict(
            initial_state=state,
            grid=grid,
            sample_times_code=np.array([0.0, 1.0e-3]),
            initial_black_hole_mass_code=0.1,
            sigma_over_m_code=1.0,
        )
        constant = evolve_mc_roe_fast(**common)
        velocity_dependent = evolve_mc_roe_fast(
            **common,
            cross_section_model="rutherford",
            cross_section_velocity_scale_code=1.0e8,
        )
        np.testing.assert_allclose(
            velocity_dependent.final_state.conservative,
            constant.final_state.conservative,
            rtol=2.0e-10,
            atol=1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()
