import unittest

import numpy as np

from sidm_bh.accretion import (
    AccretionStep,
    accretion_rate_code_to_msun_per_myr,
    accretion_rate_msun_per_myr_to_code,
    inner_boundary_sidm_accretion_rate_code,
    sidm_accretion_rate_code,
)
from sidm_bh.baryons import HernquistBaryons
from sidm_bh.fluxes import (
    ADIABATIC_INDEX,
    convective_flux_code,
    max_signal_speed_code,
    sound_speed_code,
)
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sources import (
    enclosed_baryon_mass_code,
    gravitational_acceleration_code,
    spherical_source_terms_code,
    total_enclosed_mass_code,
)
from sidm_bh.state import FluidState
from sidm_bh.units import SimulationScales


class SourceTermsTest(unittest.TestCase):
    def test_total_mass_includes_black_hole_and_baryons(self):
        dark = np.array([1.0, 2.0, 3.0])
        baryons = np.array([0.1, 0.2, 0.3])

        total = total_enclosed_mass_code(dark, black_hole_mass_code=5.0, baryon_enclosed_mass_code=baryons)

        np.testing.assert_allclose(total, [6.1, 7.2, 8.3])

    def test_baryon_mass_sampling_is_monotonic(self):
        grid = SphericalGrid.from_log_spacing(1.0e-3, 10.0, 32)
        scales = SimulationScales(radius_scale_pc=10.0, density_scale_msun_pc3=1.0)
        baryons = HernquistBaryons(total_mass_msun=1.0e5, scale_radius_pc=1.0)

        enclosed = enclosed_baryon_mass_code(baryons, grid, scales)

        self.assertTrue(np.all(np.diff(enclosed) > 0.0))

    def test_gravity_source_matches_paper_form(self):
        grid = SphericalGrid.from_log_spacing(1.0, 8.0, 3)
        state = FluidState(
            density=np.array([2.0, 2.0, 2.0]),
            radial_velocity=np.array([-1.0, 0.0, 1.0]),
            velocity_dispersion=np.array([0.5, 0.5, 0.5]),
        )
        dark_mass = np.array([1.0, 2.0, 3.0])

        source = spherical_source_terms_code(state, grid, dark_mass, black_hole_mass_code=1.0)
        acceleration = gravitational_acceleration_code(grid.centers_code, dark_mass + 1.0)

        np.testing.assert_allclose(source[0], np.zeros(3))
        np.testing.assert_allclose(source[1], state.density * acceleration + 2.0 * state.pressure / grid.centers_code)
        np.testing.assert_allclose(source[2], state.density * state.radial_velocity * acceleration)


class AccretionTest(unittest.TestCase):
    def test_inward_velocity_grows_black_hole(self):
        self.assertAlmostEqual(
            sidm_accretion_rate_code(
                inner_radius_code=0.1,
                density_code=2.0,
                radial_velocity_code=-3.0,
            ),
            0.06,
        )

    def test_outward_velocity_is_zero_when_inward_only(self):
        self.assertEqual(
            sidm_accretion_rate_code(
                inner_radius_code=0.1,
                density_code=2.0,
                radial_velocity_code=3.0,
            ),
            0.0,
        )

    def test_inner_boundary_rate_uses_inner_interface(self):
        grid = SphericalGrid.from_log_spacing(0.01, 1.0, 4)
        state = FluidState(
            density=np.ones(4) * 2.0,
            radial_velocity=np.ones(4) * -5.0,
            velocity_dispersion=np.ones(4),
        )

        self.assertAlmostEqual(
            inner_boundary_sidm_accretion_rate_code(state, grid),
            0.01**2 * 2.0 * 5.0,
        )

    def test_accretion_rate_unit_roundtrip(self):
        scales = SimulationScales(radius_scale_pc=30.0, density_scale_msun_pc3=3.7)
        rate = 1.23e-4

        converted = accretion_rate_code_to_msun_per_myr(rate, scales)

        self.assertAlmostEqual(accretion_rate_msun_per_myr_to_code(converted, scales), rate)

    def test_accretion_step_updates_mass(self):
        step = AccretionStep(initial_mass_code=2.0, rate_code=0.5, dt_code=4.0)

        self.assertEqual(step.delta_mass_code, 2.0)
        self.assertEqual(step.final_mass_code, 4.0)


class FluxTest(unittest.TestCase):
    def test_convective_flux_matches_euler_form(self):
        state = FluidState(
            density=np.array([2.0]),
            radial_velocity=np.array([3.0]),
            velocity_dispersion=np.array([4.0]),
        )

        flux = convective_flux_code(state)
        pressure = 2.0 * 4.0**2
        energy_density = 2.0 * (1.5 * 4.0**2 + 0.5 * 3.0**2)

        np.testing.assert_allclose(
            flux[:, 0],
            [
                2.0 * 3.0,
                2.0 * 3.0**2 + pressure,
                3.0 * (energy_density + pressure),
            ],
        )

    def test_sound_speed_uses_gamma_five_thirds(self):
        state = FluidState(
            density=np.ones(2),
            radial_velocity=np.zeros(2),
            velocity_dispersion=np.array([1.0, 2.0]),
        )

        np.testing.assert_allclose(
            sound_speed_code(state),
            np.sqrt(ADIABATIC_INDEX) * state.velocity_dispersion,
        )
        self.assertAlmostEqual(max_signal_speed_code(state), np.sqrt(ADIABATIC_INDEX) * 2.0)

    def test_conservative_roundtrip(self):
        state = FluidState(
            density=np.array([2.0, 3.0]),
            radial_velocity=np.array([-1.0, 4.0]),
            velocity_dispersion=np.array([5.0, 6.0]),
        )

        recovered = FluidState.from_conservative(state.conservative)

        np.testing.assert_allclose(recovered.density, state.density)
        np.testing.assert_allclose(recovered.radial_velocity, state.radial_velocity)
        np.testing.assert_allclose(recovered.velocity_dispersion, state.velocity_dispersion)


if __name__ == "__main__":
    unittest.main()
