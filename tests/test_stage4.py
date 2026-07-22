import math
import unittest

from sidm_bh.stage4 import (
    EddingtonBaryonModel,
    bondi_accretion_rate_msun_per_myr,
    effective_hernquist_binding_energy_erg,
    feedback_expanded_scale_radius,
    feedback_heated_sound_speed,
    homologous_ambient_density,
    hernquist_self_binding_energy_erg,
)
from sidm_bh.halos import NFWProfile
from sidm_bh.units import SimulationScales


class EddingtonBaryonModelTest(unittest.TestCase):
    def test_black_hole_efolding_time_uses_retained_mass(self):
        model = EddingtonBaryonModel(radiative_efficiency=0.1)
        inflow_coefficient = model.eddington_inflow_coefficient_per_myr

        self.assertAlmostEqual(
            model.black_hole_efolding_time_myr,
            1.0 / (0.9 * inflow_coefficient),
        )
        self.assertGreater(model.black_hole_efolding_time_myr, 40.0)
        self.assertLess(model.black_hole_efolding_time_myr, 60.0)

    def test_zero_duty_cycle_disables_baryon_growth(self):
        model = EddingtonBaryonModel(duty_cycle=0.0)
        self.assertEqual(model.eddington_inflow_coefficient_per_myr, 0.0)
        self.assertTrue(math.isinf(model.black_hole_efolding_time_myr))

    def test_feedback_conversion_is_dimensionless_and_positive(self):
        scales = SimulationScales(30.0, 3.7)
        model = EddingtonBaryonModel(feedback_efficiency=1.0e-4)

        ratio = model.feedback_ratio_per_consumed_mass_code(
            scales,
            initial_baryon_mass_msun=5.0e4,
            initial_scale_radius_pc=0.3,
        )

        self.assertGreater(ratio, 0.0)
        self.assertGreater(
            hernquist_self_binding_energy_erg(5.0e4, 0.3),
            0.0,
        )

    def test_feedback_expands_but_never_contracts_profile(self):
        self.assertEqual(feedback_expanded_scale_radius(0.3, 0.0, 0.5), 0.3)
        self.assertAlmostEqual(
            feedback_expanded_scale_radius(0.3, 3.0, 0.5),
            0.6,
        )

    def test_invalid_efficiencies_are_rejected(self):
        with self.assertRaises(ValueError):
            EddingtonBaryonModel(radiative_efficiency=1.0)
        with self.assertRaises(ValueError):
            EddingtonBaryonModel(eddington_ratio=1.1)
        with self.assertRaises(ValueError):
            EddingtonBaryonModel(feedback_efficiency=-0.1)

    def test_bondi_rate_has_expected_mass_density_and_speed_scaling(self):
        base = bondi_accretion_rate_msun_per_myr(100.0, 1.0e3, 10.0)
        self.assertAlmostEqual(
            bondi_accretion_rate_msun_per_myr(200.0, 1.0e3, 10.0) / base,
            4.0,
        )
        self.assertAlmostEqual(
            bondi_accretion_rate_msun_per_myr(100.0, 2.0e3, 10.0) / base,
            2.0,
        )
        self.assertAlmostEqual(
            bondi_accretion_rate_msun_per_myr(100.0, 1.0e3, 20.0) / base,
            1.0 / 8.0,
        )

    def test_bondi_code_coefficient_reproduces_physical_rate(self):
        scales = SimulationScales(30.0, 3.7)
        model = EddingtonBaryonModel(
            gas_density_msun_pc3=1.0e3,
            gas_sound_speed_km_s=10.0,
        )
        mass_code = scales.mass_to_code(100.0)
        rate_code = model.bondi_inflow_coefficient_code(scales) * mass_code**2
        rate_physical = rate_code * scales.mass_scale_msun / scales.time_scale_myr

        self.assertAlmostEqual(
            rate_physical,
            bondi_accretion_rate_msun_per_myr(100.0, 1.0e3, 10.0),
        )

    def test_effective_binding_adds_halo_and_black_hole_terms(self):
        profile = NFWProfile(3.7, 30.0)
        components = effective_hernquist_binding_energy_erg(
            total_mass_msun=5.0e4,
            scale_radius_pc=0.3,
            halo=profile,
            black_hole_mass_msun=100.0,
        )

        self.assertGreater(components.halo_erg, 0.0)
        self.assertGreater(components.black_hole_erg, 0.0)
        self.assertGreater(components.total_erg, components.self_gravity_erg)

    def test_homologous_expansion_dilutes_density_as_mass_over_radius_cubed(self):
        self.assertAlmostEqual(
            homologous_ambient_density(100.0, 0.8, 2.0),
            10.0,
        )

    def test_feedback_energy_increases_sound_speed(self):
        heated = feedback_heated_sound_speed(10.0, 90.0)
        self.assertGreater(heated, 10.0)

    def test_feedback_energy_partition_is_conservative(self):
        scales = SimulationScales(30.0, 3.7)
        model = EddingtonBaryonModel(feedback_efficiency=1.0e-4)
        full_expansion = model.feedback_ratio_per_consumed_mass_code(
            scales,
            5.0e4,
            0.3,
        )
        half_expansion = model.feedback_ratio_per_consumed_mass_code(
            scales,
            5.0e4,
            0.3,
            expansion_energy_fraction=0.5,
        )
        self.assertAlmostEqual(half_expansion, 0.5 * full_expansion)


if __name__ == "__main__":
    unittest.main()
