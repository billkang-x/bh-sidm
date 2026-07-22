import math
import unittest

from sidm_bh.constants import G_CGS, M_SUN_CGS, PC_CGS
from sidm_bh.halos import NFWProfile, SingularIsothermalSphere
from sidm_bh.sidm import conductivity_cgs, conductivity_code, knudsen_number
from sidm_bh.units import SimulationScales


class SimulationScalesTest(unittest.TestCase):
    def test_paper_scale_definitions(self):
        scales = SimulationScales(radius_scale_pc=30.0, density_scale_msun_pc3=3.7)

        self.assertAlmostEqual(
            scales.mass_scale_msun,
            4.0 * math.pi * 30.0**3 * 3.7,
        )
        self.assertAlmostEqual(
            scales.sigma_over_m_scale_cgs,
            1.0 / (scales.density_scale_cgs * scales.radius_scale_cgs),
        )

    def test_cgs_and_code_conductivity_are_consistent(self):
        scales = SimulationScales(radius_scale_pc=30.0, density_scale_msun_pc3=3.7)
        density_code = 2.5
        velocity_code = 0.4
        sigma_code = 12.0

        kappa_cgs = conductivity_cgs(
            density_cgs=density_code * scales.density_scale_cgs,
            velocity_dispersion_cms=velocity_code * scales.velocity_scale_cgs,
            sigma_over_m_cm2_g=sigma_code * scales.sigma_over_m_scale_cgs,
        )
        self.assertAlmostEqual(
            kappa_cgs / scales.conductivity_scale_cgs,
            conductivity_code(density_code, velocity_code, sigma_code),
        )


class HaloProfileTest(unittest.TestCase):
    def test_nfw_mass_concentration_constructor(self):
        profile = NFWProfile.from_mass_concentration(
            total_mass_msun=1.0e6,
            virial_radius_pc=1000.0,
            concentration=10.0,
        )

        self.assertAlmostEqual(profile.scale_radius_pc, 100.0)
        self.assertAlmostEqual(profile.enclosed_mass_msun(1000.0), 1.0e6)

    def test_nfw_recovers_concentration_for_enclosed_mass(self):
        profile = NFWProfile.from_mass_concentration(
            total_mass_msun=1.0e6,
            virial_radius_pc=400.0,
            concentration=4.0,
        )

        self.assertAlmostEqual(
            profile.concentration_for_enclosed_mass(1.0e6),
            4.0,
            places=10,
        )

    def test_nfw_self_similar_scaling_preserves_density_and_concentration(self):
        anchor = NFWProfile(3.7, 30.0)
        anchor_mass = 1.0e6
        concentration = anchor.concentration_for_enclosed_mass(anchor_mass)
        scaled = anchor.self_similar_scaled(1.0e8, anchor_mass)

        self.assertAlmostEqual(scaled.scale_density_msun_pc3, 3.7)
        self.assertAlmostEqual(
            scaled.scale_radius_pc,
            30.0 * 100.0 ** (1.0 / 3.0),
        )
        self.assertAlmostEqual(
            scaled.enclosed_mass_msun(concentration * scaled.scale_radius_pc),
            1.0e8,
            places=7,
        )

    def test_sis_enclosed_mass_matches_density_integral(self):
        sis = SingularIsothermalSphere(sound_speed_km_s=4.2)
        radius_pc = 100.0
        density = sis.density_msun_pc3(radius_pc)
        local_integral_derivative = 4.0 * math.pi * radius_pc**2 * density
        analytic_derivative = 2.0 * sis.sound_speed_cgs**2 * PC_CGS / G_CGS / M_SUN_CGS

        self.assertAlmostEqual(local_integral_derivative, analytic_derivative)

    def test_knudsen_number_positive(self):
        scales = SimulationScales(radius_scale_pc=30.0, density_scale_msun_pc3=3.7)

        self.assertGreater(
            knudsen_number(
                density_cgs=scales.density_scale_cgs,
                velocity_dispersion_cms=scales.velocity_scale_cgs,
                sigma_over_m_cm2_g=50.0,
            ),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
