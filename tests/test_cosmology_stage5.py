import unittest

from sidm_bh.cosmology import FlatLambdaCDM
from sidm_bh.halos import NFWProfile
from sidm_bh.stage4 import EddingtonBaryonModel
from sidm_bh.stage5 import (
    eddington_efolding_time_myr,
    growth_time_to_target_myr,
    required_eddington_activity,
)


class FlatLambdaCDMTest(unittest.TestCase):
    def setUp(self):
        self.cosmology = FlatLambdaCDM()

    def test_high_redshift_ages_match_planck_budget(self):
        self.assertAlmostEqual(self.cosmology.age_myr(10.0), 472.2, delta=1.0)
        self.assertAlmostEqual(self.cosmology.age_myr(4.0), 1540.0, delta=5.0)

    def test_age_redshift_roundtrip(self):
        for redshift in (0.0, 4.0, 10.0, 20.0, 30.0):
            recovered = self.cosmology.redshift_at_age_myr(
                self.cosmology.age_myr(redshift)
            )
            self.assertAlmostEqual(recovered, redshift, places=10)

    def test_elapsed_time_requires_forward_cosmic_evolution(self):
        self.assertGreater(self.cosmology.elapsed_time_myr(20.0, 6.0), 0.0)
        with self.assertRaises(ValueError):
            self.cosmology.elapsed_time_myr(6.0, 20.0)

    def test_spherical_overdensity_radius_has_mass_cube_root_scaling(self):
        small = self.cosmology.spherical_overdensity_radius_pc(1.0e8, 10.0)
        large = self.cosmology.spherical_overdensity_radius_pc(1.0e9, 10.0)
        self.assertAlmostEqual(large / small, 10.0 ** (1.0 / 3.0), places=12)

    def test_cosmological_nfw_recovers_m200(self):
        mass = 1.0e9
        radius = self.cosmology.spherical_overdensity_radius_pc(mass, 10.0)
        profile = NFWProfile.from_mass_concentration(mass, radius, 4.0)
        self.assertAlmostEqual(
            profile.enclosed_mass_msun(radius) / mass,
            1.0,
            places=12,
        )

    def test_virial_influence_radius_matches_mass_fraction(self):
        halo_mass = 1.0e9
        black_hole_mass = 1.0e5
        virial_radius = self.cosmology.spherical_overdensity_radius_pc(
            halo_mass,
            10.0,
        )
        influence_radius = self.cosmology.black_hole_influence_radius_pc(
            black_hole_mass,
            halo_mass,
            10.0,
        )
        self.assertAlmostEqual(
            influence_radius / virial_radius,
            black_hole_mass / halo_mass,
            places=14,
        )


class Stage5GrowthBudgetTest(unittest.TestCase):
    def test_efolding_time_matches_stage4_model(self):
        model = EddingtonBaryonModel(
            radiative_efficiency=0.1,
            eddington_ratio=0.7,
            duty_cycle=0.4,
        )
        self.assertAlmostEqual(
            eddington_efolding_time_myr(0.1, 0.7, 0.4),
            model.black_hole_efolding_time_myr,
            places=12,
        )

    def test_target_already_reached_requires_no_time(self):
        self.assertEqual(growth_time_to_target_myr(1.0e6, 1.0e5), 0.0)

    def test_required_activity_exactly_fills_window(self):
        duration = growth_time_to_target_myr(100.0, 1.0e6)
        self.assertAlmostEqual(
            required_eddington_activity(100.0, 1.0e6, duration),
            1.0,
            places=12,
        )


if __name__ == "__main__":
    unittest.main()
