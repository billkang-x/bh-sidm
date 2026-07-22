import math
import unittest

from sidm_bh.baryons import (
    HernquistBaryons,
    eddington_accretion_rate_msun_per_myr,
    smoothstep_mass_fraction,
)
from sidm_bh.constants import G_CGS, M_SUN_CGS, PC_CGS


class HernquistBaryonsTest(unittest.TestCase):
    def test_enclosed_mass_limits(self):
        profile = HernquistBaryons(total_mass_msun=1.0e6, scale_radius_pc=10.0)

        self.assertEqual(profile.enclosed_mass_cgs(0.0), 0.0)

        enclosed_far = profile.enclosed_mass_cgs(1.0e9)
        self.assertAlmostEqual(enclosed_far / M_SUN_CGS, 1.0e6, delta=1.0)

    def test_static_potential_matches_hernquist_formula(self):
        profile = HernquistBaryons(total_mass_msun=2.0e5, scale_radius_pc=3.0)
        radius_pc = 7.0

        expected = -G_CGS * 2.0e5 * M_SUN_CGS / ((radius_pc + 3.0) * PC_CGS)
        self.assertAlmostEqual(profile.potential_cgs(radius_pc), expected)

    def test_growing_mass_fraction(self):
        profile = HernquistBaryons(
            total_mass_msun=1.0e6,
            scale_radius_pc=10.0,
            growth_time_myr=20.0,
        )

        self.assertEqual(profile.mass_fraction(0.0), 0.0)
        self.assertAlmostEqual(profile.mass_fraction(20.0), 1.0 - math.exp(-1.0))

    def test_eddington_rate_scales_linearly(self):
        rate_100 = eddington_accretion_rate_msun_per_myr(100.0)
        rate_200 = eddington_accretion_rate_msun_per_myr(200.0)

        self.assertGreater(rate_100, 0.0)
        self.assertAlmostEqual(rate_200 / rate_100, 2.0)

    def test_smoothstep_assembly_reaches_exact_endpoints(self):
        self.assertEqual(smoothstep_mass_fraction(0.0, 2.0), 0.0)
        self.assertEqual(smoothstep_mass_fraction(2.0, 2.0), 1.0)
        self.assertAlmostEqual(smoothstep_mass_fraction(1.0, 2.0), 0.5)


if __name__ == "__main__":
    unittest.main()
