import unittest

import numpy as np

from sidm_bh.baryons import HernquistBaryons
from sidm_bh.halos import NFWProfile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.stage3 import static_baryon_equilibrium_state
from sidm_bh.units import SimulationScales


class StaticBaryonInitialConditionTest(unittest.TestCase):
    def setUp(self):
        self.profile = NFWProfile(3.7, 30.0)
        self.scales = SimulationScales(30.0, 3.7)
        self.grid = SphericalGrid.from_log_spacing(0.005 / 30.0, 20.0, 32)

    def test_no_baryon_control_returns_zero_external_mass(self):
        state, baryon_mass = static_baryon_equilibrium_state(
            self.profile,
            self.grid,
            self.scales,
            black_hole_mass_msun=100.0,
        )

        np.testing.assert_allclose(baryon_mass, 0.0)
        self.assertTrue(np.all(state.velocity_dispersion > 0.0))

    def test_baryon_equilibrium_uses_same_density_and_higher_pressure(self):
        control, _ = static_baryon_equilibrium_state(
            self.profile,
            self.grid,
            self.scales,
            black_hole_mass_msun=100.0,
        )
        baryonic, baryon_mass = static_baryon_equilibrium_state(
            self.profile,
            self.grid,
            self.scales,
            black_hole_mass_msun=100.0,
            baryons=HernquistBaryons(5.0e4, 0.3),
        )

        np.testing.assert_allclose(baryonic.density, control.density)
        self.assertTrue(np.all(np.diff(baryon_mass) > 0.0))
        self.assertTrue(np.all(baryonic.pressure > control.pressure))


if __name__ == "__main__":
    unittest.main()
