import unittest

import numpy as np

from sidm_bh.constants import M_SUN_CGS, MYR_CGS, PC_CGS
from sidm_bh.cooling import (
    cloudy_cooling_state,
    equilibrium_temperature_k,
    load_cloudy_cooling_table,
)


class CloudyCoolingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = load_cloudy_cooling_table()

    def test_pinned_table_axes_and_shape(self):
        self.assertEqual(self.table.primordial_cooling.shape, (29, 161))
        self.assertEqual(self.table.solar_metal_cooling.shape, (29, 161))
        self.assertEqual(self.table.mean_molecular_weight.shape, (29, 161))
        np.testing.assert_allclose(
            self.table.log_hydrogen_density[[0, -1]],
            [-10.0, 4.0],
        )
        np.testing.assert_allclose(
            self.table.log_temperature[[0, -1]],
            [1.0, 9.0],
        )

    def test_temperature_is_consistent_with_sound_speed_and_mmw(self):
        temperature, mmw = equilibrium_temperature_k(
            self.table,
            hydrogen_number_density_cm3=9.0e3,
            sound_speed_km_s=15.0,
        )
        self.assertAlmostEqual(temperature, 1.525e4, delta=20.0)
        self.assertAlmostEqual(mmw, 0.932, delta=0.003)

    def test_atomic_cooling_turns_on_above_hydrogen_threshold(self):
        cool = cloudy_cooling_state(self.table, 300.0, 10.0, 0.0)
        hot = cloudy_cooling_state(self.table, 300.0, 15.0, 0.0)

        self.assertGreater(cool.cooling_time_myr, 1.0e-3)
        self.assertLess(hot.cooling_time_myr, 1.0e-6)
        self.assertLess(hot.cooling_time_myr, cool.cooling_time_myr)

    def test_metallicity_shortens_cool_gas_cooling_time(self):
        primordial = cloudy_cooling_state(self.table, 300.0, 10.0, 0.0)
        solar = cloudy_cooling_state(self.table, 300.0, 10.0, 1.0)

        self.assertLess(solar.cooling_time_myr, primordial.cooling_time_myr)
        self.assertAlmostEqual(
            primordial.hydrogen_number_density_cm3,
            9.226e3,
            delta=2.0,
        )

    def test_numba_cooling_time_matches_python_diagnostic(self):
        from sidm_bh.fast_evolution import _cloudy_cooling_time_code

        expected = cloudy_cooling_state(self.table, 300.0, 10.0, 0.1)
        actual = _cloudy_cooling_time_code(
            remaining_mass=1.0,
            initial_total_mass=1.0,
            initial_scale_radius=0.1,
            current_scale_radius=0.1,
            initial_sound_speed=0.1,
            thermal_energy=0.0,
            adiabatic_index=5.0 / 3.0,
            initial_gas_density_cgs=300.0 * M_SUN_CGS / PC_CGS**3,
            hydrogen_mass_fraction=0.76,
            metallicity_solar=0.1,
            cooling_rate_multiplier=1.0,
            velocity_unit_cgs=1.0e7,
            time_unit_s=MYR_CGS,
            log_density_axis=self.table.log_hydrogen_density,
            log_temperature_axis=self.table.log_temperature,
            primordial_cooling=self.table.primordial_cooling,
            solar_metal_cooling=self.table.solar_metal_cooling,
            mean_molecular_weight=self.table.mean_molecular_weight,
        )
        self.assertLess(
            abs(actual / expected.cooling_time_myr - 1.0),
            1.0e-3,
        )


if __name__ == "__main__":
    unittest.main()
