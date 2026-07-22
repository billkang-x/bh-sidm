"""Optically thin equilibrium cooling from the Grackle Cloudy no-UVB table."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .constants import (
    K_BOLTZMANN_CGS,
    M_PROTON_CGS,
    M_SUN_CGS,
    MYR_CGS,
    PC_CGS,
)


DEFAULT_CLOUDY_NO_UVB_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "cooling"
    / "CloudyData_noUVB.h5"
)


@dataclass(frozen=True)
class CloudyCoolingTable:
    log_hydrogen_density: np.ndarray
    log_temperature: np.ndarray
    primordial_cooling: np.ndarray
    solar_metal_cooling: np.ndarray
    mean_molecular_weight: np.ndarray


@dataclass(frozen=True)
class CoolingState:
    temperature_k: float
    mean_molecular_weight: float
    hydrogen_number_density_cm3: float
    cooling_coefficient_erg_cm3_s: float
    cooling_time_myr: float


def load_cloudy_cooling_table(
    path: str | Path = DEFAULT_CLOUDY_NO_UVB_PATH,
) -> CloudyCoolingTable:
    """Load the pinned Grackle Cloudy no-UV-background cooling table."""

    source = Path(path)
    with h5py.File(source, "r") as data:
        primordial_dataset = data["CoolingRates/Primordial/Cooling"]
        metal_dataset = data["CoolingRates/Metals/Cooling"]
        mmw_dataset = data["CoolingRates/Primordial/MMW"]
        log_density = np.asarray(
            primordial_dataset.attrs["Parameter1"], dtype=float
        )
        temperature = np.asarray(
            primordial_dataset.attrs["Temperature"], dtype=float
        )
        primordial = np.asarray(primordial_dataset, dtype=float)
        metals = np.asarray(metal_dataset, dtype=float)
        mmw = np.asarray(mmw_dataset, dtype=float)
    expected_shape = (len(log_density), len(temperature))
    if primordial.shape != expected_shape:
        raise ValueError("primordial cooling table shape does not match its axes")
    if metals.shape != expected_shape or mmw.shape != expected_shape:
        raise ValueError("Cloudy cooling datasets must share one grid")
    if np.any(np.diff(log_density) <= 0.0) or np.any(np.diff(temperature) <= 0.0):
        raise ValueError("Cloudy cooling axes must be strictly increasing")
    if np.any(primordial < 0.0) or np.any(metals < 0.0):
        raise ValueError("Cloudy cooling coefficients cannot be negative")
    if np.any(mmw <= 0.0):
        raise ValueError("Cloudy mean molecular weights must be positive")
    return CloudyCoolingTable(
        log_hydrogen_density=log_density,
        log_temperature=np.log10(temperature),
        primordial_cooling=primordial,
        solar_metal_cooling=metals,
        mean_molecular_weight=mmw,
    )


def _interpolation_indices(axis: np.ndarray, value: float) -> tuple[int, float]:
    clipped = float(np.clip(value, axis[0], axis[-1]))
    upper = int(np.searchsorted(axis, clipped, side="right"))
    upper = min(max(upper, 1), len(axis) - 1)
    lower = upper - 1
    fraction = (clipped - axis[lower]) / (axis[upper] - axis[lower])
    return lower, float(fraction)


def _bilinear(
    table: np.ndarray,
    density_index: int,
    density_fraction: float,
    temperature_index: int,
    temperature_fraction: float,
) -> float:
    lower_temperature = (
        (1.0 - density_fraction) * table[density_index, temperature_index]
        + density_fraction * table[density_index + 1, temperature_index]
    )
    upper_temperature = (
        (1.0 - density_fraction) * table[density_index, temperature_index + 1]
        + density_fraction * table[density_index + 1, temperature_index + 1]
    )
    return float(
        (1.0 - temperature_fraction) * lower_temperature
        + temperature_fraction * upper_temperature
    )


def _interpolate_mmw(
    table: CloudyCoolingTable,
    log_hydrogen_density: float,
    log_temperature: float,
) -> float:
    density_index, density_fraction = _interpolation_indices(
        table.log_hydrogen_density,
        log_hydrogen_density,
    )
    temperature_index, temperature_fraction = _interpolation_indices(
        table.log_temperature,
        log_temperature,
    )
    return _bilinear(
        table.mean_molecular_weight,
        density_index,
        density_fraction,
        temperature_index,
        temperature_fraction,
    )


def _interpolate_cooling(
    values: np.ndarray,
    table: CloudyCoolingTable,
    log_hydrogen_density: float,
    log_temperature: float,
) -> float:
    density_index, density_fraction = _interpolation_indices(
        table.log_hydrogen_density,
        log_hydrogen_density,
    )
    temperature_index, temperature_fraction = _interpolation_indices(
        table.log_temperature,
        log_temperature,
    )
    log_values = np.log10(np.maximum(values, 1.0e-300))
    return 10.0 ** _bilinear(
        log_values,
        density_index,
        density_fraction,
        temperature_index,
        temperature_fraction,
    )


def equilibrium_temperature_k(
    table: CloudyCoolingTable,
    hydrogen_number_density_cm3: float,
    sound_speed_km_s: float,
    adiabatic_index: float = 5.0 / 3.0,
) -> tuple[float, float]:
    """Return temperature and table MMW consistent with adiabatic sound speed."""

    if hydrogen_number_density_cm3 <= 0.0:
        raise ValueError("hydrogen_number_density_cm3 must be positive")
    if sound_speed_km_s <= 0.0:
        raise ValueError("sound_speed_km_s must be positive")
    if adiabatic_index <= 1.0:
        raise ValueError("adiabatic_index must exceed one")
    log_density = np.log10(hydrogen_number_density_cm3)
    target_temperature = (
        (sound_speed_km_s * 1.0e5) ** 2
        * M_PROTON_CGS
        / (adiabatic_index * K_BOLTZMANN_CGS)
    )
    target_log = np.log10(target_temperature)
    lower = float(table.log_temperature[0])
    upper = float(table.log_temperature[-1])
    for _ in range(48):
        middle = 0.5 * (lower + upper)
        mmw = _interpolate_mmw(table, log_density, middle)
        residual = middle - np.log10(mmw) - target_log
        if residual > 0.0:
            upper = middle
        else:
            lower = middle
    log_temperature = 0.5 * (lower + upper)
    mmw = _interpolate_mmw(table, log_density, log_temperature)
    return 10.0**log_temperature, mmw


def cloudy_cooling_state(
    table: CloudyCoolingTable,
    gas_density_msun_pc3: float,
    sound_speed_km_s: float,
    metallicity_solar: float,
    cooling_rate_multiplier: float = 1.0,
    hydrogen_mass_fraction: float = 0.76,
    adiabatic_index: float = 5.0 / 3.0,
) -> CoolingState:
    """Return an optically thin equilibrium Cloudy cooling state."""

    if gas_density_msun_pc3 <= 0.0:
        raise ValueError("gas_density_msun_pc3 must be positive")
    if metallicity_solar < 0.0:
        raise ValueError("metallicity_solar cannot be negative")
    if cooling_rate_multiplier <= 0.0:
        raise ValueError("cooling_rate_multiplier must be positive")
    if not 0.0 < hydrogen_mass_fraction <= 1.0:
        raise ValueError("hydrogen_mass_fraction must lie in (0, 1]")
    density_cgs = gas_density_msun_pc3 * M_SUN_CGS / PC_CGS**3
    hydrogen_density = (
        hydrogen_mass_fraction * density_cgs / M_PROTON_CGS
    )
    temperature, mmw = equilibrium_temperature_k(
        table,
        hydrogen_density,
        sound_speed_km_s,
        adiabatic_index,
    )
    log_density = np.log10(hydrogen_density)
    log_temperature = np.log10(temperature)
    primordial = _interpolate_cooling(
        table.primordial_cooling,
        table,
        log_density,
        log_temperature,
    )
    metals = _interpolate_cooling(
        table.solar_metal_cooling,
        table,
        log_density,
        log_temperature,
    )
    cooling_coefficient = cooling_rate_multiplier * (
        primordial + metallicity_solar * metals
    )
    sound_speed_cgs = sound_speed_km_s * 1.0e5
    thermal_energy_density = density_cgs * sound_speed_cgs**2 / (
        adiabatic_index * (adiabatic_index - 1.0)
    )
    cooling_rate_density = hydrogen_density**2 * cooling_coefficient
    cooling_time = thermal_energy_density / cooling_rate_density / MYR_CGS
    return CoolingState(
        temperature_k=float(temperature),
        mean_molecular_weight=float(mmw),
        hydrogen_number_density_cm3=float(hydrogen_density),
        cooling_coefficient_erg_cm3_s=float(cooling_coefficient),
        cooling_time_myr=float(cooling_time),
    )
