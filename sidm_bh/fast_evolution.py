"""Numba-accelerated MC-Roe evolution for long baseline simulations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import K_BOLTZMANN_CGS, M_PROTON_CGS
from .cooling import CloudyCoolingTable
from .mesh import SphericalGrid
from .sidm import (
    CONDUCTIVITY_A,
    CONDUCTIVITY_B,
    CONDUCTIVITY_C,
    effective_cross_section_ratio_tables,
)
from .state import FluidState

try:
    from numba import njit
except ImportError as error:  # pragma: no cover - exercised only without numba
    raise ImportError("fast_evolution requires the optional numba package") from error


@dataclass(frozen=True)
class FastEvolutionResult:
    final_state: FluidState
    final_black_hole_mass_code: float
    peak_accretion_rate_code: float
    final_accretion_rate_code: float
    num_steps: int
    max_mass_budget_residual_code: float
    sample_times_code: np.ndarray
    black_hole_masses_code: np.ndarray
    accretion_rates_code: np.ndarray
    dark_matter_accreted_masses_code: np.ndarray
    dark_matter_supply_rates_code: np.ndarray
    dark_matter_supplied_masses_code: np.ndarray
    inner_dark_matter_reservoir_masses_code: np.ndarray
    baryon_accretion_rates_code: np.ndarray
    baryon_accreted_masses_code: np.ndarray
    baryon_gas_consumed_masses_code: np.ndarray
    baryon_remaining_masses_code: np.ndarray
    baryon_scale_radii_code: np.ndarray
    feedback_to_binding_ratios: np.ndarray
    baryon_thermal_energies_code: np.ndarray
    density_snapshots: np.ndarray
    radial_velocity_snapshots: np.ndarray
    velocity_dispersion_snapshots: np.ndarray


@njit(cache=True)
def _minmod3(first: float, second: float, third: float) -> float:
    if first > 0.0 and second > 0.0 and third > 0.0:
        return min(first, second, third)
    if first < 0.0 and second < 0.0 and third < 0.0:
        return max(first, second, third)
    return 0.0


@njit(cache=True)
def _effective_cross_section_ratio(
    dispersion: float,
    transition_speed: float,
    log_ratio_min: float,
    log_ratio_step: float,
    ratio_table: np.ndarray,
) -> float:
    """Interpolate the precomputed Maxwell-averaged Rutherford ratio."""

    log_ratio = np.log10(dispersion / transition_speed)
    coordinate = (log_ratio - log_ratio_min) / log_ratio_step
    if coordinate <= 0.0:
        return ratio_table[0]
    upper_index = len(ratio_table) - 1
    if coordinate >= upper_index:
        return ratio_table[upper_index]
    lower_index = int(coordinate)
    fraction = coordinate - lower_index
    return (
        (1.0 - fraction) * ratio_table[lower_index]
        + fraction * ratio_table[lower_index + 1]
    )


@njit(cache=True)
def _mc_slopes(values: np.ndarray, centers: np.ndarray, slopes: np.ndarray) -> None:
    size = len(values)
    slopes[0] = 0.0
    slopes[size - 1] = 0.0
    for i in range(1, size - 1):
        backward = (values[i] - values[i - 1]) / (centers[i] - centers[i - 1])
        forward = (values[i + 1] - values[i]) / (centers[i + 1] - centers[i])
        centered = (values[i + 1] - values[i - 1]) / (
            centers[i + 1] - centers[i - 1]
        )
        slopes[i] = _minmod3(centered, 2.0 * backward, 2.0 * forward)


@njit(cache=True)
def _physical_flux(rho: float, velocity: float, pressure: float) -> tuple:
    specific_energy = 1.5 * pressure / rho + 0.5 * velocity * velocity
    return (
        rho * velocity,
        rho * velocity * velocity + pressure,
        velocity * (rho * specific_energy + pressure),
    )


@njit(cache=True)
def _entropy_absolute(value: float, width: float) -> float:
    absolute = abs(value)
    if width > 0.0 and absolute < width:
        return 0.5 * (value * value / width + width)
    return absolute


@njit(cache=True)
def _roe_flux(
    rho_left: float,
    velocity_left: float,
    pressure_left: float,
    rho_right: float,
    velocity_right: float,
    pressure_right: float,
    entropy_fix: float,
) -> tuple:
    flux_left = _physical_flux(rho_left, velocity_left, pressure_left)
    flux_right = _physical_flux(rho_right, velocity_right, pressure_right)
    sqrt_left = np.sqrt(rho_left)
    sqrt_right = np.sqrt(rho_right)
    weight = sqrt_left + sqrt_right
    energy_left = 1.5 * pressure_left / rho_left + 0.5 * velocity_left**2
    energy_right = 1.5 * pressure_right / rho_right + 0.5 * velocity_right**2
    enthalpy_left = energy_left + pressure_left / rho_left
    enthalpy_right = energy_right + pressure_right / rho_right
    velocity = (sqrt_left * velocity_left + sqrt_right * velocity_right) / weight
    enthalpy = (sqrt_left * enthalpy_left + sqrt_right * enthalpy_right) / weight
    sound_squared = (2.0 / 3.0) * (enthalpy - 0.5 * velocity**2)
    if sound_squared <= 0.0:
        return (np.nan, np.nan, np.nan)
    sound = np.sqrt(sound_squared)
    roe_density = sqrt_left * sqrt_right
    density_jump = rho_right - rho_left
    velocity_jump = velocity_right - velocity_left
    pressure_jump = pressure_right - pressure_left
    strength_left = (
        pressure_jump - roe_density * sound * velocity_jump
    ) / (2.0 * sound_squared)
    strength_contact = density_jump - pressure_jump / sound_squared
    strength_right = (
        pressure_jump + roe_density * sound * velocity_jump
    ) / (2.0 * sound_squared)
    width = entropy_fix * sound
    lambda_left = _entropy_absolute(velocity - sound, width)
    lambda_contact = abs(velocity)
    lambda_right = _entropy_absolute(velocity + sound, width)

    dissipation_density = (
        lambda_left * strength_left
        + lambda_contact * strength_contact
        + lambda_right * strength_right
    )
    dissipation_momentum = (
        lambda_left * strength_left * (velocity - sound)
        + lambda_contact * strength_contact * velocity
        + lambda_right * strength_right * (velocity + sound)
    )
    dissipation_energy = (
        lambda_left * strength_left * (enthalpy - velocity * sound)
        + lambda_contact * strength_contact * 0.5 * velocity**2
        + lambda_right * strength_right * (enthalpy + velocity * sound)
    )
    return (
        0.5 * (flux_left[0] + flux_right[0] - dissipation_density),
        0.5 * (flux_left[1] + flux_right[1] - dissipation_momentum),
        0.5 * (flux_left[2] + flux_right[2] - dissipation_energy),
    )


@njit(cache=True)
def _record_snapshot(
    sample_index: int,
    black_hole_mass: float,
    accretion_rate: float,
    density: np.ndarray,
    velocity: np.ndarray,
    dispersion: np.ndarray,
    dark_matter_accreted_mass: float,
    dark_matter_supply_rate: float,
    dark_matter_supplied_mass: float,
    inner_dark_matter_reservoir_mass: float,
    baryon_accretion_rate: float,
    baryon_accreted_mass: float,
    baryon_gas_consumed_mass: float,
    baryon_remaining_mass: float,
    baryon_scale_radius: float,
    feedback_to_binding_ratio: float,
    baryon_thermal_energy: float,
    black_hole_history: np.ndarray,
    accretion_history: np.ndarray,
    dark_matter_accreted_history: np.ndarray,
    dark_matter_supply_history: np.ndarray,
    dark_matter_supplied_history: np.ndarray,
    inner_dark_matter_reservoir_history: np.ndarray,
    baryon_accretion_history: np.ndarray,
    baryon_accreted_history: np.ndarray,
    baryon_gas_consumed_history: np.ndarray,
    baryon_remaining_history: np.ndarray,
    baryon_scale_radius_history: np.ndarray,
    feedback_to_binding_history: np.ndarray,
    baryon_thermal_energy_history: np.ndarray,
    density_snapshots: np.ndarray,
    velocity_snapshots: np.ndarray,
    dispersion_snapshots: np.ndarray,
) -> None:
    black_hole_history[sample_index] = black_hole_mass
    accretion_history[sample_index] = accretion_rate
    dark_matter_accreted_history[sample_index] = dark_matter_accreted_mass
    dark_matter_supply_history[sample_index] = dark_matter_supply_rate
    dark_matter_supplied_history[sample_index] = dark_matter_supplied_mass
    inner_dark_matter_reservoir_history[sample_index] = (
        inner_dark_matter_reservoir_mass
    )
    baryon_accretion_history[sample_index] = baryon_accretion_rate
    baryon_accreted_history[sample_index] = baryon_accreted_mass
    baryon_gas_consumed_history[sample_index] = baryon_gas_consumed_mass
    baryon_remaining_history[sample_index] = baryon_remaining_mass
    baryon_scale_radius_history[sample_index] = baryon_scale_radius
    feedback_to_binding_history[sample_index] = feedback_to_binding_ratio
    baryon_thermal_energy_history[sample_index] = baryon_thermal_energy
    density_snapshots[sample_index] = density
    velocity_snapshots[sample_index] = velocity
    dispersion_snapshots[sample_index] = dispersion


@njit(cache=True)
def _dark_bondi_rate(
    black_hole_mass: float,
    density: float,
    velocity_dispersion: float,
    bondi_lambda: float,
) -> float:
    """Return the gamma=5/3 dark Bondi rate in the paper's code units."""

    sound_speed = np.sqrt(5.0 / 3.0) * velocity_dispersion
    if black_hole_mass <= 0.0 or density <= 0.0 or sound_speed <= 0.0:
        return 0.0
    return bondi_lambda * black_hole_mass**2 * density / sound_speed**3


@njit(cache=True)
def _current_bondi_coefficient(
    base_coefficient: float,
    evolve_ambient: int,
    remaining_mass: float,
    initial_total_mass: float,
    initial_scale_radius: float,
    current_scale_radius: float,
    initial_sound_speed: float,
    relative_velocity: float,
    thermal_energy: float,
    adiabatic_index: float,
) -> float:
    if base_coefficient < 0.0 or evolve_ambient == 0:
        return base_coefficient
    if remaining_mass <= 0.0:
        return 0.0
    density_factor = (
        remaining_mass
        / initial_total_mass
        * (initial_scale_radius / current_scale_radius) ** 3
    )
    deposited_specific_energy = thermal_energy / remaining_mass
    current_sound_squared = initial_sound_speed**2 + (
        adiabatic_index
        * (adiabatic_index - 1.0)
        * deposited_specific_energy
    )
    initial_effective_squared = initial_sound_speed**2 + relative_velocity**2
    current_effective_squared = current_sound_squared + relative_velocity**2
    speed_factor = (
        initial_effective_squared / current_effective_squared
    ) ** 1.5
    return base_coefficient * density_factor * speed_factor


@njit(cache=True)
def _axis_index(axis: np.ndarray, value: float) -> tuple:
    if value <= axis[0]:
        return 0, 0.0
    if value >= axis[-1]:
        return len(axis) - 2, 1.0
    lower = 0
    upper = len(axis) - 1
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if value >= axis[middle]:
            lower = middle
        else:
            upper = middle
    fraction = (value - axis[lower]) / (axis[upper] - axis[lower])
    return lower, fraction


@njit(cache=True)
def _density_interpolated_value(
    values: np.ndarray,
    density_index: int,
    density_fraction: float,
    temperature_index: int,
) -> float:
    return (
        (1.0 - density_fraction) * values[density_index, temperature_index]
        + density_fraction * values[density_index + 1, temperature_index]
    )


@njit(cache=True)
def _cloudy_cooling_time_code(
    remaining_mass: float,
    initial_total_mass: float,
    initial_scale_radius: float,
    current_scale_radius: float,
    initial_sound_speed: float,
    thermal_energy: float,
    adiabatic_index: float,
    initial_gas_density_cgs: float,
    hydrogen_mass_fraction: float,
    metallicity_solar: float,
    cooling_rate_multiplier: float,
    velocity_unit_cgs: float,
    time_unit_s: float,
    log_density_axis: np.ndarray,
    log_temperature_axis: np.ndarray,
    primordial_cooling: np.ndarray,
    solar_metal_cooling: np.ndarray,
    mean_molecular_weight: np.ndarray,
) -> float:
    if remaining_mass <= 0.0:
        return 1.0e300
    density_factor = (
        remaining_mass
        / initial_total_mass
        * (initial_scale_radius / current_scale_radius) ** 3
    )
    density_cgs = initial_gas_density_cgs * density_factor
    if density_cgs <= 0.0:
        return 1.0e300
    hydrogen_density = (
        hydrogen_mass_fraction * density_cgs / M_PROTON_CGS
    )
    log_density = np.log10(hydrogen_density)
    density_index, density_fraction = _axis_index(
        log_density_axis,
        log_density,
    )
    sound_squared = initial_sound_speed**2 + (
        adiabatic_index
        * (adiabatic_index - 1.0)
        * thermal_energy
        / remaining_mass
    )
    sound_squared_cgs = sound_squared * velocity_unit_cgs**2
    target_log = np.log10(
        sound_squared_cgs
        * M_PROTON_CGS
        / (adiabatic_index * K_BOLTZMANN_CGS)
    )

    lower_temperature = 0
    upper_temperature = len(log_temperature_axis) - 1
    lower_mmw = _density_interpolated_value(
        mean_molecular_weight,
        density_index,
        density_fraction,
        lower_temperature,
    )
    upper_mmw = _density_interpolated_value(
        mean_molecular_weight,
        density_index,
        density_fraction,
        upper_temperature,
    )
    lower_effective = log_temperature_axis[0] - np.log10(lower_mmw)
    upper_effective = log_temperature_axis[-1] - np.log10(upper_mmw)
    if target_log <= lower_effective:
        lower_temperature = 0
        temperature_fraction = 0.0
    elif target_log >= upper_effective:
        lower_temperature = len(log_temperature_axis) - 2
        temperature_fraction = 1.0
    else:
        while upper_temperature - lower_temperature > 1:
            middle = (lower_temperature + upper_temperature) // 2
            middle_mmw = _density_interpolated_value(
                mean_molecular_weight,
                density_index,
                density_fraction,
                middle,
            )
            middle_effective = (
                log_temperature_axis[middle] - np.log10(middle_mmw)
            )
            if target_log >= middle_effective:
                lower_temperature = middle
            else:
                upper_temperature = middle
        lower_mmw = _density_interpolated_value(
            mean_molecular_weight,
            density_index,
            density_fraction,
            lower_temperature,
        )
        upper_mmw = _density_interpolated_value(
            mean_molecular_weight,
            density_index,
            density_fraction,
            lower_temperature + 1,
        )
        lower_effective = (
            log_temperature_axis[lower_temperature] - np.log10(lower_mmw)
        )
        upper_effective = (
            log_temperature_axis[lower_temperature + 1] - np.log10(upper_mmw)
        )
        temperature_fraction = (
            (target_log - lower_effective)
            / (upper_effective - lower_effective)
        )

    primordial_lower = _density_interpolated_value(
        primordial_cooling,
        density_index,
        density_fraction,
        lower_temperature,
    )
    primordial_upper = _density_interpolated_value(
        primordial_cooling,
        density_index,
        density_fraction,
        lower_temperature + 1,
    )
    primordial_rate = 10.0 ** (
        (1.0 - temperature_fraction)
        * np.log10(max(primordial_lower, 1.0e-300))
        + temperature_fraction
        * np.log10(max(primordial_upper, 1.0e-300))
    )
    metal_lower = _density_interpolated_value(
        solar_metal_cooling,
        density_index,
        density_fraction,
        lower_temperature,
    )
    metal_upper = _density_interpolated_value(
        solar_metal_cooling,
        density_index,
        density_fraction,
        lower_temperature + 1,
    )
    metal_rate = 10.0 ** (
        (1.0 - temperature_fraction)
        * np.log10(max(metal_lower, 1.0e-300))
        + temperature_fraction
        * np.log10(max(metal_upper, 1.0e-300))
    )
    cooling_coefficient = cooling_rate_multiplier * (
        primordial_rate + metallicity_solar * metal_rate
    )
    if cooling_coefficient <= 0.0:
        return 1.0e300
    thermal_energy_density = density_cgs * sound_squared_cgs / (
        adiabatic_index * (adiabatic_index - 1.0)
    )
    cooling_rate_density = hydrogen_density**2 * cooling_coefficient
    return thermal_energy_density / cooling_rate_density / time_unit_s


@njit(cache=True)
def _evolve_kernel(
    initial_density: np.ndarray,
    initial_velocity: np.ndarray,
    initial_dispersion: np.ndarray,
    interfaces: np.ndarray,
    centers: np.ndarray,
    volumes: np.ndarray,
    widths: np.ndarray,
    sample_times: np.ndarray,
    initial_black_hole_mass: float,
    sigma_over_m: float,
    cross_section_model: int,
    cross_section_transition_speed: float,
    cross_section_log_ratio_min: float,
    cross_section_log_ratio_step: float,
    cross_section_k3_table: np.ndarray,
    cross_section_k5_table: np.ndarray,
    cfl_number: float,
    entropy_fix: float,
    inner_dark_matter_mass: float,
    dark_capture_model: int,
    dark_bondi_lambda: float,
    dark_flux_capture_mass_threshold: float,
    initial_dark_matter_reservoir: float,
    baryon_enclosed_mass: np.ndarray,
    baryon_assembly_time: float,
    evolving_baryon_total_mass: float,
    evolving_baryon_scale_radius: float,
    baryon_eddington_inflow_coefficient: float,
    baryon_bondi_inflow_coefficient: float,
    evolve_bondi_ambient: int,
    baryon_initial_sound_speed: float,
    baryon_relative_velocity: float,
    feedback_thermal_energy_per_consumed_mass: float,
    baryon_cooling_time: float,
    use_cloudy_cooling: int,
    baryon_initial_gas_density_cgs: float,
    baryon_hydrogen_mass_fraction: float,
    baryon_metallicity_solar: float,
    baryon_cooling_rate_multiplier: float,
    baryon_velocity_unit_cgs: float,
    baryon_time_unit_s: float,
    cooling_log_density_axis: np.ndarray,
    cooling_log_temperature_axis: np.ndarray,
    cooling_primordial_table: np.ndarray,
    cooling_solar_metal_table: np.ndarray,
    cooling_mmw_table: np.ndarray,
    baryon_adiabatic_index: float,
    baryon_radiative_efficiency: float,
    feedback_ratio_per_consumed_mass: float,
    feedback_expansion_exponent: float,
    source_integration: int,
    max_steps: int,
) -> tuple:
    size = len(initial_density)
    num_samples = len(sample_times)
    density = initial_density.copy()
    velocity = initial_velocity.copy()
    dispersion = initial_dispersion.copy()
    pressure = np.empty(size)
    total_mass = np.empty(size)
    slopes_density = np.empty(size)
    slopes_velocity = np.empty(size)
    slopes_pressure = np.empty(size)
    flux = np.empty((3, size + 1))
    next_density = np.empty(size)
    next_velocity = np.empty(size)
    next_dispersion = np.empty(size)
    kappa = np.empty(size)
    interface_kappa = np.empty(size + 1)
    lower = np.empty(max(size - 1, 0))
    diagonal = np.empty(size)
    upper = np.empty(max(size - 1, 0))
    rhs = np.empty(size)
    solution = np.empty(size)

    black_hole_history = np.empty(num_samples)
    accretion_history = np.empty(num_samples)
    dark_matter_accreted_history = np.empty(num_samples)
    dark_matter_supply_history = np.empty(num_samples)
    dark_matter_supplied_history = np.empty(num_samples)
    inner_dark_matter_reservoir_history = np.empty(num_samples)
    baryon_accretion_history = np.empty(num_samples)
    baryon_accreted_history = np.empty(num_samples)
    baryon_gas_consumed_history = np.empty(num_samples)
    baryon_remaining_history = np.empty(num_samples)
    baryon_scale_radius_history = np.empty(num_samples)
    feedback_to_binding_history = np.empty(num_samples)
    baryon_thermal_energy_history = np.empty(num_samples)
    density_snapshots = np.empty((num_samples, size))
    velocity_snapshots = np.empty((num_samples, size))
    dispersion_snapshots = np.empty((num_samples, size))

    black_hole_mass = initial_black_hole_mass
    dark_matter_accreted_mass = 0.0
    dark_matter_supplied_mass = 0.0
    dark_matter_reservoir = initial_dark_matter_reservoir
    baryon_accreted_mass = 0.0
    baryon_gas_consumed_mass = 0.0
    feedback_to_binding_ratio = 0.0
    baryon_thermal_energy = 0.0
    time = 0.0
    steps = 0
    cumulative_outer_flux = 0.0
    cumulative_untracked_inner_outflow = 0.0
    initial_fluid_mass = 0.0
    for i in range(size):
        initial_fluid_mass += density[i] * volumes[i]
    initial_total_mass = initial_fluid_mass + black_hole_mass + dark_matter_reservoir
    max_budget_residual = 0.0
    initial_supply_rate = max(
        -interfaces[0] ** 2 * density[0] * velocity[0], 0.0
    )
    initial_rate = initial_supply_rate
    if dark_capture_model != 0:
        initial_rate = 0.0
        if (
            dark_capture_model == 2
            and black_hole_mass >= dark_flux_capture_mass_threshold
        ):
            drain_time = np.sqrt(
                interfaces[0] ** 3
                / max(black_hole_mass + dark_matter_reservoir, 1.0e-300)
            )
            initial_rate = initial_supply_rate + dark_matter_reservoir / drain_time
        elif dark_matter_reservoir > 0.0:
            initial_rate = _dark_bondi_rate(
                black_hole_mass,
                density[0],
                dispersion[0],
                dark_bondi_lambda,
            )
    peak_rate = initial_rate
    initial_baryon_fraction = 1.0
    if baryon_assembly_time > 0.0:
        initial_baryon_fraction = 0.0
    initial_baryon_remaining = 0.0
    initial_baryon_scale_radius = 0.0
    initial_baryon_rate = 0.0
    if evolving_baryon_total_mass >= 0.0:
        initial_baryon_remaining = initial_baryon_fraction * evolving_baryon_total_mass
        initial_baryon_scale_radius = evolving_baryon_scale_radius
        if initial_baryon_remaining > 0.0:
            initial_gas_rate = (
                baryon_eddington_inflow_coefficient * black_hole_mass
            )
            initial_bondi_coefficient = _current_bondi_coefficient(
                baryon_bondi_inflow_coefficient,
                evolve_bondi_ambient,
                initial_baryon_remaining,
                evolving_baryon_total_mass,
                evolving_baryon_scale_radius,
                initial_baryon_scale_radius,
                baryon_initial_sound_speed,
                baryon_relative_velocity,
                baryon_thermal_energy,
                baryon_adiabatic_index,
            )
            if initial_bondi_coefficient >= 0.0:
                initial_gas_rate = min(
                    initial_gas_rate,
                    initial_bondi_coefficient * black_hole_mass**2,
                )
            initial_baryon_rate = (
                (1.0 - baryon_radiative_efficiency)
                * initial_gas_rate
            )
    _record_snapshot(
        0,
        black_hole_mass,
        initial_rate,
        density,
        velocity,
        dispersion,
        dark_matter_accreted_mass,
        initial_supply_rate,
        dark_matter_supplied_mass,
        dark_matter_reservoir,
        initial_baryon_rate,
        baryon_accreted_mass,
        baryon_gas_consumed_mass,
        initial_baryon_remaining,
        initial_baryon_scale_radius,
        feedback_to_binding_ratio,
        baryon_thermal_energy,
        black_hole_history,
        accretion_history,
        dark_matter_accreted_history,
        dark_matter_supply_history,
        dark_matter_supplied_history,
        inner_dark_matter_reservoir_history,
        baryon_accretion_history,
        baryon_accreted_history,
        baryon_gas_consumed_history,
        baryon_remaining_history,
        baryon_scale_radius_history,
        feedback_to_binding_history,
        baryon_thermal_energy_history,
        density_snapshots,
        velocity_snapshots,
        dispersion_snapshots,
    )
    sample_index = 1
    status = 0

    while time < sample_times[-1] and steps < max_steps:
        baryon_fraction = 1.0
        if baryon_assembly_time > 0.0:
            if time <= 0.0:
                baryon_fraction = 0.0
            elif time < baryon_assembly_time:
                phase = time / baryon_assembly_time
                baryon_fraction = phase * phase * (3.0 - 2.0 * phase)
        baryon_remaining_mass = 0.0
        current_baryon_scale_radius = 0.0
        if evolving_baryon_total_mass >= 0.0:
            baryon_remaining_mass = max(
                baryon_fraction * evolving_baryon_total_mass
                - baryon_gas_consumed_mass,
                0.0,
            )
            current_baryon_scale_radius = evolving_baryon_scale_radius * (
                1.0 + feedback_to_binding_ratio
            ) ** feedback_expansion_exponent
        running_mass = inner_dark_matter_mass + dark_matter_reservoir
        hyperbolic_dt = 1.0e300
        gravity_dt = 1.0e300
        for i in range(size):
            pressure[i] = density[i] * dispersion[i] ** 2
            center_volume = (centers[i] ** 3 - interfaces[i] ** 3) / 3.0
            dark_mass_center = running_mass + density[i] * center_volume
            baryon_mass_center = baryon_fraction * baryon_enclosed_mass[i]
            if evolving_baryon_total_mass >= 0.0:
                baryon_mass_center = (
                    baryon_remaining_mass
                    * centers[i] ** 2
                    / (centers[i] + current_baryon_scale_radius) ** 2
                )
            total_mass[i] = dark_mass_center + black_hole_mass + baryon_mass_center
            running_mass += density[i] * volumes[i]
            sound = np.sqrt(5.0 / 3.0) * dispersion[i]
            local_hyperbolic = widths[i] / (abs(velocity[i]) + sound)
            if local_hyperbolic < hyperbolic_dt:
                hyperbolic_dt = local_hyperbolic
            acceleration = total_mass[i] / centers[i] ** 2
            if acceleration > 0.0:
                local_gravity = sound / acceleration
                if local_gravity < gravity_dt:
                    gravity_dt = local_gravity
        dt = cfl_number * hyperbolic_dt
        if source_integration == 0:
            dt = min(dt, cfl_number * gravity_dt)
        if sample_index < num_samples:
            dt = min(dt, sample_times[sample_index] - time)
        if dt <= 0.0 or not np.isfinite(dt):
            status = 4
            break

        baryon_gas_inflow_rate = 0.0
        baryon_black_hole_growth_rate = 0.0
        if evolving_baryon_total_mass >= 0.0 and baryon_remaining_mass > 0.0:
            baryon_gas_inflow_rate = (
                baryon_eddington_inflow_coefficient * black_hole_mass
            )
            current_bondi_coefficient = _current_bondi_coefficient(
                baryon_bondi_inflow_coefficient,
                evolve_bondi_ambient,
                baryon_remaining_mass,
                evolving_baryon_total_mass,
                evolving_baryon_scale_radius,
                current_baryon_scale_radius,
                baryon_initial_sound_speed,
                baryon_relative_velocity,
                baryon_thermal_energy,
                baryon_adiabatic_index,
            )
            if current_bondi_coefficient >= 0.0:
                baryon_gas_inflow_rate = min(
                    baryon_gas_inflow_rate,
                    current_bondi_coefficient * black_hole_mass**2,
                )
            baryon_gas_inflow_rate = min(
                baryon_gas_inflow_rate,
                baryon_remaining_mass / dt,
            )
            baryon_black_hole_growth_rate = (
                (1.0 - baryon_radiative_efficiency)
                * baryon_gas_inflow_rate
            )

        dark_matter_supply_rate = max(
            -interfaces[0] ** 2 * density[0] * velocity[0], 0.0
        )
        dark_matter_capture_rate = dark_matter_supply_rate
        if dark_capture_model != 0:
            if (
                dark_capture_model == 2
                and black_hole_mass >= dark_flux_capture_mass_threshold
            ):
                drain_time = np.sqrt(
                    interfaces[0] ** 3
                    / max(black_hole_mass + dark_matter_reservoir, 1.0e-300)
                )
                dark_matter_capture_rate = (
                    dark_matter_supply_rate
                    + dark_matter_reservoir / max(drain_time, dt)
                )
            else:
                dark_matter_capture_rate = _dark_bondi_rate(
                    black_hole_mass,
                    density[0],
                    dispersion[0],
                    dark_bondi_lambda,
                )
        inner_flux = interfaces[0] ** 2 * density[0] * velocity[0]
        outer_flux = interfaces[-1] ** 2 * density[-1] * velocity[-1]

        _mc_slopes(density, centers, slopes_density)
        _mc_slopes(velocity, centers, slopes_velocity)
        _mc_slopes(pressure, centers, slopes_pressure)
        for interface_index in range(size + 1):
            if interface_index == 0:
                left_cell = 0
                right_cell = 0
            elif interface_index == size:
                left_cell = size - 1
                right_cell = size - 1
            else:
                left_cell = interface_index - 1
                right_cell = interface_index
            radius_interface = interfaces[interface_index]
            rho_left = density[left_cell] + slopes_density[left_cell] * (
                radius_interface - centers[left_cell]
            )
            velocity_left = velocity[left_cell] + slopes_velocity[left_cell] * (
                radius_interface - centers[left_cell]
            )
            pressure_left = pressure[left_cell] + slopes_pressure[left_cell] * (
                radius_interface - centers[left_cell]
            )
            rho_right = density[right_cell] + slopes_density[right_cell] * (
                radius_interface - centers[right_cell]
            )
            velocity_right = velocity[right_cell] + slopes_velocity[right_cell] * (
                radius_interface - centers[right_cell]
            )
            pressure_right = pressure[right_cell] + slopes_pressure[right_cell] * (
                radius_interface - centers[right_cell]
            )
            if (
                rho_left <= 0.0
                or rho_right <= 0.0
                or pressure_left <= 0.0
                or pressure_right <= 0.0
            ):
                status = 1
                break
            flux_values = _roe_flux(
                rho_left,
                velocity_left,
                pressure_left,
                rho_right,
                velocity_right,
                pressure_right,
                entropy_fix,
            )
            flux[0, interface_index] = flux_values[0]
            flux[1, interface_index] = flux_values[1]
            flux[2, interface_index] = flux_values[2]
        if status != 0:
            break

        for i in range(size):
            area_inner = interfaces[i] ** 2
            area_outer = interfaces[i + 1] ** 2
            divergence_density = (
                area_outer * flux[0, i + 1] - area_inner * flux[0, i]
            ) / volumes[i]
            divergence_momentum = (
                area_outer * flux[1, i + 1] - area_inner * flux[1, i]
            ) / volumes[i]
            divergence_energy = (
                area_outer * flux[2, i + 1] - area_inner * flux[2, i]
            ) / volumes[i]
            momentum = density[i] * velocity[i]
            energy = density[i] * (
                1.5 * dispersion[i] ** 2 + 0.5 * velocity[i] ** 2
            )
            acceleration = -total_mass[i] / centers[i] ** 2
            rho_new = density[i] - dt * divergence_density
            if source_integration == 0:
                momentum_new = momentum + dt * (
                    density[i] * acceleration
                    + 2.0 * pressure[i] / centers[i]
                    - divergence_momentum
                )
                energy_new = energy + dt * (
                    density[i] * velocity[i] * acceleration - divergence_energy
                )
            else:
                momentum_star = momentum + dt * (
                    2.0 * pressure[i] / centers[i] - divergence_momentum
                )
                energy_star = energy - dt * divergence_energy
                velocity_increment = acceleration * dt
                momentum_new = momentum_star + rho_new * velocity_increment
                energy_new = (
                    energy_star
                    + momentum_star * velocity_increment
                    + 0.5 * rho_new * velocity_increment**2
                )
            if rho_new <= 0.0:
                status = 2
                break
            velocity_new = momentum_new / rho_new
            thermal_energy = energy_new / rho_new - 0.5 * velocity_new**2
            dispersion_squared = (2.0 / 3.0) * thermal_energy
            if dispersion_squared <= 0.0 or not np.isfinite(dispersion_squared):
                status = 2
                break
            next_density[i] = rho_new
            next_velocity[i] = velocity_new
            next_dispersion[i] = np.sqrt(dispersion_squared)
        if status != 0:
            break

        density, next_density = next_density, density
        velocity, next_velocity = next_velocity, velocity
        dispersion, next_dispersion = next_dispersion, dispersion

        if sigma_over_m > 0.0:
            for i in range(size):
                lmfp_sigma_over_m = sigma_over_m
                smfp_sigma_over_m = sigma_over_m
                if cross_section_model == 1:
                    lmfp_sigma_over_m *= _effective_cross_section_ratio(
                        dispersion[i],
                        cross_section_transition_speed,
                        cross_section_log_ratio_min,
                        cross_section_log_ratio_step,
                        cross_section_k3_table,
                    )
                    smfp_sigma_over_m *= _effective_cross_section_ratio(
                        dispersion[i],
                        cross_section_transition_speed,
                        cross_section_log_ratio_min,
                        cross_section_log_ratio_step,
                        cross_section_k5_table,
                    )
                smfp_inverse = smfp_sigma_over_m / (
                    CONDUCTIVITY_B * dispersion[i]
                )
                lmfp_inverse = 1.0 / (
                    CONDUCTIVITY_A
                    * CONDUCTIVITY_C
                    * density[i]
                    * dispersion[i] ** 3
                    * lmfp_sigma_over_m
                )
                kappa[i] = 1.5 / (smfp_inverse + lmfp_inverse)
            interface_kappa[0] = kappa[0]
            interface_kappa[size] = kappa[size - 1]
            for i in range(1, size):
                interface_kappa[i] = 0.5 * (kappa[i - 1] + kappa[i])
            for i in range(size):
                factor = 2.0 * dt / (3.0 * volumes[i])
                diagonal[i] = density[i]
                rhs[i] = density[i] * dispersion[i] ** 2
                if i > 0:
                    lower[i - 1] = -factor * (
                        interfaces[i] ** 2
                        * interface_kappa[i]
                        / (centers[i] - centers[i - 1])
                    )
                    diagonal[i] -= lower[i - 1]
                if i < size - 1:
                    upper[i] = -factor * (
                        interfaces[i + 1] ** 2
                        * interface_kappa[i + 1]
                        / (centers[i + 1] - centers[i])
                    )
                    diagonal[i] -= upper[i]
            for i in range(1, size):
                multiplier = lower[i - 1] / diagonal[i - 1]
                diagonal[i] -= multiplier * upper[i - 1]
                rhs[i] -= multiplier * rhs[i - 1]
            solution[size - 1] = rhs[size - 1] / diagonal[size - 1]
            for i in range(size - 2, -1, -1):
                solution[i] = (
                    rhs[i] - upper[i] * solution[i + 1]
                ) / diagonal[i]
            for i in range(size):
                if solution[i] <= 0.0 or not np.isfinite(solution[i]):
                    status = 3
                    break
                dispersion[i] = np.sqrt(solution[i])
            if status != 0:
                break

        dark_matter_supply_delta = dark_matter_supply_rate * dt
        dark_matter_supplied_mass += dark_matter_supply_delta
        if dark_capture_model == 0:
            dark_matter_delta = dark_matter_supply_delta
            cumulative_untracked_inner_outflow += max(inner_flux, 0.0) * dt
        else:
            dark_matter_reservoir += dark_matter_supply_delta
            requested_inner_outflow = max(inner_flux, 0.0) * dt
            reservoir_outflow = min(requested_inner_outflow, dark_matter_reservoir)
            dark_matter_reservoir -= reservoir_outflow
            cumulative_untracked_inner_outflow += (
                requested_inner_outflow - reservoir_outflow
            )
            dark_matter_delta = min(
                dark_matter_capture_rate * dt,
                dark_matter_reservoir,
            )
            dark_matter_reservoir -= dark_matter_delta
            dark_matter_capture_rate = dark_matter_delta / dt
        if dark_matter_capture_rate > peak_rate:
            peak_rate = dark_matter_capture_rate
        baryon_gas_delta = baryon_gas_inflow_rate * dt
        baryon_black_hole_delta = baryon_black_hole_growth_rate * dt
        dark_matter_accreted_mass += dark_matter_delta
        baryon_gas_consumed_mass += baryon_gas_delta
        baryon_accreted_mass += baryon_black_hole_delta
        feedback_to_binding_ratio += (
            feedback_ratio_per_consumed_mass * baryon_gas_delta
        )
        if use_cloudy_cooling != 0 and baryon_remaining_mass > 0.0:
            remaining_after_inflow = max(
                baryon_remaining_mass - baryon_gas_delta,
                1.0e-300,
            )
            scale_radius_after_feedback = evolving_baryon_scale_radius * (
                1.0 + feedback_to_binding_ratio
            ) ** feedback_expansion_exponent
            source_energy_rate = (
                feedback_thermal_energy_per_consumed_mass
                * baryon_gas_inflow_rate
            )
            lower_energy = 0.0
            upper_energy = baryon_thermal_energy + source_energy_rate * dt
            trial_energy = 0.5 * (lower_energy + upper_energy)
            cooling_converged = False
            for _ in range(80):
                physical_cooling_time = _cloudy_cooling_time_code(
                    remaining_after_inflow,
                    evolving_baryon_total_mass,
                    evolving_baryon_scale_radius,
                    scale_radius_after_feedback,
                    baryon_initial_sound_speed,
                    trial_energy,
                    baryon_adiabatic_index,
                    baryon_initial_gas_density_cgs,
                    baryon_hydrogen_mass_fraction,
                    baryon_metallicity_solar,
                    baryon_cooling_rate_multiplier,
                    baryon_velocity_unit_cgs,
                    baryon_time_unit_s,
                    cooling_log_density_axis,
                    cooling_log_temperature_axis,
                    cooling_primordial_table,
                    cooling_solar_metal_table,
                    cooling_mmw_table,
                )
                cooled_fraction = -np.expm1(-dt / physical_cooling_time)
                candidate_energy = (
                    baryon_thermal_energy * (1.0 - cooled_fraction)
                    + source_energy_rate
                    * physical_cooling_time
                    * cooled_fraction
                )
                if candidate_energy > trial_energy:
                    lower_energy = trial_energy
                else:
                    upper_energy = trial_energy
                updated_trial = 0.5 * (lower_energy + upper_energy)
                if upper_energy - lower_energy <= (
                    1.0e-10 * max(updated_trial, 1.0e-300)
                ):
                    trial_energy = updated_trial
                    cooling_converged = True
                    break
                trial_energy = updated_trial
            if not cooling_converged:
                status = 6
                break
            baryon_thermal_energy = max(trial_energy, 0.0)
        elif baryon_cooling_time > 0.0:
            cooled_fraction = -np.expm1(-dt / baryon_cooling_time)
            cooling_factor = 1.0 - cooled_fraction
            baryon_thermal_energy = (
                baryon_thermal_energy * cooling_factor
                + feedback_thermal_energy_per_consumed_mass
                * baryon_gas_inflow_rate
                * baryon_cooling_time
                * cooled_fraction
            )
        else:
            baryon_thermal_energy += (
                feedback_thermal_energy_per_consumed_mass * baryon_gas_delta
            )
        black_hole_mass += dark_matter_delta + baryon_black_hole_delta
        cumulative_outer_flux += outer_flux * dt
        time += dt
        steps += 1

        fluid_mass = 0.0
        for i in range(size):
            fluid_mass += density[i] * volumes[i]
        residual = abs(
            fluid_mass
            + initial_black_hole_mass
            + dark_matter_accreted_mass
            + dark_matter_reservoir
            + cumulative_outer_flux
            - cumulative_untracked_inner_outflow
            - initial_total_mass
        )
        if residual > max_budget_residual:
            max_budget_residual = residual

        if sample_index < num_samples and time >= sample_times[sample_index] - 1.0e-14:
            current_rate = dark_matter_capture_rate
            current_supply_rate = max(
                -interfaces[0] ** 2 * density[0] * velocity[0],
                0.0,
            )
            sampled_baryon_fraction = 1.0
            if baryon_assembly_time > 0.0 and time < baryon_assembly_time:
                phase = max(time / baryon_assembly_time, 0.0)
                sampled_baryon_fraction = phase * phase * (3.0 - 2.0 * phase)
            sampled_baryon_remaining = 0.0
            sampled_baryon_scale_radius = 0.0
            sampled_baryon_rate = 0.0
            if evolving_baryon_total_mass >= 0.0:
                sampled_baryon_remaining = max(
                    sampled_baryon_fraction * evolving_baryon_total_mass
                    - baryon_gas_consumed_mass,
                    0.0,
                )
                sampled_baryon_scale_radius = evolving_baryon_scale_radius * (
                    1.0 + feedback_to_binding_ratio
                ) ** feedback_expansion_exponent
                if sampled_baryon_remaining > 0.0:
                    sampled_gas_rate = (
                        baryon_eddington_inflow_coefficient * black_hole_mass
                    )
                    sampled_bondi_coefficient = _current_bondi_coefficient(
                        baryon_bondi_inflow_coefficient,
                        evolve_bondi_ambient,
                        sampled_baryon_remaining,
                        evolving_baryon_total_mass,
                        evolving_baryon_scale_radius,
                        sampled_baryon_scale_radius,
                        baryon_initial_sound_speed,
                        baryon_relative_velocity,
                        baryon_thermal_energy,
                        baryon_adiabatic_index,
                    )
                    if sampled_bondi_coefficient >= 0.0:
                        sampled_gas_rate = min(
                            sampled_gas_rate,
                            sampled_bondi_coefficient * black_hole_mass**2,
                        )
                    sampled_baryon_rate = (
                        (1.0 - baryon_radiative_efficiency)
                        * sampled_gas_rate
                    )
            _record_snapshot(
                sample_index,
                black_hole_mass,
                current_rate,
                density,
                velocity,
                dispersion,
                dark_matter_accreted_mass,
                current_supply_rate,
                dark_matter_supplied_mass,
                dark_matter_reservoir,
                sampled_baryon_rate,
                baryon_accreted_mass,
                baryon_gas_consumed_mass,
                sampled_baryon_remaining,
                sampled_baryon_scale_radius,
                feedback_to_binding_ratio,
                baryon_thermal_energy,
                black_hole_history,
                accretion_history,
                dark_matter_accreted_history,
                dark_matter_supply_history,
                dark_matter_supplied_history,
                inner_dark_matter_reservoir_history,
                baryon_accretion_history,
                baryon_accreted_history,
                baryon_gas_consumed_history,
                baryon_remaining_history,
                baryon_scale_radius_history,
                feedback_to_binding_history,
                baryon_thermal_energy_history,
                density_snapshots,
                velocity_snapshots,
                dispersion_snapshots,
            )
            sample_index += 1

    if steps >= max_steps and time < sample_times[-1]:
        status = 5
    final_rate = accretion_history[max(sample_index - 1, 0)]
    return (
        status,
        density,
        velocity,
        dispersion,
        black_hole_mass,
        peak_rate,
        final_rate,
        steps,
        max_budget_residual,
        black_hole_history,
        accretion_history,
        dark_matter_accreted_history,
        dark_matter_supply_history,
        dark_matter_supplied_history,
        inner_dark_matter_reservoir_history,
        baryon_accretion_history,
        baryon_accreted_history,
        baryon_gas_consumed_history,
        baryon_remaining_history,
        baryon_scale_radius_history,
        feedback_to_binding_history,
        baryon_thermal_energy_history,
        density_snapshots,
        velocity_snapshots,
        dispersion_snapshots,
    )


def evolve_mc_roe_fast(
    initial_state: FluidState,
    grid: SphericalGrid,
    sample_times_code: np.ndarray,
    initial_black_hole_mass_code: float,
    sigma_over_m_code: float,
    cross_section_model: str = "constant",
    cross_section_velocity_scale_code: float | None = None,
    cfl_number: float = 0.2,
    entropy_fix: float = 0.1,
    inner_dark_matter_mass_code: float = 0.0,
    dark_capture_model: str = "boundary_flux",
    dark_bondi_lambda: float = 0.25,
    dark_flux_capture_mass_threshold_code: float = float("inf"),
    initial_dark_matter_reservoir_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
    baryon_assembly_time_code: float | None = None,
    evolving_baryon_total_mass_code: float | None = None,
    evolving_baryon_scale_radius_code: float | None = None,
    baryon_eddington_inflow_coefficient_code: float = 0.0,
    baryon_bondi_inflow_coefficient_code: float = -1.0,
    evolve_bondi_ambient: bool = False,
    baryon_initial_sound_speed_code: float = 0.0,
    baryon_relative_velocity_code: float = 0.0,
    feedback_thermal_energy_per_consumed_mass_code: float = 0.0,
    baryon_cooling_time_code: float | None = None,
    baryon_cloudy_cooling_table: CloudyCoolingTable | None = None,
    baryon_initial_gas_density_cgs: float = 0.0,
    baryon_hydrogen_mass_fraction: float = 0.76,
    baryon_metallicity_solar: float = 0.0,
    baryon_cooling_rate_multiplier: float = 1.0,
    baryon_velocity_unit_cgs: float = 0.0,
    baryon_time_unit_s: float = 0.0,
    baryon_adiabatic_index: float = 5.0 / 3.0,
    baryon_radiative_efficiency: float = 0.1,
    feedback_ratio_per_consumed_mass_code: float = 0.0,
    feedback_expansion_exponent: float = 0.5,
    source_integration: str = "euler",
    max_steps: int = 10_000_000,
) -> FastEvolutionResult:
    """Run the MC-Roe solver and retain full states at requested times."""

    sample_times = np.asarray(sample_times_code, dtype=float)
    if sample_times.ndim != 1 or len(sample_times) < 2:
        raise ValueError("sample_times_code must contain at least start and end")
    if sample_times[0] != 0.0 or np.any(np.diff(sample_times) <= 0.0):
        raise ValueError("sample_times_code must start at zero and increase")
    if len(initial_state.density) != grid.num_cells:
        raise ValueError("initial_state length must match grid")
    if initial_black_hole_mass_code < 0.0:
        raise ValueError("initial_black_hole_mass_code cannot be negative")
    if sigma_over_m_code < 0.0:
        raise ValueError("sigma_over_m_code cannot be negative")
    if cross_section_model not in ("constant", "rutherford"):
        raise ValueError("cross_section_model must be 'constant' or 'rutherford'")
    if cross_section_model == "rutherford" and (
        cross_section_velocity_scale_code is None
        or not np.isfinite(cross_section_velocity_scale_code)
        or cross_section_velocity_scale_code <= 0.0
    ):
        raise ValueError(
            "rutherford cross section requires a positive velocity scale"
        )
    if dark_capture_model not in (
        "boundary_flux",
        "bondi_reservoir",
        "influence_gated",
    ):
        raise ValueError(
            "dark_capture_model must be 'boundary_flux', 'bondi_reservoir', "
            "or 'influence_gated'"
        )
    if dark_bondi_lambda < 0.0:
        raise ValueError("dark_bondi_lambda cannot be negative")
    if initial_dark_matter_reservoir_code < 0.0:
        raise ValueError("initial_dark_matter_reservoir_code cannot be negative")
    if dark_capture_model == "influence_gated" and (
        not np.isfinite(dark_flux_capture_mass_threshold_code)
        or dark_flux_capture_mass_threshold_code <= 0.0
    ):
        raise ValueError(
            "influence_gated capture requires a positive finite mass threshold"
        )
    if (
        dark_capture_model == "boundary_flux"
        and initial_dark_matter_reservoir_code != 0.0
    ):
        raise ValueError(
            "an initial dark reservoir requires the bondi_reservoir model"
        )
    if source_integration not in ("euler", "gravity_kick"):
        raise ValueError("source_integration must be 'euler' or 'gravity_kick'")
    if baryon_assembly_time_code is not None and baryon_assembly_time_code <= 0.0:
        raise ValueError("baryon_assembly_time_code must be positive when supplied")
    if (evolving_baryon_total_mass_code is None) != (
        evolving_baryon_scale_radius_code is None
    ):
        raise ValueError(
            "evolving baryon total mass and scale radius must be supplied together"
        )
    if evolving_baryon_total_mass_code is not None:
        if evolving_baryon_total_mass_code <= 0.0:
            raise ValueError("evolving_baryon_total_mass_code must be positive")
        if evolving_baryon_scale_radius_code <= 0.0:
            raise ValueError("evolving_baryon_scale_radius_code must be positive")
        if baryon_enclosed_mass_code is not None:
            raise ValueError(
                "fixed and evolving baryon mass profiles cannot be combined"
            )
    elif (
        baryon_eddington_inflow_coefficient_code != 0.0
        or baryon_bondi_inflow_coefficient_code >= 0.0
    ):
        raise ValueError("baryon accretion requires an evolving baryon reservoir")
    if baryon_eddington_inflow_coefficient_code < 0.0:
        raise ValueError("baryon Eddington coefficient cannot be negative")
    if baryon_bondi_inflow_coefficient_code < 0.0 and (
        baryon_bondi_inflow_coefficient_code != -1.0
    ):
        raise ValueError("baryon Bondi coefficient must be non-negative or -1")
    if evolve_bondi_ambient:
        if evolving_baryon_total_mass_code is None:
            raise ValueError("evolving Bondi ambient requires a baryon reservoir")
        if baryon_bondi_inflow_coefficient_code < 0.0:
            raise ValueError("evolving Bondi ambient requires Bondi accretion")
        if baryon_initial_sound_speed_code <= 0.0:
            raise ValueError("baryon_initial_sound_speed_code must be positive")
    if baryon_relative_velocity_code < 0.0:
        raise ValueError("baryon_relative_velocity_code cannot be negative")
    if feedback_thermal_energy_per_consumed_mass_code < 0.0:
        raise ValueError(
            "feedback thermal energy per consumed mass cannot be negative"
        )
    if baryon_cooling_time_code is not None and (
        np.isnan(baryon_cooling_time_code) or baryon_cooling_time_code <= 0.0
    ):
        raise ValueError("baryon_cooling_time_code must be positive when supplied")
    if baryon_cloudy_cooling_table is not None:
        if baryon_cooling_time_code is not None:
            raise ValueError("fixed and Cloudy cooling cannot be combined")
        if not evolve_bondi_ambient:
            raise ValueError("Cloudy cooling requires an evolving Bondi ambient")
        if baryon_initial_gas_density_cgs <= 0.0:
            raise ValueError("Cloudy cooling requires positive initial gas density")
        if not 0.0 < baryon_hydrogen_mass_fraction <= 1.0:
            raise ValueError("baryon_hydrogen_mass_fraction must lie in (0, 1]")
        if baryon_metallicity_solar < 0.0:
            raise ValueError("baryon_metallicity_solar cannot be negative")
        if baryon_cooling_rate_multiplier <= 0.0:
            raise ValueError("baryon_cooling_rate_multiplier must be positive")
        if baryon_velocity_unit_cgs <= 0.0 or baryon_time_unit_s <= 0.0:
            raise ValueError("Cloudy cooling requires positive code-unit scales")
    if baryon_adiabatic_index <= 1.0:
        raise ValueError("baryon_adiabatic_index must exceed one")
    if not 0.0 < baryon_radiative_efficiency < 1.0:
        raise ValueError("baryon_radiative_efficiency must lie in (0, 1)")
    if feedback_ratio_per_consumed_mass_code < 0.0:
        raise ValueError("feedback ratio per consumed mass cannot be negative")
    if feedback_expansion_exponent < 0.0:
        raise ValueError("feedback_expansion_exponent cannot be negative")
    baryons = (
        np.zeros(grid.num_cells)
        if baryon_enclosed_mass_code is None
        else np.asarray(baryon_enclosed_mass_code, dtype=float)
    )
    if baryons.shape != (grid.num_cells,):
        raise ValueError("baryon_enclosed_mass_code shape must match grid")
    if baryon_cloudy_cooling_table is None:
        cooling_log_density_axis = np.array([0.0, 1.0])
        cooling_log_temperature_axis = np.array([0.0, 1.0])
        cooling_primordial_table = np.ones((2, 2))
        cooling_solar_metal_table = np.ones((2, 2))
        cooling_mmw_table = np.ones((2, 2))
    else:
        cooling_log_density_axis = np.ascontiguousarray(
            baryon_cloudy_cooling_table.log_hydrogen_density,
            dtype=float,
        )
        cooling_log_temperature_axis = np.ascontiguousarray(
            baryon_cloudy_cooling_table.log_temperature,
            dtype=float,
        )
        cooling_primordial_table = np.ascontiguousarray(
            baryon_cloudy_cooling_table.primordial_cooling,
            dtype=float,
        )
        cooling_solar_metal_table = np.ascontiguousarray(
            baryon_cloudy_cooling_table.solar_metal_cooling,
            dtype=float,
        )
        cooling_mmw_table = np.ascontiguousarray(
            baryon_cloudy_cooling_table.mean_molecular_weight,
            dtype=float,
        )

    if cross_section_model == "constant":
        cross_section_log_ratio_min = 0.0
        cross_section_log_ratio_step = 1.0
        cross_section_k3_table = np.ones(2, dtype=float)
        cross_section_k5_table = np.ones(2, dtype=float)
        cross_section_transition_speed = 1.0
    else:
        log_ratio, cross_section_k3_table, cross_section_k5_table = (
            effective_cross_section_ratio_tables()
        )
        cross_section_log_ratio_min = float(log_ratio[0])
        cross_section_log_ratio_step = float(log_ratio[1] - log_ratio[0])
        cross_section_transition_speed = float(
            cross_section_velocity_scale_code
        )

    output = _evolve_kernel(
        initial_state.density,
        initial_state.radial_velocity,
        initial_state.velocity_dispersion,
        grid.interfaces_code,
        grid.centers_code,
        grid.cell_volumes_code,
        grid.widths_code,
        sample_times,
        initial_black_hole_mass_code,
        sigma_over_m_code,
        0 if cross_section_model == "constant" else 1,
        cross_section_transition_speed,
        cross_section_log_ratio_min,
        cross_section_log_ratio_step,
        np.ascontiguousarray(cross_section_k3_table, dtype=float),
        np.ascontiguousarray(cross_section_k5_table, dtype=float),
        cfl_number,
        entropy_fix,
        inner_dark_matter_mass_code,
        {
            "boundary_flux": 0,
            "bondi_reservoir": 1,
            "influence_gated": 2,
        }[dark_capture_model],
        dark_bondi_lambda,
        dark_flux_capture_mass_threshold_code,
        initial_dark_matter_reservoir_code,
        baryons,
        -1.0 if baryon_assembly_time_code is None else baryon_assembly_time_code,
        (
            -1.0
            if evolving_baryon_total_mass_code is None
            else evolving_baryon_total_mass_code
        ),
        (
            0.0
            if evolving_baryon_scale_radius_code is None
            else evolving_baryon_scale_radius_code
        ),
        baryon_eddington_inflow_coefficient_code,
        baryon_bondi_inflow_coefficient_code,
        1 if evolve_bondi_ambient else 0,
        baryon_initial_sound_speed_code,
        baryon_relative_velocity_code,
        feedback_thermal_energy_per_consumed_mass_code,
        (
            -1.0
            if baryon_cooling_time_code is None
            or np.isinf(baryon_cooling_time_code)
            else baryon_cooling_time_code
        ),
        1 if baryon_cloudy_cooling_table is not None else 0,
        baryon_initial_gas_density_cgs,
        baryon_hydrogen_mass_fraction,
        baryon_metallicity_solar,
        baryon_cooling_rate_multiplier,
        baryon_velocity_unit_cgs,
        baryon_time_unit_s,
        cooling_log_density_axis,
        cooling_log_temperature_axis,
        cooling_primordial_table,
        cooling_solar_metal_table,
        cooling_mmw_table,
        baryon_adiabatic_index,
        baryon_radiative_efficiency,
        feedback_ratio_per_consumed_mass_code,
        feedback_expansion_exponent,
        0 if source_integration == "euler" else 1,
        max_steps,
    )
    status = output[0]
    if status != 0:
        messages = {
            1: "MC reconstruction or Roe average became non-physical",
            2: "hyperbolic update produced a non-positive state",
            3: "conduction update produced a non-positive state",
            4: "time integration produced a non-positive timestep",
            5: "maximum step count reached before the final sample",
            6: "Cloudy cooling fixed-point iteration did not converge",
        }
        raise RuntimeError(messages.get(status, f"unknown fast solver status {status}"))
    return FastEvolutionResult(
        final_state=FluidState(output[1], output[2], output[3]),
        final_black_hole_mass_code=output[4],
        peak_accretion_rate_code=output[5],
        final_accretion_rate_code=output[6],
        num_steps=output[7],
        max_mass_budget_residual_code=output[8],
        sample_times_code=sample_times,
        black_hole_masses_code=output[9],
        accretion_rates_code=output[10],
        dark_matter_accreted_masses_code=output[11],
        dark_matter_supply_rates_code=output[12],
        dark_matter_supplied_masses_code=output[13],
        inner_dark_matter_reservoir_masses_code=output[14],
        baryon_accretion_rates_code=output[15],
        baryon_accreted_masses_code=output[16],
        baryon_gas_consumed_masses_code=output[17],
        baryon_remaining_masses_code=output[18],
        baryon_scale_radii_code=output[19],
        feedback_to_binding_ratios=output[20],
        baryon_thermal_energies_code=output[21],
        density_snapshots=output[22],
        radial_velocity_snapshots=output[23],
        velocity_dispersion_snapshots=output[24],
    )
