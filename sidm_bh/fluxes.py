"""Convective Euler fluxes for the dimensionless SIDM fluid equations."""

from __future__ import annotations

import numpy as np

from .state import FluidState

ADIABATIC_INDEX = 5.0 / 3.0


def sound_speed_code(state: FluidState) -> np.ndarray:
    """Return adiabatic sound speed sqrt(gamma) * v."""

    return np.sqrt(ADIABATIC_INDEX) * state.velocity_dispersion


def max_signal_speed_code(state: FluidState) -> float:
    """Return max |u| + c_s for CFL estimates."""

    return float(np.max(np.abs(state.radial_velocity) + sound_speed_code(state)))


def convective_flux_code(state: FluidState) -> np.ndarray:
    """Return the non-diffusive flux vector F_c(U) in code units."""

    density = state.density
    velocity = state.radial_velocity
    pressure = state.pressure
    energy_density = density * state.specific_energy
    return np.vstack(
        [
            density * velocity,
            density * velocity**2 + pressure,
            velocity * (energy_density + pressure),
        ]
    )


def rusanov_flux_code(left: FluidState, right: FluidState) -> np.ndarray:
    """Return the local Lax-Friedrichs flux at one or more interfaces.

    ``left`` and ``right`` contain the reconstructed primitive states on
    either side of each interface. The current baseline solver uses
    piecewise-constant reconstruction; keeping the Riemann flux independent
    of that choice makes it straightforward to add MC reconstruction later.
    """

    if left.density.shape != right.density.shape:
        raise ValueError("left and right states must have matching shapes")

    left_flux = convective_flux_code(left)
    right_flux = convective_flux_code(right)
    signal_speed = np.maximum(
        np.abs(left.radial_velocity) + sound_speed_code(left),
        np.abs(right.radial_velocity) + sound_speed_code(right),
    )
    return 0.5 * (left_flux + right_flux) - 0.5 * signal_speed * (
        right.conservative - left.conservative
    )


def roe_flux_code(
    left: FluidState,
    right: FluidState,
    entropy_fix: float = 0.1,
) -> np.ndarray:
    """Return the Roe flux in Eqs. (A7)-(A16).

    A Harten entropy fix is applied to the two acoustic eigenvalues. The
    contact wave is left unchanged so stationary contacts remain exact.
    """

    if left.density.shape != right.density.shape:
        raise ValueError("left and right states must have matching shapes")
    if entropy_fix < 0.0:
        raise ValueError("entropy_fix cannot be negative")

    sqrt_left_density = np.sqrt(left.density)
    sqrt_right_density = np.sqrt(right.density)
    density_weight = sqrt_left_density + sqrt_right_density
    left_enthalpy = left.specific_energy + left.velocity_dispersion**2
    right_enthalpy = right.specific_energy + right.velocity_dispersion**2
    roe_velocity = (
        sqrt_left_density * left.radial_velocity
        + sqrt_right_density * right.radial_velocity
    ) / density_weight
    roe_enthalpy = (
        sqrt_left_density * left_enthalpy
        + sqrt_right_density * right_enthalpy
    ) / density_weight
    sound_speed_squared = (ADIABATIC_INDEX - 1.0) * (
        roe_enthalpy - 0.5 * roe_velocity**2
    )
    if np.any(sound_speed_squared <= 0.0):
        return rusanov_flux_code(left, right)
    roe_sound_speed = np.sqrt(sound_speed_squared)
    roe_density = sqrt_left_density * sqrt_right_density

    density_jump = right.density - left.density
    velocity_jump = right.radial_velocity - left.radial_velocity
    pressure_jump = right.pressure - left.pressure
    acoustic_left_strength = (
        pressure_jump - roe_density * roe_sound_speed * velocity_jump
    ) / (2.0 * sound_speed_squared)
    contact_strength = density_jump - pressure_jump / sound_speed_squared
    acoustic_right_strength = (
        pressure_jump + roe_density * roe_sound_speed * velocity_jump
    ) / (2.0 * sound_speed_squared)

    eigenvalues = np.vstack(
        [
            roe_velocity - roe_sound_speed,
            roe_velocity,
            roe_velocity + roe_sound_speed,
        ]
    )
    absolute_eigenvalues = np.abs(eigenvalues)
    if entropy_fix > 0.0:
        entropy_width = entropy_fix * roe_sound_speed
        for wave in (0, 2):
            inside = absolute_eigenvalues[wave] < entropy_width
            absolute_eigenvalues[wave] = np.where(
                inside,
                0.5
                * (
                    eigenvalues[wave] ** 2 / entropy_width
                    + entropy_width
                ),
                absolute_eigenvalues[wave],
            )

    eigenvectors = np.empty((3, 3, left.density.size), dtype=float)
    eigenvectors[:, 0] = np.vstack(
        [
            np.ones_like(roe_velocity),
            roe_velocity - roe_sound_speed,
            roe_enthalpy - roe_velocity * roe_sound_speed,
        ]
    )
    eigenvectors[:, 1] = np.vstack(
        [
            np.ones_like(roe_velocity),
            roe_velocity,
            0.5 * roe_velocity**2,
        ]
    )
    eigenvectors[:, 2] = np.vstack(
        [
            np.ones_like(roe_velocity),
            roe_velocity + roe_sound_speed,
            roe_enthalpy + roe_velocity * roe_sound_speed,
        ]
    )
    strengths = np.vstack(
        [acoustic_left_strength, contact_strength, acoustic_right_strength]
    )
    dissipation = np.sum(
        eigenvectors * (absolute_eigenvalues * strengths)[None, :, :],
        axis=1,
    )
    return 0.5 * (
        convective_flux_code(left) + convective_flux_code(right) - dissipation
    )
