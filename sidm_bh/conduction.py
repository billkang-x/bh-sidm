"""Implicit SIDM heat conduction from Appendix B of the baseline paper."""

from __future__ import annotations

import numpy as np

from .mesh import SphericalGrid
from .sidm import (
    CONDUCTIVITY_A,
    CONDUCTIVITY_B,
    CONDUCTIVITY_C,
)
from .state import FluidState


def solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    right_hand_side: np.ndarray,
) -> np.ndarray:
    """Solve a tridiagonal system with the Thomas algorithm."""

    lower_array = np.asarray(lower, dtype=float)
    diagonal_array = np.asarray(diagonal, dtype=float)
    upper_array = np.asarray(upper, dtype=float)
    rhs = np.asarray(right_hand_side, dtype=float)
    size = len(diagonal_array)
    if diagonal_array.ndim != 1 or rhs.shape != diagonal_array.shape:
        raise ValueError("diagonal and right_hand_side must be matching 1D arrays")
    if lower_array.shape != (max(size - 1, 0),):
        raise ValueError("lower must have length n - 1")
    if upper_array.shape != (max(size - 1, 0),):
        raise ValueError("upper must have length n - 1")
    if size == 0:
        raise ValueError("tridiagonal system cannot be empty")

    modified_diagonal = diagonal_array.copy()
    modified_rhs = rhs.copy()
    if modified_diagonal[0] == 0.0:
        raise ValueError("tridiagonal system has a zero pivot")
    for i in range(1, size):
        multiplier = lower_array[i - 1] / modified_diagonal[i - 1]
        modified_diagonal[i] -= multiplier * upper_array[i - 1]
        modified_rhs[i] -= multiplier * modified_rhs[i - 1]
        if modified_diagonal[i] == 0.0:
            raise ValueError("tridiagonal system has a zero pivot")

    solution = np.empty(size, dtype=float)
    solution[-1] = modified_rhs[-1] / modified_diagonal[-1]
    for i in range(size - 2, -1, -1):
        solution[i] = (
            modified_rhs[i] - upper_array[i] * solution[i + 1]
        ) / modified_diagonal[i]
    return solution


def cell_conductivity_code(
    state: FluidState,
    sigma_over_m_code: float | np.ndarray,
    calibration_c: float = CONDUCTIVITY_C,
) -> np.ndarray:
    """Evaluate the frozen Eq. (9) conductivity at cell centers."""

    cross_section = np.asarray(sigma_over_m_code, dtype=float)
    if cross_section.ndim > 1 or (
        cross_section.ndim == 1 and cross_section.shape != state.density.shape
    ):
        raise ValueError("sigma_over_m_code must be scalar or match the state")
    if np.any(cross_section <= 0.0):
        raise ValueError("sigma_over_m_code must be positive")
    if calibration_c <= 0.0:
        raise ValueError("calibration_c must be positive")
    smfp_inverse = cross_section / (
        CONDUCTIVITY_B * state.velocity_dispersion
    )
    lmfp_inverse = 1.0 / (
        CONDUCTIVITY_A
        * calibration_c
        * state.density
        * state.velocity_dispersion**3
        * cross_section
    )
    return 1.5 / (smfp_inverse + lmfp_inverse)


def interface_conductivity_code(cell_conductivity: np.ndarray) -> np.ndarray:
    """Interpolate conductivity to interfaces with arithmetic averaging."""

    cell_values = np.asarray(cell_conductivity, dtype=float)
    if cell_values.ndim != 1 or len(cell_values) == 0:
        raise ValueError("cell_conductivity must be a non-empty 1D array")
    if np.any(cell_values <= 0.0):
        raise ValueError("cell_conductivity values must be positive")

    interfaces = np.empty(len(cell_values) + 1, dtype=float)
    interfaces[0] = cell_values[0]
    interfaces[-1] = cell_values[-1]
    interfaces[1:-1] = 0.5 * (cell_values[:-1] + cell_values[1:])
    return interfaces


def conduction_system_code(
    state: FluidState,
    grid: SphericalGrid,
    dt_code: float,
    conductivity_at_interfaces_code: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble the tridiagonal coefficients in Eqs. (B5a)-(B5d)."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if dt_code <= 0.0:
        raise ValueError("dt_code must be positive")
    interface_kappa = np.asarray(conductivity_at_interfaces_code, dtype=float)
    if interface_kappa.shape != (grid.num_cells + 1,):
        raise ValueError("interface conductivity must have length num_cells + 1")
    if np.any(interface_kappa <= 0.0):
        raise ValueError("interface conductivity values must be positive")

    radius = grid.centers_code
    factor = 2.0 * dt_code / (3.0 * grid.cell_volumes_code)
    lower = -factor[1:] * (
        grid.interface_areas_code[1:-1]
        * interface_kappa[1:-1]
        / (radius[1:] - radius[:-1])
    )
    upper = -factor[:-1] * (
        grid.interface_areas_code[1:-1]
        * interface_kappa[1:-1]
        / (radius[1:] - radius[:-1])
    )
    diagonal = state.density.copy()
    diagonal[1:] -= lower
    diagonal[:-1] -= upper
    rhs = state.density * state.velocity_dispersion**2
    return lower, diagonal, upper, rhs


def implicit_conduction_step(
    state: FluidState,
    grid: SphericalGrid,
    dt_code: float,
    sigma_over_m_code: float,
    calibration_c: float = CONDUCTIVITY_C,
) -> FluidState:
    """Advance only thermal conduction using frozen coefficients.

    Density and radial velocity remain fixed. Zero-gradient boundary
    conditions make the conductive flux vanish at both domain edges.
    """

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    cell_kappa = cell_conductivity_code(
        state,
        sigma_over_m_code,
        calibration_c=calibration_c,
    )
    interface_kappa = interface_conductivity_code(cell_kappa)
    lower, diagonal, upper, rhs = conduction_system_code(
        state,
        grid,
        dt_code,
        interface_kappa,
    )
    dispersion_squared = solve_tridiagonal(lower, diagonal, upper, rhs)
    if np.any(dispersion_squared <= 0.0):
        raise ValueError("conduction update produced non-positive temperature")
    return FluidState(
        density=state.density,
        radial_velocity=state.radial_velocity,
        velocity_dispersion=np.sqrt(dispersion_squared),
    )
