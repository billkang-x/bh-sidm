"""Slope-limited spatial reconstruction for spherical finite volumes."""

from __future__ import annotations

import numpy as np

from .mesh import SphericalGrid
from .state import FluidState


def minmod(*values: np.ndarray) -> np.ndarray:
    """Return the elementwise minmod of two or more arrays."""

    if len(values) < 2:
        raise ValueError("minmod requires at least two arguments")
    arrays = [np.asarray(value, dtype=float) for value in values]
    shape = arrays[0].shape
    if any(array.shape != shape for array in arrays[1:]):
        raise ValueError("minmod arguments must have matching shapes")
    stacked = np.stack(arrays)
    all_positive = np.all(stacked > 0.0, axis=0)
    all_negative = np.all(stacked < 0.0, axis=0)
    return np.where(
        all_positive,
        np.min(stacked, axis=0),
        np.where(all_negative, np.max(stacked, axis=0), 0.0),
    )


def mc_slopes(values: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Return nonuniform-grid monotonized-central slopes.

    Boundary slopes are zero to implement the zero-gradient boundary
    prescription used by the baseline paper.
    """

    values_array = np.asarray(values, dtype=float)
    centers_array = np.asarray(centers, dtype=float)
    if values_array.ndim != 2:
        raise ValueError("values must have shape (num_variables, num_cells)")
    if centers_array.ndim != 1 or values_array.shape[1] != len(centers_array):
        raise ValueError("centers length must match the number of cells")
    if len(centers_array) < 2 or np.any(np.diff(centers_array) <= 0.0):
        raise ValueError("centers must contain at least two increasing values")

    slopes = np.zeros_like(values_array)
    if len(centers_array) == 2:
        return slopes
    backward = (values_array[:, 1:-1] - values_array[:, :-2]) / (
        centers_array[1:-1] - centers_array[:-2]
    )
    forward = (values_array[:, 2:] - values_array[:, 1:-1]) / (
        centers_array[2:] - centers_array[1:-1]
    )
    centered = (values_array[:, 2:] - values_array[:, :-2]) / (
        centers_array[2:] - centers_array[:-2]
    )
    slopes[:, 1:-1] = minmod(centered, 2.0 * backward, 2.0 * forward)
    return slopes


def mc_reconstruct_primitive(
    state: FluidState,
    grid: SphericalGrid,
) -> tuple[FluidState, FluidState]:
    """Reconstruct ``rho``, ``u``, and ``p`` at every grid interface."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    primitive = np.vstack(
        [state.density, state.radial_velocity, state.pressure]
    )
    slopes = mc_slopes(primitive, grid.centers_code)
    num_interfaces = grid.num_cells + 1
    left_values = np.empty((3, num_interfaces), dtype=float)
    right_values = np.empty((3, num_interfaces), dtype=float)

    left_values[:, 0] = primitive[:, 0]
    right_values[:, 0] = primitive[:, 0]
    left_values[:, -1] = primitive[:, -1]
    right_values[:, -1] = primitive[:, -1]
    if grid.num_cells > 1:
        interfaces = grid.interfaces_code[1:-1]
        left_values[:, 1:-1] = primitive[:, :-1] + slopes[:, :-1] * (
            interfaces - grid.centers_code[:-1]
        )
        right_values[:, 1:-1] = primitive[:, 1:] + slopes[:, 1:] * (
            interfaces - grid.centers_code[1:]
        )

    if np.any(left_values[0] <= 0.0) or np.any(right_values[0] <= 0.0):
        raise ValueError("MC reconstruction produced non-positive density")
    if np.any(left_values[2] <= 0.0) or np.any(right_values[2] <= 0.0):
        raise ValueError("MC reconstruction produced non-positive pressure")
    left = FluidState.from_pressure(
        left_values[0],
        left_values[1],
        left_values[2],
    )
    right = FluidState.from_pressure(
        right_values[0],
        right_values[1],
        right_values[2],
    )
    return left, right
