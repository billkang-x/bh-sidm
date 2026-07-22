"""Source terms for the spherical SIDM fluid equations."""

from __future__ import annotations

import numpy as np

from .baryons import HernquistBaryons
from .mesh import SphericalGrid
from .state import FluidState
from .units import SimulationScales


def enclosed_baryon_mass_code(
    baryons: HernquistBaryons,
    grid: SphericalGrid,
    scales: SimulationScales,
    time_myr: float = 0.0,
) -> np.ndarray:
    """Sample baryonic enclosed mass at cell centers in code units."""

    radii_pc = scales.radius_from_code(grid.centers_code)
    masses_cgs = np.array(
        [baryons.enclosed_mass_cgs(radius, time_myr=time_myr) for radius in radii_pc]
    )
    return masses_cgs / scales.mass_scale_cgs


def total_enclosed_mass_code(
    dark_matter_enclosed_mass_code: np.ndarray,
    black_hole_mass_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
) -> np.ndarray:
    """Combine dark matter, black-hole, and optional baryonic masses."""

    dark_matter = np.asarray(dark_matter_enclosed_mass_code, dtype=float)
    if dark_matter.ndim != 1:
        raise ValueError("dark_matter_enclosed_mass_code must be one-dimensional")
    if np.any(dark_matter < 0.0):
        raise ValueError("dark matter enclosed mass cannot be negative")
    if black_hole_mass_code < 0.0:
        raise ValueError("black_hole_mass_code cannot be negative")

    total = dark_matter + black_hole_mass_code
    if baryon_enclosed_mass_code is not None:
        baryons = np.asarray(baryon_enclosed_mass_code, dtype=float)
        if baryons.shape != dark_matter.shape:
            raise ValueError("baryon_enclosed_mass_code shape must match dark matter")
        if np.any(baryons < 0.0):
            raise ValueError("baryon enclosed mass cannot be negative")
        total = total + baryons
    return total


def gravitational_acceleration_code(
    radius_code: np.ndarray,
    total_mass_code: np.ndarray,
) -> np.ndarray:
    """Outward-positive gravitational acceleration in code units."""

    radius = np.asarray(radius_code, dtype=float)
    total_mass = np.asarray(total_mass_code, dtype=float)
    if radius.shape != total_mass.shape:
        raise ValueError("radius_code shape must match total_mass_code")
    if np.any(radius <= 0.0):
        raise ValueError("radius_code values must be positive")
    if np.any(total_mass < 0.0):
        raise ValueError("total_mass_code cannot be negative")
    return -total_mass / radius**2


def spherical_source_terms_code(
    state: FluidState,
    grid: SphericalGrid,
    dark_matter_enclosed_mass_code: np.ndarray,
    black_hole_mass_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
) -> np.ndarray:
    """Return dimensionless source vector for the spherical Euler system.

    The form matches the baseline paper's conservation-law source:

    S = [0,
         -rho M(<r)/r^2 + 2p/r,
         -rho M(<r) u/r^2].
    """

    radius = grid.centers_code
    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")

    total_mass = total_enclosed_mass_code(
        dark_matter_enclosed_mass_code,
        black_hole_mass_code=black_hole_mass_code,
        baryon_enclosed_mass_code=baryon_enclosed_mass_code,
    )
    if total_mass.shape != state.density.shape:
        raise ValueError("enclosed mass shape must match state")

    acceleration = gravitational_acceleration_code(radius, total_mass)
    source = np.zeros((3, grid.num_cells), dtype=float)
    source[1] = state.density * acceleration + 2.0 * state.pressure / radius
    source[2] = state.density * state.radial_velocity * acceleration
    return source

