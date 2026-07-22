"""Initial-condition helpers for SIDM baseline and baryonic-potential runs."""

from __future__ import annotations

from typing import Protocol

import numpy as np

from .mesh import SphericalGrid
from .state import FluidState
from .units import SimulationScales


class SphericalMassProfile(Protocol):
    def density_msun_pc3(self, radius_pc: float) -> float:
        ...

    def enclosed_mass_msun(self, radius_pc: float) -> float:
        ...


def sample_density_code(
    profile: SphericalMassProfile,
    grid: SphericalGrid,
    scales: SimulationScales,
) -> np.ndarray:
    radii_pc = scales.radius_from_code(grid.centers_code)
    densities = np.array([profile.density_msun_pc3(r) for r in radii_pc])
    return densities / scales.density_scale_msun_pc3


def sample_enclosed_mass_code(
    profile: SphericalMassProfile,
    grid: SphericalGrid,
    scales: SimulationScales,
) -> np.ndarray:
    radii_pc = scales.radius_from_code(grid.centers_code)
    masses = np.array([profile.enclosed_mass_msun(r) for r in radii_pc])
    return masses / scales.mass_scale_msun


def hydrostatic_pressure_code(
    grid: SphericalGrid,
    density_code: np.ndarray,
    enclosed_mass_code: np.ndarray,
    point_mass_code: float = 0.0,
    extra_enclosed_mass_code: np.ndarray | None = None,
    outer_pressure_code: float = 0.0,
) -> np.ndarray:
    """Integrate dp/dr = -rho M(<r) / r^2 inward on cell centers.

    outer_pressure_code is interpreted at the outer interface, not at the
    outermost cell center. This keeps the outermost cell pressure positive
    even when the boundary pressure itself is set to zero.
    """

    radius = grid.centers_code
    density = np.asarray(density_code, dtype=float)
    enclosed = np.asarray(enclosed_mass_code, dtype=float)
    if density.shape != radius.shape:
        raise ValueError("density_code shape must match grid centers")
    if enclosed.shape != radius.shape:
        raise ValueError("enclosed_mass_code shape must match grid centers")
    if np.any(density <= 0.0):
        raise ValueError("density_code must be positive")
    if point_mass_code < 0.0:
        raise ValueError("point_mass_code cannot be negative")

    total_mass = enclosed + point_mass_code
    if extra_enclosed_mass_code is not None:
        extra = np.asarray(extra_enclosed_mass_code, dtype=float)
        if extra.shape != radius.shape:
            raise ValueError("extra_enclosed_mass_code shape must match grid centers")
        total_mass = total_mass + extra

    integrand = density * total_mass / radius**2
    pressure = np.empty_like(density)
    pressure[-1] = outer_pressure_code + integrand[-1] * (
        grid.outer_edges_code[-1] - radius[-1]
    )
    for i in range(len(density) - 2, -1, -1):
        dr = radius[i + 1] - radius[i]
        pressure[i] = pressure[i + 1] + 0.5 * (integrand[i] + integrand[i + 1]) * dr
    return pressure


def hydrostatic_state_from_profile(
    profile: SphericalMassProfile,
    grid: SphericalGrid,
    scales: SimulationScales,
    point_mass_msun: float = 0.0,
    extra_enclosed_mass_msun: np.ndarray | None = None,
    outer_pressure_code: float = 0.0,
) -> FluidState:
    density = sample_density_code(profile, grid, scales)
    enclosed_mass = sample_enclosed_mass_code(profile, grid, scales)
    point_mass = scales.mass_to_code(point_mass_msun)
    extra_enclosed = None
    if extra_enclosed_mass_msun is not None:
        extra_enclosed = np.asarray(extra_enclosed_mass_msun, dtype=float) / scales.mass_scale_msun
    pressure = hydrostatic_pressure_code(
        grid,
        density,
        enclosed_mass,
        point_mass_code=point_mass,
        extra_enclosed_mass_code=extra_enclosed,
        outer_pressure_code=outer_pressure_code,
    )
    return FluidState.from_pressure(
        density=density,
        radial_velocity=np.zeros_like(density),
        pressure=pressure,
    )


def isothermal_state_from_profile(
    profile: SphericalMassProfile,
    grid: SphericalGrid,
    scales: SimulationScales,
    velocity_dispersion_km_s: float,
) -> FluidState:
    """Sample a profile with spatially constant physical velocity dispersion."""

    if velocity_dispersion_km_s <= 0.0:
        raise ValueError("velocity_dispersion_km_s must be positive")
    density = sample_density_code(profile, grid, scales)
    dispersion_code = scales.velocity_to_code(velocity_dispersion_km_s)
    return FluidState(
        density=density,
        radial_velocity=np.zeros_like(density),
        velocity_dispersion=np.full_like(density, dispersion_code),
    )
