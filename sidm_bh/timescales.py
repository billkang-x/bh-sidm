"""Local physical timescales for stage-3 SIDM flow diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .conduction import cell_conductivity_code
from .mesh import SphericalGrid
from .sidm import CONDUCTIVITY_A
from .solver import enclosed_dark_matter_mass_code
from .sources import total_enclosed_mass_code
from .state import FluidState


@dataclass(frozen=True)
class TimescaleProfiles:
    dynamical_code: np.ndarray
    collision_code: np.ndarray
    conduction_radius_code: np.ndarray
    conduction_gradient_code: np.ndarray
    inflow_code: np.ndarray
    thermal_length_code: np.ndarray
    knudsen_number: np.ndarray


def local_timescale_profiles_code(
    state: FluidState,
    grid: SphericalGrid,
    sigma_over_m_code: float,
    black_hole_mass_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
    inner_dark_matter_mass_code: float = 0.0,
) -> TimescaleProfiles:
    """Return timescales using the same transport closure as the solver.

    The collision time is ``1 / (a rho sigma v)``. Conductive diffusion of
    ``v^2`` obeys ``t_cond = 3 rho L^2 / (2 kappa)``. The gradient-scale
    estimate limits L to one cell width through the local radius, while the
    radius-scale estimate always uses L = r.
    """

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if sigma_over_m_code <= 0.0:
        raise ValueError("sigma_over_m_code must be positive")
    if black_hole_mass_code < 0.0:
        raise ValueError("black_hole_mass_code cannot be negative")

    dark_mass = enclosed_dark_matter_mass_code(
        state,
        grid,
        inner_enclosed_mass_code=inner_dark_matter_mass_code,
    )
    total_mass = total_enclosed_mass_code(
        dark_mass,
        black_hole_mass_code=black_hole_mass_code,
        baryon_enclosed_mass_code=baryon_enclosed_mass_code,
    )
    radius = grid.centers_code
    dynamical = np.sqrt(radius**3 / total_mass)
    collision = 1.0 / (
        CONDUCTIVITY_A
        * state.density
        * sigma_over_m_code
        * state.velocity_dispersion
    )

    conductivity = cell_conductivity_code(state, sigma_over_m_code)
    temperature = state.velocity_dispersion**2
    edge_order = 2 if grid.num_cells >= 3 else 1
    temperature_gradient = np.gradient(
        temperature,
        radius,
        edge_order=edge_order,
    )
    raw_thermal_length = np.full_like(radius, np.inf)
    nonzero_gradient = np.abs(temperature_gradient) > 0.0
    raw_thermal_length[nonzero_gradient] = (
        temperature[nonzero_gradient]
        / np.abs(temperature_gradient[nonzero_gradient])
    )
    thermal_length = np.clip(raw_thermal_length, grid.widths_code, radius)
    conduction_radius = 1.5 * state.density * radius**2 / conductivity
    conduction_gradient = (
        1.5 * state.density * thermal_length**2 / conductivity
    )

    inflow = np.full_like(radius, np.inf)
    inward = state.radial_velocity < 0.0
    inflow[inward] = radius[inward] / -state.radial_velocity[inward]
    knudsen = 1.0 / (
        sigma_over_m_code
        * state.velocity_dispersion
        * np.sqrt(state.density)
    )
    return TimescaleProfiles(
        dynamical_code=dynamical,
        collision_code=collision,
        conduction_radius_code=conduction_radius,
        conduction_gradient_code=conduction_gradient,
        inflow_code=inflow,
        thermal_length_code=thermal_length,
        knudsen_number=knudsen,
    )


def inward_flux_median_radius_code(
    state: FluidState,
    grid: SphericalGrid,
    maximum_radius_code: float | None = None,
) -> float:
    """Return the log-radius median of inward mass flux inside a radius."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if maximum_radius_code is not None and maximum_radius_code <= 0.0:
        raise ValueError("maximum_radius_code must be positive when supplied")
    selected = np.ones(grid.num_cells, dtype=bool)
    if maximum_radius_code is not None:
        selected &= grid.centers_code <= maximum_radius_code
    inward_flux = np.maximum(
        -grid.centers_code**2 * state.density * state.radial_velocity,
        0.0,
    )
    log_width = np.log(grid.outer_edges_code / grid.inner_edges_code)
    weights = np.where(selected, inward_flux * log_width, 0.0)
    total = float(np.sum(weights))
    if total <= 0.0:
        return float("nan")
    index = int(np.searchsorted(np.cumsum(weights), 0.5 * total))
    return float(grid.centers_code[index])
