"""First-order hyperbolic time integration for the 1D SIDM equations."""

from __future__ import annotations

import numpy as np

from .fluxes import roe_flux_code, rusanov_flux_code, sound_speed_code
from .mesh import SphericalGrid
from .reconstruction import mc_reconstruct_primitive
from .sources import spherical_source_terms_code, total_enclosed_mass_code
from .state import FluidState


def cfl_timestep_code(
    state: FluidState,
    grid: SphericalGrid,
    cfl_number: float = 0.4,
) -> float:
    """Return a conservative explicit timestep in dimensionless units."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if not 0.0 < cfl_number <= 1.0:
        raise ValueError("cfl_number must lie in (0, 1]")

    signal_speed = np.abs(state.radial_velocity) + sound_speed_code(state)
    return float(cfl_number * np.min(grid.widths_code / signal_speed))


def gravity_timestep_code(
    state: FluidState,
    grid: SphericalGrid,
    total_mass_code: np.ndarray,
    cfl_number: float = 0.4,
) -> float:
    """Limit the fractional velocity change from explicit gravity."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if not 0.0 < cfl_number <= 1.0:
        raise ValueError("cfl_number must lie in (0, 1]")
    total_mass = np.asarray(total_mass_code, dtype=float)
    if total_mass.shape != state.density.shape:
        raise ValueError("total_mass_code shape must match state")
    if np.any(total_mass < 0.0):
        raise ValueError("total_mass_code cannot be negative")

    acceleration = total_mass / grid.centers_code**2
    active = acceleration > 0.0
    if not np.any(active):
        return float("inf")
    characteristic_speed = sound_speed_code(state)
    return float(
        cfl_number * np.min(characteristic_speed[active] / acceleration[active])
    )


def enclosed_dark_matter_mass_code(
    state: FluidState,
    grid: SphericalGrid,
    inner_enclosed_mass_code: float = 0.0,
) -> np.ndarray:
    """Return cell-center enclosed SIDM mass from finite-volume averages.

    Dimensionless halo mass obeys dM/dr = r^2 rho. Mass inside the excised
    region may be supplied separately through ``inner_enclosed_mass_code``.
    """

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if inner_enclosed_mass_code < 0.0:
        raise ValueError("inner_enclosed_mass_code cannot be negative")

    shell_mass = state.density * grid.cell_volumes_code
    mass_at_inner_edge = inner_enclosed_mass_code + np.concatenate(
        ([0.0], np.cumsum(shell_mass[:-1]))
    )
    volume_to_center = (
        grid.centers_code**3 - grid.inner_edges_code**3
    ) / 3.0
    return mass_at_inner_edge + state.density * volume_to_center


def stable_timestep_code(
    state: FluidState,
    grid: SphericalGrid,
    cfl_number: float = 0.4,
    black_hole_mass_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
    dark_matter_enclosed_mass_code: np.ndarray | None = None,
    inner_dark_matter_mass_code: float = 0.0,
    include_gravity_limit: bool = True,
) -> float:
    """Return the minimum hyperbolic and explicit-gravity timestep."""

    if dark_matter_enclosed_mass_code is None:
        dark_mass = enclosed_dark_matter_mass_code(
            state,
            grid,
            inner_enclosed_mass_code=inner_dark_matter_mass_code,
        )
    else:
        dark_mass = np.asarray(dark_matter_enclosed_mass_code, dtype=float)
    total_mass = total_enclosed_mass_code(
        dark_mass,
        black_hole_mass_code=black_hole_mass_code,
        baryon_enclosed_mass_code=baryon_enclosed_mass_code,
    )
    hyperbolic_dt = cfl_timestep_code(state, grid, cfl_number)
    if not include_gravity_limit:
        return hyperbolic_dt
    return min(
        hyperbolic_dt,
        gravity_timestep_code(state, grid, total_mass, cfl_number),
    )


def zero_gradient_interface_states(state: FluidState) -> tuple[FluidState, FluidState]:
    """Build piecewise-constant interface states with transmissive edges."""

    def extend(values: np.ndarray, side: str) -> np.ndarray:
        if side == "left":
            return np.concatenate(([values[0]], values))
        return np.concatenate((values, [values[-1]]))

    left = FluidState(
        density=extend(state.density, "left"),
        radial_velocity=extend(state.radial_velocity, "left"),
        velocity_dispersion=extend(state.velocity_dispersion, "left"),
    )
    right = FluidState(
        density=extend(state.density, "right"),
        radial_velocity=extend(state.radial_velocity, "right"),
        velocity_dispersion=extend(state.velocity_dispersion, "right"),
    )
    return left, right


def reconstructed_interface_states(
    state: FluidState,
    grid: SphericalGrid,
    reconstruction: str = "constant",
) -> tuple[FluidState, FluidState]:
    """Return interface states for the selected reconstruction method."""

    if reconstruction == "constant":
        return zero_gradient_interface_states(state)
    if reconstruction == "mc":
        return mc_reconstruct_primitive(state, grid)
    raise ValueError("reconstruction must be 'constant' or 'mc'")


def interface_flux_code(
    left: FluidState,
    right: FluidState,
    riemann_solver: str = "rusanov",
    entropy_fix: float = 0.1,
) -> np.ndarray:
    """Return interface fluxes for the selected approximate solver."""

    if riemann_solver == "rusanov":
        return rusanov_flux_code(left, right)
    if riemann_solver == "roe":
        return roe_flux_code(left, right, entropy_fix=entropy_fix)
    raise ValueError("riemann_solver must be 'rusanov' or 'roe'")


def apply_gravity_kick(
    state: FluidState,
    acceleration_code: np.ndarray,
    dt_code: float,
) -> FluidState:
    """Apply a frozen gravitational acceleration without changing heat."""

    acceleration = np.asarray(acceleration_code, dtype=float)
    if acceleration.shape != state.density.shape:
        raise ValueError("acceleration_code shape must match state")
    if dt_code <= 0.0:
        raise ValueError("dt_code must be positive")
    velocity_increment = acceleration * dt_code
    kicked = state.conservative.copy()
    kicked[1] += state.density * velocity_increment
    kicked[2] += (
        state.density * state.radial_velocity * velocity_increment
        + 0.5 * state.density * velocity_increment**2
    )
    return FluidState.from_conservative(kicked)


def advance_hyperbolic_step(
    state: FluidState,
    grid: SphericalGrid,
    dt_code: float,
    black_hole_mass_code: float = 0.0,
    baryon_enclosed_mass_code: np.ndarray | None = None,
    dark_matter_enclosed_mass_code: np.ndarray | None = None,
    inner_dark_matter_mass_code: float = 0.0,
    reconstruction: str = "constant",
    riemann_solver: str = "rusanov",
    entropy_fix: float = 0.1,
    positivity_fallback: bool = True,
    source_integration: str = "euler",
) -> FluidState:
    """Advance Eq. (14) by one first-order explicit Euler step.

    The update includes spherical geometry and gravitational source terms,
    but not thermal conduction. Boundary states use the paper's zero-gradient
    prescription. A non-positive state raises ``ValueError`` through
    ``FluidState.from_conservative`` instead of silently applying floors.
    """

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    if dt_code <= 0.0:
        raise ValueError("dt_code must be positive")

    if dark_matter_enclosed_mass_code is None:
        dark_mass = enclosed_dark_matter_mass_code(
            state,
            grid,
            inner_enclosed_mass_code=inner_dark_matter_mass_code,
        )
    else:
        dark_mass = np.asarray(dark_matter_enclosed_mass_code, dtype=float)
        if dark_mass.shape != state.density.shape:
            raise ValueError("dark_matter_enclosed_mass_code shape must match state")

    total_mass = total_enclosed_mass_code(
        dark_mass,
        black_hole_mass_code=black_hole_mass_code,
        baryon_enclosed_mass_code=baryon_enclosed_mass_code,
    )
    source = spherical_source_terms_code(
        state,
        grid,
        dark_mass,
        black_hole_mass_code=black_hole_mass_code,
        baryon_enclosed_mass_code=baryon_enclosed_mass_code,
    )
    left, right = reconstructed_interface_states(state, grid, reconstruction)

    if source_integration not in ("euler", "gravity_kick"):
        raise ValueError("source_integration must be 'euler' or 'gravity_kick'")

    def update_with_flux(flux: np.ndarray) -> FluidState:
        area_weighted_flux = flux * grid.interface_areas_code
        flux_divergence = (
            area_weighted_flux[:, 1:] - area_weighted_flux[:, :-1]
        ) / grid.cell_volumes_code
        if source_integration == "euler":
            updated = state.conservative + dt_code * (source - flux_divergence)
            return FluidState.from_conservative(updated)

        geometric_source = np.zeros_like(source)
        geometric_source[1] = 2.0 * state.pressure / grid.centers_code
        intermediate_conservative = state.conservative + dt_code * (
            geometric_source - flux_divergence
        )
        intermediate = FluidState.from_conservative(intermediate_conservative)
        acceleration = -total_mass / grid.centers_code**2
        return apply_gravity_kick(intermediate, acceleration, dt_code)

    primary_flux = interface_flux_code(
        left,
        right,
        riemann_solver=riemann_solver,
        entropy_fix=entropy_fix,
    )
    try:
        return update_with_flux(primary_flux)
    except ValueError:
        if not positivity_fallback or riemann_solver != "roe":
            raise
        return update_with_flux(rusanov_flux_code(left, right))
