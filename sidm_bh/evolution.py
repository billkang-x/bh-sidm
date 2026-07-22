"""Time-evolution driver and conservation diagnostics for SIDM accretion."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .accretion import inner_boundary_sidm_accretion_rate_code
from .baryons import smoothstep_mass_fraction
from .conduction import implicit_conduction_step
from .mesh import SphericalGrid
from .solver import advance_hyperbolic_step, stable_timestep_code
from .state import FluidState


def fluid_mass_code(state: FluidState, grid: SphericalGrid) -> float:
    """Return SIDM mass currently resolved on the finite-volume grid."""

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    return float(np.sum(state.density * grid.cell_volumes_code))


def boundary_mass_fluxes_code(
    state: FluidState,
    grid: SphericalGrid,
) -> tuple[float, float]:
    """Return signed outward-oriented mass fluxes at inner and outer edges.

    At either boundary the flux is ``A rho u`` with radial velocity positive
    toward increasing radius. Thus a negative inner flux feeds the black hole,
    while a positive outer flux removes mass from the computational domain.
    """

    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    inner_flux = (
        grid.interface_areas_code[0]
        * state.density[0]
        * state.radial_velocity[0]
    )
    outer_flux = (
        grid.interface_areas_code[-1]
        * state.density[-1]
        * state.radial_velocity[-1]
    )
    return float(inner_flux), float(outer_flux)


@dataclass(frozen=True)
class EvolutionHistory:
    """Scalar histories recorded at every accepted timestep."""

    times_code: np.ndarray
    timesteps_code: np.ndarray
    black_hole_masses_code: np.ndarray
    sidm_accretion_rates_code: np.ndarray
    fluid_masses_code: np.ndarray
    cumulative_outer_flux_code: np.ndarray
    cumulative_inner_inflow_code: np.ndarray
    mass_budget_residuals_code: np.ndarray

    @property
    def num_steps(self) -> int:
        return len(self.timesteps_code)

    @property
    def max_absolute_mass_budget_residual_code(self) -> float:
        return float(np.max(np.abs(self.mass_budget_residuals_code)))


@dataclass(frozen=True)
class EvolutionResult:
    """Final fluid state, black-hole mass, and scalar evolution history."""

    final_state: FluidState
    final_black_hole_mass_code: float
    history: EvolutionHistory


def evolve_sidm(
    initial_state: FluidState,
    grid: SphericalGrid,
    end_time_code: float,
    initial_black_hole_mass_code: float = 0.0,
    sigma_over_m_code: float | None = None,
    cfl_number: float = 0.2,
    baryon_enclosed_mass_code: np.ndarray | None = None,
    baryon_assembly_time_code: float | None = None,
    inner_dark_matter_mass_code: float = 0.0,
    max_steps: int = 1_000_000,
    reconstruction: str = "constant",
    riemann_solver: str = "rusanov",
    entropy_fix: float = 0.1,
    positivity_fallback: bool = True,
    source_integration: str = "euler",
) -> EvolutionResult:
    """Evolve the coupled SIDM fluid and central black-hole mass.

    Each accepted step applies the selected reconstruction/Riemann solver,
    followed by the optional implicit conduction update. The black hole gains
    exactly the inward SIDM mass flux used by the zero-gradient inner boundary
    state.
    """

    if len(initial_state.density) != grid.num_cells:
        raise ValueError("initial_state length must match grid")
    if end_time_code <= 0.0:
        raise ValueError("end_time_code must be positive")
    if initial_black_hole_mass_code < 0.0:
        raise ValueError("initial_black_hole_mass_code cannot be negative")
    if sigma_over_m_code is not None and sigma_over_m_code <= 0.0:
        raise ValueError("sigma_over_m_code must be positive when conduction is enabled")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive")
    if baryon_assembly_time_code is not None and baryon_assembly_time_code <= 0.0:
        raise ValueError("baryon_assembly_time_code must be positive when supplied")

    full_baryon_mass = (
        None
        if baryon_enclosed_mass_code is None
        else np.asarray(baryon_enclosed_mass_code, dtype=float)
    )
    if full_baryon_mass is not None and full_baryon_mass.shape != (grid.num_cells,):
        raise ValueError("baryon_enclosed_mass_code shape must match grid")

    state = initial_state
    black_hole_mass = float(initial_black_hole_mass_code)
    time = 0.0
    cumulative_outer_flux = 0.0
    cumulative_inner_inflow = 0.0
    initial_total_mass = fluid_mass_code(state, grid) + black_hole_mass

    times = [time]
    timesteps: list[float] = []
    black_hole_masses = [black_hole_mass]
    accretion_rates = [inner_boundary_sidm_accretion_rate_code(state, grid)]
    fluid_masses = [fluid_mass_code(state, grid)]
    cumulative_outer_fluxes = [cumulative_outer_flux]
    cumulative_inner_inflows = [cumulative_inner_inflow]
    mass_budget_residuals = [0.0]

    for _ in range(max_steps):
        if time >= end_time_code:
            break
        if full_baryon_mass is None:
            current_baryon_mass = None
        elif baryon_assembly_time_code is None:
            current_baryon_mass = full_baryon_mass
        else:
            current_baryon_mass = full_baryon_mass * smoothstep_mass_fraction(
                time,
                baryon_assembly_time_code,
            )
        dt = min(
            stable_timestep_code(
                state,
                grid,
                cfl_number=cfl_number,
                black_hole_mass_code=black_hole_mass,
                baryon_enclosed_mass_code=current_baryon_mass,
                inner_dark_matter_mass_code=inner_dark_matter_mass_code,
                include_gravity_limit=source_integration != "gravity_kick",
            ),
            end_time_code - time,
        )
        if not np.isfinite(dt) or dt <= 0.0:
            raise RuntimeError("time integration produced a non-positive timestep")

        accretion_rate = inner_boundary_sidm_accretion_rate_code(state, grid)
        inner_flux, outer_flux = boundary_mass_fluxes_code(state, grid)
        next_state = advance_hyperbolic_step(
            state,
            grid,
            dt,
            black_hole_mass_code=black_hole_mass,
            baryon_enclosed_mass_code=current_baryon_mass,
            inner_dark_matter_mass_code=inner_dark_matter_mass_code,
            reconstruction=reconstruction,
            riemann_solver=riemann_solver,
            entropy_fix=entropy_fix,
            positivity_fallback=positivity_fallback,
            source_integration=source_integration,
        )
        if sigma_over_m_code is not None:
            next_state = implicit_conduction_step(
                next_state,
                grid,
                dt,
                sigma_over_m_code,
            )

        black_hole_mass += accretion_rate * dt
        cumulative_outer_flux += outer_flux * dt
        cumulative_inner_inflow += max(inner_flux, 0.0) * dt
        time += dt
        state = next_state

        current_fluid_mass = fluid_mass_code(state, grid)
        budget_residual = (
            current_fluid_mass
            + black_hole_mass
            + cumulative_outer_flux
            - cumulative_inner_inflow
            - initial_total_mass
        )
        timesteps.append(dt)
        times.append(time)
        black_hole_masses.append(black_hole_mass)
        accretion_rates.append(inner_boundary_sidm_accretion_rate_code(state, grid))
        fluid_masses.append(current_fluid_mass)
        cumulative_outer_fluxes.append(cumulative_outer_flux)
        cumulative_inner_inflows.append(cumulative_inner_inflow)
        mass_budget_residuals.append(budget_residual)
    if time < end_time_code:
        raise RuntimeError(
            f"evolution did not reach end_time_code={end_time_code} "
            f"within max_steps={max_steps}"
        )

    history = EvolutionHistory(
        times_code=np.asarray(times),
        timesteps_code=np.asarray(timesteps),
        black_hole_masses_code=np.asarray(black_hole_masses),
        sidm_accretion_rates_code=np.asarray(accretion_rates),
        fluid_masses_code=np.asarray(fluid_masses),
        cumulative_outer_flux_code=np.asarray(cumulative_outer_fluxes),
        cumulative_inner_inflow_code=np.asarray(cumulative_inner_inflows),
        mass_budget_residuals_code=np.asarray(mass_budget_residuals),
    )
    return EvolutionResult(
        final_state=state,
        final_black_hole_mass_code=black_hole_mass,
        history=history,
    )
