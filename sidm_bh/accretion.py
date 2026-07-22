"""Black-hole accretion diagnostics and updates."""

from __future__ import annotations

from dataclasses import dataclass

from .mesh import SphericalGrid
from .state import FluidState
from .units import SimulationScales


def sidm_accretion_rate_code(
    inner_radius_code: float,
    density_code: float,
    radial_velocity_code: float,
    inward_only: bool = True,
) -> float:
    """Return dM_BH/dt from the dimensionless inner-boundary SIDM flux.

    The baseline paper's normalization absorbs the common 4 pi factor into
    M0, giving dM_tilde/dt_tilde = -r_min^2 rho_min u_min.
    """

    if inner_radius_code <= 0.0:
        raise ValueError("inner_radius_code must be positive")
    if density_code <= 0.0:
        raise ValueError("density_code must be positive")

    rate = -(inner_radius_code**2) * density_code * radial_velocity_code
    if inward_only:
        return max(rate, 0.0)
    return rate


def inner_boundary_sidm_accretion_rate_code(
    state: FluidState,
    grid: SphericalGrid,
    inward_only: bool = True,
) -> float:
    if len(state.density) != grid.num_cells:
        raise ValueError("state length must match grid")
    return sidm_accretion_rate_code(
        inner_radius_code=grid.interfaces_code[0],
        density_code=float(state.density[0]),
        radial_velocity_code=float(state.radial_velocity[0]),
        inward_only=inward_only,
    )


def accretion_rate_code_to_msun_per_myr(
    rate_code: float,
    scales: SimulationScales,
) -> float:
    return rate_code * scales.mass_scale_msun / scales.time_scale_myr


def accretion_rate_msun_per_myr_to_code(
    rate_msun_per_myr: float,
    scales: SimulationScales,
) -> float:
    return rate_msun_per_myr * scales.time_scale_myr / scales.mass_scale_msun


@dataclass(frozen=True)
class AccretionStep:
    """Mass update from a constant accretion rate over one timestep."""

    initial_mass_code: float
    rate_code: float
    dt_code: float

    def __post_init__(self) -> None:
        if self.initial_mass_code < 0.0:
            raise ValueError("initial_mass_code cannot be negative")
        if self.dt_code < 0.0:
            raise ValueError("dt_code cannot be negative")

    @property
    def delta_mass_code(self) -> float:
        return self.rate_code * self.dt_code

    @property
    def final_mass_code(self) -> float:
        return max(self.initial_mass_code + self.delta_mass_code, 0.0)

