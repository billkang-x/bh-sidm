"""Matched initial conditions for static-baryon stage-3 experiments."""

from __future__ import annotations

import numpy as np

from .baryons import HernquistBaryons
from .initial_conditions import hydrostatic_state_from_profile, SphericalMassProfile
from .mesh import SphericalGrid
from .sources import enclosed_baryon_mass_code
from .state import FluidState
from .units import SimulationScales


def static_baryon_equilibrium_state(
    dark_matter_profile: SphericalMassProfile,
    grid: SphericalGrid,
    scales: SimulationScales,
    black_hole_mass_msun: float,
    baryons: HernquistBaryons | None = None,
) -> tuple[FluidState, np.ndarray]:
    """Return a hydrostatic state and the matching fixed baryon mass profile."""

    if black_hole_mass_msun < 0.0:
        raise ValueError("black_hole_mass_msun cannot be negative")
    if baryons is None or baryons.total_mass_msun == 0.0:
        baryon_mass_code = np.zeros(grid.num_cells)
        extra_mass_msun = None
    else:
        baryon_mass_code = enclosed_baryon_mass_code(baryons, grid, scales)
        extra_mass_msun = baryon_mass_code * scales.mass_scale_msun
    state = hydrostatic_state_from_profile(
        dark_matter_profile,
        grid,
        scales,
        point_mass_msun=black_hole_mass_msun,
        extra_enclosed_mass_msun=extra_mass_msun,
    )
    return state, baryon_mass_code
