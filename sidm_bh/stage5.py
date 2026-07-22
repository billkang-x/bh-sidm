"""Stage-5 observable targets and analytic growth-budget filters."""

from __future__ import annotations

import math

from .baryons import eddington_accretion_rate_msun_per_myr


LRD_TARGET_MASSES_MSUN = (1.0e5, 1.0e6, 1.0e7)


def eddington_efolding_time_myr(
    radiative_efficiency: float = 0.1,
    eddington_ratio: float = 1.0,
    duty_cycle: float = 1.0,
) -> float:
    """Return the retained black-hole mass e-folding time."""

    if not 0.0 < radiative_efficiency < 1.0:
        raise ValueError("radiative_efficiency must lie in (0, 1)")
    if eddington_ratio < 0.0:
        raise ValueError("eddington_ratio cannot be negative")
    if not 0.0 <= duty_cycle <= 1.0:
        raise ValueError("duty_cycle must lie in [0, 1]")
    activity = eddington_ratio * duty_cycle
    if activity == 0.0:
        return float("inf")
    retained_rate_per_mass = (
        (1.0 - radiative_efficiency)
        * activity
        * eddington_accretion_rate_msun_per_myr(
            1.0,
            radiative_efficiency,
        )
    )
    return 1.0 / retained_rate_per_mass


def growth_time_to_target_myr(
    seed_mass_msun: float,
    target_mass_msun: float,
    radiative_efficiency: float = 0.1,
    eddington_ratio: float = 1.0,
    duty_cycle: float = 1.0,
) -> float:
    if seed_mass_msun <= 0.0 or target_mass_msun <= 0.0:
        raise ValueError("seed and target masses must be positive")
    if target_mass_msun <= seed_mass_msun:
        return 0.0
    return eddington_efolding_time_myr(
        radiative_efficiency,
        eddington_ratio,
        duty_cycle,
    ) * math.log(target_mass_msun / seed_mass_msun)


def required_eddington_activity(
    seed_mass_msun: float,
    target_mass_msun: float,
    available_time_myr: float,
    radiative_efficiency: float = 0.1,
) -> float:
    """Return the minimum product `(Eddington ratio) * duty cycle`."""

    if available_time_myr <= 0.0:
        raise ValueError("available_time_myr must be positive")
    unit_activity_time = growth_time_to_target_myr(
        seed_mass_msun,
        target_mass_msun,
        radiative_efficiency=radiative_efficiency,
    )
    return unit_activity_time / available_time_myr
