"""Spherical baryonic mass profiles for the SIDM accretion prototype."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import G_CGS, MYR_CGS, M_SUN_CGS, PC_CGS


def smoothstep_mass_fraction(time: float, assembly_time: float) -> float:
    """Return a finite-time, zero-slope baryonic assembly fraction."""

    if assembly_time <= 0.0:
        raise ValueError("assembly_time must be positive")
    if time <= 0.0:
        return 0.0
    if time >= assembly_time:
        return 1.0
    phase = time / assembly_time
    return phase * phase * (3.0 - 2.0 * phase)


@dataclass(frozen=True)
class HernquistBaryons:
    """Hernquist baryonic component with optional exponential mass growth.

    Parameters are expressed in astrophysical units for readability:
    total_mass_msun is the asymptotic baryonic mass, scale_radius_pc is the
    Hernquist scale radius, and growth_time_myr controls mass assembly.
    If growth_time_myr is None or non-positive, the profile is static.
    """

    total_mass_msun: float
    scale_radius_pc: float
    growth_time_myr: float | None = None

    def mass_fraction(self, time_myr: float = 0.0) -> float:
        if self.growth_time_myr is None or self.growth_time_myr <= 0.0:
            return 1.0
        if time_myr <= 0.0:
            return 0.0
        return 1.0 - math.exp(-time_myr / self.growth_time_myr)

    def total_mass_cgs(self, time_myr: float = 0.0) -> float:
        return self.total_mass_msun * M_SUN_CGS * self.mass_fraction(time_myr)

    def enclosed_mass_cgs(self, radius_pc: float, time_myr: float = 0.0) -> float:
        if radius_pc <= 0.0:
            return 0.0
        radius_cm = radius_pc * PC_CGS
        scale_cm = self.scale_radius_pc * PC_CGS
        return self.total_mass_cgs(time_myr) * radius_cm**2 / (radius_cm + scale_cm) ** 2

    def density_cgs(self, radius_pc: float, time_myr: float = 0.0) -> float:
        if radius_pc <= 0.0:
            return math.inf
        radius_cm = radius_pc * PC_CGS
        scale_cm = self.scale_radius_pc * PC_CGS
        mass_cgs = self.total_mass_cgs(time_myr)
        return mass_cgs * scale_cm / (2.0 * math.pi * radius_cm * (radius_cm + scale_cm) ** 3)

    def potential_cgs(self, radius_pc: float, time_myr: float = 0.0) -> float:
        radius_cm = max(radius_pc, 0.0) * PC_CGS
        scale_cm = self.scale_radius_pc * PC_CGS
        return -G_CGS * self.total_mass_cgs(time_myr) / (radius_cm + scale_cm)

    def acceleration_cgs(self, radius_pc: float, time_myr: float = 0.0) -> float:
        """Return inward radial acceleration magnitude in cm/s^2."""

        if radius_pc <= 0.0:
            return 0.0
        radius_cm = radius_pc * PC_CGS
        enclosed = self.enclosed_mass_cgs(radius_pc, time_myr)
        return G_CGS * enclosed / radius_cm**2


def eddington_accretion_rate_msun_per_myr(
    black_hole_mass_msun: float,
    radiative_efficiency: float = 0.1,
) -> float:
    """Eddington-limited baryonic accretion rate in Msun/Myr."""

    if black_hole_mass_msun <= 0.0:
        return 0.0
    if radiative_efficiency <= 0.0:
        raise ValueError("radiative_efficiency must be positive")

    from .constants import C_CGS, SIGMA_T_CGS, M_PROTON_CGS

    mdot_g_s = (
        4.0
        * math.pi
        * G_CGS
        * black_hole_mass_msun
        * M_SUN_CGS
        * M_PROTON_CGS
        / (radiative_efficiency * SIGMA_T_CGS * C_CGS)
    )
    return mdot_g_s * MYR_CGS / M_SUN_CGS
