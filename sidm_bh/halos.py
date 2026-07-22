"""Spherical dark-matter halo profiles used for baseline experiments."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import G_CGS, M_SUN_CGS, PC_CGS


@dataclass(frozen=True)
class NFWProfile:
    """Navarro-Frenk-White profile.

    Parameters use pc and Msun/pc^3. The density is
    rho = rho_s / [x (1 + x)^2], with x = r / r_s.
    """

    scale_density_msun_pc3: float
    scale_radius_pc: float

    def __post_init__(self) -> None:
        if self.scale_density_msun_pc3 <= 0.0:
            raise ValueError("scale_density_msun_pc3 must be positive")
        if self.scale_radius_pc <= 0.0:
            raise ValueError("scale_radius_pc must be positive")

    @classmethod
    def from_mass_concentration(
        cls,
        total_mass_msun: float,
        virial_radius_pc: float,
        concentration: float,
    ) -> "NFWProfile":
        if total_mass_msun <= 0.0:
            raise ValueError("total_mass_msun must be positive")
        if virial_radius_pc <= 0.0:
            raise ValueError("virial_radius_pc must be positive")
        if concentration <= 0.0:
            raise ValueError("concentration must be positive")

        scale_radius_pc = virial_radius_pc / concentration
        shape = math.log(1.0 + concentration) - concentration / (1.0 + concentration)
        scale_density = total_mass_msun / (
            4.0 * math.pi * scale_radius_pc**3 * shape
        )
        return cls(scale_density, scale_radius_pc)

    def density_msun_pc3(self, radius_pc: float) -> float:
        if radius_pc <= 0.0:
            return math.inf
        x = radius_pc / self.scale_radius_pc
        return self.scale_density_msun_pc3 / (x * (1.0 + x) ** 2)

    def enclosed_mass_msun(self, radius_pc: float) -> float:
        if radius_pc <= 0.0:
            return 0.0
        x = radius_pc / self.scale_radius_pc
        shape = math.log(1.0 + x) - x / (1.0 + x)
        return 4.0 * math.pi * self.scale_density_msun_pc3 * self.scale_radius_pc**3 * shape

    def concentration_for_enclosed_mass(self, total_mass_msun: float) -> float:
        """Return ``r / r_s`` where the profile encloses ``total_mass_msun``."""

        if total_mass_msun <= 0.0:
            raise ValueError("total_mass_msun must be positive")
        lower = 0.0
        upper = 1.0
        while self.enclosed_mass_msun(upper * self.scale_radius_pc) < total_mass_msun:
            upper *= 2.0
        for _ in range(80):
            middle = 0.5 * (lower + upper)
            enclosed = self.enclosed_mass_msun(middle * self.scale_radius_pc)
            if enclosed < total_mass_msun:
                lower = middle
            else:
                upper = middle
        return 0.5 * (lower + upper)

    def self_similar_scaled(
        self,
        total_mass_msun: float,
        anchor_mass_msun: float,
    ) -> "NFWProfile":
        """Scale ``r_s`` at fixed scale density and concentration.

        The anchor mass is understood to be enclosed at one fixed multiple of
        ``r_s``. Consequently ``r_s`` scales as mass to the one-third power.
        """

        if total_mass_msun <= 0.0:
            raise ValueError("total_mass_msun must be positive")
        if anchor_mass_msun <= 0.0:
            raise ValueError("anchor_mass_msun must be positive")
        radius_factor = (total_mass_msun / anchor_mass_msun) ** (1.0 / 3.0)
        return NFWProfile(
            self.scale_density_msun_pc3,
            self.scale_radius_pc * radius_factor,
        )


@dataclass(frozen=True)
class SingularIsothermalSphere:
    """Singular isothermal sphere with rho = c_s^2 / (2 pi G r^2)."""

    sound_speed_km_s: float

    def __post_init__(self) -> None:
        if self.sound_speed_km_s <= 0.0:
            raise ValueError("sound_speed_km_s must be positive")

    @property
    def sound_speed_cgs(self) -> float:
        return self.sound_speed_km_s * 1.0e5

    def density_msun_pc3(self, radius_pc: float) -> float:
        if radius_pc <= 0.0:
            return math.inf
        radius_cm = radius_pc * PC_CGS
        density_cgs = self.sound_speed_cgs**2 / (2.0 * math.pi * G_CGS * radius_cm**2)
        return density_cgs * PC_CGS**3 / M_SUN_CGS

    def enclosed_mass_msun(self, radius_pc: float) -> float:
        if radius_pc <= 0.0:
            return 0.0
        radius_cm = radius_pc * PC_CGS
        mass_cgs = 2.0 * self.sound_speed_cgs**2 * radius_cm / G_CGS
        return mass_cgs / M_SUN_CGS
