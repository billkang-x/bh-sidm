"""Code-unit scales for the 1D SIDM fluid equations."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import G_CGS, M_SUN_CGS, PC_CGS, MYR_CGS


@dataclass(frozen=True)
class SimulationScales:
    """Characteristic scales used by the dimensionless SIDM equations.

    The baseline paper defines radius and density scales r0 and rho0, then
    derives M0, v0, t0, p0, and the cross-section-per-mass scale from them.
    Inputs are in pc and Msun/pc^3 to match the paper's halo parameters.
    """

    radius_scale_pc: float
    density_scale_msun_pc3: float

    def __post_init__(self) -> None:
        if self.radius_scale_pc <= 0.0:
            raise ValueError("radius_scale_pc must be positive")
        if self.density_scale_msun_pc3 <= 0.0:
            raise ValueError("density_scale_msun_pc3 must be positive")

    @classmethod
    def for_singular_isothermal_sphere(
        cls,
        sound_speed_km_s: float,
        radius_scale_pc: float = 1.0,
    ) -> "SimulationScales":
        """Return natural scales with rho0 equal to the SIS density at r0."""

        if sound_speed_km_s <= 0.0:
            raise ValueError("sound_speed_km_s must be positive")
        if radius_scale_pc <= 0.0:
            raise ValueError("radius_scale_pc must be positive")
        sound_speed_cgs = sound_speed_km_s * 1.0e5
        radius_cgs = radius_scale_pc * PC_CGS
        density_cgs = sound_speed_cgs**2 / (
            2.0 * math.pi * G_CGS * radius_cgs**2
        )
        density_msun_pc3 = density_cgs * PC_CGS**3 / M_SUN_CGS
        return cls(radius_scale_pc, density_msun_pc3)

    @property
    def radius_scale_cgs(self) -> float:
        return self.radius_scale_pc * PC_CGS

    @property
    def density_scale_cgs(self) -> float:
        return self.density_scale_msun_pc3 * M_SUN_CGS / PC_CGS**3

    @property
    def mass_scale_msun(self) -> float:
        return 4.0 * math.pi * self.radius_scale_pc**3 * self.density_scale_msun_pc3

    @property
    def mass_scale_cgs(self) -> float:
        return self.mass_scale_msun * M_SUN_CGS

    @property
    def velocity_scale_cgs(self) -> float:
        return math.sqrt(G_CGS * self.mass_scale_cgs / self.radius_scale_cgs)

    @property
    def velocity_scale_km_s(self) -> float:
        return self.velocity_scale_cgs / 1.0e5

    @property
    def time_scale_s(self) -> float:
        return 1.0 / math.sqrt(4.0 * math.pi * G_CGS * self.density_scale_cgs)

    @property
    def time_scale_myr(self) -> float:
        return self.time_scale_s / MYR_CGS

    @property
    def pressure_scale_cgs(self) -> float:
        return self.density_scale_cgs * self.velocity_scale_cgs**2

    @property
    def sigma_over_m_scale_cgs(self) -> float:
        return 1.0 / (self.density_scale_cgs * self.radius_scale_cgs)

    @property
    def conductivity_scale_cgs(self) -> float:
        return self.density_scale_cgs * self.radius_scale_cgs**2 / self.time_scale_s

    def radius_to_code(self, radius_pc: float) -> float:
        return radius_pc / self.radius_scale_pc

    def radius_from_code(self, radius_code: float) -> float:
        return radius_code * self.radius_scale_pc

    def density_to_code(self, density_msun_pc3: float) -> float:
        return density_msun_pc3 / self.density_scale_msun_pc3

    def density_from_code(self, density_code: float) -> float:
        return density_code * self.density_scale_msun_pc3

    def mass_to_code(self, mass_msun: float) -> float:
        return mass_msun / self.mass_scale_msun

    def mass_from_code(self, mass_code: float) -> float:
        return mass_code * self.mass_scale_msun

    def velocity_to_code(self, velocity_km_s: float) -> float:
        return velocity_km_s / self.velocity_scale_km_s

    def velocity_from_code(self, velocity_code: float) -> float:
        return velocity_code * self.velocity_scale_km_s

    def time_to_code(self, time_myr: float) -> float:
        return time_myr / self.time_scale_myr

    def time_from_code(self, time_code: float) -> float:
        return time_code * self.time_scale_myr

    def sigma_over_m_to_code(self, sigma_over_m_cm2_g: float) -> float:
        return sigma_over_m_cm2_g / self.sigma_over_m_scale_cgs

    def sigma_over_m_from_code(self, sigma_over_m_code: float) -> float:
        return sigma_over_m_code * self.sigma_over_m_scale_cgs
