"""Minimal flat-LambdaCDM utilities for the stage-5 redshift budget."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import G_CGS, M_SUN_CGS, MPC_CGS, MYR_CGS, PC_CGS


@dataclass(frozen=True)
class FlatLambdaCDM:
    """Flat matter-plus-Lambda cosmology.

    The defaults are the Planck 2018 base-LambdaCDM values. Radiation is
    omitted; over the stage-5 range `4 <= z <= 30` this is a percent-level
    timing approximation rather than a recombination-era cosmology.
    """

    hubble_km_s_mpc: float = 67.4
    omega_matter: float = 0.315
    omega_lambda: float = 0.685

    def __post_init__(self) -> None:
        if self.hubble_km_s_mpc <= 0.0:
            raise ValueError("hubble_km_s_mpc must be positive")
        if self.omega_matter <= 0.0 or self.omega_lambda <= 0.0:
            raise ValueError("density parameters must be positive")
        if not math.isclose(
            self.omega_matter + self.omega_lambda,
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-10,
        ):
            raise ValueError("FlatLambdaCDM requires omega_matter + omega_lambda = 1")

    @property
    def hubble_s(self) -> float:
        return self.hubble_km_s_mpc * 1.0e5 / MPC_CGS

    def expansion_rate_s(self, redshift: float) -> float:
        if redshift < 0.0:
            raise ValueError("redshift cannot be negative")
        return self.hubble_s * math.sqrt(
            self.omega_matter * (1.0 + redshift) ** 3
            + self.omega_lambda
        )

    def age_myr(self, redshift: float) -> float:
        """Return the analytic matter-plus-Lambda cosmic age."""

        if redshift < 0.0:
            raise ValueError("redshift cannot be negative")
        argument = math.sqrt(self.omega_lambda / self.omega_matter) / (
            1.0 + redshift
        ) ** 1.5
        age_s = (
            2.0
            * math.asinh(argument)
            / (3.0 * self.hubble_s * math.sqrt(self.omega_lambda))
        )
        return age_s / MYR_CGS

    def redshift_at_age_myr(self, age_myr: float) -> float:
        """Invert :meth:`age_myr` within the physical age of this cosmology."""

        if age_myr <= 0.0:
            raise ValueError("age_myr must be positive")
        if age_myr > self.age_myr(0.0):
            raise ValueError("age_myr exceeds the present cosmic age")
        argument = (
            1.5
            * self.hubble_s
            * math.sqrt(self.omega_lambda)
            * age_myr
            * MYR_CGS
        )
        one_plus_redshift = (
            math.sqrt(self.omega_lambda / self.omega_matter)
            / math.sinh(argument)
        ) ** (2.0 / 3.0)
        return max(0.0, one_plus_redshift - 1.0)

    def elapsed_time_myr(
        self,
        formation_redshift: float,
        observation_redshift: float,
    ) -> float:
        if observation_redshift > formation_redshift:
            raise ValueError("observation redshift must not exceed formation redshift")
        return self.age_myr(observation_redshift) - self.age_myr(
            formation_redshift
        )

    def critical_density_msun_pc3(self, redshift: float) -> float:
        hubble = self.expansion_rate_s(redshift)
        density_cgs = 3.0 * hubble**2 / (8.0 * math.pi * G_CGS)
        return density_cgs * PC_CGS**3 / M_SUN_CGS

    def spherical_overdensity_radius_pc(
        self,
        mass_msun: float,
        redshift: float,
        overdensity: float = 200.0,
    ) -> float:
        """Return `r_Delta` for an overdensity relative to critical density."""

        if mass_msun <= 0.0:
            raise ValueError("mass_msun must be positive")
        if overdensity <= 0.0:
            raise ValueError("overdensity must be positive")
        density = overdensity * self.critical_density_msun_pc3(redshift)
        return (3.0 * mass_msun / (4.0 * math.pi * density)) ** (1.0 / 3.0)

    def spherical_overdensity_velocity_km_s(
        self,
        mass_msun: float,
        redshift: float,
        overdensity: float = 200.0,
    ) -> float:
        radius_pc = self.spherical_overdensity_radius_pc(
            mass_msun,
            redshift,
            overdensity,
        )
        return math.sqrt(
            G_CGS * mass_msun * M_SUN_CGS / (radius_pc * PC_CGS)
        ) / 1.0e5

    def black_hole_influence_radius_pc(
        self,
        black_hole_mass_msun: float,
        halo_mass_msun: float,
        redshift: float,
        overdensity: float = 200.0,
    ) -> float:
        """Return `G M_BH / V_Delta^2` for a halo virial velocity."""

        if black_hole_mass_msun <= 0.0:
            raise ValueError("black_hole_mass_msun must be positive")
        radius_pc = self.spherical_overdensity_radius_pc(
            halo_mass_msun,
            redshift,
            overdensity,
        )
        return black_hole_mass_msun / halo_mass_msun * radius_pc
