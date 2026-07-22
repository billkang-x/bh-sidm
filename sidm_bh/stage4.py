"""Stage-4 baryonic accretion and parametric feedback utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .baryons import eddington_accretion_rate_msun_per_myr
from .constants import C_CGS, G_CGS, M_SUN_CGS, MYR_CGS, PC_CGS
from .halos import NFWProfile
from .units import SimulationScales


@dataclass(frozen=True)
class EddingtonBaryonModel:
    """Finite-reservoir Eddington accretion with Hernquist expansion.

    ``eddington_ratio`` and ``duty_cycle`` multiply the gas inflow rate.
    A fraction ``1-radiative_efficiency`` of that inflow is retained by the
    black hole. Feedback couples a fraction ``feedback_efficiency`` of the
    radiated rest-mass energy to the baryonic component.
    """

    radiative_efficiency: float = 0.1
    eddington_ratio: float = 1.0
    duty_cycle: float = 1.0
    feedback_efficiency: float = 0.0
    feedback_expansion_exponent: float = 0.5
    gas_density_msun_pc3: float | None = None
    gas_sound_speed_km_s: float = 10.0
    gas_relative_velocity_km_s: float = 0.0
    bondi_alpha: float = 1.0

    def __post_init__(self) -> None:
        if not 0.0 < self.radiative_efficiency < 1.0:
            raise ValueError("radiative_efficiency must lie in (0, 1)")
        if not 0.0 <= self.eddington_ratio <= 1.0:
            raise ValueError("eddington_ratio must lie in [0, 1]")
        if not 0.0 <= self.duty_cycle <= 1.0:
            raise ValueError("duty_cycle must lie in [0, 1]")
        if not 0.0 <= self.feedback_efficiency <= 1.0:
            raise ValueError("feedback_efficiency must lie in [0, 1]")
        if self.feedback_expansion_exponent < 0.0:
            raise ValueError("feedback_expansion_exponent cannot be negative")
        if self.gas_density_msun_pc3 is not None and self.gas_density_msun_pc3 <= 0.0:
            raise ValueError("gas_density_msun_pc3 must be positive when supplied")
        if self.gas_sound_speed_km_s <= 0.0:
            raise ValueError("gas_sound_speed_km_s must be positive")
        if self.gas_relative_velocity_km_s < 0.0:
            raise ValueError("gas_relative_velocity_km_s cannot be negative")
        if self.bondi_alpha <= 0.0:
            raise ValueError("bondi_alpha must be positive")

    @property
    def retained_fraction(self) -> float:
        return 1.0 - self.radiative_efficiency

    @property
    def eddington_inflow_coefficient_per_myr(self) -> float:
        return (
            self.eddington_ratio
            * self.duty_cycle
            * eddington_accretion_rate_msun_per_myr(
                1.0,
                self.radiative_efficiency,
            )
        )

    @property
    def black_hole_efolding_time_myr(self) -> float:
        coefficient = (
            self.retained_fraction
            * self.eddington_inflow_coefficient_per_myr
        )
        if coefficient == 0.0:
            return float("inf")
        return 1.0 / coefficient

    def eddington_inflow_coefficient_code(
        self,
        scales: SimulationScales,
    ) -> float:
        return self.eddington_inflow_coefficient_per_myr * scales.time_scale_myr

    @property
    def bondi_inflow_coefficient_msun_inv_myr(self) -> float | None:
        if self.gas_density_msun_pc3 is None:
            return None
        return bondi_accretion_rate_msun_per_myr(
            black_hole_mass_msun=1.0,
            gas_density_msun_pc3=self.gas_density_msun_pc3,
            sound_speed_km_s=self.gas_sound_speed_km_s,
            relative_velocity_km_s=self.gas_relative_velocity_km_s,
            alpha=self.bondi_alpha,
        )

    def bondi_inflow_coefficient_code(
        self,
        scales: SimulationScales,
    ) -> float:
        coefficient = self.bondi_inflow_coefficient_msun_inv_myr
        if coefficient is None:
            return -1.0
        return coefficient * scales.mass_scale_msun * scales.time_scale_myr

    def feedback_ratio_per_consumed_mass_code(
        self,
        scales: SimulationScales,
        initial_baryon_mass_msun: float,
        initial_scale_radius_pc: float,
        binding_energy_erg: float | None = None,
        expansion_energy_fraction: float = 1.0,
    ) -> float:
        """Return ``d(E_fb/E_bind)`` per code-unit gas mass consumed."""

        binding_energy = binding_energy_erg
        if binding_energy is None:
            binding_energy = hernquist_self_binding_energy_erg(
                initial_baryon_mass_msun,
                initial_scale_radius_pc,
            )
        if binding_energy <= 0.0:
            raise ValueError("binding_energy_erg must be positive")
        if not 0.0 <= expansion_energy_fraction <= 1.0:
            raise ValueError("expansion_energy_fraction must lie in [0, 1]")
        coupled_energy_per_code_mass = (
            expansion_energy_fraction * self.feedback_efficiency
            * self.radiative_efficiency
            * scales.mass_scale_cgs
            * C_CGS**2
        )
        return coupled_energy_per_code_mass / binding_energy

    def feedback_thermal_energy_per_consumed_mass_code(
        self,
        scales: SimulationScales,
        heating_energy_fraction: float,
    ) -> float:
        """Return deposited thermal energy per gas mass in code velocity^2."""

        if not 0.0 <= heating_energy_fraction <= 1.0:
            raise ValueError("heating_energy_fraction must lie in [0, 1]")
        light_speed_code = C_CGS / scales.velocity_scale_cgs
        return (
            heating_energy_fraction
            * self.feedback_efficiency
            * self.radiative_efficiency
            * light_speed_code**2
        )


@dataclass(frozen=True)
class BindingEnergyComponents:
    self_gravity_erg: float
    halo_erg: float
    black_hole_erg: float

    @property
    def total_erg(self) -> float:
        return self.self_gravity_erg + self.halo_erg + self.black_hole_erg


def bondi_accretion_rate_msun_per_myr(
    black_hole_mass_msun: float,
    gas_density_msun_pc3: float,
    sound_speed_km_s: float,
    relative_velocity_km_s: float = 0.0,
    alpha: float = 1.0,
) -> float:
    """Return the spherical Bondi-Hoyle gas inflow rate."""

    if black_hole_mass_msun <= 0.0:
        return 0.0
    if gas_density_msun_pc3 <= 0.0:
        raise ValueError("gas_density_msun_pc3 must be positive")
    if sound_speed_km_s <= 0.0:
        raise ValueError("sound_speed_km_s must be positive")
    if relative_velocity_km_s < 0.0:
        raise ValueError("relative_velocity_km_s cannot be negative")
    if alpha <= 0.0:
        raise ValueError("alpha must be positive")
    mass_cgs = black_hole_mass_msun * M_SUN_CGS
    density_cgs = gas_density_msun_pc3 * M_SUN_CGS / PC_CGS**3
    effective_speed_cgs = (
        np.sqrt(sound_speed_km_s**2 + relative_velocity_km_s**2) * 1.0e5
    )
    rate_cgs = (
        4.0
        * np.pi
        * alpha
        * G_CGS**2
        * mass_cgs**2
        * density_cgs
        / effective_speed_cgs**3
    )
    return float(rate_cgs * MYR_CGS / M_SUN_CGS)


def effective_hernquist_binding_energy_erg(
    total_mass_msun: float,
    scale_radius_pc: float,
    halo: NFWProfile,
    black_hole_mass_msun: float,
    integration_points: int = 4096,
) -> BindingEnergyComponents:
    """Return initial self, NFW, and black-hole binding components."""

    if black_hole_mass_msun < 0.0:
        raise ValueError("black_hole_mass_msun cannot be negative")
    if integration_points < 128:
        raise ValueError("integration_points must be at least 128")
    self_energy = hernquist_self_binding_energy_erg(
        total_mass_msun,
        scale_radius_pc,
    )
    baryon_mass_cgs = total_mass_msun * M_SUN_CGS
    scale_radius_cgs = scale_radius_pc * PC_CGS
    black_hole_energy = (
        G_CGS
        * black_hole_mass_msun
        * M_SUN_CGS
        * baryon_mass_cgs
        / scale_radius_cgs
    )

    x = np.geomspace(1.0e-8, 1.0e8, integration_points)
    radius_cgs = x * scale_radius_cgs
    halo_scale_radius_cgs = halo.scale_radius_pc * PC_CGS
    halo_density_cgs = (
        halo.scale_density_msun_pc3 * M_SUN_CGS / PC_CGS**3
    )
    halo_potential = -(
        4.0
        * np.pi
        * G_CGS
        * halo_density_cgs
        * halo_scale_radius_cgs**3
        * np.log1p(radius_cgs / halo_scale_radius_cgs)
        / radius_cgs
    )
    mass_derivative = 2.0 * baryon_mass_cgs * x / (1.0 + x) ** 3
    integrand = halo_potential * mass_derivative
    halo_energy = -float(
        np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(x))
    )
    return BindingEnergyComponents(
        self_gravity_erg=self_energy,
        halo_erg=halo_energy,
        black_hole_erg=black_hole_energy,
    )


def hernquist_self_binding_energy_erg(
    total_mass_msun: float,
    scale_radius_pc: float,
) -> float:
    """Return the magnitude of the Hernquist self-potential energy."""

    if total_mass_msun <= 0.0:
        raise ValueError("total_mass_msun must be positive")
    if scale_radius_pc <= 0.0:
        raise ValueError("scale_radius_pc must be positive")
    mass_cgs = total_mass_msun * M_SUN_CGS
    radius_cgs = scale_radius_pc * PC_CGS
    return G_CGS * mass_cgs**2 / (6.0 * radius_cgs)


def feedback_expanded_scale_radius(
    initial_scale_radius: float,
    feedback_to_binding_ratio: float,
    expansion_exponent: float,
) -> float:
    if initial_scale_radius <= 0.0:
        raise ValueError("initial_scale_radius must be positive")
    if feedback_to_binding_ratio < 0.0:
        raise ValueError("feedback_to_binding_ratio cannot be negative")
    if expansion_exponent < 0.0:
        raise ValueError("expansion_exponent cannot be negative")
    return initial_scale_radius * (
        1.0 + feedback_to_binding_ratio
    ) ** expansion_exponent


def homologous_ambient_density(
    initial_density: float,
    remaining_mass_fraction: float,
    scale_radius_expansion_factor: float,
) -> float:
    """Return ambient density for a homologously expanding gas reservoir."""

    if initial_density <= 0.0:
        raise ValueError("initial_density must be positive")
    if not 0.0 <= remaining_mass_fraction <= 1.0:
        raise ValueError("remaining_mass_fraction must lie in [0, 1]")
    if scale_radius_expansion_factor <= 0.0:
        raise ValueError("scale_radius_expansion_factor must be positive")
    return (
        initial_density
        * remaining_mass_fraction
        / scale_radius_expansion_factor**3
    )


def feedback_heated_sound_speed(
    initial_sound_speed: float,
    deposited_energy_per_mass: float,
    adiabatic_index: float = 5.0 / 3.0,
) -> float:
    """Return sound speed after adding specific thermal energy.

    The sound speed and deposited specific energy must use consistent units,
    with the latter expressed as velocity squared.
    """

    if initial_sound_speed <= 0.0:
        raise ValueError("initial_sound_speed must be positive")
    if deposited_energy_per_mass < 0.0:
        raise ValueError("deposited_energy_per_mass cannot be negative")
    if adiabatic_index <= 1.0:
        raise ValueError("adiabatic_index must exceed one")
    sound_squared = initial_sound_speed**2 + (
        adiabatic_index
        * (adiabatic_index - 1.0)
        * deposited_energy_per_mass
    )
    return float(np.sqrt(sound_squared))
