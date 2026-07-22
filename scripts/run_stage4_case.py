"""Run one stage-4 Eddington baryon-accretion and feedback experiment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sidm_bh.baryons import HernquistBaryons, smoothstep_mass_fraction
from sidm_bh.constants import G_CGS, M_SUN_CGS, PC_CGS
from sidm_bh.cosmology import FlatLambdaCDM
from sidm_bh.cooling import cloudy_cooling_state, load_cloudy_cooling_table
from sidm_bh.fast_evolution import evolve_mc_roe_fast
from sidm_bh.halos import NFWProfile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sidm import maxwellian_viscosity_cross_section_ratio
from sidm_bh.stage3 import static_baryon_equilibrium_state
from sidm_bh.stage4 import (
    EddingtonBaryonModel,
    effective_hernquist_binding_energy_erg,
)
from sidm_bh.units import SimulationScales


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--halo-mass-msun", type=float, default=1.0e6)
    parser.add_argument("--halo-redshift", type=float)
    parser.add_argument("--halo-concentration", type=float)
    parser.add_argument("--black-hole-mass-msun", type=float, default=100.0)
    parser.add_argument("--baryon-fraction", type=float, default=0.05)
    parser.add_argument("--scale-radius-over-rs", type=float, default=0.01)
    parser.add_argument("--assembly-time-myr", type=float, default=0.65)
    parser.add_argument("--eddington-ratio", type=float, default=1.0)
    parser.add_argument("--duty-cycle", type=float, default=1.0)
    parser.add_argument("--radiative-efficiency", type=float, default=0.1)
    parser.add_argument("--feedback-efficiency", type=float, default=0.0)
    parser.add_argument("--feedback-eta", type=float, default=0.5)
    parser.add_argument("--feedback-heating-fraction", type=float, default=0.0)
    parser.add_argument("--feedback-binding-model", choices=("self", "effective_initial"), default="self")
    parser.add_argument("--gas-density-msun-pc3", type=float)
    parser.add_argument("--gas-sound-speed-km-s", type=float, default=10.0)
    parser.add_argument("--gas-relative-velocity-km-s", type=float, default=0.0)
    parser.add_argument("--gas-adiabatic-index", type=float, default=5.0 / 3.0)
    parser.add_argument("--gas-cooling-time-myr", type=float)
    parser.add_argument("--cloudy-cooling", action="store_true")
    parser.add_argument("--gas-metallicity-solar", type=float, default=0.0)
    parser.add_argument("--gas-hydrogen-mass-fraction", type=float, default=0.76)
    parser.add_argument("--cooling-rate-multiplier", type=float, default=1.0)
    parser.add_argument("--evolve-bondi-ambient", action="store_true")
    parser.add_argument("--bondi-alpha", type=float, default=1.0)
    parser.add_argument("--sigma-over-m-cm2-g", type=float, default=50.0)
    parser.add_argument(
        "--cross-section-model",
        choices=("constant", "rutherford"),
        default="constant",
    )
    parser.add_argument("--cross-section-velocity-km-s", type=float)
    parser.add_argument(
        "--dark-capture-model",
        choices=("boundary_flux", "bondi_reservoir", "influence_gated"),
        default="boundary_flux",
    )
    parser.add_argument("--dark-bondi-lambda", type=float, default=0.25)
    parser.add_argument(
        "--dark-flux-capture-rmin-over-influence",
        type=float,
        default=0.10416666666666667,
    )
    parser.add_argument("--r-min-pc", type=float)
    parser.add_argument("--r-min-over-rs", type=float)
    parser.add_argument("--r-max-pc", type=float)
    parser.add_argument("--r-max-over-rs", type=float)
    parser.add_argument("--cells", type=int)
    parser.add_argument("--duration-myr", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=201)
    parser.add_argument("--cfl", type=float, default=0.2)
    parser.add_argument("--entropy-fix", type=float, default=0.1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=20_000_000)
    return parser.parse_args()


def sustained_dark_dominance_onset_myr(
    times_myr: np.ndarray,
    dark_rate: np.ndarray,
    baryon_rate: np.ndarray,
    minimum_time_myr: float,
) -> float:
    eligible = times_myr >= minimum_time_myr
    for index in np.flatnonzero(eligible):
        if np.all(dark_rate[index:] >= baryon_rate[index:]):
            return float(times_myr[index])
    return float("nan")


def first_target_time_myr(
    times_myr: np.ndarray,
    masses_msun: np.ndarray,
    target_msun: float,
) -> float:
    reached = np.flatnonzero(masses_msun >= target_msun)
    return float(times_myr[reached[0]]) if len(reached) else float("nan")


def main() -> None:
    args = parse_args()
    if args.halo_mass_msun <= 0.0:
        raise ValueError("halo_mass_msun must be positive")
    if args.black_hole_mass_msun <= 0.0:
        raise ValueError("black_hole_mass_msun must be positive")
    if (args.halo_redshift is None) != (args.halo_concentration is None):
        raise ValueError("halo redshift and concentration must be supplied together")
    if args.halo_redshift is not None and args.halo_redshift < 0.0:
        raise ValueError("halo_redshift cannot be negative")
    if args.halo_concentration is not None and args.halo_concentration <= 0.0:
        raise ValueError("halo_concentration must be positive")
    if not 0.0 < args.baryon_fraction <= 1.0:
        raise ValueError("baryon_fraction must lie in (0, 1]")
    if args.scale_radius_over_rs <= 0.0:
        raise ValueError("scale_radius_over_rs must be positive")
    if args.assembly_time_myr < 0.0:
        raise ValueError("assembly_time_myr cannot be negative")
    if args.duration_myr <= 0.0:
        raise ValueError("duration_myr must be positive")
    if args.samples < 2:
        raise ValueError("samples must be at least two")
    if not 0.0 <= args.feedback_heating_fraction <= 1.0:
        raise ValueError("feedback_heating_fraction must lie in [0, 1]")
    if args.gas_adiabatic_index <= 1.0:
        raise ValueError("gas_adiabatic_index must exceed one")
    if args.gas_cooling_time_myr is not None and (
        np.isnan(args.gas_cooling_time_myr) or args.gas_cooling_time_myr <= 0.0
    ):
        raise ValueError("gas_cooling_time_myr must be positive when supplied")
    if args.cloudy_cooling and args.gas_cooling_time_myr is not None:
        raise ValueError("Cloudy and fixed cooling cannot be combined")
    if args.gas_metallicity_solar < 0.0:
        raise ValueError("gas_metallicity_solar cannot be negative")
    if not 0.0 < args.gas_hydrogen_mass_fraction <= 1.0:
        raise ValueError("gas_hydrogen_mass_fraction must lie in (0, 1]")
    if args.cooling_rate_multiplier <= 0.0:
        raise ValueError("cooling_rate_multiplier must be positive")
    if args.dark_bondi_lambda < 0.0:
        raise ValueError("dark_bondi_lambda cannot be negative")
    if args.sigma_over_m_cm2_g < 0.0:
        raise ValueError("sigma_over_m_cm2_g cannot be negative")
    if args.cross_section_model == "rutherford" and (
        args.cross_section_velocity_km_s is None
        or args.cross_section_velocity_km_s <= 0.0
    ):
        raise ValueError(
            "rutherford cross section requires a positive velocity scale"
        )
    if args.dark_flux_capture_rmin_over_influence <= 0.0:
        raise ValueError(
            "dark_flux_capture_rmin_over_influence must be positive"
        )
    if args.evolve_bondi_ambient and args.gas_density_msun_pc3 is None:
        raise ValueError("evolving Bondi ambient requires gas density")
    if args.r_min_pc is not None and args.r_min_over_rs is not None:
        raise ValueError("r_min_pc and r_min_over_rs are mutually exclusive")
    if args.r_max_pc is not None and args.r_max_over_rs is not None:
        raise ValueError("r_max_pc and r_max_over_rs are mutually exclusive")
    for name, value in (
        ("r_min_pc", args.r_min_pc),
        ("r_min_over_rs", args.r_min_over_rs),
        ("r_max_pc", args.r_max_pc),
        ("r_max_over_rs", args.r_max_over_rs),
    ):
        if value is not None and value <= 0.0:
            raise ValueError(f"{name} must be positive when supplied")

    anchor_mass_msun = 1.0e6
    anchor_profile = NFWProfile(3.7, 30.0)
    if args.halo_redshift is None:
        profile = anchor_profile.self_similar_scaled(
            args.halo_mass_msun,
            anchor_mass_msun,
        )
        halo_profile_mode = "fixed_density_self_similar"
        halo_concentration = profile.concentration_for_enclosed_mass(
            args.halo_mass_msun
        )
        halo_virial_radius_pc = halo_concentration * profile.scale_radius_pc
        halo_virial_velocity_km_s = np.sqrt(
            G_CGS
            * args.halo_mass_msun
            * M_SUN_CGS
            / (halo_virial_radius_pc * PC_CGS)
        ) / 1.0e5
        black_hole_influence_radius_pc = (
            args.black_hole_mass_msun
            / args.halo_mass_msun
            * halo_virial_radius_pc
        )
    else:
        cosmology = FlatLambdaCDM()
        halo_virial_radius_pc = cosmology.spherical_overdensity_radius_pc(
            args.halo_mass_msun,
            args.halo_redshift,
        )
        halo_concentration = args.halo_concentration
        halo_virial_velocity_km_s = (
            cosmology.spherical_overdensity_velocity_km_s(
                args.halo_mass_msun,
                args.halo_redshift,
            )
        )
        black_hole_influence_radius_pc = cosmology.black_hole_influence_radius_pc(
            args.black_hole_mass_msun,
            args.halo_mass_msun,
            args.halo_redshift,
        )
        profile = NFWProfile.from_mass_concentration(
            args.halo_mass_msun,
            halo_virial_radius_pc,
            halo_concentration,
        )
        halo_profile_mode = "m200c_at_redshift"
    scales = SimulationScales(
        profile.scale_radius_pc,
        profile.scale_density_msun_pc3,
    )
    radius_factor = profile.scale_radius_pc / anchor_profile.scale_radius_pc
    r_min_pc = (
        args.r_min_pc
        if args.r_min_pc is not None
        else (
            args.r_min_over_rs * profile.scale_radius_pc
            if args.r_min_over_rs is not None
            else 0.005
        )
    )
    r_max_pc = args.r_max_pc
    if r_max_pc is None:
        r_max_pc = (
            args.r_max_over_rs * profile.scale_radius_pc
            if args.r_max_over_rs is not None
            else 5000.0 * radius_factor
        )
    if r_max_pc <= r_min_pc:
        raise ValueError("outer radius must exceed inner radius")
    baseline_log_width = np.log(5000.0 / 0.005) / 256.0
    cells = (
        args.cells
        if args.cells is not None
        else int(round(np.log(r_max_pc / r_min_pc) / baseline_log_width))
    )
    grid = SphericalGrid.from_log_spacing(
        scales.radius_to_code(r_min_pc),
        scales.radius_to_code(r_max_pc),
        cells,
    )
    initial_dark_reservoir_msun = (
        profile.enclosed_mass_msun(r_min_pc)
        if args.dark_capture_model in ("bondi_reservoir", "influence_gated")
        else 0.0
    )
    dark_flux_capture_mass_threshold_msun = (
        r_min_pc
        / halo_virial_radius_pc
        * args.halo_mass_msun
        / args.dark_flux_capture_rmin_over_influence
    )
    initial_state, _ = static_baryon_equilibrium_state(
        profile,
        grid,
        scales,
        args.black_hole_mass_msun,
        baryons=None,
    )
    baryon_mass_msun = args.baryon_fraction * args.halo_mass_msun
    baryon_scale_radius_pc = (
        args.scale_radius_over_rs * profile.scale_radius_pc
    )
    baryons = HernquistBaryons(
        baryon_mass_msun,
        baryon_scale_radius_pc,
    )
    model = EddingtonBaryonModel(
        radiative_efficiency=args.radiative_efficiency,
        eddington_ratio=args.eddington_ratio,
        duty_cycle=args.duty_cycle,
        feedback_efficiency=args.feedback_efficiency,
        feedback_expansion_exponent=args.feedback_eta,
        gas_density_msun_pc3=args.gas_density_msun_pc3,
        gas_sound_speed_km_s=args.gas_sound_speed_km_s,
        gas_relative_velocity_km_s=args.gas_relative_velocity_km_s,
        bondi_alpha=args.bondi_alpha,
    )
    binding_components = effective_hernquist_binding_energy_erg(
        baryon_mass_msun,
        baryon_scale_radius_pc,
        profile,
        args.black_hole_mass_msun,
    )
    feedback_binding_energy_erg = (
        binding_components.self_gravity_erg
        if args.feedback_binding_model == "self"
        else binding_components.total_erg
    )
    times_myr = np.linspace(0.0, args.duration_myr, args.samples)
    assembly_time_code = (
        None
        if args.assembly_time_myr == 0.0
        else scales.time_to_code(args.assembly_time_myr)
    )
    expansion_energy_fraction = 1.0 - args.feedback_heating_fraction
    expansion_feedback_per_consumed_mass_code = (
        model.feedback_ratio_per_consumed_mass_code(
            scales,
            baryon_mass_msun,
            baryon_scale_radius_pc,
            binding_energy_erg=feedback_binding_energy_erg,
            expansion_energy_fraction=expansion_energy_fraction,
        )
    )
    total_feedback_per_consumed_mass_code = (
        model.feedback_ratio_per_consumed_mass_code(
            scales,
            baryon_mass_msun,
            baryon_scale_radius_pc,
            binding_energy_erg=feedback_binding_energy_erg,
        )
    )
    thermal_energy_per_consumed_mass_code = (
        model.feedback_thermal_energy_per_consumed_mass_code(
            scales,
            args.feedback_heating_fraction,
        )
    )
    cooling_time_myr = (
        None
        if args.gas_cooling_time_myr is None
        or np.isinf(args.gas_cooling_time_myr)
        else args.gas_cooling_time_myr
    )
    cloudy_table = load_cloudy_cooling_table() if args.cloudy_cooling else None

    started = perf_counter()
    result = evolve_mc_roe_fast(
        initial_state,
        grid,
        scales.time_to_code(times_myr),
        scales.mass_to_code(args.black_hole_mass_msun),
        scales.sigma_over_m_to_code(args.sigma_over_m_cm2_g),
        cross_section_model=args.cross_section_model,
        cross_section_velocity_scale_code=(
            None
            if args.cross_section_velocity_km_s is None
            else scales.velocity_to_code(args.cross_section_velocity_km_s)
        ),
        cfl_number=args.cfl,
        entropy_fix=args.entropy_fix,
        dark_capture_model=args.dark_capture_model,
        dark_bondi_lambda=args.dark_bondi_lambda,
        dark_flux_capture_mass_threshold_code=scales.mass_to_code(
            dark_flux_capture_mass_threshold_msun
        ),
        initial_dark_matter_reservoir_code=scales.mass_to_code(
            initial_dark_reservoir_msun
        ),
        baryon_assembly_time_code=assembly_time_code,
        evolving_baryon_total_mass_code=scales.mass_to_code(
            baryon_mass_msun
        ),
        evolving_baryon_scale_radius_code=scales.radius_to_code(
            baryon_scale_radius_pc
        ),
        baryon_eddington_inflow_coefficient_code=(
            model.eddington_inflow_coefficient_code(scales)
        ),
        baryon_bondi_inflow_coefficient_code=(
            model.bondi_inflow_coefficient_code(scales)
        ),
        evolve_bondi_ambient=args.evolve_bondi_ambient,
        baryon_initial_sound_speed_code=scales.velocity_to_code(
            args.gas_sound_speed_km_s
        ),
        baryon_relative_velocity_code=scales.velocity_to_code(
            args.gas_relative_velocity_km_s
        ),
        feedback_thermal_energy_per_consumed_mass_code=(
            thermal_energy_per_consumed_mass_code
        ),
        baryon_cooling_time_code=(
            None
            if cooling_time_myr is None
            else scales.time_to_code(cooling_time_myr)
        ),
        baryon_cloudy_cooling_table=cloudy_table,
        baryon_initial_gas_density_cgs=(
            0.0
            if args.gas_density_msun_pc3 is None
            else args.gas_density_msun_pc3 * M_SUN_CGS / PC_CGS**3
        ),
        baryon_hydrogen_mass_fraction=args.gas_hydrogen_mass_fraction,
        baryon_metallicity_solar=args.gas_metallicity_solar,
        baryon_cooling_rate_multiplier=args.cooling_rate_multiplier,
        baryon_velocity_unit_cgs=scales.velocity_scale_cgs,
        baryon_time_unit_s=scales.time_scale_s,
        baryon_adiabatic_index=args.gas_adiabatic_index,
        baryon_radiative_efficiency=args.radiative_efficiency,
        feedback_ratio_per_consumed_mass_code=(
            expansion_feedback_per_consumed_mass_code
        ),
        feedback_expansion_exponent=args.feedback_eta,
        source_integration="euler",
        max_steps=args.max_steps,
    )
    elapsed = perf_counter() - started

    black_hole_mass_msun = result.black_hole_masses_code * scales.mass_scale_msun
    dark_accreted_msun = (
        result.dark_matter_accreted_masses_code * scales.mass_scale_msun
    )
    baryon_accreted_msun = (
        result.baryon_accreted_masses_code * scales.mass_scale_msun
    )
    baryon_gas_consumed_msun = (
        result.baryon_gas_consumed_masses_code * scales.mass_scale_msun
    )
    energy_scale_erg = scales.mass_scale_cgs * scales.velocity_scale_cgs**2
    baryon_thermal_energy_erg = (
        result.baryon_thermal_energies_code * energy_scale_erg
    )
    injected_thermal_energy_erg = (
        thermal_energy_per_consumed_mass_code
        * result.baryon_gas_consumed_masses_code
        * energy_scale_erg
    )
    cooling_loss_energy_erg = np.maximum(
        injected_thermal_energy_erg - baryon_thermal_energy_erg,
        0.0,
    )
    dark_rate_msun_myr = (
        result.accretion_rates_code
        * scales.mass_scale_msun
        / scales.time_scale_myr
    )
    dark_supply_rate_msun_myr = (
        result.dark_matter_supply_rates_code
        * scales.mass_scale_msun
        / scales.time_scale_myr
    )
    dark_supplied_msun = (
        result.dark_matter_supplied_masses_code * scales.mass_scale_msun
    )
    inner_dark_reservoir_msun = (
        result.inner_dark_matter_reservoir_masses_code
        * scales.mass_scale_msun
    )
    baryon_rate_msun_myr = (
        result.baryon_accretion_rates_code
        * scales.mass_scale_msun
        / scales.time_scale_myr
    )
    eddington_gas_limit_msun_myr = (
        model.eddington_inflow_coefficient_per_myr * black_hole_mass_msun
    )
    bondi_coefficient = model.bondi_inflow_coefficient_msun_inv_myr
    remaining_fraction = (
        result.baryon_remaining_masses_code
        / scales.mass_to_code(baryon_mass_msun)
    )
    scale_radius_expansion = (
        result.baryon_scale_radii_code
        / scales.radius_to_code(baryon_scale_radius_pc)
    )
    ambient_density_factor = np.ones_like(times_myr)
    ambient_speed_factor = np.ones_like(times_myr)
    ambient_sound_speed_code = np.full_like(
        times_myr,
        scales.velocity_to_code(args.gas_sound_speed_km_s),
    )
    if args.evolve_bondi_ambient:
        ambient_density_factor = (
            remaining_fraction / scale_radius_expansion**3
        )
        populated = result.baryon_remaining_masses_code > 0.0
        deposited_specific_energy_code = np.zeros_like(times_myr)
        deposited_specific_energy_code[populated] = (
            result.baryon_thermal_energies_code[populated]
            / result.baryon_remaining_masses_code[populated]
        )
        ambient_sound_speed_code = np.sqrt(
            scales.velocity_to_code(args.gas_sound_speed_km_s) ** 2
            + args.gas_adiabatic_index
            * (args.gas_adiabatic_index - 1.0)
            * deposited_specific_energy_code
        )
        relative_velocity_code = scales.velocity_to_code(
            args.gas_relative_velocity_km_s
        )
        initial_effective_speed_squared = (
            scales.velocity_to_code(args.gas_sound_speed_km_s) ** 2
            + relative_velocity_code**2
        )
        current_effective_speed_squared = (
            ambient_sound_speed_code**2 + relative_velocity_code**2
        )
        ambient_speed_factor = (
            initial_effective_speed_squared / current_effective_speed_squared
        ) ** 1.5
    ambient_gas_density_msun_pc3 = np.full_like(times_myr, np.nan)
    if args.gas_density_msun_pc3 is not None:
        ambient_gas_density_msun_pc3 = (
            args.gas_density_msun_pc3 * ambient_density_factor
        )
    ambient_gas_sound_speed_km_s = scales.velocity_from_code(
        ambient_sound_speed_code
    )
    if args.cross_section_model == "constant":
        lmfp_effective_sigma_over_m_cm2_g = np.full_like(
            result.velocity_dispersion_snapshots,
            args.sigma_over_m_cm2_g,
        )
        smfp_effective_sigma_over_m_cm2_g = (
            lmfp_effective_sigma_over_m_cm2_g.copy()
        )
    else:
        dispersion_ratio = (
            scales.velocity_from_code(result.velocity_dispersion_snapshots)
            / args.cross_section_velocity_km_s
        )
        lmfp_effective_sigma_over_m_cm2_g = (
            args.sigma_over_m_cm2_g
            * maxwellian_viscosity_cross_section_ratio(
                dispersion_ratio,
                velocity_power=3,
            )
        )
        smfp_effective_sigma_over_m_cm2_g = (
            args.sigma_over_m_cm2_g
            * maxwellian_viscosity_cross_section_ratio(
                dispersion_ratio,
                velocity_power=5,
            )
        )
    ambient_gas_temperature_k = np.full_like(times_myr, np.nan)
    ambient_mean_molecular_weight = np.full_like(times_myr, np.nan)
    physical_cooling_time_myr = np.full_like(times_myr, np.inf)
    cooling_coefficient_erg_cm3_s = np.zeros_like(times_myr)
    if cloudy_table is not None:
        for index in np.flatnonzero(ambient_gas_density_msun_pc3 > 0.0):
            cooling_state = cloudy_cooling_state(
                cloudy_table,
                ambient_gas_density_msun_pc3[index],
                ambient_gas_sound_speed_km_s[index],
                args.gas_metallicity_solar,
                cooling_rate_multiplier=args.cooling_rate_multiplier,
                hydrogen_mass_fraction=args.gas_hydrogen_mass_fraction,
                adiabatic_index=args.gas_adiabatic_index,
            )
            ambient_gas_temperature_k[index] = cooling_state.temperature_k
            ambient_mean_molecular_weight[index] = (
                cooling_state.mean_molecular_weight
            )
            physical_cooling_time_myr[index] = cooling_state.cooling_time_myr
            cooling_coefficient_erg_cm3_s[index] = (
                cooling_state.cooling_coefficient_erg_cm3_s
            )
    bondi_gas_limit_msun_myr = np.full_like(times_myr, np.inf)
    if bondi_coefficient is not None:
        bondi_gas_limit_msun_myr = (
            bondi_coefficient
            * ambient_density_factor
            * ambient_speed_factor
            * black_hole_mass_msun**2
        )
    retained_eddington_limit_msun_myr = (
        model.retained_fraction * eddington_gas_limit_msun_myr
    )
    retained_bondi_limit_msun_myr = (
        model.retained_fraction * bondi_gas_limit_msun_myr
    )
    bondi_limited = retained_bondi_limit_msun_myr < retained_eddington_limit_msun_myr
    assembled_fraction = np.array(
        [
            1.0
            if args.assembly_time_myr == 0.0
            else smoothstep_mass_fraction(time, args.assembly_time_myr)
            for time in times_myr
        ]
    )
    onset = sustained_dark_dominance_onset_myr(
        times_myr,
        dark_rate_msun_myr,
        baryon_rate_msun_myr,
        args.assembly_time_myr,
    )
    assembled_samples = times_myr >= args.assembly_time_myr
    bondi_limited_fraction = (
        float(np.mean(bondi_limited[assembled_samples]))
        if bondi_coefficient is not None and np.any(assembled_samples)
        else 0.0
    )
    eddington_transition_time = float("nan")
    if bondi_coefficient is not None:
        transitioned = np.flatnonzero(assembled_samples & ~bondi_limited)
        if len(transitioned):
            eddington_transition_time = float(times_myr[transitioned[0]])
    total_growth = black_hole_mass_msun[-1] - args.black_hole_mass_msun
    metadata = {
        "stage": 4,
        "profile": "nfw_plus_evolving_hernquist",
        "initial_condition": "hydrostatic_dm_bh_no_baryon",
        "halo_mass_msun": args.halo_mass_msun,
        "halo_profile_mode": halo_profile_mode,
        "halo_redshift": args.halo_redshift,
        "halo_concentration": halo_concentration,
        "halo_virial_radius_pc": halo_virial_radius_pc,
        "halo_virial_velocity_km_s": halo_virial_velocity_km_s,
        "black_hole_influence_radius_pc": black_hole_influence_radius_pc,
        "halo_scale_radius_pc": profile.scale_radius_pc,
        "halo_scale_density_msun_pc3": profile.scale_density_msun_pc3,
        "black_hole_seed_msun": args.black_hole_mass_msun,
        "baryon_fraction": args.baryon_fraction,
        "initial_baryon_reservoir_msun": baryon_mass_msun,
        "initial_baryon_scale_radius_pc": baryon_scale_radius_pc,
        "scale_radius_over_rs": args.scale_radius_over_rs,
        "assembly_time_myr": args.assembly_time_myr,
        "duration_myr": args.duration_myr,
        "eddington_ratio": args.eddington_ratio,
        "duty_cycle": args.duty_cycle,
        "radiative_efficiency": args.radiative_efficiency,
        "eddington_black_hole_efolding_time_myr": (
            model.black_hole_efolding_time_myr
        ),
        "gas_density_msun_pc3": args.gas_density_msun_pc3,
        "gas_sound_speed_km_s": args.gas_sound_speed_km_s,
        "gas_relative_velocity_km_s": args.gas_relative_velocity_km_s,
        "gas_adiabatic_index": args.gas_adiabatic_index,
        "gas_cooling_time_myr": cooling_time_myr,
        "gas_cooling_model": (
            "cloudy_no_uvb"
            if args.cloudy_cooling
            else ("fixed_exponential" if cooling_time_myr is not None else "none")
        ),
        "gas_metallicity_solar": args.gas_metallicity_solar,
        "gas_hydrogen_mass_fraction": args.gas_hydrogen_mass_fraction,
        "cooling_rate_multiplier": args.cooling_rate_multiplier,
        "bondi_alpha": args.bondi_alpha,
        "bondi_enabled": bondi_coefficient is not None,
        "evolve_bondi_ambient": args.evolve_bondi_ambient,
        "post_assembly_bondi_limited_sample_fraction": bondi_limited_fraction,
        "bondi_to_eddington_transition_myr": eddington_transition_time,
        "feedback_model": "hernquist_scale_radius_expansion",
        "feedback_efficiency": args.feedback_efficiency,
        "feedback_heating_fraction": args.feedback_heating_fraction,
        "feedback_expansion_fraction": expansion_energy_fraction,
        "feedback_eta": args.feedback_eta,
        "feedback_binding_model": args.feedback_binding_model,
        "baryon_self_binding_energy_erg": binding_components.self_gravity_erg,
        "baryon_halo_binding_energy_erg": binding_components.halo_erg,
        "baryon_black_hole_binding_energy_erg": binding_components.black_hole_erg,
        "feedback_binding_energy_erg": feedback_binding_energy_erg,
        "sigma_over_m_cm2_g": args.sigma_over_m_cm2_g,
        "cross_section_model": args.cross_section_model,
        "cross_section_velocity_scale_km_s": (
            args.cross_section_velocity_km_s
        ),
        "initial_inner_effective_sigma_over_m_cm2_g": float(
            lmfp_effective_sigma_over_m_cm2_g[0, 0]
        ),
        "final_inner_effective_sigma_over_m_cm2_g": float(
            lmfp_effective_sigma_over_m_cm2_g[-1, 0]
        ),
        "initial_inner_lmfp_effective_sigma_over_m_cm2_g": float(
            lmfp_effective_sigma_over_m_cm2_g[0, 0]
        ),
        "final_inner_lmfp_effective_sigma_over_m_cm2_g": float(
            lmfp_effective_sigma_over_m_cm2_g[-1, 0]
        ),
        "initial_inner_smfp_effective_sigma_over_m_cm2_g": float(
            smfp_effective_sigma_over_m_cm2_g[0, 0]
        ),
        "final_inner_smfp_effective_sigma_over_m_cm2_g": float(
            smfp_effective_sigma_over_m_cm2_g[-1, 0]
        ),
        "dark_capture_model": args.dark_capture_model,
        "dark_bondi_lambda": args.dark_bondi_lambda,
        "dark_flux_capture_rmin_over_influence": (
            args.dark_flux_capture_rmin_over_influence
        ),
        "dark_flux_capture_mass_threshold_msun": (
            dark_flux_capture_mass_threshold_msun
        ),
        "initial_inner_dark_matter_reservoir_msun": initial_dark_reservoir_msun,
        "r_min_pc": r_min_pc,
        "r_min_over_rs": r_min_pc / profile.scale_radius_pc,
        "r_min_over_black_hole_influence_radius": (
            r_min_pc / black_hole_influence_radius_pc
        ),
        "r_max_pc": r_max_pc,
        "r_max_over_rs": r_max_pc / profile.scale_radius_pc,
        "cells": cells,
        "reconstruction": "mc",
        "riemann_solver": "roe",
        "cfl": args.cfl,
        "entropy_fix": args.entropy_fix,
        "steps": result.num_steps,
        "elapsed_seconds": elapsed,
        "final_black_hole_mass_msun": float(black_hole_mass_msun[-1]),
        "dark_matter_accreted_msun": float(dark_accreted_msun[-1]),
        "dark_matter_supplied_to_inner_boundary_msun": float(
            dark_supplied_msun[-1]
        ),
        "final_inner_dark_matter_reservoir_msun": float(
            inner_dark_reservoir_msun[-1]
        ),
        "dark_capture_fraction_of_available_supply": float(
            dark_accreted_msun[-1]
            / (initial_dark_reservoir_msun + dark_supplied_msun[-1])
            if initial_dark_reservoir_msun + dark_supplied_msun[-1] > 0.0
            else 0.0
        ),
        "baryon_accreted_onto_bh_msun": float(baryon_accreted_msun[-1]),
        "baryon_gas_consumed_msun": float(baryon_gas_consumed_msun[-1]),
        "radiated_mass_equivalent_msun": float(
            args.radiative_efficiency * baryon_gas_consumed_msun[-1]
        ),
        "dark_fraction_of_black_hole_growth": float(
            dark_accreted_msun[-1] / total_growth
            if total_growth > 0.0
            else 0.0
        ),
        "dark_dominated_onset_myr": onset,
        "final_baryon_scale_radius_pc": float(
            scales.radius_from_code(result.baryon_scale_radii_code[-1])
        ),
        "final_feedback_to_binding_ratio": float(
            result.feedback_to_binding_ratios[-1]
        ),
        "final_total_feedback_to_binding_ratio": float(
            total_feedback_per_consumed_mass_code
            * result.baryon_gas_consumed_masses_code[-1]
        ),
        "final_ambient_gas_density_msun_pc3": float(
            ambient_gas_density_msun_pc3[-1]
        ),
        "final_ambient_gas_sound_speed_km_s": float(
            ambient_gas_sound_speed_km_s[-1]
        ),
        "final_baryon_thermal_energy_erg": float(
            baryon_thermal_energy_erg[-1]
        ),
        "final_ambient_gas_temperature_k": float(
            ambient_gas_temperature_k[-1]
        ),
        "final_physical_cooling_time_myr": float(
            physical_cooling_time_myr[-1]
        ),
        "minimum_physical_cooling_time_myr": float(
            np.min(physical_cooling_time_myr)
        ),
        "cumulative_injected_thermal_energy_erg": float(
            injected_thermal_energy_erg[-1]
        ),
        "cumulative_cooling_loss_energy_erg": float(
            cooling_loss_energy_erg[-1]
        ),
        "retained_thermal_energy_fraction": float(
            baryon_thermal_energy_erg[-1] / injected_thermal_energy_erg[-1]
            if injected_thermal_energy_erg[-1] > 0.0
            else 0.0
        ),
        "mass_budget_residual_code": result.max_mass_budget_residual_code,
        "time_to_1e4_msun_myr": first_target_time_myr(
            times_myr, black_hole_mass_msun, 1.0e4
        ),
        "time_to_1e5_msun_myr": first_target_time_myr(
            times_myr, black_hole_mass_msun, 1.0e5
        ),
        "time_to_1e6_msun_myr": first_target_time_myr(
            times_myr, black_hole_mass_msun, 1.0e6
        ),
        "time_to_1e7_msun_myr": first_target_time_myr(
            times_myr, black_hole_mass_msun, 1.0e7
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata_json=json.dumps(metadata, sort_keys=True),
        times_myr=times_myr,
        assembled_baryon_fraction=assembled_fraction,
        radii_pc=scales.radius_from_code(grid.centers_code),
        black_hole_mass_msun=black_hole_mass_msun,
        dark_matter_accreted_msun=dark_accreted_msun,
        dark_matter_supplied_to_inner_boundary_msun=dark_supplied_msun,
        inner_dark_matter_reservoir_msun=inner_dark_reservoir_msun,
        baryon_accreted_onto_bh_msun=baryon_accreted_msun,
        baryon_gas_consumed_msun=baryon_gas_consumed_msun,
        dark_matter_accretion_rate_msun_myr=dark_rate_msun_myr,
        dark_matter_supply_rate_msun_myr=dark_supply_rate_msun_myr,
        baryon_accretion_rate_msun_myr=baryon_rate_msun_myr,
        bondi_baryon_growth_limit_msun_myr=retained_bondi_limit_msun_myr,
        eddington_baryon_growth_limit_msun_myr=retained_eddington_limit_msun_myr,
        bondi_limited=bondi_limited,
        baryon_remaining_msun=(
            result.baryon_remaining_masses_code * scales.mass_scale_msun
        ),
        baryon_scale_radius_pc=scales.radius_from_code(
            result.baryon_scale_radii_code
        ),
        feedback_to_binding_ratio=result.feedback_to_binding_ratios,
        ambient_gas_density_msun_pc3=ambient_gas_density_msun_pc3,
        ambient_gas_sound_speed_km_s=ambient_gas_sound_speed_km_s,
        ambient_density_factor=ambient_density_factor,
        ambient_speed_factor=ambient_speed_factor,
        baryon_thermal_energy_erg=baryon_thermal_energy_erg,
        injected_thermal_energy_erg=injected_thermal_energy_erg,
        cooling_loss_energy_erg=cooling_loss_energy_erg,
        ambient_gas_temperature_k=ambient_gas_temperature_k,
        ambient_mean_molecular_weight=ambient_mean_molecular_weight,
        physical_cooling_time_myr=physical_cooling_time_myr,
        cooling_coefficient_erg_cm3_s=cooling_coefficient_erg_cm3_s,
        density_msun_pc3=(
            result.density_snapshots * scales.density_scale_msun_pc3
        ),
        radial_velocity_km_s=(
            result.radial_velocity_snapshots * scales.velocity_scale_km_s
        ),
        velocity_dispersion_km_s=(
            result.velocity_dispersion_snapshots * scales.velocity_scale_km_s
        ),
        effective_sigma_over_m_cm2_g=lmfp_effective_sigma_over_m_cm2_g,
        lmfp_effective_sigma_over_m_cm2_g=(
            lmfp_effective_sigma_over_m_cm2_g
        ),
        smfp_effective_sigma_over_m_cm2_g=(
            smfp_effective_sigma_over_m_cm2_g
        ),
    )
    print(json.dumps(metadata, sort_keys=True))
    print(f"saved={args.output.resolve()}")


if __name__ == "__main__":
    main()
