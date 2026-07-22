"""Run one Hernquist assembly-history experiment for stage 3."""

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
from sidm_bh.fast_evolution import evolve_mc_roe_fast
from sidm_bh.halos import NFWProfile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.sources import enclosed_baryon_mass_code
from sidm_bh.stage3 import static_baryon_equilibrium_state
from sidm_bh.units import SimulationScales


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baryon-fraction", type=float, required=True)
    parser.add_argument("--scale-radius-over-rs", type=float, required=True)
    parser.add_argument("--assembly-time-myr", type=float, required=True)
    parser.add_argument("--halo-mass-msun", type=float, default=1.0e6)
    parser.add_argument("--black-hole-mass-msun", type=float, default=100.0)
    parser.add_argument("--r-min-pc", type=float, default=0.005)
    parser.add_argument("--r-max-pc", type=float)
    parser.add_argument("--cells", type=int)
    parser.add_argument("--duration-myr", type=float, default=2.0)
    parser.add_argument("--cfl", type=float, default=0.2)
    parser.add_argument("--entropy-fix", type=float, default=0.1)
    parser.add_argument("--sigma-over-m-cm2-g", type=float, default=50.0)
    parser.add_argument("--samples", type=int, default=101)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=20_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.baryon_fraction <= 1.0:
        raise ValueError("baryon_fraction must lie in [0, 1]")
    if args.scale_radius_over_rs <= 0.0:
        raise ValueError("scale_radius_over_rs must be positive")
    if args.assembly_time_myr < 0.0:
        raise ValueError("assembly_time_myr cannot be negative")
    if args.duration_myr <= 0.0:
        raise ValueError("duration_myr must be positive")
    if args.halo_mass_msun <= 0.0:
        raise ValueError("halo_mass_msun must be positive")
    if args.black_hole_mass_msun < 0.0:
        raise ValueError("black_hole_mass_msun cannot be negative")
    if args.sigma_over_m_cm2_g < 0.0:
        raise ValueError("sigma_over_m_cm2_g cannot be negative")
    if args.samples < 2:
        raise ValueError("samples must be at least two")

    anchor_halo_mass_msun = 1.0e6
    anchor_profile = NFWProfile(3.7, 30.0)
    profile = anchor_profile.self_similar_scaled(
        args.halo_mass_msun,
        anchor_halo_mass_msun,
    )
    scales = SimulationScales(
        profile.scale_radius_pc,
        profile.scale_density_msun_pc3,
    )
    radius_factor = profile.scale_radius_pc / anchor_profile.scale_radius_pc
    r_max_pc = (
        args.r_max_pc if args.r_max_pc is not None else 5000.0 * radius_factor
    )
    if r_max_pc <= args.r_min_pc:
        raise ValueError("r_max_pc must exceed r_min_pc")
    baseline_log_width = np.log(5000.0 / 0.005) / 256.0
    cells = (
        args.cells
        if args.cells is not None
        else int(round(np.log(r_max_pc / args.r_min_pc) / baseline_log_width))
    )
    if cells < 2:
        raise ValueError("cells must be at least two")
    grid = SphericalGrid.from_log_spacing(
        scales.radius_to_code(args.r_min_pc),
        scales.radius_to_code(r_max_pc),
        cells,
    )
    initial_state, _ = static_baryon_equilibrium_state(
        profile,
        grid,
        scales,
        args.black_hole_mass_msun,
        baryons=None,
    )
    baryons = HernquistBaryons(
        total_mass_msun=args.baryon_fraction * args.halo_mass_msun,
        scale_radius_pc=args.scale_radius_over_rs * profile.scale_radius_pc,
    )
    full_baryon_mass_code = enclosed_baryon_mass_code(baryons, grid, scales)
    times_myr = np.linspace(0.0, args.duration_myr, args.samples)
    assembly_time_code = (
        None
        if args.assembly_time_myr == 0.0
        else scales.time_to_code(args.assembly_time_myr)
    )

    started = perf_counter()
    result = evolve_mc_roe_fast(
        initial_state,
        grid,
        scales.time_to_code(times_myr),
        scales.mass_to_code(args.black_hole_mass_msun),
        scales.sigma_over_m_to_code(args.sigma_over_m_cm2_g),
        cfl_number=args.cfl,
        entropy_fix=args.entropy_fix,
        baryon_enclosed_mass_code=full_baryon_mass_code,
        baryon_assembly_time_code=assembly_time_code,
        source_integration="euler",
        max_steps=args.max_steps,
    )
    elapsed = perf_counter() - started
    final_mass_msun = scales.mass_from_code(result.final_black_hole_mass_code)
    if args.baryon_fraction == 0.0:
        mass_fractions = np.zeros_like(times_myr)
        protocol = "no_baryon_control"
    elif args.assembly_time_myr == 0.0:
        mass_fractions = np.ones_like(times_myr)
        protocol = "instantaneous_turn_on"
    else:
        mass_fractions = np.array(
            [
                smoothstep_mass_fraction(time, args.assembly_time_myr)
                for time in times_myr
            ]
        )
        protocol = "finite_smoothstep_assembly"
    metadata = {
        "profile": "nfw",
        "initial_condition": "hydrostatic_dm_bh_no_baryon",
        "baryon_protocol": protocol,
        "halo_mass_msun": args.halo_mass_msun,
        "halo_scaling": "self_similar_fixed_scale_density_and_concentration",
        "nfw_scale_density_msun_pc3": profile.scale_density_msun_pc3,
        "nfw_scale_radius_pc": profile.scale_radius_pc,
        "nfw_concentration": anchor_profile.concentration_for_enclosed_mass(
            anchor_halo_mass_msun
        ),
        "nfw_virial_radius_pc": (
            anchor_profile.concentration_for_enclosed_mass(anchor_halo_mass_msun)
            * profile.scale_radius_pc
        ),
        "black_hole_seed_msun": args.black_hole_mass_msun,
        "baryon_fraction": args.baryon_fraction,
        "baryon_mass_msun": args.baryon_fraction * args.halo_mass_msun,
        "scale_radius_over_rs": args.scale_radius_over_rs,
        "scale_radius_pc": args.scale_radius_over_rs * profile.scale_radius_pc,
        "assembly_time_myr": args.assembly_time_myr,
        "duration_myr": args.duration_myr,
        "r_min_pc": args.r_min_pc,
        "r_max_pc": r_max_pc,
        "cells": cells,
        "reconstruction": "mc",
        "riemann_solver": "roe",
        "cfl": args.cfl,
        "entropy_fix": args.entropy_fix,
        "conduction": (
            "implicit" if args.sigma_over_m_cm2_g > 0.0 else "disabled"
        ),
        "sigma_over_m_cm2_g": args.sigma_over_m_cm2_g,
        "steps": result.num_steps,
        "elapsed_seconds": elapsed,
        "final_black_hole_mass_msun": final_mass_msun,
        "accreted_dark_matter_msun": (
            final_mass_msun - args.black_hole_mass_msun
        ),
        "peak_accretion_rate_msun_myr": (
            result.peak_accretion_rate_code
            * scales.mass_scale_msun
            / scales.time_scale_myr
        ),
        "final_accretion_rate_msun_myr": (
            result.final_accretion_rate_code
            * scales.mass_scale_msun
            / scales.time_scale_myr
        ),
        "mass_budget_residual_code": result.max_mass_budget_residual_code,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata_json=json.dumps(metadata, sort_keys=True),
        times_myr=times_myr,
        baryon_mass_fraction=mass_fractions,
        radii_pc=scales.radius_from_code(grid.centers_code),
        full_baryon_enclosed_mass_msun=(
            full_baryon_mass_code * scales.mass_scale_msun
        ),
        black_hole_mass_msun=(
            result.black_hole_masses_code * scales.mass_scale_msun
        ),
        accretion_rate_msun_myr=(
            result.accretion_rates_code
            * scales.mass_scale_msun
            / scales.time_scale_myr
        ),
        density_msun_pc3=(
            result.density_snapshots * scales.density_scale_msun_pc3
        ),
        radial_velocity_km_s=(
            result.radial_velocity_snapshots * scales.velocity_scale_km_s
        ),
        velocity_dispersion_km_s=(
            result.velocity_dispersion_snapshots * scales.velocity_scale_km_s
        ),
    )
    print(json.dumps(metadata, sort_keys=True))
    print(f"saved={args.output.resolve()}")


if __name__ == "__main__":
    main()
