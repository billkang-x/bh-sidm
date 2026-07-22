"""Run one matched static-Hernquist stage-3 experiment."""

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

from sidm_bh.baryons import HernquistBaryons
from sidm_bh.fast_evolution import evolve_mc_roe_fast
from sidm_bh.halos import NFWProfile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.stage3 import static_baryon_equilibrium_state
from sidm_bh.units import SimulationScales


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baryon-fraction", type=float, required=True)
    parser.add_argument("--scale-radius-over-rs", type=float, required=True)
    parser.add_argument("--r-min-pc", type=float, default=0.005)
    parser.add_argument("--r-max-pc", type=float, default=5000.0)
    parser.add_argument("--cells", type=int, default=256)
    parser.add_argument("--cfl", type=float, default=0.2)
    parser.add_argument("--entropy-fix", type=float, default=0.1)
    parser.add_argument("--samples", type=int, default=101)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-steps", type=int, default=20_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.baryon_fraction < 0.0 or args.baryon_fraction > 1.0:
        raise ValueError("baryon_fraction must lie in [0, 1]")
    if args.scale_radius_over_rs <= 0.0:
        raise ValueError("scale_radius_over_rs must be positive")
    if args.samples < 2:
        raise ValueError("samples must be at least two")

    halo_mass_msun = 1.0e6
    black_hole_mass_msun = 100.0
    scale_radius_pc = 30.0
    profile = NFWProfile(3.7, scale_radius_pc)
    scales = SimulationScales(scale_radius_pc, 3.7)
    grid = SphericalGrid.from_log_spacing(
        scales.radius_to_code(args.r_min_pc),
        scales.radius_to_code(args.r_max_pc),
        args.cells,
    )
    baryons = None
    if args.baryon_fraction > 0.0:
        baryons = HernquistBaryons(
            total_mass_msun=args.baryon_fraction * halo_mass_msun,
            scale_radius_pc=args.scale_radius_over_rs * scale_radius_pc,
        )
    state, baryon_mass_code = static_baryon_equilibrium_state(
        profile,
        grid,
        scales,
        black_hole_mass_msun,
        baryons=baryons,
    )
    times_myr = np.linspace(0.0, 2.0, args.samples)
    started = perf_counter()
    result = evolve_mc_roe_fast(
        state,
        grid,
        scales.time_to_code(times_myr),
        scales.mass_to_code(black_hole_mass_msun),
        scales.sigma_over_m_to_code(50.0),
        cfl_number=args.cfl,
        entropy_fix=args.entropy_fix,
        baryon_enclosed_mass_code=baryon_mass_code,
        source_integration="euler",
        max_steps=args.max_steps,
    )
    elapsed = perf_counter() - started
    final_mass_msun = scales.mass_from_code(result.final_black_hole_mass_code)
    metadata = {
        "profile": "nfw",
        "initial_condition": "hydrostatic_dm_bh_static_baryon",
        "halo_mass_msun": halo_mass_msun,
        "black_hole_seed_msun": black_hole_mass_msun,
        "baryon_fraction": args.baryon_fraction,
        "baryon_mass_msun": args.baryon_fraction * halo_mass_msun,
        "scale_radius_over_rs": args.scale_radius_over_rs,
        "scale_radius_pc": args.scale_radius_over_rs * scale_radius_pc,
        "r_min_pc": args.r_min_pc,
        "r_max_pc": args.r_max_pc,
        "cells": args.cells,
        "cfl": args.cfl,
        "entropy_fix": args.entropy_fix,
        "steps": result.num_steps,
        "elapsed_seconds": elapsed,
        "final_black_hole_mass_msun": final_mass_msun,
        "accreted_dark_matter_msun": final_mass_msun - black_hole_mass_msun,
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
        radii_pc=scales.radius_from_code(grid.centers_code),
        baryon_enclosed_mass_msun=baryon_mass_code * scales.mass_scale_msun,
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
