"""Run and save one MC-Roe stage-2 baseline case."""

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

from sidm_bh.fast_evolution import evolve_mc_roe_fast
from sidm_bh.halos import NFWProfile, SingularIsothermalSphere
from sidm_bh.initial_conditions import (
    hydrostatic_state_from_profile,
    isothermal_state_from_profile,
)
from sidm_bh.mesh import SphericalGrid
from sidm_bh.units import SimulationScales


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("nfw", "sis"), required=True)
    parser.add_argument("--heat", choices=("on", "off"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--r-min-pc", type=float)
    parser.add_argument("--r-max-pc", type=float)
    parser.add_argument("--cells", type=int)
    parser.add_argument("--cfl", type=float)
    parser.add_argument("--entropy-fix", type=float, default=0.1)
    parser.add_argument("--samples", type=int, default=201)
    parser.add_argument("--max-steps", type=int, default=50_000_000)
    return parser.parse_args()


def build_case(args: argparse.Namespace):
    if args.profile == "nfw":
        profile = NFWProfile(3.7, 30.0)
        scales = SimulationScales(30.0, 3.7)
        r_min_pc = args.r_min_pc if args.r_min_pc is not None else 0.005
        r_max_pc = args.r_max_pc if args.r_max_pc is not None else 5000.0
        cells = args.cells if args.cells is not None else 256
        cfl = args.cfl if args.cfl is not None else 0.2
        grid = SphericalGrid.from_log_spacing(
            scales.radius_to_code(r_min_pc),
            scales.radius_to_code(r_max_pc),
            cells,
        )
        state = hydrostatic_state_from_profile(profile, grid, scales)
    else:
        profile = SingularIsothermalSphere(4.2)
        scales = SimulationScales.for_singular_isothermal_sphere(4.2)
        r_min_pc = args.r_min_pc if args.r_min_pc is not None else 0.001
        r_max_pc = args.r_max_pc if args.r_max_pc is not None else 1000.0
        cells = args.cells if args.cells is not None else 128
        cfl = args.cfl if args.cfl is not None else 0.8
        grid = SphericalGrid.from_log_spacing(
            scales.radius_to_code(r_min_pc),
            scales.radius_to_code(r_max_pc),
            cells,
        )
        state = isothermal_state_from_profile(profile, grid, scales, 4.2)
    return state, grid, scales, r_min_pc, r_max_pc, cells, cfl


def main() -> None:
    args = parse_args()
    if args.samples < 2:
        raise ValueError("samples must be at least two")
    state, grid, scales, r_min_pc, r_max_pc, cells, cfl = build_case(args)
    times_myr = np.linspace(0.0, 2.0, args.samples)
    sigma_code = scales.sigma_over_m_to_code(50.0) if args.heat == "on" else 0.0
    started = perf_counter()
    result = evolve_mc_roe_fast(
        state,
        grid,
        scales.time_to_code(times_myr),
        scales.mass_to_code(100.0),
        sigma_code,
        cfl_number=cfl,
        entropy_fix=args.entropy_fix,
        source_integration="euler",
        max_steps=args.max_steps,
    )
    elapsed = perf_counter() - started
    metadata = {
        "profile": args.profile,
        "heat": args.heat,
        "r_min_pc": r_min_pc,
        "r_max_pc": r_max_pc,
        "cells": cells,
        "cfl": cfl,
        "entropy_fix": args.entropy_fix,
        "sigma_over_m_cm2_g": 50.0 if args.heat == "on" else 0.0,
        "steps": result.num_steps,
        "elapsed_seconds": elapsed,
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
