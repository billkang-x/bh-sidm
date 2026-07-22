"""Build the stage-3 seed, boundary, and transport similarity check."""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sidm_bh.halos import NFWProfile
from sidm_bh.mesh import SphericalGrid
from sidm_bh.stage3 import static_baryon_equilibrium_state
from sidm_bh.timescales import local_timescale_profiles_code
from sidm_bh.units import SimulationScales


OUTPUT = ROOT / "hpc" / "stage3_similarity_matrix.tsv"
ANCHOR_MASS_MSUN = 1.0e6
ANCHOR_PROFILE = NFWProfile(3.7, 30.0)
HALO_MASSES_MSUN = [1.0e6, 1.0e7, 1.0e8]
PROTOCOLS = [
    "scaled_seed_boundary_fixed_sigma",
    "fully_dimensionless_self_similar",
]
BARYON_FRACTION = 0.05
SCALE_RADIUS_OVER_RS = 0.01
ASSEMBLY_MULTIPLIER = 1.25


def build_rows() -> list[dict[str, float | int | str]]:
    rows = []
    for protocol in PROTOCOLS:
        for halo_mass_msun in HALO_MASSES_MSUN:
            mass_ratio = halo_mass_msun / ANCHOR_MASS_MSUN
            profile = ANCHOR_PROFILE.self_similar_scaled(
                halo_mass_msun,
                ANCHOR_MASS_MSUN,
            )
            radius_factor = profile.scale_radius_pc / ANCHOR_PROFILE.scale_radius_pc
            black_hole_mass_msun = 100.0 * mass_ratio
            r_min_pc = 0.005 * radius_factor
            r_max_pc = 5000.0 * radius_factor
            cells = 256
            sigma = 50.0
            if protocol == "fully_dimensionless_self_similar":
                sigma /= radius_factor
            scales = SimulationScales(
                profile.scale_radius_pc,
                profile.scale_density_msun_pc3,
            )
            grid = SphericalGrid.from_log_spacing(
                scales.radius_to_code(r_min_pc),
                scales.radius_to_code(r_max_pc),
                cells,
            )
            state, _ = static_baryon_equilibrium_state(
                profile,
                grid,
                scales,
                black_hole_mass_msun,
                baryons=None,
            )
            supply_radius_pc = 3.5 * SCALE_RADIUS_OVER_RS * profile.scale_radius_pc
            supply_index = int(
                np.argmin(
                    np.abs(
                        scales.radius_from_code(grid.centers_code)
                        - supply_radius_pc
                    )
                )
            )
            timescales = local_timescale_profiles_code(
                state,
                grid,
                scales.sigma_over_m_to_code(sigma),
                black_hole_mass_code=scales.mass_to_code(
                    black_hole_mass_msun
                ),
            )
            predicted_time_myr = float(
                scales.time_from_code(timescales.dynamical_code[supply_index])
            )
            for case_type, baryon_fraction in (
                ("control", 0.0),
                ("baryon", BARYON_FRACTION),
            ):
                assembly_time_myr = (
                    0.0
                    if case_type == "control"
                    else ASSEMBLY_MULTIPLIER * predicted_time_myr
                )
                rows.append(
                    {
                        "task_id": len(rows),
                        "protocol": protocol,
                        "case_type": case_type,
                        "halo_mass_msun": halo_mass_msun,
                        "black_hole_mass_msun": black_hole_mass_msun,
                        "baryon_fraction": baryon_fraction,
                        "scale_radius_over_rs": SCALE_RADIUS_OVER_RS,
                        "predicted_time_myr": predicted_time_myr,
                        "assembly_time_myr": assembly_time_myr,
                        "sigma_over_m_cm2_g": sigma,
                        "r_min_pc": r_min_pc,
                        "r_max_pc": r_max_pc,
                        "cells": cells,
                    }
                )
    return rows


def main() -> None:
    rows = build_rows()
    if len(rows) != 12:
        raise RuntimeError(f"expected 12 tasks, generated {len(rows)}")
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
