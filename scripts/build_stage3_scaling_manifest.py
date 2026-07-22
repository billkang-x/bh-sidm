"""Build the stage-3 baryon-fraction and halo-mass scaling experiment."""

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


OUTPUT = ROOT / "hpc" / "stage3_scaling_matrix.tsv"
ANCHOR_HALO_MASS_MSUN = 1.0e6
ANCHOR_PROFILE = NFWProfile(3.7, 30.0)
HALO_MASSES_MSUN = [1.0e6, 1.0e7, 1.0e8]
BARYON_FRACTIONS = [0.01, 0.05, 0.16]
ASSEMBLY_MULTIPLIERS = [0.0, 0.5, 1.0, 2.0]
SCALE_RADIUS_OVER_RS = 0.01
SUPPLY_RADIUS_OVER_AB = 3.5
CROSS_SECTION_CM2_G = 50.0
BLACK_HOLE_MASS_MSUN = 100.0
R_MIN_PC = 0.005


def halo_numerics(halo_mass_msun: float) -> tuple[NFWProfile, float, int]:
    profile = ANCHOR_PROFILE.self_similar_scaled(
        halo_mass_msun,
        ANCHOR_HALO_MASS_MSUN,
    )
    radius_factor = profile.scale_radius_pc / ANCHOR_PROFILE.scale_radius_pc
    r_max_pc = 5000.0 * radius_factor
    baseline_log_width = np.log(5000.0 / R_MIN_PC) / 256.0
    cells = int(round(np.log(r_max_pc / R_MIN_PC) / baseline_log_width))
    return profile, r_max_pc, cells


def predicted_assembly_time_myr(
    halo_mass_msun: float,
) -> tuple[float, float, float, int]:
    profile, r_max_pc, cells = halo_numerics(halo_mass_msun)
    scales = SimulationScales(
        profile.scale_radius_pc,
        profile.scale_density_msun_pc3,
    )
    grid = SphericalGrid.from_log_spacing(
        scales.radius_to_code(R_MIN_PC),
        scales.radius_to_code(r_max_pc),
        cells,
    )
    state, _ = static_baryon_equilibrium_state(
        profile,
        grid,
        scales,
        BLACK_HOLE_MASS_MSUN,
        baryons=None,
    )
    supply_radius_pc = (
        SUPPLY_RADIUS_OVER_AB
        * SCALE_RADIUS_OVER_RS
        * profile.scale_radius_pc
    )
    supply_index = int(
        np.argmin(
            np.abs(
                scales.radius_from_code(grid.centers_code) - supply_radius_pc
            )
        )
    )
    timescales = local_timescale_profiles_code(
        state,
        grid,
        scales.sigma_over_m_to_code(CROSS_SECTION_CM2_G),
        black_hole_mass_code=scales.mass_to_code(BLACK_HOLE_MASS_MSUN),
    )
    sampled_radius_pc = float(
        scales.radius_from_code(grid.centers_code[supply_index])
    )
    predicted_time_myr = float(
        scales.time_from_code(timescales.dynamical_code[supply_index])
    )
    return predicted_time_myr, sampled_radius_pc, r_max_pc, cells


def build_rows() -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for halo_mass_msun in HALO_MASSES_MSUN:
        predicted_time, supply_radius, r_max_pc, cells = (
            predicted_assembly_time_myr(halo_mass_msun)
        )
        rows.append(
            {
                "task_id": len(rows),
                "case_type": "control",
                "halo_mass_msun": halo_mass_msun,
                "baryon_fraction": 0.0,
                "scale_radius_over_rs": SCALE_RADIUS_OVER_RS,
                "supply_radius_pc": supply_radius,
                "predicted_time_myr": predicted_time,
                "assembly_multiplier": 0.0,
                "assembly_time_myr": 0.0,
                "sigma_over_m_cm2_g": CROSS_SECTION_CM2_G,
                "r_max_pc": r_max_pc,
                "cells": cells,
            }
        )
        for baryon_fraction in BARYON_FRACTIONS:
            for multiplier in ASSEMBLY_MULTIPLIERS:
                rows.append(
                    {
                        "task_id": len(rows),
                        "case_type": "baryon",
                        "halo_mass_msun": halo_mass_msun,
                        "baryon_fraction": baryon_fraction,
                        "scale_radius_over_rs": SCALE_RADIUS_OVER_RS,
                        "supply_radius_pc": supply_radius,
                        "predicted_time_myr": predicted_time,
                        "assembly_multiplier": multiplier,
                        "assembly_time_myr": multiplier * predicted_time,
                        "sigma_over_m_cm2_g": CROSS_SECTION_CM2_G,
                        "r_max_pc": r_max_pc,
                        "cells": cells,
                    }
                )
    return rows


def main() -> None:
    rows = build_rows()
    if len(rows) != 39:
        raise RuntimeError(f"expected 39 tasks, generated {len(rows)}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    for halo_mass_msun in HALO_MASSES_MSUN:
        row = next(row for row in rows if row["halo_mass_msun"] == halo_mass_msun)
        print(
            f"M={halo_mass_msun:.0e} Msun: "
            f"r_supply={row['supply_radius_pc']:.6g} pc, "
            f"T_pred={row['predicted_time_myr']:.6g} Myr, "
            f"r_max={row['r_max_pc']:.6g} pc, cells={row['cells']}"
        )
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
