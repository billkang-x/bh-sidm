"""Build the refinement around stage-3 scaling-matrix optima."""

from __future__ import annotations

import csv
from pathlib import Path

from build_stage3_scaling_manifest import (
    BARYON_FRACTIONS,
    CROSS_SECTION_CM2_G,
    HALO_MASSES_MSUN,
    SCALE_RADIUS_OVER_RS,
    halo_numerics,
    predicted_assembly_time_myr,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage3_scaling_refinement.tsv"
MULTIPLIERS = {
    0.01: [0.25, 0.75],
    0.05: [0.75, 1.25],
    0.16: [1.5, 2.5, 3.0],
}


def main() -> None:
    rows = []
    for halo_mass_msun in HALO_MASSES_MSUN:
        predicted_time, supply_radius, r_max_pc, cells = (
            predicted_assembly_time_myr(halo_mass_msun)
        )
        _, expected_r_max, expected_cells = halo_numerics(halo_mass_msun)
        if r_max_pc != expected_r_max or cells != expected_cells:
            raise RuntimeError("inconsistent halo numerics")
        for baryon_fraction in BARYON_FRACTIONS:
            for multiplier in MULTIPLIERS[baryon_fraction]:
                rows.append(
                    {
                        "task_id": len(rows),
                        "case_type": "refinement",
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
    if len(rows) != 21:
        raise RuntimeError(f"expected 21 tasks, generated {len(rows)}")
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
