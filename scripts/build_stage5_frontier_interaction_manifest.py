"""Build the targeted stage-5 concentration/transport/compactness matrix."""

from __future__ import annotations

import csv
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_frontier_interaction.tsv"
CONCENTRATIONS = (6.0, 8.0, 10.0, 12.0)
LOW_CROSS_SECTIONS = (0.3, 1.0, 3.0)
HIGH_CROSS_SECTIONS = (10.0, 30.0, 100.0)
SCALE_RADII = (0.0005, 0.001, 0.003)


def main() -> None:
    rows = []
    # Keep the original low-cross-section task IDs stable, then append the
    # high-cross-section closure points.
    for cross_sections in (LOW_CROSS_SECTIONS, HIGH_CROSS_SECTIONS):
        for concentration, cross_section, scale_radius in product(
            CONCENTRATIONS, cross_sections, SCALE_RADII
        ):
            rows.append(
                {
                    "task_id": len(rows),
                    "halo_concentration": concentration,
                    "sigma_over_m_cm2_g": cross_section,
                    "scale_radius_over_rs": scale_radius,
                    "r_min_over_rs": concentration / 48_000.0,
                    "r_max_over_rs": concentration * (125.0 / 6.0),
                    "cells": 256,
                }
            )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 72:
        raise RuntimeError(f"expected 72 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
