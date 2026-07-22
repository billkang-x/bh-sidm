"""Build the resolved-seed calibration matrix for the dark Bondi reservoir."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_light_seed_validation.tsv"
CONCENTRATION = 8.0
REFERENCE_SEED_TO_HALO = 1.0e-4
R_MAX_OVER_RS = 166.66666666666667
BOUNDARIES = (0.05208333333333334, 0.10416666666666667, 5.0 / 24.0)


def logarithmic_cells(r_min_over_rs: float) -> int:
    return int(round(256 * log(R_MAX_OVER_RS / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    rows = []

    def add(reference_boundary: float, bondi_lambda: float) -> None:
        r_min = reference_boundary * REFERENCE_SEED_TO_HALO * CONCENTRATION
        rows.append(
            {
                "task_id": len(rows),
                "reference_r_min_over_influence": reference_boundary,
                "dark_bondi_lambda": bondi_lambda,
                "r_min_over_rs": r_min,
                "r_max_over_rs": R_MAX_OVER_RS,
                "cells": logarithmic_cells(r_min),
            }
        )

    for boundary in BOUNDARIES:
        add(boundary, 0.25)
    for bondi_lambda in (0.20, 0.30):
        add(0.10416666666666667, bondi_lambda)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 5:
        raise RuntimeError(f"expected 5 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
