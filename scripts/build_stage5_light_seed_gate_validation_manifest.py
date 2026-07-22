"""Build validation cases for the influence-gated dark capture closure."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_light_seed_gate_validation.tsv"
CONCENTRATION = 8.0
REFERENCE_SEED_TO_HALO = 1.0e-4
R_MAX_OVER_RS = 166.66666666666667


def logarithmic_cells(r_min_over_rs: float) -> int:
    return int(round(256 * log(R_MAX_OVER_RS / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    rows = []

    def add(feeding_boundary: float, flux_capture_boundary: float) -> None:
        r_min = feeding_boundary * REFERENCE_SEED_TO_HALO * CONCENTRATION
        rows.append(
            {
                "task_id": len(rows),
                "reference_r_min_over_influence": feeding_boundary,
                "flux_capture_r_min_over_influence": flux_capture_boundary,
                "dark_bondi_lambda": 0.25,
                "r_min_over_rs": r_min,
                "r_max_over_rs": R_MAX_OVER_RS,
                "cells": logarithmic_cells(r_min),
            }
        )

    for boundary in (0.05208333333333334, 0.10416666666666667, 5.0 / 24.0):
        add(boundary, 0.10416666666666667)
    for switch in (0.05208333333333334, 5.0 / 24.0):
        add(0.10416666666666667, switch)

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
