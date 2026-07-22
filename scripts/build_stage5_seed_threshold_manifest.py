"""Build a seed-mass refinement around the gated-capture transition."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_seed_threshold.tsv"
SEEDS_MSUN = (2.0e4, 3.0e4, 4.0e4, 5.0e4, 6.0e4, 7.0e4, 8.0e4)
BOUNDARIES = (0.05208333333333334, 0.10416666666666667)
CONCENTRATION = 8.0
REFERENCE_SEED_TO_HALO = 1.0e-4
R_MAX_OVER_RS = 166.66666666666667


def logarithmic_cells(r_min_over_rs: float) -> int:
    return int(round(256 * log(R_MAX_OVER_RS / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    rows = []
    for seed in SEEDS_MSUN:
        for boundary in BOUNDARIES:
            r_min = boundary * REFERENCE_SEED_TO_HALO * CONCENTRATION
            rows.append(
                {
                    "task_id": len(rows),
                    "black_hole_seed_msun": seed,
                    "reference_r_min_over_influence": boundary,
                    "seed_r_min_over_influence": boundary * 1.0e5 / seed,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": R_MAX_OVER_RS,
                    "cells": logarithmic_cells(r_min),
                    "dark_bondi_lambda": 0.25,
                    "flux_capture_r_min_over_influence": 0.10416666666666667,
                }
            )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 14:
        raise RuntimeError(f"expected 14 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
