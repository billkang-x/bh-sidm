"""Build the final seed-threshold and grid-convergence refinement."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_seed_threshold_final.tsv"
SEEDS_MSUN = (3.25e4, 3.5e4, 3.75e4)
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
                    "axis": "inner_boundary",
                    "reference_r_min_over_influence": boundary,
                    "seed_r_min_over_influence": boundary * 1.0e5 / seed,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": R_MAX_OVER_RS,
                    "cells": logarithmic_cells(r_min),
                    "dark_bondi_lambda": 0.25,
                    "flux_capture_r_min_over_influence": 0.10416666666666667,
                }
            )
        boundary = 0.10416666666666667
        r_min = boundary * REFERENCE_SEED_TO_HALO * CONCENTRATION
        rows.append(
            {
                "task_id": len(rows),
                "black_hole_seed_msun": seed,
                "axis": "grid",
                "reference_r_min_over_influence": boundary,
                "seed_r_min_over_influence": boundary * 1.0e5 / seed,
                "r_min_over_rs": r_min,
                "r_max_over_rs": R_MAX_OVER_RS,
                "cells": 2 * logarithmic_cells(r_min),
                "dark_bondi_lambda": 0.25,
                "flux_capture_r_min_over_influence": 0.10416666666666667,
            }
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 9:
        raise RuntimeError(f"expected 9 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
