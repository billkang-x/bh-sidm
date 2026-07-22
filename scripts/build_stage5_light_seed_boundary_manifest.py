"""Build the fixed-halo light-seed and feeding-boundary closure matrix."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_light_seed_boundary.tsv"
SEEDS_MSUN = (1.0e2, 1.0e3, 1.0e4, 1.0e5)
CONCENTRATION = 8.0
REFERENCE_SEED_MSUN = 1.0e5
HALO_MASS_MSUN = 1.0e9
R_MAX_OVER_RS = 166.66666666666667
BOUNDARIES = (0.05208333333333334, 0.10416666666666667, 5.0 / 24.0)


def logarithmic_cells(r_min_over_rs: float) -> int:
    return int(round(256 * log(R_MAX_OVER_RS / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    rows = []
    for seed in SEEDS_MSUN:
        for boundary in BOUNDARIES:
            r_min = boundary * REFERENCE_SEED_MSUN / HALO_MASS_MSUN * CONCENTRATION
            rows.append(
                {
                    "task_id": len(rows),
                    "black_hole_seed_msun": seed,
                    "axis": "inner_boundary",
                    "reference_r_min_over_influence": boundary,
                    "seed_r_min_over_influence": boundary * REFERENCE_SEED_MSUN / seed,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": R_MAX_OVER_RS,
                    "cells": logarithmic_cells(r_min),
                    "dark_bondi_lambda": 0.25,
                    "flux_capture_r_min_over_influence": 0.10416666666666667,
                }
            )
        boundary = 0.10416666666666667
        r_min = boundary * REFERENCE_SEED_MSUN / HALO_MASS_MSUN * CONCENTRATION
        rows.append(
            {
                "task_id": len(rows),
                "black_hole_seed_msun": seed,
                "axis": "grid",
                "reference_r_min_over_influence": boundary,
                "seed_r_min_over_influence": boundary * REFERENCE_SEED_MSUN / seed,
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
    if len(rows) != 16:
        raise RuntimeError(f"expected 16 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
