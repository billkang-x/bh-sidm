"""Build targeted convergence checks for moderate-sigma high-c cases."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_moderate_convergence.tsv"


def main() -> None:
    rows = []
    for concentration in (10.0, 12.0):
        r_max = concentration * (125.0 / 6.0)
        for ratio in (0.05208333333333334, 0.10416666666666667):
            r_min = ratio * 1.0e-4 * concentration
            cells = int(round(256 * log(r_max / r_min) / log(1.0e6)))
            rows.append(
                {
                    "task_id": len(rows),
                    "halo_concentration": concentration,
                    "axis": "inner_boundary",
                    "axis_value": ratio,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": r_max,
                    "cells": cells,
                }
            )
        rows.append(
            {
                "task_id": len(rows),
                "halo_concentration": concentration,
                "axis": "grid",
                "axis_value": 512.0,
                "r_min_over_rs": (5.0 / 24.0) * 1.0e-4 * concentration,
                "r_max_over_rs": r_max,
                "cells": 512,
            }
        )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 6:
        raise RuntimeError(f"expected 6 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
