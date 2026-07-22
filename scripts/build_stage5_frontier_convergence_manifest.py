"""Build numerical convergence matrices for three stage-5 frontier cases."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_frontier_convergence.tsv"
MODELS = {
    "representative": {
        "halo_concentration": 8.0,
        "sigma_over_m_cm2_g": 1.0,
        "scale_radius_over_rs": 0.003,
    },
    "resolved_high": {
        "halo_concentration": 12.0,
        "sigma_over_m_cm2_g": 100.0,
        "scale_radius_over_rs": 0.003,
    },
    "extreme": {
        "halo_concentration": 12.0,
        "sigma_over_m_cm2_g": 100.0,
        "scale_radius_over_rs": 0.0005,
    },
}
BASE_RMIN_OVER_INFLUENCE = 5.0 / 24.0
RMIN_OVER_INFLUENCE = (0.02604166666666667, 0.05208333333333334, 0.10416666666666667, BASE_RMIN_OVER_INFLUENCE, 0.4166666666666667)


def logarithmic_cells(r_min_over_rs: float, r_max_over_rs: float) -> int:
    return int(round(256 * log(r_max_over_rs / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    rows = []

    def add(
        model: str,
        axis: str,
        axis_value: float,
        *,
        cells: int = 256,
        r_min_over_influence: float = BASE_RMIN_OVER_INFLUENCE,
        cfl: float = 0.2,
        entropy_fix: float = 0.1,
    ) -> None:
        parameters = MODELS[model]
        concentration = parameters["halo_concentration"]
        r_min_over_rs = r_min_over_influence * 1.0e-4 * concentration
        r_max_over_rs = concentration * (125.0 / 6.0)
        rows.append(
            {
                "task_id": len(rows),
                "model": model,
                "axis": axis,
                "axis_value": axis_value,
                **parameters,
                "r_min_over_influence": r_min_over_influence,
                "r_min_over_rs": r_min_over_rs,
                "r_max_over_rs": r_max_over_rs,
                "cells": cells,
                "cfl": cfl,
                "entropy_fix": entropy_fix,
            }
        )

    for model in MODELS:
        add(model, "baseline", 0.0)
        for cells in (128, 192, 384, 512):
            add(model, "grid", float(cells), cells=cells)
        for ratio in RMIN_OVER_INFLUENCE:
            if ratio == BASE_RMIN_OVER_INFLUENCE:
                continue
            concentration = MODELS[model]["halo_concentration"]
            r_min = ratio * 1.0e-4 * concentration
            r_max = concentration * (125.0 / 6.0)
            add(
                model,
                "inner_boundary",
                ratio,
                cells=logarithmic_cells(r_min, r_max),
                r_min_over_influence=ratio,
            )
        for cfl in (0.1, 0.3):
            add(model, "cfl", cfl, cfl=cfl)
        for entropy_fix in (0.0, 0.05, 0.2):
            add(model, "entropy_fix", entropy_fix, entropy_fix=entropy_fix)

    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 42:
        raise RuntimeError(f"expected 42 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
