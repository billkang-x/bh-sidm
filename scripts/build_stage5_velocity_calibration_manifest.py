"""Build the velocity-dependent SIDM microphysics calibration matrix."""

from __future__ import annotations

import csv
from math import log
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_velocity_calibration.tsv"
HALO_MASS_MSUN = 1.0e9
CONCENTRATION = 8.0
REFERENCE_SEED_MSUN = 1.0e5
R_MAX_OVER_RS = CONCENTRATION * (125.0 / 6.0)
BOUNDARIES = (5.0 / 24.0, 5.0 / 48.0)
SIGMA0_VALUES = (1.0, 3.0, 10.0, 30.0, 100.0)
VELOCITY_SCALES_KM_S = (10.0, 30.0, 100.0, 300.0)


def logarithmic_cells(r_min_over_rs: float) -> int:
    return int(round(256 * log(R_MAX_OVER_RS / r_min_over_rs) / log(1.0e6)))


def main() -> None:
    models = [
        ("constant_sigma1", "constant", 1.0, 0.0),
        ("constant_sigma10", "constant", 10.0, 0.0),
        ("constant_sigma100", "constant", 100.0, 0.0),
    ]
    models.extend(
        (
            f"rutherford_sigma{sigma0:g}_w{velocity:g}",
            "rutherford",
            sigma0,
            velocity,
        )
        for sigma0 in SIGMA0_VALUES
        for velocity in VELOCITY_SCALES_KM_S
    )

    rows = []
    for label, model, sigma0, velocity in models:
        for boundary in BOUNDARIES:
            r_min = (
                boundary
                * REFERENCE_SEED_MSUN
                / HALO_MASS_MSUN
                * CONCENTRATION
            )
            rows.append(
                {
                    "task_id": len(rows),
                    "model_label": label,
                    "cross_section_model": model,
                    "sigma0_over_m_cm2_g": sigma0,
                    "velocity_scale_km_s": velocity,
                    "r_min_over_reference_influence": boundary,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": R_MAX_OVER_RS,
                    "cells": logarithmic_cells(r_min),
                }
            )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 46:
        raise RuntimeError(f"expected 46 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
