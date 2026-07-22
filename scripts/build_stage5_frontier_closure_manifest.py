"""Build the stage-5 one-axis physical-frontier closure matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_frontier_closure.tsv"
BASE = {
    "halo_concentration": 8.0,
    "sigma_over_m_cm2_g": 10.0,
    "baryon_fraction": 0.16,
    "scale_radius_over_rs": 0.003,
}
CONCENTRATIONS = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
CROSS_SECTIONS = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
SCALE_RADII = (0.0005, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.03)
BARYON_FRACTIONS = (0.01, 0.03, 0.05, 0.10, 0.16)


def boundary_ratios(concentration: float) -> tuple[float, float]:
    # Keep the physical domain fixed when concentration changes.
    return concentration / 48_000.0, concentration * (125.0 / 6.0)


def main() -> None:
    rows: list[dict] = []
    seen: set[tuple[float, float, float, float]] = set()

    def add(axis: str, axis_value: float, **updates: float) -> None:
        parameters = {**BASE, **updates}
        key = tuple(parameters.values())
        if key in seen:
            return
        seen.add(key)
        r_min, r_max = boundary_ratios(parameters["halo_concentration"])
        rows.append(
            {
                "task_id": len(rows),
                "axis": axis,
                "axis_value": axis_value,
                **parameters,
                "r_min_over_rs": r_min,
                "r_max_over_rs": r_max,
                "cells": 256,
            }
        )

    add("baseline", 0.0)
    for value in CONCENTRATIONS:
        add("concentration", value, halo_concentration=value)
    for value in CROSS_SECTIONS:
        add("cross_section", value, sigma_over_m_cm2_g=value)
    for value in SCALE_RADII:
        add("scale_radius", value, scale_radius_over_rs=value)
    for value in BARYON_FRACTIONS:
        add("baryon_fraction", value, baryon_fraction=value)

    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 26:
        raise RuntimeError(f"expected 26 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
