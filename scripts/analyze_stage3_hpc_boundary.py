"""Summarize boundary checks for selected stage-3 HPC matrix points."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage3"
MAIN = RESULTS / "hpc_matrix"
BOUNDARY = RESULTS / "hpc_boundary"
RADII = [0.0025, 0.005, 0.01]
CASES = [
    {
        "label": "Global maximum",
        "sigma_over_m_cm2_g": 10.0,
        "scale_radius_over_rs": 0.01,
        "assembly_time_myr": 0.5,
        "main_control": 0,
        "main_science": 7,
        "boundary_controls": [0, 1],
        "boundary_science": [4, 5],
    },
    {
        "label": "High-sigma compact",
        "sigma_over_m_cm2_g": 100.0,
        "scale_radius_over_rs": 0.01,
        "assembly_time_myr": 1.0,
        "main_control": 3,
        "main_science": 128,
        "boundary_controls": [2, 3],
        "boundary_science": [6, 7],
    },
    {
        "label": "Weak extended",
        "sigma_over_m_cm2_g": 100.0,
        "scale_radius_over_rs": 0.3,
        "assembly_time_myr": 1.0,
        "main_control": 3,
        "main_science": 163,
        "boundary_controls": [2, 3],
        "boundary_science": [8, 9],
    },
]


def growth(path: Path) -> float:
    with np.load(path, allow_pickle=False) as data:
        return float(data["black_hole_mass_msun"][-1] - 100.0)


def main() -> None:
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for axis, case in zip(axes, CASES, strict=True):
        factors = [
            growth(BOUNDARY / f"task_b{case['boundary_science'][0]:02d}.npz")
            / growth(BOUNDARY / f"task_b{case['boundary_controls'][0]:02d}.npz"),
            growth(MAIN / f"task_{case['main_science']:03d}.npz")
            / growth(MAIN / f"task_{case['main_control']:03d}.npz"),
            growth(BOUNDARY / f"task_b{case['boundary_science'][1]:02d}.npz")
            / growth(BOUNDARY / f"task_b{case['boundary_controls'][1]:02d}.npz"),
        ]
        for radius, factor in zip(RADII, factors, strict=True):
            rows.append(
                {
                    "case": case["label"],
                    "sigma_over_m_cm2_g": case["sigma_over_m_cm2_g"],
                    "scale_radius_over_rs": case["scale_radius_over_rs"],
                    "assembly_time_myr": case["assembly_time_myr"],
                    "r_min_pc": radius,
                    "enhancement_over_matched_control": factor,
                }
            )
        axis.plot(RADII, factors, marker="o", color="#d62828")
        axis.axhline(1.0, color="#6c757d", linestyle=":")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(case["label"])
        axis.set_xlabel("Inner boundary [pc]")
        axis.set_ylabel("Dark-growth enhancement")
        axis.grid(alpha=0.25)
    fig.savefig(RESULTS / "figures" / "stage3_hpc_boundary_followup.png", dpi=180)
    plt.close(fig)

    with (RESULTS / "hpc_boundary_followup_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(RESULTS / "hpc_boundary_followup_summary.csv")


if __name__ == "__main__":
    main()
