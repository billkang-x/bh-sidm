"""Analyze final convergence for the trusted stage-5 compactness peak."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_trusted_peak_convergence.tsv"
RESULTS = ROOT / "results" / "stage5" / "trusted_peak_convergence"
BASELINE = ROOT / "results" / "stage5" / "compactness_peak" / "task_002.npz"
SUMMARY = ROOT / "results" / "stage5" / "trusted_peak_convergence_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "trusted_peak_convergence_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_trusted_peak_convergence.png"


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with np.load(BASELINE, allow_pickle=False) as data:
        baseline_metadata = json.loads(str(data["metadata_json"]))
    baseline_mass = float(baseline_metadata["final_black_hole_mass_msun"])
    baseline_time_to_1e7 = float(baseline_metadata["time_to_1e7_msun_myr"])
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        cases.append(
            {
                "task_id": task_id,
                "axis": row["axis"],
                "axis_value": float(row["axis_value"]),
                "r_min_over_influence": float(
                    metadata["r_min_over_black_hole_influence_radius"]
                ),
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "time_to_1e7_msun_myr": float(
                    metadata["time_to_1e7_msun_myr"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    boundary = sorted(
        (case for case in cases if case["axis"] == "inner_boundary"),
        key=lambda case: case["axis_value"],
    )
    grid = next(case for case in cases if case["axis"] == "grid")
    boundary_difference = relative_difference(
        boundary[0]["final_black_hole_mass_msun"],
        boundary[1]["final_black_hole_mass_msun"],
    )
    grid_difference = relative_difference(
        baseline_mass, grid["final_black_hole_mass_msun"]
    )
    statistics = {
        "case_count": len(cases),
        "baseline_mass_msun": baseline_mass,
        "baseline_time_to_1e7_myr": baseline_time_to_1e7,
        "smallest_boundary_mass_msun": boundary[0]["final_black_hole_mass_msun"],
        "smallest_boundary_time_to_1e7_myr": boundary[0]["time_to_1e7_msun_myr"],
        "smallest_boundary_pair_relative_difference": boundary_difference,
        "grid_256_to_512_relative_difference": grid_difference,
        "passes_five_percent_convergence": bool(
            boundary_difference < 0.05 and grid_difference < 0.05
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    axes[0].plot(
        [case["axis_value"] for case in boundary] + [5.0 / 24.0],
        [case["final_black_hole_mass_msun"] for case in boundary] + [baseline_mass],
        color="#1b9e77", marker="o",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"$r_{\min}/r_{\rm infl}$")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].set_title("Fiducial model: inner-boundary sensitivity")
    axes[1].plot(
        [256, 512],
        [baseline_mass, grid["final_black_hole_mass_msun"]],
        color="#1b9e77", marker="o",
    )
    axes[1].set_xlabel("Cells")
    axes[1].set_ylabel("Final black-hole mass [M_sun]")
    axes[1].set_title("Fiducial model: grid sensitivity")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
