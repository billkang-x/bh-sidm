"""Analyze targeted convergence for moderate-sigma high-c cases."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_moderate_convergence.tsv"
RESULTS = ROOT / "results" / "stage5" / "moderate_convergence"
INTERACTION_RESULTS = ROOT / "results" / "stage5" / "frontier_interaction"
SUMMARY = ROOT / "results" / "stage5" / "moderate_convergence_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "moderate_convergence_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_moderate_convergence.png"
BASELINE_TASKS = {10.0: 26, 12.0: 35}


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
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
                "halo_concentration": float(row["halo_concentration"]),
                "axis": row["axis"],
                "axis_value": float(row["axis_value"]),
                "r_min_over_influence": float(
                    metadata["r_min_over_black_hole_influence_radius"]
                ),
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )
    baselines = {}
    for concentration, task_id in BASELINE_TASKS.items():
        with np.load(
            INTERACTION_RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False
        ) as data:
            metadata = json.loads(str(data["metadata_json"]))
        baselines[concentration] = float(metadata["final_black_hole_mass_msun"])
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    statistics = {"case_count": len(cases), "convergence_by_concentration": {}}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    colors = {10.0: "#1b9e77", 12.0: "#b2182b"}
    for index, concentration in enumerate((10.0, 12.0)):
        boundary = sorted(
            (
                case
                for case in cases
                if case["halo_concentration"] == concentration
                and case["axis"] == "inner_boundary"
            ),
            key=lambda case: case["axis_value"],
        )
        grid = next(
            case
            for case in cases
            if case["halo_concentration"] == concentration
            and case["axis"] == "grid"
        )
        baseline = baselines[concentration]
        boundary_difference = relative_difference(
            boundary[0]["final_black_hole_mass_msun"],
            boundary[1]["final_black_hole_mass_msun"],
        )
        grid_difference = relative_difference(
            baseline, grid["final_black_hole_mass_msun"]
        )
        statistics["convergence_by_concentration"][f"c{concentration:g}"] = {
            "baseline_mass_msun": baseline,
            "smallest_boundary_mass_msun": boundary[0]["final_black_hole_mass_msun"],
            "smallest_boundary_pair_relative_difference": boundary_difference,
            "grid_256_to_512_relative_difference": grid_difference,
            "passes_five_percent_convergence": bool(
                boundary_difference < 0.05 and grid_difference < 0.05
            ),
        }
        axes[0].plot(
            [case["axis_value"] for case in boundary] + [5.0 / 24.0],
            [case["final_black_hole_mass_msun"] for case in boundary] + [baseline],
            color=colors[concentration], marker="o", label=f"c={concentration:g}",
        )
        axes[1].plot(
            [256, 512],
            [baseline, grid["final_black_hole_mass_msun"]],
            color=colors[concentration], marker="o", label=f"c={concentration:g}",
        )
    statistics["maximum_mass_budget_residual_code"] = max(
        case["mass_budget_residual_code"] for case in cases
    )
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("r_min/r_influence")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].set_title("Inner-boundary response")
    axes[1].set_xlabel("Cells")
    axes[1].set_ylabel("Final black-hole mass [M_sun]")
    axes[1].set_title("Grid response")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
