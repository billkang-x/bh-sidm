"""Analyze the final seed threshold with boundary and grid checks."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_seed_threshold_final.tsv"
RESULTS = ROOT / "results" / "stage5" / "seed_threshold_final"
SUMMARY = ROOT / "results" / "stage5" / "seed_threshold_final_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "seed_threshold_final_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_seed_threshold_final.png"


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    histories = {}
    for row in manifest:
        task_id = int(row["task_id"])
        seed = float(row["black_hole_seed_msun"])
        axis = row["axis"]
        boundary = float(row["reference_r_min_over_influence"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            if axis == "inner_boundary":
                histories[(seed, boundary)] = (
                    data["times_myr"].copy(),
                    data["black_hole_mass_msun"].copy(),
                )
        cases.append(
            {
                "task_id": task_id,
                "black_hole_seed_msun": seed,
                "axis": axis,
                "reference_r_min_over_influence": boundary,
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
                "dark_capture_fraction": float(metadata["dark_capture_fraction_of_available_supply"]),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    by_seed = {}
    robust_success = []
    robust_failure = []
    ambiguous = []
    for seed in sorted({case["black_hole_seed_msun"] for case in cases}):
        boundary = sorted(
            (
                case
                for case in cases
                if case["black_hole_seed_msun"] == seed
                and case["axis"] == "inner_boundary"
            ),
            key=lambda case: case["reference_r_min_over_influence"],
        )
        grid = next(
            case
            for case in cases
            if case["black_hole_seed_msun"] == seed and case["axis"] == "grid"
        )
        reaches = [np.isfinite(case["time_to_1e7_msun_myr"]) for case in boundary]
        grid_reaches = np.isfinite(grid["time_to_1e7_msun_myr"])
        classification = "ambiguous"
        if all(reaches) and grid_reaches:
            classification = "robust_success"
            robust_success.append(seed)
        elif not any(reaches) and not grid_reaches:
            classification = "robust_failure"
            robust_failure.append(seed)
        else:
            ambiguous.append(seed)
        grid_base = boundary[1]
        by_seed[f"{seed:.0f}"] = {
            "classification": classification,
            "small_boundary_final_mass_msun": boundary[0]["final_black_hole_mass_msun"],
            "large_boundary_final_mass_msun": boundary[1]["final_black_hole_mass_msun"],
            "boundary_relative_difference": relative_difference(
                boundary[0]["final_black_hole_mass_msun"],
                boundary[1]["final_black_hole_mass_msun"],
            ),
            "grid_relative_difference": relative_difference(
                grid_base["final_black_hole_mass_msun"],
                grid["final_black_hole_mass_msun"],
            ),
            "small_boundary_time_to_1e7_myr": (
                boundary[0]["time_to_1e7_msun_myr"] if reaches[0] else None
            ),
            "large_boundary_time_to_1e7_myr": (
                boundary[1]["time_to_1e7_msun_myr"] if reaches[1] else None
            ),
            "fine_grid_time_to_1e7_myr": (
                grid["time_to_1e7_msun_myr"] if grid_reaches else None
            ),
        }
    statistics = {
        "case_count": len(cases),
        "classification_by_seed": by_seed,
        "highest_robust_failure_seed_msun": max(robust_failure) if robust_failure else None,
        "lowest_robust_success_seed_msun": min(robust_success) if robust_success else None,
        "ambiguous_seed_masses_msun": ambiguous,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for boundary, marker in zip((0.05208333333333334, 0.10416666666666667), ("o", "s")):
        selected = sorted(
            (
                case
                for case in cases
                if case["axis"] == "inner_boundary"
                and np.isclose(case["reference_r_min_over_influence"], boundary)
            ),
            key=lambda case: case["black_hole_seed_msun"],
        )
        axes[0].plot(
            [case["black_hole_seed_msun"] for case in selected],
            [case["final_black_hole_mass_msun"] for case in selected],
            marker=marker,
            label=f"feed q={boundary:.3f}",
        )
    axes[0].axhline(1.0e7, color="black", linestyle="--", alpha=0.5)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Seed mass [M_sun]")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].legend()

    for seed in sorted({case["black_hole_seed_msun"] for case in cases}):
        time, mass = histories[(seed, 0.10416666666666667)]
        axes[1].plot(time, mass, label=f"{seed:.2e}")
    axes[1].axhline(1.0e7, color="black", linestyle="--", alpha=0.5)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time [Myr]")
    axes[1].set_ylabel("Black-hole mass [M_sun]")
    axes[1].legend(title="Seed")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
