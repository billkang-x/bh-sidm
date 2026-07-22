"""Close the numerical boundary of the stage-5 applicability map."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCREEN_SUMMARY = ROOT / "results" / "stage5" / "applicability_map_summary.csv"
MANIFEST = ROOT / "hpc" / "stage5_applicability_refinement.tsv"
RESULTS = ROOT / "results" / "stage5" / "applicability_refinement"
SUMMARY = ROOT / "results" / "stage5" / "applicability_refinement_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "applicability_refinement_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_applicability_refinement.png"
TARGET_MASS_MSUN = 1.0e7


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with SCREEN_SUMMARY.open(newline="", encoding="ascii") as stream:
        screen_rows = list(csv.DictReader(stream))
    screen = {int(row["task_id"]): row for row in screen_rows}
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))

    refinement = []
    for row in manifest:
        task_id = int(row["task_id"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        refinement.append(
            {
                "task_id": task_id,
                "source_task_id": int(row["source_task_id"]),
                "variant": row["variant"],
                "model_label": row["model_label"],
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )

    refinement_lookup = {
        (row["source_task_id"], row["variant"]): row for row in refinement
    }
    closed = []
    for source_id in sorted(
        {row["source_task_id"] for row in refinement if row["variant"] == "inner_boundary"}
    ):
        base = screen[source_id]
        boundary = refinement_lookup[(source_id, "inner_boundary")]
        grid = refinement_lookup.get((source_id, "grid"))
        screen_mass = float(base["final_black_hole_mass_msun"])
        boundary_mass = boundary["final_black_hole_mass_msun"]
        screen_success = screen_mass >= TARGET_MASS_MSUN
        boundary_success = boundary_mass >= TARGET_MASS_MSUN
        if screen_success and boundary_success:
            classification = "robust_success"
        elif not screen_success and not boundary_success:
            classification = "robust_failure"
        else:
            classification = "boundary_ambiguous"
        boundary_difference = relative_difference(screen_mass, boundary_mass)
        grid_difference = (
            relative_difference(screen_mass, grid["final_black_hole_mass_msun"])
            if grid is not None
            else float("nan")
        )
        closed.append(
            {
                "source_task_id": source_id,
                "design": base["design"],
                "axis": base["axis"],
                "model_label": base["model_label"],
                "halo_mass_msun": float(base["halo_mass_msun"]),
                "halo_redshift": float(base["halo_redshift"]),
                "halo_concentration": float(base["halo_concentration"]),
                "black_hole_seed_msun": float(base["black_hole_seed_msun"]),
                "scale_radius_over_rs": float(base["scale_radius_over_rs"]),
                "assembly_time_myr": float(base["assembly_time_myr"]),
                "assembly_over_dynamical_time": float(
                    base["assembly_over_dynamical_time"]
                ),
                "screen_final_mass_msun": screen_mass,
                "small_boundary_final_mass_msun": boundary_mass,
                "grid_final_mass_msun": (
                    grid["final_black_hole_mass_msun"] if grid is not None else float("nan")
                ),
                "boundary_relative_difference": boundary_difference,
                "grid_relative_difference": grid_difference,
                "terminal_mass_passes_five_percent": bool(
                    boundary_difference < 0.05
                    and (not np.isfinite(grid_difference) or grid_difference < 0.05)
                ),
                "threshold_classification": classification,
                "small_boundary_time_to_1e7_msun_myr": boundary[
                    "time_to_1e7_msun_myr"
                ],
            }
        )

    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(closed[0]))
        writer.writeheader()
        writer.writerows(closed)

    by_model = {}
    for label in sorted({row["model_label"] for row in closed}):
        selected = [row for row in closed if row["model_label"] == label]
        by_model[label] = {
            "audited_count": len(selected),
            "robust_success_count": sum(
                row["threshold_classification"] == "robust_success"
                for row in selected
            ),
            "robust_failure_count": sum(
                row["threshold_classification"] == "robust_failure"
                for row in selected
            ),
            "boundary_ambiguous_count": sum(
                row["threshold_classification"] == "boundary_ambiguous"
                for row in selected
            ),
            "terminal_mass_five_percent_count": sum(
                row["terminal_mass_passes_five_percent"] for row in selected
            ),
            "maximum_boundary_relative_difference": max(
                row["boundary_relative_difference"] for row in selected
            ),
        }
    statistics = {
        "refinement_case_count": len(refinement),
        "audited_screen_point_count": len(closed),
        "target_mass_msun": TARGET_MASS_MSUN,
        "models": by_model,
        "maximum_mass_budget_residual_code": max(
            row["mass_budget_residual_code"] for row in refinement
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    colors = {
        "robust_success": "#2ca02c",
        "robust_failure": "#7f7f7f",
        "boundary_ambiguous": "#d62728",
    }
    display_labels = {
        "robust_success": "Robust success",
        "robust_failure": "Robust failure",
        "boundary_ambiguous": "Boundary ambiguous",
    }
    for classification, color in colors.items():
        selected = [
            row for row in closed if row["threshold_classification"] == classification
        ]
        axes[0].scatter(
            [row["screen_final_mass_msun"] for row in selected],
            [row["small_boundary_final_mass_msun"] for row in selected],
            color=color,
            label=display_labels[classification],
            alpha=0.8,
        )
        axes[1].scatter(
            [row["screen_final_mass_msun"] for row in selected],
            [100.0 * row["boundary_relative_difference"] for row in selected],
            color=color,
            label=display_labels[classification],
            alpha=0.8,
        )
    limits = [1.0e2, max(row["screen_final_mass_msun"] for row in closed) * 1.2]
    axes[0].plot(limits, limits, color="black", linestyle="--")
    axes[0].axvline(TARGET_MASS_MSUN, color="red", alpha=0.5)
    axes[0].axhline(TARGET_MASS_MSUN, color="red", alpha=0.5)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Initial-survey final mass [M_sun]")
    axes[0].set_ylabel("Reduced-boundary final mass [M_sun]")
    axes[0].legend(fontsize=8)
    axes[1].set_xscale("log")
    axes[1].axhline(5.0, color="black", linestyle="--")
    axes[1].set_xlabel("Initial-survey final mass [M_sun]")
    axes[1].set_ylabel("Boundary difference [%]")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
