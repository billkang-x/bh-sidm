"""Analyze numerical convergence for three stage-5 frontier cases."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_frontier_convergence.tsv"
RESULTS = ROOT / "results" / "stage5" / "frontier_convergence"
SUMMARY = ROOT / "results" / "stage5" / "frontier_convergence_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "frontier_convergence_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_frontier_convergence.png"
MODELS = ("representative", "resolved_high", "extreme")
AXES = ("grid", "inner_boundary", "cfl", "entropy_fix")


def selected(cases: list[dict], model: str, axis: str) -> list[dict]:
    result = [
        case for case in cases if case["model"] == model and case["axis"] == axis
    ]
    baseline = next(
        case for case in cases if case["model"] == model and case["axis"] == "baseline"
    )
    result.append({**baseline, "axis_value": baseline[axis]})
    return sorted(result, key=lambda case: case["axis_value"])


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    missing_task_ids = []
    for row in manifest:
        task_id = int(row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        if not path.exists():
            missing_task_ids.append(task_id)
            continue
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        cases.append(
            {
                "task_id": task_id,
                "model": row["model"],
                "axis": row["axis"],
                "axis_value": float(row["axis_value"]),
                "halo_concentration": float(row["halo_concentration"]),
                "sigma_over_m_cm2_g": float(row["sigma_over_m_cm2_g"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "inner_boundary": float(row["r_min_over_influence"]),
                "grid": int(row["cells"]),
                "cfl": float(row["cfl"]),
                "entropy_fix": float(row["entropy_fix"]),
                "baryon_scale_radius_over_r_min": float(
                    metadata["initial_baryon_scale_radius_pc"] / metadata["r_min_pc"]
                ),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    convergence = {}
    for model in MODELS:
        grid = selected(cases, model, "grid")
        boundary = selected(cases, model, "inner_boundary")
        cfl = selected(cases, model, "cfl")
        entropy = selected(cases, model, "entropy_fix")
        grid_difference = relative_difference(
            grid[-2]["final_black_hole_mass_msun"],
            grid[-1]["final_black_hole_mass_msun"],
        )
        boundary_difference = relative_difference(
            boundary[0]["final_black_hole_mass_msun"],
            boundary[1]["final_black_hole_mass_msun"],
        )
        cfl_spread = (
            max(case["final_black_hole_mass_msun"] for case in cfl)
            - min(case["final_black_hole_mass_msun"] for case in cfl)
        ) / np.mean([case["final_black_hole_mass_msun"] for case in cfl])
        entropy_spread = (
            max(case["final_black_hole_mass_msun"] for case in entropy)
            - min(case["final_black_hole_mass_msun"] for case in entropy)
        ) / np.mean([case["final_black_hole_mass_msun"] for case in entropy])
        convergence[model] = {
            "baseline": next(
                case
                for case in cases
                if case["model"] == model and case["axis"] == "baseline"
            ),
            "grid_384_to_512_relative_difference": grid_difference,
            "smallest_boundary_pair_relative_difference": boundary_difference,
            "cfl_mass_spread_fraction": cfl_spread,
            "entropy_mass_spread_fraction": entropy_spread,
            "smallest_boundary_mass_msun": boundary[0]["final_black_hole_mass_msun"],
            "smallest_boundary_a_over_r_min": boundary[0]["baryon_scale_radius_over_r_min"],
            "passes_five_percent_convergence": bool(
                grid_difference < 0.05
                and boundary_difference < 0.05
                and cfl_spread < 0.05
                and entropy_spread < 0.05
            ),
        }
    statistics = {
        "case_count": len(cases),
        "expected_case_count": len(manifest),
        "missing_task_ids": missing_task_ids,
        "convergence_by_model": convergence,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
        "maximum_steps": max(case["steps"] for case in cases),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )

    labels = {
        "grid": "Cells",
        "inner_boundary": "r_min/r_influence",
        "cfl": "CFL",
        "entropy_fix": "Entropy fix",
    }
    fig, axes = plt.subplots(3, 4, figsize=(14, 10), constrained_layout=True)
    colors = {"representative": "#1b9e77", "resolved_high": "#2166ac", "extreme": "#b2182b"}
    for row_index, model in enumerate(MODELS):
        for column_index, axis_name in enumerate(AXES):
            axis = axes[row_index, column_index]
            values = selected(cases, model, axis_name)
            axis.plot(
                [case["axis_value"] for case in values],
                [case["final_black_hole_mass_msun"] for case in values],
                color=colors[model], marker="o",
            )
            if axis_name == "inner_boundary":
                axis.set_xscale("log")
            axis.set_xlabel(labels[axis_name])
            axis.set_ylabel("Final black-hole mass [M_sun]")
            axis.set_title(f"{model}: {axis_name.replace('_', ' ')}")
            axis.grid(alpha=0.25)
            axis.tick_params(labelsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
