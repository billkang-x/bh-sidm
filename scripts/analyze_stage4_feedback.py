"""Analyze parametric Hernquist-expansion feedback in stage 4."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_feedback.tsv"
RESULTS = ROOT / "results" / "stage4" / "feedback"
REFINEMENT_MANIFEST = ROOT / "hpc" / "stage4_feedback_refinement.tsv"
REFINEMENT_RESULTS = ROOT / "results" / "stage4" / "feedback_refinement"
NO_FEEDBACK_SUMMARY = ROOT / "results" / "stage4" / "no_feedback_summary.csv"
SUMMARY = ROOT / "results" / "stage4" / "feedback_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "feedback_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_feedback.png"


def load_cases() -> tuple[list[dict], float]:
    with NO_FEEDBACK_SUMMARY.open(newline="", encoding="utf-8") as stream:
        baseline_rows = list(csv.DictReader(stream))
    no_eddington_control = next(
        float(row["dark_matter_accreted_msun"])
        for row in baseline_rows
        if float(row["seed_mass_msun"]) == 100.0
        and float(row["eddington_ratio"]) == 0.0
    )
    matrices = [("base", MANIFEST, RESULTS)]
    if REFINEMENT_MANIFEST.exists() and REFINEMENT_RESULTS.exists():
        matrices.append(
            ("refinement", REFINEMENT_MANIFEST, REFINEMENT_RESULTS)
        )
    cases = []
    for matrix_name, manifest_path, result_directory in matrices:
        with manifest_path.open(newline="", encoding="ascii") as stream:
            manifest = list(csv.DictReader(stream, delimiter="\t"))
        for row in manifest:
            task_id = int(row["task_id"])
            path = result_directory / f"task_{task_id:03d}.npz"
            if not path.exists():
                raise FileNotFoundError(path)
            with np.load(path, allow_pickle=False) as data:
                metadata = json.loads(str(data["metadata_json"]))
            case = {
                "matrix": matrix_name,
                "task_id": task_id,
                "feedback_efficiency": float(row["feedback_efficiency"]),
                "feedback_eta": float(row["feedback_eta"]),
            "final_black_hole_mass_msun": float(
                metadata["final_black_hole_mass_msun"]
            ),
            "dark_matter_accreted_msun": float(
                metadata["dark_matter_accreted_msun"]
            ),
            "baryon_accreted_onto_bh_msun": float(
                metadata["baryon_accreted_onto_bh_msun"]
            ),
            "dark_fraction_of_black_hole_growth": float(
                metadata["dark_fraction_of_black_hole_growth"]
            ),
            "final_baryon_scale_radius_pc": float(
                metadata["final_baryon_scale_radius_pc"]
            ),
            "final_feedback_to_binding_ratio": float(
                metadata["final_feedback_to_binding_ratio"]
            ),
            "mass_budget_residual_code": float(
                metadata["mass_budget_residual_code"]
            ),
            "path": path,
            }
            case["scale_radius_expansion_factor"] = (
                case["final_baryon_scale_radius_pc"] / 0.3
            )
            case["extra_dark_matter_vs_no_eddington_msun"] = (
                case["dark_matter_accreted_msun"] - no_eddington_control
            )
            case["catalytic_dark_mass_per_baryon_mass"] = (
                case["extra_dark_matter_vs_no_eddington_msun"]
                / case["baryon_accreted_onto_bh_msun"]
            )
            cases.append(case)
    reference = next(case for case in cases if case["feedback_efficiency"] == 0.0)
    for case in cases:
        case["dark_growth_fraction_of_no_feedback"] = (
            case["dark_matter_accreted_msun"]
            / reference["dark_matter_accreted_msun"]
        )
        case["total_growth_fraction_of_no_feedback"] = (
            (case["final_black_hole_mass_msun"] - 100.0)
            / (reference["final_black_hole_mass_msun"] - 100.0)
        )
    return cases, no_eddington_control


def save_summary(cases: list[dict]) -> None:
    rows = [{key: value for key, value in case.items() if key != "path"} for case in cases]
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(cases: list[dict], no_eddington_control: float) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {0.5: "#3b4cc0", 1.0: "#d62828"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    reference = next(case for case in cases if case["feedback_efficiency"] == 0.0)
    for eta in (0.5, 1.0):
        selected = [reference] + sorted(
            [
                case
                for case in cases
                if case["feedback_efficiency"] > 0.0
                and case["feedback_eta"] == eta
            ],
            key=lambda case: case["feedback_efficiency"],
        )
        efficiencies = [
            max(case["feedback_efficiency"], 3.0e-8)
            for case in selected
        ]
        axes[0, 0].plot(
            efficiencies,
            [case["scale_radius_expansion_factor"] for case in selected],
            marker="o",
            color=colors[eta],
            label=f"eta={eta:g}",
        )
        axes[0, 1].plot(
            efficiencies,
            [case["dark_matter_accreted_msun"] for case in selected],
            marker="o",
            color=colors[eta],
            label=f"eta={eta:g}",
        )
        axes[1, 0].plot(
            efficiencies,
            [case["catalytic_dark_mass_per_baryon_mass"] for case in selected],
            marker="o",
            color=colors[eta],
            label=f"eta={eta:g}",
        )
    for axis in axes[0, :]:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Feedback efficiency")
        axis.legend(frameon=False)
        axis.grid(alpha=0.25)
    axes[0, 0].set_ylabel("Final Hernquist scale radius / initial")
    axes[0, 0].set_title("Feedback-driven potential expansion")
    axes[0, 1].axhline(
        no_eddington_control,
        color="black",
        linestyle=":",
        label="No-Eddington control",
    )
    axes[0, 1].set_ylabel("Dark matter accreted [M_sun]")
    axes[0, 1].set_title("Dark-growth suppression")

    axes[1, 0].axhline(0.0, color="black", linestyle=":")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xlabel("Feedback efficiency")
    axes[1, 0].set_ylabel("Extra DM / baryon BH mass")
    axes[1, 0].set_title("Feedback reverses baryonic catalysis")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(frameon=False)

    representative = [
        next(case for case in cases if case["feedback_efficiency"] == 0.0),
        next(
            case
            for case in cases
            if case["feedback_efficiency"] == 1.0e-4
            and case["feedback_eta"] == 0.5
        ),
        next(
            case
            for case in cases
            if case["feedback_efficiency"] == 1.0e-2
            and case["feedback_eta"] == 0.5
        ),
    ]
    for case in representative:
        with np.load(case["path"], allow_pickle=False) as data:
            label = f"eps_f={case['feedback_efficiency']:g}"
            axes[1, 1].plot(
                data["times_myr"],
                data["baryon_scale_radius_pc"] / 0.3,
                label=label,
            )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Time [Myr]")
    axes[1, 1].set_ylabel("Hernquist scale radius / initial")
    axes[1, 1].set_title("Expansion history, eta=0.5")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases, no_eddington_control = load_cases()
    save_summary(cases)
    plot(cases, no_eddington_control)
    active = [case for case in cases if case["feedback_efficiency"] > 0.0]
    reversed_cases = [
        case
        for case in active
        if case["extra_dark_matter_vs_no_eddington_msun"] < 0.0
    ]
    thresholds = {}
    estimated_thresholds = {}
    threshold_energy_ratios = {}
    for eta in (0.5, 1.0):
        selected = sorted(
            [case for case in reversed_cases if case["feedback_eta"] == eta],
            key=lambda case: case["feedback_efficiency"],
        )
        thresholds[str(eta)] = (
            selected[0]["feedback_efficiency"] if selected else None
        )
        threshold_energy_ratios[str(eta)] = (
            selected[0]["final_feedback_to_binding_ratio"] if selected else None
        )
        ordered = sorted(
            [
                case
                for case in active
                if case["feedback_eta"] == eta
            ],
            key=lambda case: case["feedback_efficiency"],
        )
        estimate = None
        for lower, upper in zip(ordered[:-1], ordered[1:]):
            y_lower = lower["extra_dark_matter_vs_no_eddington_msun"]
            y_upper = upper["extra_dark_matter_vs_no_eddington_msun"]
            if y_lower >= 0.0 and y_upper < 0.0:
                log_lower = np.log10(lower["feedback_efficiency"])
                log_upper = np.log10(upper["feedback_efficiency"])
                fraction = y_lower / (y_lower - y_upper)
                estimate = float(10.0 ** (
                    log_lower + fraction * (log_upper - log_lower)
                ))
                break
        estimated_thresholds[str(eta)] = estimate
    statistics = {
        "all_cases_present": len(cases) in (9, 17),
        "maximum_mass_budget_residual_code": max(
            abs(case["mass_budget_residual_code"]) for case in cases
        ),
        "feedback_reversal_sampled_threshold_by_eta": thresholds,
        "feedback_reversal_log_interpolated_threshold_by_eta": estimated_thresholds,
        "feedback_to_binding_ratio_at_sampled_reversal": threshold_energy_ratios,
        "minimum_dark_growth_fraction_of_no_feedback": min(
            case["dark_growth_fraction_of_no_feedback"] for case in cases
        ),
        "minimum_total_growth_fraction_of_no_feedback": min(
            case["total_growth_fraction_of_no_feedback"] for case in cases
        ),
        "maximum_scale_radius_expansion_factor": max(
            case["scale_radius_expansion_factor"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(SUMMARY)
    print(STATISTICS)
    print(FIGURE)
    print(json.dumps(statistics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
