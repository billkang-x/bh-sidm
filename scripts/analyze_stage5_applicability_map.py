"""Analyze the controlled six-parameter SIDM applicability map."""

from __future__ import annotations

import csv
import json
from math import log
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_applicability_map.tsv"
RESULTS = ROOT / "results" / "stage5" / "applicability_map"
SUMMARY = ROOT / "results" / "stage5" / "applicability_map_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "applicability_map_statistics.json"
REFINEMENT = ROOT / "hpc" / "stage5_applicability_refinement.tsv"
MAIN_FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_applicability_main_effects.png"
JOINT_FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_applicability_joint_map.png"
TARGET_MASS_MSUN = 1.0e7
FEATURE_NAMES = (
    "log10_halo_mass",
    "halo_redshift",
    "halo_concentration",
    "log10_seed_mass",
    "log10_scale_radius_over_rs",
    "assembly_over_dynamical_time",
)


def sigmoid(value: np.ndarray) -> np.ndarray:
    positive = value >= 0.0
    result = np.empty_like(value)
    result[positive] = 1.0 / (1.0 + np.exp(-value[positive]))
    exponential = np.exp(value[~positive])
    result[~positive] = exponential / (1.0 + exponential)
    return result


def fit_logistic(
    features: np.ndarray,
    target: np.ndarray,
    regularization: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(features, axis=0)
    scale = np.std(features, axis=0)
    scale[scale == 0.0] = 1.0
    standardized = (features - mean) / scale
    design = np.column_stack((np.ones(len(features)), standardized))
    coefficients = np.zeros(design.shape[1])
    penalty = np.ones_like(coefficients)
    penalty[0] = 0.0
    for _ in range(100):
        probability = sigmoid(design @ coefficients)
        weight = np.maximum(probability * (1.0 - probability), 1.0e-8)
        hessian = design.T @ (weight[:, None] * design)
        hessian += regularization * np.diag(penalty)
        gradient = design.T @ (target - probability)
        gradient -= regularization * penalty * coefficients
        update = np.linalg.solve(hessian, gradient)
        coefficients += update
        if np.max(np.abs(update)) < 1.0e-9:
            break
    return coefficients, mean, scale


def predict_logistic(
    features: np.ndarray,
    coefficients: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
) -> np.ndarray:
    design = np.column_stack((np.ones(len(features)), (features - mean) / scale))
    return sigmoid(design @ coefficients)


def binary_auc(target: np.ndarray, score: np.ndarray) -> float:
    positive = score[target == 1]
    negative = score[target == 0]
    if len(positive) == 0 or len(negative) == 0:
        return float("nan")
    comparisons = positive[:, None] - negative[None, :]
    return float(
        np.mean(
            (comparisons > 0.0).astype(float)
            + 0.5 * (comparisons == 0.0).astype(float)
        )
    )


def feature_matrix(cases: list[dict]) -> np.ndarray:
    return np.array(
        [
            [
                np.log10(case["halo_mass_msun"]),
                case["halo_redshift"],
                case["halo_concentration"],
                np.log10(case["black_hole_seed_msun"]),
                np.log10(case["scale_radius_over_rs"]),
                case["assembly_over_dynamical_time"],
            ]
            for case in cases
        ]
    )


def model_classifier(cases: list[dict]) -> dict:
    features = feature_matrix(cases)
    target = np.array(
        [case["final_black_hole_mass_msun"] >= TARGET_MASS_MSUN for case in cases],
        dtype=int,
    )
    if np.all(target == target[0]):
        return {
            "class_count": {"failure": int(np.sum(target == 0)), "success": int(np.sum(target == 1))},
            "fit_available": False,
        }
    coefficients, mean, scale = fit_logistic(features, target)
    fold_score = np.empty(len(cases))
    for fold in range(5):
        validation = np.arange(len(cases)) % 5 == fold
        train = ~validation
        fold_coefficients, fold_mean, fold_scale = fit_logistic(
            features[train], target[train]
        )
        fold_score[validation] = predict_logistic(
            features[validation], fold_coefficients, fold_mean, fold_scale
        )
    return {
        "class_count": {
            "failure": int(np.sum(target == 0)),
            "success": int(np.sum(target == 1)),
        },
        "fit_available": True,
        "standardized_coefficients": {
            "intercept": float(coefficients[0]),
            **{
                name: float(value)
                for name, value in zip(FEATURE_NAMES, coefficients[1:], strict=True)
            },
        },
        "five_fold_auc": binary_auc(target, fold_score),
        "five_fold_accuracy": float(np.mean((fold_score >= 0.5) == target)),
    }


def build_refinement(
    cases: list[dict],
    manifest_lookup: dict[int, dict[str, str]],
) -> list[dict]:
    rows = []
    for label in sorted({case["model_label"] for case in cases}):
        selected_model = [case for case in cases if case["model_label"] == label]
        proximity = lambda case: abs(
            np.log10(case["final_black_hole_mass_msun"] / TARGET_MASS_MSUN)
        )
        main_cases = sorted(
            (case for case in selected_model if case["design"] == "main_effect"),
            key=proximity,
        )[:8]
        joint_cases = sorted(
            (case for case in selected_model if case["design"] == "latin_hypercube"),
            key=proximity,
        )[:8]
        audit_cases = main_cases + joint_cases
        grid_cases = sorted(audit_cases, key=proximity)[:4]
        for case in audit_cases:
            source = manifest_lookup[case["task_id"]]
            screen_r_min = float(source["r_min_over_rs"])
            r_max = float(source["r_max_over_rs"])
            screen_cells = int(source["cells"])
            log_width = log(r_max / screen_r_min) / screen_cells
            small_r_min = 0.5 * screen_r_min
            rows.append(
                {
                    "task_id": len(rows),
                    "source_task_id": case["task_id"],
                    "variant": "inner_boundary",
                    "model_label": source["model_label"],
                    "cross_section_model": source["cross_section_model"],
                    "sigma0_over_m_cm2_g": source["sigma0_over_m_cm2_g"],
                    "velocity_scale_km_s": source["velocity_scale_km_s"],
                    "halo_mass_msun": source["halo_mass_msun"],
                    "halo_redshift": source["halo_redshift"],
                    "halo_concentration": source["halo_concentration"],
                    "black_hole_seed_msun": source["black_hole_seed_msun"],
                    "scale_radius_over_rs": source["scale_radius_over_rs"],
                    "assembly_time_myr": source["assembly_time_myr"],
                    "r_min_over_reference_influence": 5.0 / 48.0,
                    "r_min_over_rs": small_r_min,
                    "r_max_over_rs": r_max,
                    "cells": int(round(log(r_max / small_r_min) / log_width)),
                }
            )
        for case in grid_cases:
            source = manifest_lookup[case["task_id"]]
            rows.append(
                {
                    "task_id": len(rows),
                    "source_task_id": case["task_id"],
                    "variant": "grid",
                    "model_label": source["model_label"],
                    "cross_section_model": source["cross_section_model"],
                    "sigma0_over_m_cm2_g": source["sigma0_over_m_cm2_g"],
                    "velocity_scale_km_s": source["velocity_scale_km_s"],
                    "halo_mass_msun": source["halo_mass_msun"],
                    "halo_redshift": source["halo_redshift"],
                    "halo_concentration": source["halo_concentration"],
                    "black_hole_seed_msun": source["black_hole_seed_msun"],
                    "scale_radius_over_rs": source["scale_radius_over_rs"],
                    "assembly_time_myr": source["assembly_time_myr"],
                    "r_min_over_reference_influence": 5.0 / 24.0,
                    "r_min_over_rs": source["r_min_over_rs"],
                    "r_max_over_rs": source["r_max_over_rs"],
                    "cells": 2 * int(source["cells"]),
                }
            )
    if len(rows) != 80:
        raise RuntimeError(f"expected 80 refinement cases, generated {len(rows)}")
    with REFINEMENT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    manifest_lookup = {int(row["task_id"]): row for row in manifest}
    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        final_mass = float(metadata["final_black_hole_mass_msun"])
        cases.append(
            {
                "task_id": task_id,
                "design": row["design"],
                "axis": row["axis"],
                "axis_value": float(row["axis_value"]),
                "model_label": row["model_label"],
                "cross_section_model": row["cross_section_model"],
                "sigma0_over_m_cm2_g": float(row["sigma0_over_m_cm2_g"]),
                "velocity_scale_km_s": float(row["velocity_scale_km_s"]),
                "halo_mass_msun": float(row["halo_mass_msun"]),
                "halo_redshift": float(row["halo_redshift"]),
                "halo_concentration": float(row["halo_concentration"]),
                "black_hole_seed_msun": float(row["black_hole_seed_msun"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "assembly_time_myr": float(row["assembly_time_myr"]),
                "initial_dynamical_time_myr": float(row["initial_dynamical_time_myr"]),
                "assembly_over_dynamical_time": float(row["assembly_over_dynamical_time"]),
                "r_min_over_reference_influence": float(
                    row["r_min_over_reference_influence"]
                ),
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": final_mass,
                "growth_factor": final_mass / float(row["black_hole_seed_msun"]),
                "reaches_1e7": bool(final_mass >= TARGET_MASS_MSUN),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "dark_fraction_of_growth": float(metadata["dark_fraction_of_black_hole_growth"]),
                "initial_inner_lmfp_effective_sigma_cm2_g": float(
                    metadata["initial_inner_lmfp_effective_sigma_over_m_cm2_g"]
                ),
                "final_inner_lmfp_effective_sigma_cm2_g": float(
                    metadata["final_inner_lmfp_effective_sigma_over_m_cm2_g"]
                ),
                "initial_inner_smfp_effective_sigma_cm2_g": float(
                    metadata["initial_inner_smfp_effective_sigma_over_m_cm2_g"]
                ),
                "final_inner_smfp_effective_sigma_cm2_g": float(
                    metadata["final_inner_smfp_effective_sigma_over_m_cm2_g"]
                ),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )

    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    refinement = build_refinement(cases, manifest_lookup)
    model_statistics = {}
    for label in sorted({case["model_label"] for case in cases}):
        selected = [case for case in cases if case["model_label"] == label]
        joint = [case for case in selected if case["design"] == "latin_hypercube"]
        model_statistics[label] = {
            "screen_success_count": sum(case["reaches_1e7"] for case in selected),
            "screen_success_fraction": float(np.mean([case["reaches_1e7"] for case in selected])),
            "joint_design_classifier": model_classifier(joint),
            "maximum_final_mass_msun": max(case["final_black_hole_mass_msun"] for case in selected),
            "minimum_successful_seed_msun": min(
                (case["black_hole_seed_msun"] for case in selected if case["reaches_1e7"]),
                default=None,
            ),
        }
    statistics = {
        "case_count": len(cases),
        "refinement_case_count": len(refinement),
        "target_mass_msun": TARGET_MASS_MSUN,
        "models": model_statistics,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    labels = sorted({case["model_label"] for case in cases})
    colors = dict(zip(labels, plt.cm.tab10(np.linspace(0.0, 0.7, len(labels))), strict=True))
    display_labels = {
        "constant_sigma1": r"Constant $\sigma/m=1$",
        "vd_high_transport": "Rutherford: high transport",
        "vd_low_transport": "Rutherford: low transport",
        "vd_matched_transport": "Rutherford: virial matched",
    }
    axes_names = (
        "halo_mass_msun",
        "halo_redshift",
        "halo_concentration",
        "black_hole_seed_msun",
        "scale_radius_over_rs",
        "assembly_time_myr",
    )
    axis_labels = (
        "M200 [M_sun]",
        "Redshift",
        "Concentration",
        "Seed mass [M_sun]",
        "a_b/r_s",
        "Assembly time [Myr]",
    )
    MAIN_FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 2, figsize=(12, 13), constrained_layout=True)
    for axis_plot, axis_name, axis_label in zip(
        axes.flat, axes_names, axis_labels, strict=True
    ):
        for label in labels:
            selected = sorted(
                (
                    case
                    for case in cases
                    if case["model_label"] == label
                    and case["design"] == "main_effect"
                    and case["axis"] in (axis_name, "baseline")
                ),
                key=lambda case: case[axis_name],
            )
            axis_plot.plot(
                [case[axis_name] for case in selected],
                [case["final_black_hole_mass_msun"] for case in selected],
                marker="o",
                color=colors[label],
                label=display_labels.get(label, label),
            )
        if axis_name in (
            "halo_mass_msun",
            "black_hole_seed_msun",
            "scale_radius_over_rs",
        ):
            axis_plot.set_xscale("log")
        axis_plot.set_yscale("log")
        axis_plot.axhline(TARGET_MASS_MSUN, color="red", linestyle="--", alpha=0.7)
        axis_plot.set_xlabel(axis_label)
        axis_plot.set_ylabel("Final black-hole mass [M_sun]")
        axis_plot.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    fig.savefig(MAIN_FIGURE, dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(len(labels), 3, figsize=(14, 3.4 * len(labels)), constrained_layout=True)
    for row_index, label in enumerate(labels):
        selected = [
            case
            for case in cases
            if case["model_label"] == label and case["design"] == "latin_hypercube"
        ]
        color_values = np.log10(
            [case["final_black_hole_mass_msun"] for case in selected]
        )
        projections = (
            (
                [case["halo_mass_msun"] for case in selected],
                [case["black_hole_seed_msun"] for case in selected],
                "M200 [M_sun]",
                "Seed mass [M_sun]",
                True,
                True,
            ),
            (
                [case["scale_radius_over_rs"] for case in selected],
                [case["assembly_over_dynamical_time"] for case in selected],
                "a_b/r_s",
                "T_asm/t_dyn",
                True,
                False,
            ),
            (
                [case["halo_redshift"] for case in selected],
                [case["halo_concentration"] for case in selected],
                "Redshift",
                "Concentration",
                False,
                False,
            ),
        )
        for column, (x, y, xlabel, ylabel, logx, logy) in enumerate(projections):
            axis_plot = axes[row_index, column]
            scatter = axis_plot.scatter(
                x,
                y,
                c=color_values,
                vmin=2.0,
                vmax=8.0,
                cmap="viridis",
                edgecolor="black",
                linewidth=0.25,
            )
            if logx:
                axis_plot.set_xscale("log")
            if logy:
                axis_plot.set_yscale("log")
            axis_plot.set_xlabel(xlabel)
            axis_plot.set_ylabel(ylabel)
            axis_plot.grid(alpha=0.2)
            if column == 0:
                axis_plot.set_title(display_labels.get(label, label))
        fig.colorbar(scatter, ax=axes[row_index, :], label="log10 final mass [M_sun]")
    fig.savefig(JOINT_FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
