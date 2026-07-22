"""Render reader-facing applicability figures from archived CSV summaries."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MAP_SUMMARY = ROOT / "results" / "stage5" / "applicability_map_summary.csv"
REFINEMENT_SUMMARY = (
    ROOT / "results" / "stage5" / "applicability_refinement_summary.csv"
)
OUTPUT_DIR = ROOT / "paper" / "figures"
TARGET_MASS_MSUN = 1.0e7

DISPLAY_LABELS = {
    "constant_sigma1": r"Constant $\sigma/m=1$",
    "vd_high_transport": "Rutherford: high transport",
    "vd_low_transport": "Rutherford: low transport",
    "vd_matched_transport": "Rutherford: virial matched",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def plot_main_effects(rows: list[dict[str, str]]) -> None:
    labels = sorted({row["model_label"] for row in rows})
    colors = dict(
        zip(labels, plt.cm.tab10(np.linspace(0.0, 0.7, len(labels))), strict=True)
    )
    axes_names = (
        "halo_mass_msun",
        "halo_redshift",
        "halo_concentration",
        "black_hole_seed_msun",
        "scale_radius_over_rs",
        "assembly_time_myr",
    )
    axis_labels = (
        r"$M_{200}$ [$M_\odot$]",
        "Redshift",
        "Concentration",
        r"Seed mass [$M_\odot$]",
        r"$a_b/r_s$",
        "Assembly time [Myr]",
    )
    fig, axes = plt.subplots(3, 2, figsize=(12, 13), constrained_layout=True)
    for axis_plot, axis_name, axis_label in zip(
        axes.flat, axes_names, axis_labels, strict=True
    ):
        for label in labels:
            selected = sorted(
                (
                    row
                    for row in rows
                    if row["model_label"] == label
                    and row["design"] == "main_effect"
                    and row["axis"] in (axis_name, "baseline")
                ),
                key=lambda row: float(row[axis_name]),
            )
            axis_plot.plot(
                [float(row[axis_name]) for row in selected],
                [float(row["final_black_hole_mass_msun"]) for row in selected],
                marker="o",
                color=colors[label],
                label=DISPLAY_LABELS.get(label, label),
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
        axis_plot.set_ylabel(r"Final black-hole mass [$M_\odot$]")
        axis_plot.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    fig.savefig(OUTPUT_DIR / "stage5_applicability_main_effects.png", dpi=180)
    plt.close(fig)


def plot_refinement(rows: list[dict[str, str]]) -> None:
    colors = {
        "robust_success": "#2ca02c",
        "robust_failure": "#7f7f7f",
        "boundary_ambiguous": "#d62728",
    }
    labels = {
        "robust_success": "Robust success",
        "robust_failure": "Robust failure",
        "boundary_ambiguous": "Boundary ambiguous",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for classification, color in colors.items():
        selected = [
            row for row in rows if row["threshold_classification"] == classification
        ]
        initial_mass = [float(row["screen_final_mass_msun"]) for row in selected]
        reduced_mass = [
            float(row["small_boundary_final_mass_msun"]) for row in selected
        ]
        difference = [
            100.0 * float(row["boundary_relative_difference"]) for row in selected
        ]
        axes[0].scatter(
            initial_mass,
            reduced_mass,
            color=color,
            label=labels[classification],
            alpha=0.8,
        )
        axes[1].scatter(
            initial_mass,
            difference,
            color=color,
            label=labels[classification],
            alpha=0.8,
        )
    upper = max(float(row["screen_final_mass_msun"]) for row in rows) * 1.2
    limits = [1.0e2, upper]
    axes[0].plot(limits, limits, color="black", linestyle="--")
    axes[0].axvline(TARGET_MASS_MSUN, color="red", alpha=0.5)
    axes[0].axhline(TARGET_MASS_MSUN, color="red", alpha=0.5)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"Initial-survey final mass [$M_\odot$]")
    axes[0].set_ylabel(r"Reduced-boundary final mass [$M_\odot$]")
    axes[0].legend(fontsize=8)
    axes[1].set_xscale("log")
    axes[1].axhline(5.0, color="black", linestyle="--")
    axes[1].set_xlabel(r"Initial-survey final mass [$M_\odot$]")
    axes[1].set_ylabel("Boundary difference [%]")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(OUTPUT_DIR / "stage5_applicability_refinement.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_main_effects(load_rows(MAP_SUMMARY))
    plot_refinement(load_rows(REFINEMENT_SUMMARY))


if __name__ == "__main__":
    main()
