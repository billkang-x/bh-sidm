"""Analyze the stage-3 baryon-fraction and halo-mass scaling matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage3_scaling_matrix.tsv"
RESULTS = ROOT / "results" / "stage3" / "scaling"
REFINEMENT_MANIFEST = ROOT / "hpc" / "stage3_scaling_refinement.tsv"
REFINEMENT_RESULTS = ROOT / "results" / "stage3" / "scaling_refinement"
SUMMARY = ROOT / "results" / "stage3" / "scaling_summary.csv"
OPTIMA = ROOT / "results" / "stage3" / "scaling_optima.csv"
STATISTICS = ROOT / "results" / "stage3" / "scaling_statistics.json"
FIGURE = ROOT / "results" / "stage3" / "figures" / "stage3_scaling.png"


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_rows() -> list[dict[str, float | int | str]]:
    rows = []
    matrices = [("base", MANIFEST, RESULTS)]
    if REFINEMENT_MANIFEST.exists() and REFINEMENT_RESULTS.exists():
        matrices.append(
            ("refinement", REFINEMENT_MANIFEST, REFINEMENT_RESULTS)
        )
    for matrix_name, manifest, result_directory in matrices:
        for manifest_row in load_manifest(manifest):
            task_id = int(manifest_row["task_id"])
            path = result_directory / f"task_{task_id:03d}.npz"
            if not path.exists():
                raise FileNotFoundError(path)
            with np.load(path, allow_pickle=False) as data:
                metadata = json.loads(str(data["metadata_json"]))
                seed_mass = float(metadata["black_hole_seed_msun"])
                growth = float(data["black_hole_mass_msun"][-1] - seed_mass)
            rows.append(
                {
                    "task_id": task_id,
                    "matrix": matrix_name,
                    "case_type": manifest_row["case_type"],
                    "halo_mass_msun": float(manifest_row["halo_mass_msun"]),
                    "baryon_fraction": float(manifest_row["baryon_fraction"]),
                    "supply_radius_pc": float(manifest_row["supply_radius_pc"]),
                    "predicted_time_myr": float(manifest_row["predicted_time_myr"]),
                    "assembly_multiplier": float(
                        manifest_row["assembly_multiplier"]
                    ),
                    "assembly_time_myr": float(manifest_row["assembly_time_myr"]),
                    "accreted_dark_matter_msun": growth,
                    "final_black_hole_mass_msun": float(
                        metadata["final_black_hole_mass_msun"]
                    ),
                    "peak_accretion_rate_msun_myr": float(
                        metadata["peak_accretion_rate_msun_myr"]
                    ),
                    "steps": int(metadata["steps"]),
                    "elapsed_seconds": float(metadata["elapsed_seconds"]),
                    "mass_budget_residual_code": float(
                        metadata["mass_budget_residual_code"]
                    ),
                    "cells": int(metadata["cells"]),
                    "r_max_pc": float(metadata["r_max_pc"]),
                    "nfw_scale_radius_pc": float(metadata["nfw_scale_radius_pc"]),
                    "nfw_concentration": float(metadata["nfw_concentration"]),
                }
            )
    controls = {
        row["halo_mass_msun"]: row["accreted_dark_matter_msun"]
        for row in rows
        if row["case_type"] == "control"
    }
    for row in rows:
        row["enhancement_over_matched_control"] = (
            row["accreted_dark_matter_msun"] / controls[row["halo_mass_msun"]]
        )
    return rows


def best_rows(rows: list[dict[str, float | int | str]]) -> list[dict]:
    optima = []
    halo_masses = sorted(
        {float(row["halo_mass_msun"]) for row in rows}
    )
    baryon_fractions = sorted(
        {
            float(row["baryon_fraction"])
            for row in rows
            if row["case_type"] != "control"
        }
    )
    for halo_mass in halo_masses:
        for baryon_fraction in baryon_fractions:
            selected = [
                row
                for row in rows
                if row["case_type"] != "control"
                and row["halo_mass_msun"] == halo_mass
                and row["baryon_fraction"] == baryon_fraction
            ]
            selected.sort(key=lambda row: row["assembly_multiplier"])
            best_index = int(
                np.argmax(
                    [row["accreted_dark_matter_msun"] for row in selected]
                )
            )
            best = dict(selected[best_index])
            bracketed = 0 < best_index < len(selected) - 1
            peak_multiplier = float(best["assembly_multiplier"])
            if bracketed:
                local = selected[best_index - 1 : best_index + 2]
                x = np.array([row["assembly_multiplier"] for row in local])
                y = np.log(
                    [row["accreted_dark_matter_msun"] for row in local]
                )
                quadratic, linear, _ = np.polyfit(x, y, 2)
                candidate = -linear / (2.0 * quadratic)
                if quadratic < 0.0 and x[0] <= candidate <= x[-1]:
                    peak_multiplier = float(candidate)
            best["quadratic_peak_bracketed"] = bracketed
            best["quadratic_peak_multiplier"] = peak_multiplier
            best["quadratic_peak_time_myr"] = (
                peak_multiplier * best["predicted_time_myr"]
            )
            optima.append(best)
    return optima


def power_law_slope(x: list[float], y: list[float]) -> float:
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


def analyze_statistics(rows: list[dict], optima: list[dict]) -> dict:
    halo_masses = sorted({row["halo_mass_msun"] for row in optima})
    baryon_fractions = sorted({row["baryon_fraction"] for row in optima})
    controls = [row for row in rows if row["case_type"] == "control"]
    controls.sort(key=lambda row: row["halo_mass_msun"])
    median_peaks = {
        baryon_fraction: float(
            np.median(
                [
                    row["quadratic_peak_multiplier"]
                    for row in optima
                    if row["baryon_fraction"] == baryon_fraction
                ]
            )
        )
        for baryon_fraction in baryon_fractions
    }
    clock_fit = np.polyfit(
        np.log10(np.array(baryon_fractions) / 0.05),
        np.log10([median_peaks[value] for value in baryon_fractions]),
        1,
    )
    return {
        "all_tasks_present": len(rows) in (39, 60),
        "maximum_mass_budget_residual_code": max(
            abs(row["mass_budget_residual_code"]) for row in rows
        ),
        "maximum_nfw_concentration_spread": (
            max(row["nfw_concentration"] for row in rows)
            - min(row["nfw_concentration"] for row in rows)
        ),
        "quadratic_peak_bracketed_fraction": float(
            np.mean([row["quadratic_peak_bracketed"] for row in optima])
        ),
        "baseline_predictor_within_25_percent_fraction": float(
            np.mean(
                [
                    0.75 <= row["quadratic_peak_multiplier"] <= 1.25
                    for row in optima
                ]
            )
        ),
        "baseline_predictor_within_25_percent_by_baryon_fraction": {
            f"{baryon_fraction:g}": float(
                np.mean(
                    [
                        0.75 <= row["quadratic_peak_multiplier"] <= 1.25
                        for row in optima
                        if row["baryon_fraction"] == baryon_fraction
                    ]
                )
            )
            for baryon_fraction in baryon_fractions
        },
        "median_quadratic_peak_multiplier_by_baryon_fraction": {
            f"{key:g}": value for key, value in median_peaks.items()
        },
        "assembly_clock_baryon_fraction_exponent": float(clock_fit[0]),
        "assembly_clock_normalization_at_fb_0p05": float(10.0 ** clock_fit[1]),
        "optimal_assembly_multipliers": [
            row["assembly_multiplier"] for row in optima
        ],
        "baryon_fraction_growth_slopes_by_halo_mass": {
            f"{halo_mass:.0e}": power_law_slope(
                baryon_fractions,
                [
                    row["accreted_dark_matter_msun"]
                    for row in optima
                    if row["halo_mass_msun"] == halo_mass
                ],
            )
            for halo_mass in halo_masses
        },
        "halo_mass_growth_slopes_by_baryon_fraction": {
            f"{baryon_fraction:g}": power_law_slope(
                halo_masses,
                [
                    row["accreted_dark_matter_msun"]
                    for row in optima
                    if row["baryon_fraction"] == baryon_fraction
                ],
            )
            for baryon_fraction in baryon_fractions
        },
        "control_halo_mass_growth_slope": power_law_slope(
            [row["halo_mass_msun"] for row in controls],
            [row["accreted_dark_matter_msun"] for row in controls],
        ),
    }


def save_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict], optima: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    halo_masses = sorted({row["halo_mass_msun"] for row in optima})
    baryon_fractions = sorted({row["baryon_fraction"] for row in optima})
    colors = {0.01: "#3b4cc0", 0.05: "#d62828", 0.16: "#2a9d8f"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2), constrained_layout=True)
    for baryon_fraction in baryon_fractions:
        selected = [
            row
            for row in rows
            if row["case_type"] != "control"
            and row["baryon_fraction"] == baryon_fraction
        ]
        multipliers = sorted({row["assembly_multiplier"] for row in selected})
        values = [
            [
                row["enhancement_over_matched_control"]
                for row in selected
                if row["assembly_multiplier"] == multiplier
            ]
            for multiplier in multipliers
        ]
        lower = np.array([min(local) for local in values])
        upper = np.array([max(local) for local in values])
        median = np.array([np.median(local) for local in values])
        axes[0, 0].fill_between(
            multipliers,
            lower,
            upper,
            color=colors[baryon_fraction],
            alpha=0.14,
        )
        axes[0, 0].plot(
            multipliers,
            median,
            marker="o",
            color=colors[baryon_fraction],
            label=f"f_b={baryon_fraction:g}",
        )
    axes[0, 0].axvline(1.0, color="black", linestyle=":", linewidth=1)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("Assembly time / predicted dynamical time")
    axes[0, 0].set_ylabel("Enhancement over matched control")
    axes[0, 0].set_title("Prediction test across the full matrix")
    axes[0, 0].legend(frameon=False, fontsize=8)

    for halo_mass in halo_masses:
        selected = [row for row in optima if row["halo_mass_msun"] == halo_mass]
        selected.sort(key=lambda row: row["baryon_fraction"])
        axes[0, 1].plot(
            [row["baryon_fraction"] for row in selected],
            [row["accreted_dark_matter_msun"] for row in selected],
            marker="o",
            label=f"M={halo_mass:.0e} M_sun",
        )
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel("Baryon fraction")
    axes[0, 1].set_ylabel("Optimal accreted dark mass [M_sun]")
    axes[0, 1].set_title("Baryon-fraction scaling")
    axes[0, 1].legend(frameon=False)

    for baryon_fraction in baryon_fractions:
        selected = [
            row for row in optima if row["baryon_fraction"] == baryon_fraction
        ]
        selected.sort(key=lambda row: row["halo_mass_msun"])
        axes[1, 0].plot(
            [row["halo_mass_msun"] for row in selected],
            [row["accreted_dark_matter_msun"] for row in selected],
            marker="o",
            color=colors[baryon_fraction],
            label=f"f_b={baryon_fraction:g}",
        )
        axes[1, 1].plot(
            [row["halo_mass_msun"] for row in selected],
            [
                row["accreted_dark_matter_msun"] / row["halo_mass_msun"]
                for row in selected
            ],
            marker="o",
            color=colors[baryon_fraction],
            label=f"f_b={baryon_fraction:g}",
        )
    for axis in axes[1]:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Halo mass [M_sun]")
        axis.legend(frameon=False)
    axes[1, 0].set_ylabel("Optimal accreted dark mass [M_sun]")
    axes[1, 0].set_title("Halo-mass scaling")
    axes[1, 1].set_ylabel("Accreted dark mass / halo mass")
    axes[1, 1].set_title("Growth efficiency")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    rows = load_rows()
    optima = best_rows(rows)
    statistics = analyze_statistics(rows, optima)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    save_csv(SUMMARY, rows)
    save_csv(OPTIMA, optima)
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    plot(rows, optima)
    print(SUMMARY)
    print(OPTIMA)
    print(STATISTICS)
    print(FIGURE)
    print(json.dumps(statistics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
