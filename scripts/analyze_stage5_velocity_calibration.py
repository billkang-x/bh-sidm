"""Analyze the velocity-dependent SIDM calibration and boundary audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sidm_bh.sidm import maxwellian_viscosity_cross_section_ratio


MANIFEST = ROOT / "hpc" / "stage5_velocity_calibration.tsv"
RESULTS = ROOT / "results" / "stage5" / "velocity_calibration"
SUMMARY = ROOT / "results" / "stage5" / "velocity_calibration_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "velocity_calibration_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_velocity_calibration.png"
TARGET_MASS_MSUN = 1.0e7
VIRIAL_VELOCITY_KM_S = 65.48612036054037


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    histories = {}
    for row in manifest:
        task_id = int(row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            if np.isclose(
                float(row["r_min_over_reference_influence"]),
                5.0 / 48.0,
            ):
                histories[row["model_label"]] = {
                    "time": data["times_myr"].copy(),
                    "inner_sigma": data[
                        "lmfp_effective_sigma_over_m_cm2_g"
                    ][:, 0].copy(),
                }
        velocity = float(row["velocity_scale_km_s"])
        sigma0 = float(row["sigma0_over_m_cm2_g"])
        effective_virial = sigma0
        effective_virial_smfp = sigma0
        if row["cross_section_model"] == "rutherford":
            dispersion_ratio = VIRIAL_VELOCITY_KM_S / np.sqrt(2.0) / velocity
            effective_virial *= maxwellian_viscosity_cross_section_ratio(
                dispersion_ratio,
                velocity_power=3,
            )
            effective_virial_smfp *= maxwellian_viscosity_cross_section_ratio(
                dispersion_ratio,
                velocity_power=5,
            )
        cases.append(
            {
                "task_id": task_id,
                "model_label": row["model_label"],
                "cross_section_model": row["cross_section_model"],
                "sigma0_over_m_cm2_g": sigma0,
                "velocity_scale_km_s": velocity,
                "effective_sigma_at_virial_dispersion_cm2_g": effective_virial,
                "smfp_effective_sigma_at_virial_dispersion_cm2_g": (
                    effective_virial_smfp
                ),
                "r_min_over_reference_influence": float(
                    row["r_min_over_reference_influence"]
                ),
                "cells": int(row["cells"]),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "initial_inner_effective_sigma_cm2_g": float(
                    metadata["initial_inner_lmfp_effective_sigma_over_m_cm2_g"]
                ),
                "final_inner_effective_sigma_cm2_g": float(
                    metadata["final_inner_lmfp_effective_sigma_over_m_cm2_g"]
                ),
                "initial_inner_smfp_effective_sigma_cm2_g": float(
                    metadata["initial_inner_smfp_effective_sigma_over_m_cm2_g"]
                ),
                "final_inner_smfp_effective_sigma_cm2_g": float(
                    metadata["final_inner_smfp_effective_sigma_over_m_cm2_g"]
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

    by_model = {}
    for label in sorted({case["model_label"] for case in cases}):
        pair = sorted(
            (case for case in cases if case["model_label"] == label),
            key=lambda case: case["r_min_over_reference_influence"],
        )
        conservative, screen = pair
        difference = relative_difference(
            conservative["final_black_hole_mass_msun"],
            screen["final_black_hole_mass_msun"],
        )
        by_model[label] = {
            "cross_section_model": conservative["cross_section_model"],
            "sigma0_over_m_cm2_g": conservative["sigma0_over_m_cm2_g"],
            "velocity_scale_km_s": conservative["velocity_scale_km_s"],
            "effective_sigma_at_virial_dispersion_cm2_g": conservative[
                "effective_sigma_at_virial_dispersion_cm2_g"
            ],
            "smfp_effective_sigma_at_virial_dispersion_cm2_g": conservative[
                "smfp_effective_sigma_at_virial_dispersion_cm2_g"
            ],
            "conservative_final_mass_msun": conservative[
                "final_black_hole_mass_msun"
            ],
            "conservative_time_to_1e7_msun_myr": (
                conservative["time_to_1e7_msun_myr"]
                if np.isfinite(conservative["time_to_1e7_msun_myr"])
                else None
            ),
            "boundary_relative_difference": difference,
            "passes_five_percent_boundary_test": bool(difference < 0.05),
            "robustly_crosses_1e7": bool(
                conservative["final_black_hole_mass_msun"] >= TARGET_MASS_MSUN
                and screen["final_black_hole_mass_msun"] >= TARGET_MASS_MSUN
            ),
        }

    statistics = {
        "case_count": len(cases),
        "model_count": len(by_model),
        "target_mass_msun": TARGET_MASS_MSUN,
        "models": by_model,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )

    sigma_values = (1.0, 3.0, 10.0, 30.0, 100.0)
    velocity_values = (10.0, 30.0, 100.0, 300.0)
    mass_grid = np.empty((len(sigma_values), len(velocity_values)))
    boundary_grid = np.empty_like(mass_grid)
    for i, sigma0 in enumerate(sigma_values):
        for j, velocity in enumerate(velocity_values):
            item = by_model[f"rutherford_sigma{sigma0:g}_w{velocity:g}"]
            mass_grid[i, j] = item["conservative_final_mass_msun"]
            boundary_grid[i, j] = item["boundary_relative_difference"]

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    image = axes[0, 0].imshow(
        np.log10(mass_grid), origin="lower", aspect="auto", cmap="viridis"
    )
    axes[0, 0].set_xticks(range(len(velocity_values)), velocity_values)
    axes[0, 0].set_yticks(range(len(sigma_values)), sigma_values)
    axes[0, 0].set_xlabel("Rutherford transition speed [km/s]")
    axes[0, 0].set_ylabel("Low-speed sigma0/m [cm2/g]")
    axes[0, 0].set_title("Conservative final mass")
    fig.colorbar(image, ax=axes[0, 0], label="log10 final mass [M_sun]")

    image = axes[0, 1].imshow(
        100.0 * boundary_grid,
        origin="lower",
        aspect="auto",
        cmap="magma",
    )
    axes[0, 1].set_xticks(range(len(velocity_values)), velocity_values)
    axes[0, 1].set_yticks(range(len(sigma_values)), sigma_values)
    axes[0, 1].set_xlabel("Rutherford transition speed [km/s]")
    axes[0, 1].set_ylabel("Low-speed sigma0/m [cm2/g]")
    axes[0, 1].set_title("Inner-boundary sensitivity")
    fig.colorbar(image, ax=axes[0, 1], label="Relative difference [%]")

    conservative_cases = [
        case
        for case in cases
        if np.isclose(case["r_min_over_reference_influence"], 5.0 / 48.0)
    ]
    for model, marker, color in (
        ("constant", "s", "black"),
        ("rutherford", "o", "#1f77b4"),
    ):
        selected = [
            case for case in conservative_cases if case["cross_section_model"] == model
        ]
        axes[1, 0].scatter(
            [case["effective_sigma_at_virial_dispersion_cm2_g"] for case in selected],
            [case["final_black_hole_mass_msun"] for case in selected],
            marker=marker,
            color=color,
            alpha=0.8,
            label=model,
        )
    axes[1, 0].axhline(TARGET_MASS_MSUN, color="red", linestyle="--")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel("LMFP sigma0 K3/m at Vvir/sqrt(2) [cm2/g]")
    axes[1, 0].set_ylabel("Final black-hole mass [M_sun]")
    axes[1, 0].set_title("Virial-scale matching is not the full closure")
    axes[1, 0].legend()

    labels = (
        ("constant_sigma1", r"Constant $\sigma/m=1$"),
        ("rutherford_sigma30_w10", r"Rutherford $w=10$ km/s"),
        ("rutherford_sigma30_w30", r"Rutherford $w=30$ km/s"),
        ("rutherford_sigma30_w100", r"Rutherford $w=100$ km/s"),
    )
    colors = ("black", "#d62728", "#2ca02c", "#1f77b4")
    for (label, display_label), color in zip(labels, colors, strict=True):
        history = histories[label]
        axes[1, 1].plot(
            history["time"],
            history["inner_sigma"],
            color=color,
            label=display_label,
        )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Time [Myr]")
    axes[1, 1].set_ylabel("Inner LMFP sigma0 K3/m [cm2/g]")
    axes[1, 1].set_title("Local velocity dependence during growth")
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
