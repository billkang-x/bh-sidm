"""Summarize and plot the matched static-Hernquist experiment matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage3"
FIGURES = RESULTS / "figures"
FRACTIONS = [0.01, 0.05, 0.16]
SCALE_RATIOS = [0.001, 0.01, 0.1]


def load(path: Path) -> dict:
    with np.load(path) as data:
        result = {name: data[name].copy() for name in data.files}
    result["metadata"] = json.loads(str(result.pop("metadata_json")))
    return result


def matrix_path(fraction: float, scale_ratio: float) -> Path:
    return RESULTS / f"fb{fraction}_a{scale_ratio}_rmin0.005.npz"


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    control = load(RESULTS / "control_fb0_rmin0.005.npz")
    control_delta = control["black_hole_mass_msun"][-1] - 100.0
    enhancement = np.empty((len(FRACTIONS), len(SCALE_RATIOS)))
    rows = []
    matrix = {}
    for row_index, fraction in enumerate(FRACTIONS):
        for column_index, scale_ratio in enumerate(SCALE_RATIOS):
            data = load(matrix_path(fraction, scale_ratio))
            delta_mass = data["black_hole_mass_msun"][-1] - 100.0
            factor = delta_mass / control_delta
            enhancement[row_index, column_index] = factor
            matrix[(fraction, scale_ratio)] = data
            rows.append(
                {
                    "baryon_fraction": fraction,
                    "scale_radius_over_rs": scale_ratio,
                    "final_black_hole_mass_msun": data["black_hole_mass_msun"][-1],
                    "accreted_dark_matter_msun": delta_mass,
                    "dark_growth_enhancement": factor,
                    "peak_accretion_rate_msun_myr": data["metadata"][
                        "peak_accretion_rate_msun_myr"
                    ],
                    "final_accretion_rate_msun_myr": data["metadata"][
                        "final_accretion_rate_msun_myr"
                    ],
                }
            )

    with (RESULTS / "static_matrix_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axis = plt.subplots(figsize=(7, 4.8), constrained_layout=True)
    image = axis.imshow(np.log10(enhancement), cmap="RdYlBu_r", aspect="auto")
    for i in range(enhancement.shape[0]):
        for j in range(enhancement.shape[1]):
            axis.text(j, i, f"{enhancement[i, j]:.2f}x", ha="center", va="center")
    axis.set_xticks(range(len(SCALE_RATIOS)), [str(value) for value in SCALE_RATIOS])
    axis.set_yticks(range(len(FRACTIONS)), [str(value) for value in FRACTIONS])
    axis.set_xlabel("Hernquist scale radius / NFW scale radius")
    axis.set_ylabel("Baryon fraction")
    axis.set_title("Static-baryon enhancement of accreted dark mass")
    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label("log10 enhancement")
    fig.savefig(FIGURES / "stage3_static_enhancement.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(FRACTIONS)))
    for axis, scale_ratio in zip(axes, SCALE_RATIOS, strict=True):
        axis.plot(
            control["times_myr"],
            control["black_hole_mass_msun"],
            "k--",
            label="No baryons",
        )
        for fraction, color in zip(FRACTIONS, colors, strict=True):
            data = matrix[(fraction, scale_ratio)]
            axis.semilogy(
                data["times_myr"],
                data["black_hole_mass_msun"],
                color=color,
                label=f"f_b={fraction}",
            )
        axis.set_title(f"a_b/r_s={scale_ratio}")
        axis.set_xlabel("Time [Myr]")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Black-hole mass [M_sun]")
    axes[-1].legend(frameon=False)
    fig.savefig(FIGURES / "stage3_static_growth.png", dpi=180)
    plt.close(fig)

    boundary_data = {
        "Enhanced: f_b=0.05, a/r_s=0.01": [],
        "Suppressed: f_b=0.05, a/r_s=0.1": [],
    }
    for r_min in (0.0025, 0.005, 0.01):
        if r_min == 0.005:
            control_case = control
            enhanced_case = matrix[(0.05, 0.01)]
            suppressed_case = matrix[(0.05, 0.1)]
        else:
            suffix = str(r_min)
            control_case = load(RESULTS / f"control_fb0_a0.01_rmin{suffix}.npz")
            enhanced_case = load(
                RESULTS / f"enhance_fb0.05_a0.01_rmin{suffix}.npz"
            )
            suppressed_case = load(
                RESULTS / f"suppress_fb0.05_a0.1_rmin{suffix}.npz"
            )
        control_growth = control_case["black_hole_mass_msun"][-1] - 100.0
        boundary_data["Enhanced: f_b=0.05, a/r_s=0.01"].append(
            (enhanced_case["black_hole_mass_msun"][-1] - 100.0) / control_growth
        )
        boundary_data["Suppressed: f_b=0.05, a/r_s=0.1"].append(
            (suppressed_case["black_hole_mass_msun"][-1] - 100.0) / control_growth
        )
    fig, axis = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    for label, values in boundary_data.items():
        axis.plot([0.0025, 0.005, 0.01], values, marker="o", label=label)
    axis.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Inner boundary [pc]")
    axis.set_ylabel("Dark-growth enhancement")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    fig.savefig(FIGURES / "stage3_boundary_robustness.png", dpi=180)
    plt.close(fig)
    print(RESULTS.resolve())


if __name__ == "__main__":
    main()
