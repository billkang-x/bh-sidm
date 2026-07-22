"""Analyze the stage-3 HPC concentration/assembly/cross-section matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage3"
MATRIX = RESULTS / "hpc_matrix"
FIGURES = RESULTS / "figures"
SCALE_RATIOS = [0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3]
ASSEMBLY_TIMES = [0.0, 0.05, 0.2, 0.5, 1.0]
CROSS_SECTIONS = [10.0, 30.0, 50.0, 100.0]


def load_records() -> list[dict]:
    records = []
    for path in sorted(MATRIX.glob("task_*.npz")):
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            records.append(
                {
                    **metadata,
                    "task_id": int(path.stem.split("_")[1]),
                    "accreted_dark_matter_msun": float(
                        data["black_hole_mass_msun"][-1] - 100.0
                    ),
                }
            )
    return records


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    records = load_records()
    controls = {
        record["sigma_over_m_cm2_g"]: record["accreted_dark_matter_msun"]
        for record in records
        if record["baryon_fraction"] == 0.0
    }
    science = [record for record in records if record["baryon_fraction"] > 0.0]
    for record in science:
        record["enhancement_over_matched_control"] = (
            record["accreted_dark_matter_msun"]
            / controls[record["sigma_over_m_cm2_g"]]
        )

    summary_fields = [
        "task_id",
        "baryon_fraction",
        "scale_radius_over_rs",
        "assembly_time_myr",
        "sigma_over_m_cm2_g",
        "final_black_hole_mass_msun",
        "accreted_dark_matter_msun",
        "enhancement_over_matched_control",
        "peak_accretion_rate_msun_myr",
        "final_accretion_rate_msun_myr",
        "steps",
        "mass_budget_residual_code",
    ]
    with (RESULTS / "hpc_matrix_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(
            {field: record[field] for field in summary_fields}
            for record in science
        )

    lookup = {
        (
            record["sigma_over_m_cm2_g"],
            record["scale_radius_over_rs"],
            record["assembly_time_myr"],
        ): record
        for record in science
    }
    optima = []
    for cross_section in CROSS_SECTIONS:
        for scale_ratio in SCALE_RATIOS:
            candidates = [
                lookup[(cross_section, scale_ratio, time)]
                for time in ASSEMBLY_TIMES
            ]
            best = max(
                candidates,
                key=lambda record: record["accreted_dark_matter_msun"],
            )
            optima.append(
                {
                    "sigma_over_m_cm2_g": cross_section,
                    "scale_radius_over_rs": scale_ratio,
                    "optimal_assembly_time_myr": best["assembly_time_myr"],
                    "maximum_enhancement": best[
                        "enhancement_over_matched_control"
                    ],
                    "maximum_accreted_dark_matter_msun": best[
                        "accreted_dark_matter_msun"
                    ],
                    "task_id": best["task_id"],
                }
            )
    with (RESULTS / "hpc_optima_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(optima[0]))
        writer.writeheader()
        writer.writerows(optima)

    all_enhancements = np.array(
        [record["enhancement_over_matched_control"] for record in science]
    )
    color_limits = (float(all_enhancements.min()), float(all_enhancements.max()))
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for axis, cross_section in zip(axes.flat, CROSS_SECTIONS, strict=True):
        matrix = np.array(
            [
                [
                    lookup[(cross_section, scale_ratio, time)][
                        "enhancement_over_matched_control"
                    ]
                    for scale_ratio in SCALE_RATIOS
                ]
                for time in ASSEMBLY_TIMES
            ]
        )
        image = axis.imshow(
            matrix,
            aspect="auto",
            origin="lower",
            cmap="viridis",
            norm=LogNorm(vmin=color_limits[0], vmax=color_limits[1]),
        )
        axis.set_xticks(range(len(SCALE_RATIOS)), [str(v) for v in SCALE_RATIOS])
        axis.set_yticks(range(len(ASSEMBLY_TIMES)), [str(v) for v in ASSEMBLY_TIMES])
        axis.set_xlabel("Hernquist scale radius / NFW scale radius")
        axis.set_ylabel("Assembly time [Myr]")
        axis.set_title(f"sigma/m = {cross_section:g} cm2/g")
        maximum = np.unravel_index(np.argmax(matrix), matrix.shape)
        axis.plot(maximum[1], maximum[0], marker="*", color="white", markersize=11)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.9)
    colorbar.set_label("Dark-growth enhancement")
    fig.savefig(FIGURES / "stage3_hpc_enhancement_surface.png", dpi=180)
    plt.close(fig)

    peak_matrix = np.empty((len(CROSS_SECTIONS), len(SCALE_RATIOS)))
    time_matrix = np.empty_like(peak_matrix)
    for i, cross_section in enumerate(CROSS_SECTIONS):
        for j, scale_ratio in enumerate(SCALE_RATIOS):
            row = next(
                item
                for item in optima
                if item["sigma_over_m_cm2_g"] == cross_section
                and item["scale_radius_over_rs"] == scale_ratio
            )
            peak_matrix[i, j] = row["maximum_enhancement"]
            time_matrix[i, j] = row["optimal_assembly_time_myr"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    peak_image = axes[0].imshow(
        peak_matrix,
        aspect="auto",
        cmap="magma",
        norm=LogNorm(vmin=peak_matrix.min(), vmax=peak_matrix.max()),
    )
    time_image = axes[1].imshow(
        time_matrix,
        aspect="auto",
        cmap="YlGnBu",
        vmin=0.0,
        vmax=1.0,
    )
    for axis in axes:
        axis.set_xticks(range(len(SCALE_RATIOS)), [str(v) for v in SCALE_RATIOS])
        axis.set_yticks(range(len(CROSS_SECTIONS)), [str(v) for v in CROSS_SECTIONS])
        axis.set_xlabel("Hernquist scale radius / NFW scale radius")
        axis.set_ylabel("sigma/m [cm2/g]")
    axes[0].set_title("Maximum enhancement over assembly time")
    axes[1].set_title("Assembly time that maximizes 2 Myr growth")
    fig.colorbar(peak_image, ax=axes[0], label="Maximum enhancement")
    fig.colorbar(time_image, ax=axes[1], label="Optimal assembly time [Myr]")
    fig.savefig(FIGURES / "stage3_hpc_optimum_map.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    axes[0].plot(
        CROSS_SECTIONS,
        [controls[value] for value in CROSS_SECTIONS],
        marker="o",
        color="#264653",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("sigma/m [cm2/g]")
    axes[0].set_ylabel("No-baryon accreted dark mass [M_sun]")
    axes[0].set_title("Matched controls are non-monotonic")
    for scale_ratio, color in zip(
        [0.01, 0.02, 0.1, 0.3],
        ["#d62828", "#f4a261", "#2a9d8f", "#3b4cc0"],
        strict=True,
    ):
        axes[1].plot(
            CROSS_SECTIONS,
            [
                max(
                    lookup[(cross_section, scale_ratio, time)][
                        "enhancement_over_matched_control"
                    ]
                    for time in ASSEMBLY_TIMES
                )
                for cross_section in CROSS_SECTIONS
            ],
            marker="o",
            color=color,
            label=f"a_b/r_s={scale_ratio}",
        )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("sigma/m [cm2/g]")
    axes[1].set_ylabel("Maximum dark-growth enhancement")
    axes[1].set_title("Cross-section response depends on concentration")
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURES / "stage3_hpc_cross_section_response.png", dpi=180)
    plt.close(fig)

    global_maximum = max(
        science,
        key=lambda record: record["enhancement_over_matched_control"],
    )
    global_minimum = min(
        science,
        key=lambda record: record["enhancement_over_matched_control"],
    )
    findings = {
        "controls_accreted_dark_matter_msun": controls,
        "global_maximum": {
            key: global_maximum[key]
            for key in (
                "task_id",
                "scale_radius_over_rs",
                "assembly_time_myr",
                "sigma_over_m_cm2_g",
                "accreted_dark_matter_msun",
                "enhancement_over_matched_control",
            )
        },
        "global_minimum": {
            key: global_minimum[key]
            for key in (
                "task_id",
                "scale_radius_over_rs",
                "assembly_time_myr",
                "sigma_over_m_cm2_g",
                "accreted_dark_matter_msun",
                "enhancement_over_matched_control",
            )
        },
    }
    (RESULTS / "hpc_key_findings.json").write_text(
        json.dumps(findings, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(RESULTS / "hpc_matrix_summary.csv")
    print(RESULTS / "hpc_optima_summary.csv")


if __name__ == "__main__":
    main()
