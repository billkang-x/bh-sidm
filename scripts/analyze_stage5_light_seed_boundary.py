"""Analyze the fixed-halo light-seed and feeding-boundary closure matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_light_seed_boundary.tsv"
RESULTS = ROOT / "results" / "stage5" / "light_seed_boundary"
SUMMARY = ROOT / "results" / "stage5" / "light_seed_boundary_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "light_seed_boundary_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_light_seed_boundary.png"
TARGETS = (1.0e5, 1.0e6, 1.0e7)


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def finite_or_none(value: float) -> float | None:
    return value if np.isfinite(value) else None


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    histories = {}
    for row in manifest:
        task_id = int(row["task_id"])
        seed = float(row["black_hole_seed_msun"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            if row["axis"] == "inner_boundary" and np.isclose(
                float(row["reference_r_min_over_influence"]),
                0.05208333333333334,
            ):
                histories[seed] = {
                    "time": data["times_myr"].copy(),
                    "mass": data["black_hole_mass_msun"].copy(),
                }
        case = {
            "task_id": task_id,
            "black_hole_seed_msun": seed,
            "axis": row["axis"],
            "reference_r_min_over_influence": float(row["reference_r_min_over_influence"]),
            "seed_r_min_over_influence": float(row["seed_r_min_over_influence"]),
            "cells": int(row["cells"]),
            "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
            "growth_factor": float(metadata["final_black_hole_mass_msun"]) / seed,
            "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
            "dark_matter_supplied_msun": float(metadata["dark_matter_supplied_to_inner_boundary_msun"]),
            "final_dark_reservoir_msun": float(metadata["final_inner_dark_matter_reservoir_msun"]),
            "dark_capture_fraction": float(metadata["dark_capture_fraction_of_available_supply"]),
            "baryon_accreted_msun": float(metadata["baryon_accreted_onto_bh_msun"]),
            "dark_fraction_of_growth": float(metadata["dark_fraction_of_black_hole_growth"]),
            "time_to_1e5_msun_myr": float(metadata["time_to_1e5_msun_myr"]),
            "time_to_1e6_msun_myr": float(metadata["time_to_1e6_msun_myr"]),
            "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
            "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
            "steps": int(metadata["steps"]),
        }
        cases.append(case)

    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    by_seed = {}
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
        grid_base = next(
            case
            for case in boundary
            if np.isclose(case["reference_r_min_over_influence"], 0.10416666666666667)
        )
        boundary_difference = relative_difference(
            boundary[0]["final_black_hole_mass_msun"],
            boundary[1]["final_black_hole_mass_msun"],
        )
        grid_difference = relative_difference(
            grid_base["final_black_hole_mass_msun"],
            grid["final_black_hole_mass_msun"],
        )
        conservative = boundary[0]
        by_seed[f"{seed:.0f}"] = {
            "smallest_boundary_final_mass_msun": conservative["final_black_hole_mass_msun"],
            "smallest_boundary_growth_factor": conservative["growth_factor"],
            "smallest_boundary_pair_relative_difference": boundary_difference,
            "grid_relative_difference": grid_difference,
            "passes_five_percent_convergence": bool(
                boundary_difference < 0.05 and grid_difference < 0.05
            ),
            "dark_capture_fraction": conservative["dark_capture_fraction"],
            "dark_fraction_of_growth": conservative["dark_fraction_of_growth"],
            "time_to_1e5_msun_myr": finite_or_none(conservative["time_to_1e5_msun_myr"]),
            "time_to_1e6_msun_myr": finite_or_none(conservative["time_to_1e6_msun_myr"]),
            "time_to_1e7_msun_myr": finite_or_none(conservative["time_to_1e7_msun_myr"]),
        }

    statistics = {
        "case_count": len(cases),
        "convergence_by_seed": by_seed,
        "all_seed_masses_pass_five_percent_convergence": bool(
            all(item["passes_five_percent_convergence"] for item in by_seed.values())
        ),
        "light_seed_reaches_1e7_within_2_myr": bool(
            any(
                float(seed) <= 1.0e3 and item["time_to_1e7_msun_myr"] is not None
                for seed, item in by_seed.items()
            )
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(by_seed)))
    for color, seed in zip(colors, sorted(float(value) for value in by_seed)):
        boundary = sorted(
            (
                case
                for case in cases
                if case["black_hole_seed_msun"] == seed
                and case["axis"] == "inner_boundary"
            ),
            key=lambda case: case["reference_r_min_over_influence"],
        )
        axes[0, 0].plot(
            [case["reference_r_min_over_influence"] for case in boundary],
            [case["final_black_hole_mass_msun"] for case in boundary],
            marker="o", color=color, label=f"seed={seed:.0e}",
        )
        conservative = boundary[0]
        axes[0, 1].scatter(seed, conservative["growth_factor"], color=color, s=45)
        axes[1, 0].scatter(
            seed,
            conservative["dark_capture_fraction"],
            color=color,
            s=45,
        )
        history = histories[seed]
        axes[1, 1].plot(history["time"], history["mass"], color=color, label=f"seed={seed:.0e}")
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("Feeding radius / reference influence radius")
    axes[0, 0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0, 0].legend()
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel("Seed mass [M_sun]")
    axes[0, 1].set_ylabel("Conservative growth factor")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_xlabel("Seed mass [M_sun]")
    axes[1, 0].set_ylabel("Captured / available dark supply")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Time [Myr]")
    axes[1, 1].set_ylabel("Black-hole mass [M_sun]")
    for target in TARGETS:
        axes[1, 1].axhline(target, color="black", alpha=0.2, linestyle="--")
    axes[1, 1].legend()
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
