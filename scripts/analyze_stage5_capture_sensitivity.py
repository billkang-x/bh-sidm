"""Analyze seed growth sensitivity to the dark Bondi capture coefficient."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_capture_sensitivity.tsv"
RESULTS = ROOT / "results" / "stage5" / "capture_sensitivity"
SUMMARY = ROOT / "results" / "stage5" / "capture_sensitivity_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "capture_sensitivity_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_capture_sensitivity.png"


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    histories = {}
    for row in manifest:
        task_id = int(row["task_id"])
        seed = float(row["black_hole_seed_msun"])
        bondi_lambda = float(row["dark_bondi_lambda"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            histories[(seed, bondi_lambda)] = (
                data["times_myr"].copy(),
                data["black_hole_mass_msun"].copy(),
            )
        cases.append(
            {
                "task_id": task_id,
                "black_hole_seed_msun": seed,
                "reference_r_min_over_influence": float(row["reference_r_min_over_influence"]),
                "dark_bondi_lambda": bondi_lambda,
                "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "growth_factor": float(metadata["final_black_hole_mass_msun"]) / seed,
                "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
                "dark_capture_fraction": float(metadata["dark_capture_fraction_of_available_supply"]),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    by_seed = {}
    for seed in sorted({case["black_hole_seed_msun"] for case in cases}):
        pair = sorted(
            (case for case in cases if case["black_hole_seed_msun"] == seed),
            key=lambda case: case["dark_bondi_lambda"],
        )
        by_seed[f"{seed:.0f}"] = {
            "mass_at_lambda_020_msun": pair[0]["final_black_hole_mass_msun"],
            "mass_at_lambda_030_msun": pair[1]["final_black_hole_mass_msun"],
            "reaches_1e7_at_lambda_020": bool(np.isfinite(pair[0]["time_to_1e7_msun_myr"])),
            "reaches_1e7_at_lambda_030": bool(np.isfinite(pair[1]["time_to_1e7_msun_myr"])),
        }
    statistics = {
        "case_count": len(cases),
        "sensitivity_by_seed": by_seed,
        "stellar_light_seeds_fail_across_lambda_range": bool(
            all(
                not item["reaches_1e7_at_lambda_020"]
                and not item["reaches_1e7_at_lambda_030"]
                for seed, item in by_seed.items()
                if float(seed) <= 1.0e3
            )
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for bondi_lambda, marker in ((0.20, "o"), (0.30, "s")):
        selected = sorted(
            (case for case in cases if np.isclose(case["dark_bondi_lambda"], bondi_lambda)),
            key=lambda case: case["black_hole_seed_msun"],
        )
        axes[0].plot(
            [case["black_hole_seed_msun"] for case in selected],
            [case["final_black_hole_mass_msun"] for case in selected],
            marker=marker,
            label=f"lambda={bondi_lambda:.2f}",
        )
    axes[0].axhline(1.0e7, color="black", linestyle="--", alpha=0.5)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Seed mass [M_sun]")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].legend()

    for seed in (1.0e3, 3.0e4, 4.0e4):
        for bondi_lambda, linestyle in ((0.20, "--"), (0.30, "-")):
            time, mass = histories[(seed, bondi_lambda)]
            axes[1].plot(
                time,
                mass,
                linestyle=linestyle,
                label=f"seed={seed:.0e}, lambda={bondi_lambda:.2f}",
            )
    axes[1].axhline(1.0e7, color="black", linestyle=":", alpha=0.5)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time [Myr]")
    axes[1].set_ylabel("Black-hole mass [M_sun]")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
