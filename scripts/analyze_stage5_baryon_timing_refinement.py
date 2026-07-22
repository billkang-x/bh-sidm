"""Analyze the stage-5 high-baryon assembly-time refinement."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_baryon_timing_refinement.tsv"
RESULTS = ROOT / "results" / "stage5" / "baryon_timing_refinement"
SUMMARY = ROOT / "results" / "stage5" / "baryon_timing_refinement_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "baryon_timing_refinement_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_baryon_timing_refinement.png"


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    rows = []
    histories = []
    for manifest_row in manifest:
        task_id = int(manifest_row["task_id"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            row = {
                "task_id": task_id,
                "halo_redshift": float(manifest_row["halo_redshift"]),
                "assembly_time_myr": float(manifest_row["assembly_time_myr"]),
                "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
                "baryon_accreted_onto_bh_msun": float(metadata["baryon_accreted_onto_bh_msun"]),
                "dark_fraction_of_growth": float(metadata["dark_fraction_of_black_hole_growth"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
            rows.append(row)
            histories.append(
                (
                    row,
                    data["times_myr"].copy(),
                    data["black_hole_mass_msun"].copy(),
                )
            )
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    best = {
        f"z{redshift:g}": max(
            (row for row in rows if row["halo_redshift"] == redshift),
            key=lambda row: row["final_black_hole_mass_msun"],
        )
        for redshift in (20.0, 30.0)
    }
    statistics = {
        "case_count": len(rows),
        "best_by_redshift": best,
        "target_1e7_reached_count": sum(
            row["final_black_hole_mass_msun"] >= 1.0e7 for row in rows
        ),
        "maximum_mass_budget_residual_code": max(
            row["mass_budget_residual_code"] for row in rows
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {20.0: "#1b9e77", 30.0: "#b2182b"}
    for redshift in (20.0, 30.0):
        selected = sorted(
            (row for row in rows if row["halo_redshift"] == redshift),
            key=lambda row: row["assembly_time_myr"],
        )
        axes[0].plot(
            [row["assembly_time_myr"] for row in selected],
            [row["final_black_hole_mass_msun"] for row in selected],
            color=colors[redshift],
            marker="o",
            label=f"z={redshift:g}",
        )
    for row, times, masses in histories:
        if row in best.values():
            axes[1].plot(
                times,
                masses,
                color=colors[row["halo_redshift"]],
                label=(
                    f"z={row['halo_redshift']:g}, "
                    f"Tasm={row['assembly_time_myr']:g} Myr"
                ),
            )
    axes[0].axhline(1.0e7, color="black", linestyle=":")
    axes[0].set_xlabel("Assembly time [Myr]")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].set_title("Assembly-time refinement")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time [Myr]")
    axes[1].set_ylabel("Black-hole mass [M_sun]")
    axes[1].set_title("Best refined histories")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
