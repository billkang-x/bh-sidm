"""Analyze the baryon screen on the stage-5 transport frontier."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_baryon_frontier.tsv"
RESULTS = ROOT / "results" / "stage5" / "baryon_frontier"
SUMMARY = ROOT / "results" / "stage5" / "baryon_frontier_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "baryon_frontier_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_baryon_frontier.png"


def load_cases() -> list[dict]:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        final_mass = float(metadata["final_black_hole_mass_msun"])
        cases.append(
            {
                "task_id": task_id,
                "halo_redshift": float(row["halo_redshift"]),
                "baryon_fraction": float(row["baryon_fraction"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "assembly_time_myr": float(row["assembly_time_myr"]),
                "final_black_hole_mass_msun": final_mass,
                "black_hole_growth_factor": final_mass / 1.0e5,
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_onto_bh_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "dark_fraction_of_growth": float(
                    metadata["dark_fraction_of_black_hole_growth"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )
    return cases


def save(cases: list[dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    best_by_redshift = {
        f"z{redshift:g}": max(
            (case for case in cases if case["halo_redshift"] == redshift),
            key=lambda case: case["final_black_hole_mass_msun"],
        )
        for redshift in sorted({case["halo_redshift"] for case in cases})
    }
    statistics = {
        "case_count": len(cases),
        "best_by_redshift": best_by_redshift,
        "global_best_case": max(
            cases, key=lambda case: case["final_black_hole_mass_msun"]
        ),
        "target_reached_counts": {
            f"{target:.0e}": sum(
                case["final_black_hole_mass_msun"] >= target
                for case in cases
            )
            for target in (1.0e6, 1.0e7)
        },
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    redshifts = sorted({case["halo_redshift"] for case in cases})
    baryon_fractions = sorted({case["baryon_fraction"] for case in cases})
    colors = {0.003: "#2166ac", 0.01: "#1b9e77", 0.03: "#b2182b"}
    fig, axes = plt.subplots(3, 3, figsize=(12, 11), constrained_layout=True)
    for row_index, redshift in enumerate(redshifts):
        for column_index, baryon_fraction in enumerate(baryon_fractions):
            axis = axes[row_index, column_index]
            for scale_radius in sorted(
                {case["scale_radius_over_rs"] for case in cases}
            ):
                selected = sorted(
                    (
                        case
                        for case in cases
                        if case["halo_redshift"] == redshift
                        and case["baryon_fraction"] == baryon_fraction
                        and case["scale_radius_over_rs"] == scale_radius
                    ),
                    key=lambda case: case["assembly_time_myr"],
                )
                axis.plot(
                    [case["assembly_time_myr"] for case in selected],
                    [case["black_hole_growth_factor"] for case in selected],
                    color=colors[scale_radius],
                    marker="o",
                    label=f"a_b/r_s={scale_radius:g}",
                )
            axis.set_yscale("log")
            axis.grid(alpha=0.25)
            axis.set_title(f"z={redshift:g}, f_b={baryon_fraction:g}")
            if row_index == 2:
                axis.set_xlabel("Assembly time [Myr]")
            if column_index == 0:
                axis.set_ylabel("2 Myr growth factor")
            if row_index == 0 and column_index == 0:
                axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save(cases)
    plot(cases)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
