"""Analyze the targeted stage-5 physical-frontier interaction matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_frontier_interaction.tsv"
RESULTS = ROOT / "results" / "stage5" / "frontier_interaction"
SUMMARY = ROOT / "results" / "stage5" / "frontier_interaction_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "frontier_interaction_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_frontier_interaction.png"


def main() -> None:
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
        cases.append(
            {
                "task_id": task_id,
                "halo_concentration": float(row["halo_concentration"]),
                "sigma_over_m_cm2_g": float(row["sigma_over_m_cm2_g"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "r_min_over_influence": float(
                    metadata["r_min_over_black_hole_influence_radius"]
                ),
                "baryon_scale_radius_over_r_min": float(
                    metadata["initial_baryon_scale_radius_pc"]
                    / metadata["r_min_pc"]
                ),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "scale_radius_expansion_factor": float(
                    metadata["final_baryon_scale_radius_pc"]
                    / metadata["initial_baryon_scale_radius_pc"]
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
    best_by_concentration = {
        f"c{concentration:g}": max(
            (case for case in cases if case["halo_concentration"] == concentration),
            key=lambda case: case["final_black_hole_mass_msun"],
        )
        for concentration in sorted({case["halo_concentration"] for case in cases})
    }
    statistics = {
        "case_count": len(cases),
        "global_best_case": max(
            cases, key=lambda case: case["final_black_hole_mass_msun"]
        ),
        "best_by_concentration": best_by_concentration,
        "target_reached_counts": {
            f"{target:.0e}": sum(
                case["final_black_hole_mass_msun"] >= target for case in cases
            )
            for target in (1.0e7, 3.0e7, 1.0e8)
        },
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
        "maximum_steps": max(case["steps"] for case in cases),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )

    concentrations = sorted({case["halo_concentration"] for case in cases})
    radii = sorted({case["scale_radius_over_rs"] for case in cases})
    colors = {0.0005: "#b2182b", 0.001: "#1b9e77", 0.003: "#2166ac"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)
    for axis, concentration in zip(axes.flat, concentrations):
        for radius in radii:
            selected = sorted(
                (
                    case
                    for case in cases
                    if case["halo_concentration"] == concentration
                    and case["scale_radius_over_rs"] == radius
                ),
                key=lambda case: case["sigma_over_m_cm2_g"],
            )
            axis.plot(
                [case["sigma_over_m_cm2_g"] for case in selected],
                [case["final_black_hole_mass_msun"] for case in selected],
                color=colors[radius], marker="o", label=f"a_b/r_s={radius:g}",
            )
        axis.set_xscale("log")
        axis.axhline(1.0e7, color="black", linestyle=":")
        axis.set_xlabel("sigma/m [cm2/g]")
        axis.set_ylabel("Final black-hole mass [M_sun]")
        axis.set_title(f"c={concentration:g}")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
