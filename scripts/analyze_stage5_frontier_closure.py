"""Analyze the stage-5 one-axis physical-frontier closure matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_frontier_closure.tsv"
RESULTS = ROOT / "results" / "stage5" / "frontier_closure"
SUMMARY = ROOT / "results" / "stage5" / "frontier_closure_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "frontier_closure_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_frontier_closure.png"


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
            eddington = np.flatnonzero(~data["bondi_limited"])
            first_eddington = (
                float(data["times_myr"][eddington[0]])
                if len(eddington)
                else float("nan")
            )
            final_bondi_to_eddington = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
        cases.append(
            {
                "task_id": task_id,
                "axis": row["axis"],
                "axis_value": float(row["axis_value"]),
                "halo_concentration": float(row["halo_concentration"]),
                "sigma_over_m_cm2_g": float(row["sigma_over_m_cm2_g"]),
                "baryon_fraction": float(row["baryon_fraction"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "r_min_over_rs": float(row["r_min_over_rs"]),
                "r_min_over_influence": float(
                    metadata["r_min_over_black_hole_influence_radius"]
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
                "first_eddington_time_myr": first_eddington,
                "final_bondi_to_eddington_ratio": final_bondi_to_eddington,
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
    return cases


def main() -> None:
    cases = load_cases()
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    axis_parameters = {
        "concentration": "halo_concentration",
        "cross_section": "sigma_over_m_cm2_g",
        "scale_radius": "scale_radius_over_rs",
        "baryon_fraction": "baryon_fraction",
    }
    best_by_axis = {
        axis: max(
            (case for case in cases if case["axis"] in (axis, "baseline")),
            key=lambda case: case["final_black_hole_mass_msun"],
        )
        for axis in axis_parameters
    }
    statistics = {
        "case_count": len(cases),
        "best_by_axis": best_by_axis,
        "global_best_case": max(
            cases, key=lambda case: case["final_black_hole_mass_msun"]
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
        "maximum_steps": max(case["steps"] for case in cases),
        "r_min_over_influence_range": [
            min(case["r_min_over_influence"] for case in cases),
            max(case["r_min_over_influence"] for case in cases),
        ],
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )

    axes_order = (
        ("concentration", "halo_concentration", "Halo concentration"),
        ("cross_section", "sigma_over_m_cm2_g", "sigma/m [cm2/g]"),
        ("scale_radius", "scale_radius_over_rs", "a_b/r_s"),
        ("baryon_fraction", "baryon_fraction", "Baryon fraction"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)
    colors = ("#b2182b", "#2166ac", "#1b9e77", "#6a3d9a")
    for axis, (name, parameter, label), color in zip(axes.flat, axes_order, colors):
        selected = sorted(
            (case for case in cases if case["axis"] in (name, "baseline")),
            key=lambda case: case[parameter],
        )
        axis.plot(
            [case[parameter] for case in selected],
            [case["final_black_hole_mass_msun"] for case in selected],
            color=color,
            marker="o",
        )
        if name in ("cross_section", "scale_radius"):
            axis.set_xscale("log")
        axis.axhline(1.0e7, color="black", linestyle=":")
        axis.set_xlabel(label)
        axis.set_ylabel("Final black-hole mass [M_sun]")
        axis.set_title(name.replace("_", " ").title())
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
