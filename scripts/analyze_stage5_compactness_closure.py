"""Analyze the final compactness closure around the trusted frontier."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_compactness_closure.tsv"
RESULTS = ROOT / "results" / "stage5" / "compactness_closure"
BASELINE = ROOT / "results" / "stage5" / "frontier_interaction" / "task_014.npz"
SUMMARY = ROOT / "results" / "stage5" / "compactness_closure_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "compactness_closure_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_compactness_closure.png"


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    with np.load(BASELINE, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
    cases = [
        {
            "task_id": -1,
            "scale_radius_over_rs": 0.003,
            "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
            "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
            "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
            "steps": int(metadata["steps"]),
        }
    ]
    for row in manifest:
        task_id = int(row["task_id"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        cases.append(
            {
                "task_id": task_id,
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    cases.sort(key=lambda case: case["scale_radius_over_rs"])
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    best = max(cases, key=lambda case: case["final_black_hole_mass_msun"])
    statistics = {
        "case_count": len(cases),
        "best_case": best,
        "compactness_peak_is_internal": bool(
            best is not cases[0] and best is not cases[-1]
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )
    fig, axis = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    axis.plot(
        [case["scale_radius_over_rs"] for case in cases],
        [case["final_black_hole_mass_msun"] for case in cases],
        color="#1b9e77", marker="o",
    )
    axis.set_xscale("log")
    axis.set_xlabel("a_b/r_s")
    axis.set_ylabel("Final black-hole mass [M_sun]")
    axis.set_title("Trusted-frontier compactness closure")
    axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
