"""Combine the trusted-frontier compactness closure and peak refinement."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_SUMMARY = ROOT / "results" / "stage5" / "compactness_closure_summary.csv"
MANIFEST = ROOT / "hpc" / "stage5_compactness_peak.tsv"
RESULTS = ROOT / "results" / "stage5" / "compactness_peak"
SUMMARY = ROOT / "results" / "stage5" / "compactness_peak_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "compactness_peak_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_compactness_peak.png"


def main() -> None:
    with BASE_SUMMARY.open(newline="", encoding="ascii") as stream:
        cases = [
            {
                "source": "closure",
                "task_id": int(row["task_id"]),
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "final_black_hole_mass_msun": float(row["final_black_hole_mass_msun"]),
                "dark_matter_accreted_msun": float(row["dark_matter_accreted_msun"]),
                "mass_budget_residual_code": float(row["mass_budget_residual_code"]),
                "steps": int(row["steps"]),
            }
            for row in csv.DictReader(stream)
        ]
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    for row in manifest:
        task_id = int(row["task_id"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        cases.append(
            {
                "source": "peak_refinement",
                "task_id": task_id,
                "scale_radius_over_rs": float(row["scale_radius_over_rs"]),
                "final_black_hole_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "dark_matter_accreted_msun": float(metadata["dark_matter_accreted_msun"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    cases.sort(key=lambda case: case["scale_radius_over_rs"])
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
    fig, axis = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    axis.plot(
        [case["scale_radius_over_rs"] for case in cases],
        [case["final_black_hole_mass_msun"] for case in cases],
        color="#1b9e77", marker="o",
    )
    axis.set_xscale("log")
    axis.set_xlabel("a_b/r_s")
    axis.set_ylabel("Final black-hole mass [M_sun]")
    axis.set_title("Fiducial-model compactness scan")
    axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
