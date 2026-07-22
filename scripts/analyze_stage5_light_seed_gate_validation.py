"""Validate the influence-gated closure against resolved direct-flux runs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_light_seed_gate_validation.tsv"
RESULTS = ROOT / "results" / "stage5" / "light_seed_gate_validation"
SUMMARY = ROOT / "results" / "stage5" / "light_seed_gate_validation_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "light_seed_gate_validation_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_light_seed_gate_validation.png"
DIRECT = {
    0.05208333333333334: ROOT / "results" / "stage5" / "trusted_peak_convergence" / "task_000.npz",
    0.10416666666666667: ROOT / "results" / "stage5" / "trusted_peak_convergence" / "task_001.npz",
    5.0 / 24.0: ROOT / "results" / "stage5" / "compactness_peak" / "task_002.npz",
}


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    direct = {}
    for boundary, path in DIRECT.items():
        with np.load(path, allow_pickle=False) as data:
            direct[boundary] = json.loads(str(data["metadata_json"]))

    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        feeding_boundary = float(row["reference_r_min_over_influence"])
        flux_boundary = float(row["flux_capture_r_min_over_influence"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            masses = data["black_hole_mass_msun"]
            times = data["times_myr"]
        threshold = float(metadata["dark_flux_capture_mass_threshold_msun"])
        reached = np.flatnonzero(masses >= threshold)
        cases.append(
            {
                "task_id": task_id,
                "reference_r_min_over_influence": feeding_boundary,
                "flux_capture_r_min_over_influence": flux_boundary,
                "flux_capture_mass_threshold_msun": threshold,
                "flux_capture_onset_myr": float(times[reached[0]]) if len(reached) else float("nan"),
                "direct_final_mass_msun": float(direct[feeding_boundary]["final_black_hole_mass_msun"]),
                "gated_final_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "final_mass_relative_difference": relative_difference(
                    float(direct[feeding_boundary]["final_black_hole_mass_msun"]),
                    float(metadata["final_black_hole_mass_msun"]),
                ),
                "final_dark_reservoir_msun": float(metadata["final_inner_dark_matter_reservoir_msun"]),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    nominal = sorted(
        (case for case in cases if np.isclose(case["flux_capture_r_min_over_influence"], 0.10416666666666667)),
        key=lambda case: case["reference_r_min_over_influence"],
    )
    resolved_nominal = nominal[:2]
    sensitivity = [
        case
        for case in cases
        if np.isclose(case["reference_r_min_over_influence"], 0.10416666666666667)
    ]
    statistics = {
        "case_count": len(cases),
        "resolved_smallest_two_maximum_relative_difference": max(
            case["final_mass_relative_difference"] for case in resolved_nominal
        ),
        "passes_five_percent_resolved_validation": bool(
            all(case["final_mass_relative_difference"] < 0.05 for case in resolved_nominal)
        ),
        "switch_threshold_mass_spread_fraction": (
            max(case["gated_final_mass_msun"] for case in sensitivity)
            - min(case["gated_final_mass_msun"] for case in sensitivity)
        )
        / np.mean([case["gated_final_mass_msun"] for case in sensitivity]),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    axes[0].plot(
        [case["reference_r_min_over_influence"] for case in nominal],
        [case["direct_final_mass_msun"] for case in nominal],
        marker="o", label="resolved boundary flux",
    )
    axes[0].plot(
        [case["reference_r_min_over_influence"] for case in nominal],
        [case["gated_final_mass_msun"] for case in nominal],
        marker="s", label="influence-gated",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Feeding radius / reference influence radius")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].legend()
    sensitivity.sort(key=lambda case: case["flux_capture_r_min_over_influence"])
    axes[1].plot(
        [case["flux_capture_r_min_over_influence"] for case in sensitivity],
        [case["gated_final_mass_msun"] for case in sensitivity],
        marker="o",
    )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Flux-capture switch r_min/r_influence")
    axes[1].set_ylabel("Final black-hole mass [M_sun]")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
