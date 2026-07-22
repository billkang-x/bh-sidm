"""Analyze the 10 Myr stage-5 no-feedback cosmological anchor runs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_long_pilot.tsv"
RESULTS = ROOT / "results" / "stage5" / "long_pilot"
SUMMARY = ROOT / "results" / "stage5" / "long_pilot_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "long_pilot_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_long_pilot.png"


def load_cases() -> list[dict]:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        if float(row["halo_redshift"]) == 30.0:
            refinement = RESULTS / "z30_12myr.npz"
            if refinement.exists():
                path = refinement
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            cases.append(
                {
                    "task_id": task_id,
                    "halo_redshift": float(row["halo_redshift"]),
                    "duration_myr": float(metadata["duration_myr"]),
                    "final_black_hole_mass_msun": float(
                        metadata["final_black_hole_mass_msun"]
                    ),
                    "dark_matter_accreted_msun": float(
                        metadata["dark_matter_accreted_msun"]
                    ),
                    "baryon_accreted_onto_bh_msun": float(
                        metadata["baryon_accreted_onto_bh_msun"]
                    ),
                    "dark_fraction_of_growth": float(
                        metadata["dark_fraction_of_black_hole_growth"]
                    ),
                    "time_to_1e6_msun_myr": float(
                        metadata["time_to_1e6_msun_myr"]
                    ),
                    "time_to_1e7_msun_myr": float(
                        metadata["time_to_1e7_msun_myr"]
                    ),
                    "mass_budget_residual_code": float(
                        metadata["mass_budget_residual_code"]
                    ),
                    "steps": int(metadata["steps"]),
                    "times_myr": data["times_myr"].copy(),
                    "black_hole_mass_msun": data["black_hole_mass_msun"].copy(),
                    "dark_matter_accreted_history_msun": data[
                        "dark_matter_accreted_msun"
                    ].copy(),
                    "baryon_accreted_history_msun": data[
                        "baryon_accreted_onto_bh_msun"
                    ].copy(),
                }
            )
    return cases


def save(cases: list[dict]) -> None:
    scalar_cases = [
        {key: value for key, value in case.items() if not isinstance(value, np.ndarray)}
        for case in cases
    ]
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(scalar_cases[0]))
        writer.writeheader()
        writer.writerows(scalar_cases)
    statistics = {
        "case_count": len(cases),
        "target_reached_count": {
            "1e6": int(
                sum(np.isfinite(case["time_to_1e6_msun_myr"]) for case in cases)
            ),
            "1e7": int(
                sum(np.isfinite(case["time_to_1e7_msun_myr"]) for case in cases)
            ),
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
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {10.0: "#2166ac", 20.0: "#1b9e77", 30.0: "#b2182b"}
    for case in cases:
        redshift = case["halo_redshift"]
        axes[0].plot(
            case["times_myr"],
            case["black_hole_mass_msun"],
            color=colors[redshift],
            label=f"z={redshift:g}",
        )
        axes[1].plot(
            case["times_myr"],
            case["dark_matter_accreted_history_msun"],
            color=colors[redshift],
            label=f"DM, z={redshift:g}",
        )
        axes[1].plot(
            case["times_myr"],
            case["baryon_accreted_history_msun"],
            color=colors[redshift],
            linestyle="--",
            label=f"baryon, z={redshift:g}",
        )
    axes[0].axhline(1.0e6, color="black", linestyle=":")
    axes[0].axhline(1.0e7, color="black", linestyle="--")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Time [Myr]")
    axes[0].set_ylabel("Black-hole mass [M_sun]")
    axes[0].set_title("Resolved 1e9 M_sun halo anchors")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time [Myr]")
    axes[1].set_ylabel("Accreted mass [M_sun]")
    axes[1].set_title("Growth-channel decomposition")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save(cases)
    plot(cases)
    print(SUMMARY.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
