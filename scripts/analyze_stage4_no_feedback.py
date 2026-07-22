"""Analyze stage-4 Eddington-limited baryon accretion without feedback."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_no_feedback.tsv"
RESULTS = ROOT / "results" / "stage4" / "no_feedback"
SUMMARY = ROOT / "results" / "stage4" / "no_feedback_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "no_feedback_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_no_feedback.png"


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
            case = {
                "task_id": task_id,
                "seed_mass_msun": float(row["black_hole_mass_msun"]),
                "eddington_ratio": float(row["eddington_ratio"]),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_onto_bh_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "baryon_gas_consumed_msun": float(
                    metadata["baryon_gas_consumed_msun"]
                ),
                "dark_fraction_of_black_hole_growth": float(
                    metadata["dark_fraction_of_black_hole_growth"]
                ),
                "dark_dominated_onset_myr": float(
                    metadata["dark_dominated_onset_myr"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
                "elapsed_seconds": float(metadata["elapsed_seconds"]),
                "path": path,
                "eddington_efolding_time_myr": float(
                    metadata["eddington_black_hole_efolding_time_myr"]
                ),
                "duration_myr": float(metadata["duration_myr"]),
            }
            cases.append(case)
    controls = {
        case["seed_mass_msun"]: case
        for case in cases
        if case["eddington_ratio"] == 0.0
    }
    for case in cases:
        control = controls[case["seed_mass_msun"]]
        case["extra_dark_matter_vs_no_eddington_msun"] = (
            case["dark_matter_accreted_msun"]
            - control["dark_matter_accreted_msun"]
        )
        baryon_growth = case["baryon_accreted_onto_bh_msun"]
        case["catalytic_dark_mass_per_baryon_mass"] = (
            case["extra_dark_matter_vs_no_eddington_msun"] / baryon_growth
            if baryon_growth > 0.0
            else float("nan")
        )
        if case["eddington_ratio"] > 0.0:
            isolated_efolding = case["eddington_efolding_time_myr"]
            isolated_growth = case["seed_mass_msun"] * (
                np.exp(case["duration_myr"] / isolated_efolding) - 1.0
            )
            case["baryon_growth_boost_over_isolated_eddington"] = (
                baryon_growth / isolated_growth
            )
        else:
            case["baryon_growth_boost_over_isolated_eddington"] = float("nan")
    return cases


def save_summary(cases: list[dict]) -> None:
    rows = [{key: value for key, value in case.items() if key != "path"} for case in cases]
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {10.0: "#3b4cc0", 100.0: "#d62828", 1000.0: "#2a9d8f"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for case in cases:
        if case["eddington_ratio"] != 1.0:
            continue
        with np.load(case["path"], allow_pickle=False) as data:
            axes[0, 0].plot(
                data["times_myr"],
                data["black_hole_mass_msun"],
                color=colors[case["seed_mass_msun"]],
                label=f"seed={case['seed_mass_msun']:g} M_sun",
            )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("Time [Myr]")
    axes[0, 0].set_ylabel("Black-hole mass [M_sun]")
    axes[0, 0].set_title("Full-Eddington total growth")
    axes[0, 0].legend(frameon=False)

    seed_cases = sorted(
        [case for case in cases if case["seed_mass_msun"] == 100.0],
        key=lambda case: case["eddington_ratio"],
    )
    for case in seed_cases:
        with np.load(case["path"], allow_pickle=False) as data:
            axes[0, 1].plot(
                data["times_myr"],
                data["dark_matter_accreted_msun"],
                label=f"DM, f_Edd={case['eddington_ratio']:g}",
            )
            if case["eddington_ratio"] > 0.0:
                axes[0, 1].plot(
                    data["times_myr"],
                    data["baryon_accreted_onto_bh_msun"],
                    linestyle="--",
                    label=f"baryon, f_Edd={case['eddington_ratio']:g}",
                )
    axes[0, 1].set_yscale("symlog", linthresh=0.1)
    axes[0, 1].set_xlabel("Time [Myr]")
    axes[0, 1].set_ylabel("Cumulative accreted mass [M_sun]")
    axes[0, 1].set_title("Separated channels for 100 M_sun seed")
    axes[0, 1].legend(frameon=False, fontsize=8)

    eddington_ratios = sorted({case["eddington_ratio"] for case in cases})
    for seed in sorted(colors):
        selected = sorted(
            [case for case in cases if case["seed_mass_msun"] == seed],
            key=lambda case: case["eddington_ratio"],
        )
        axes[1, 0].plot(
            eddington_ratios,
            [
                case["extra_dark_matter_vs_no_eddington_msun"]
                for case in selected
            ],
            marker="o",
            color=colors[seed],
            label=f"seed={seed:g} M_sun",
        )
    axes[1, 0].set_yscale("symlog", linthresh=0.1)
    axes[1, 0].set_xlabel("Eddington ratio")
    axes[1, 0].set_ylabel("Extra dark matter vs no-Eddington [M_sun]")
    axes[1, 0].set_title("Baryonic growth catalyzes dark inflow")
    axes[1, 0].legend(frameon=False)

    active = [case for case in cases if case["eddington_ratio"] > 0.0]
    positions = np.arange(len(active))
    axes[1, 1].bar(
        positions - 0.18,
        [case["catalytic_dark_mass_per_baryon_mass"] for case in active],
        width=0.36,
        label="extra DM / baryon BH mass",
        color="#3b4cc0",
    )
    axes[1, 1].bar(
        positions + 0.18,
        [case["baryon_growth_boost_over_isolated_eddington"] for case in active],
        width=0.36,
        label="baryon / isolated Eddington",
        color="#f4a261",
    )
    axes[1, 1].set_xticks(
        positions,
        [
            f"{case['seed_mass_msun']:g}\n{case['eddington_ratio']:g}"
            for case in active
        ],
    )
    axes[1, 1].set_xlabel("Seed mass / Eddington ratio")
    axes[1, 1].set_ylabel("Coupling factor")
    axes[1, 1].set_title("Two-way growth coupling")
    axes[1, 1].legend(frameon=False, fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save_summary(cases)
    plot(cases)
    active = [case for case in cases if case["eddington_ratio"] > 0.0]
    statistics = {
        "all_cases_present": len(cases) == 9,
        "maximum_mass_budget_residual_code": max(
            abs(case["mass_budget_residual_code"]) for case in cases
        ),
        "dark_fraction_of_growth_range_full_eddington": [
            min(
                case["dark_fraction_of_black_hole_growth"]
                for case in cases
                if case["eddington_ratio"] == 1.0
            ),
            max(
                case["dark_fraction_of_black_hole_growth"]
                for case in cases
                if case["eddington_ratio"] == 1.0
            ),
        ],
        "catalytic_dark_mass_per_baryon_mass_range": [
            min(case["catalytic_dark_mass_per_baryon_mass"] for case in active),
            max(case["catalytic_dark_mass_per_baryon_mass"] for case in active),
        ],
        "baryon_growth_boost_over_isolated_eddington_range": [
            min(
                case["baryon_growth_boost_over_isolated_eddington"]
                for case in active
            ),
            max(
                case["baryon_growth_boost_over_isolated_eddington"]
                for case in active
            ),
        ],
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(SUMMARY)
    print(STATISTICS)
    print(FIGURE)
    print(json.dumps(statistics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
