"""Analyze Bondi supply with effective-binding expansion feedback."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_effective_feedback.tsv"
RESULTS = ROOT / "results" / "stage4" / "effective_feedback"
NO_FEEDBACK_SUMMARY = ROOT / "results" / "stage4" / "no_feedback_summary.csv"
SUMMARY = ROOT / "results" / "stage4" / "effective_feedback_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "effective_feedback_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_effective_feedback.png"


def load_cases() -> tuple[list[dict], float]:
    with NO_FEEDBACK_SUMMARY.open(newline="", encoding="utf-8") as stream:
        no_feedback = list(csv.DictReader(stream))
    no_eddington_control = next(
        float(row["dark_matter_accreted_msun"])
        for row in no_feedback
        if float(row["seed_mass_msun"]) == 100.0
        and float(row["eddington_ratio"]) == 0.0
    )
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
            final_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
        case = {
            "task_id": task_id,
            "gas_regime": row["gas_regime"],
            "gas_density_msun_pc3": float(row["gas_density_msun_pc3"]),
            "feedback_efficiency": float(row["feedback_efficiency"]),
            "feedback_eta": float(row["feedback_eta"]),
            "final_black_hole_mass_msun": float(
                metadata["final_black_hole_mass_msun"]
            ),
            "dark_matter_accreted_msun": float(
                metadata["dark_matter_accreted_msun"]
            ),
            "baryon_accreted_onto_bh_msun": float(
                metadata["baryon_accreted_onto_bh_msun"]
            ),
            "scale_radius_expansion_factor": float(
                metadata["final_baryon_scale_radius_pc"] / 0.3
            ),
            "feedback_to_binding_ratio": float(
                metadata["final_feedback_to_binding_ratio"]
            ),
            "feedback_binding_energy_erg": float(
                metadata["feedback_binding_energy_erg"]
            ),
            "bondi_limited_sample_fraction": float(
                metadata["post_assembly_bondi_limited_sample_fraction"]
            ),
            "transition_time_myr": float(
                metadata["bondi_to_eddington_transition_myr"]
            ),
            "final_bondi_to_eddington_ratio": final_ratio,
            "mass_budget_residual_code": float(
                metadata["mass_budget_residual_code"]
            ),
            "path": path,
        }
        case["extra_dark_matter_vs_no_eddington_msun"] = (
            case["dark_matter_accreted_msun"] - no_eddington_control
        )
        cases.append(case)
    return cases, no_eddington_control


def save_summary(cases: list[dict]) -> None:
    rows = [{key: value for key, value in case.items() if key != "path"} for case in cases]
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        ("transition", 0.5): "#3b4cc0",
        ("transition", 1.0): "#7b2cbf",
        ("eddington_saturated", 0.5): "#2a9d8f",
        ("eddington_saturated", 1.0): "#d62828",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for regime in ("transition", "eddington_saturated"):
        reference = next(
            case
            for case in cases
            if case["gas_regime"] == regime
            and case["feedback_efficiency"] == 0.0
        )
        for eta in (0.5, 1.0):
            selected = [reference] + sorted(
                [
                    case
                    for case in cases
                    if case["gas_regime"] == regime
                    and case["feedback_efficiency"] > 0.0
                    and case["feedback_eta"] == eta
                ],
                key=lambda case: case["feedback_efficiency"],
            )
            x = [max(case["feedback_efficiency"], 3.0e-7) for case in selected]
            label = f"{regime}, eta={eta:g}"
            axes[0, 0].plot(
                x,
                [case["extra_dark_matter_vs_no_eddington_msun"] for case in selected],
                marker="o",
                color=colors[(regime, eta)],
                label=label,
            )
            axes[0, 1].plot(
                x,
                [case["scale_radius_expansion_factor"] for case in selected],
                marker="o",
                color=colors[(regime, eta)],
                label=label,
            )
            axes[1, 0].plot(
                x,
                [case["bondi_limited_sample_fraction"] for case in selected],
                marker="o",
                color=colors[(regime, eta)],
                label=label,
            )
    axes[0, 0].axhline(0.0, color="black", linestyle=":")
    axes[0, 0].set_ylabel("Extra DM vs no-Eddington [M_sun]")
    axes[0, 0].set_title("Effective binding shifts feedback reversal")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_ylabel("Final Hernquist radius / initial")
    axes[0, 1].set_title("Potential expansion")
    axes[1, 0].set_ylabel("Post-assembly Bondi-limited fraction")
    axes[1, 0].set_title("Feedback delays the limiter transition")
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xscale("log")
        axis.set_xlabel("Feedback efficiency")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)

    representative = [
        next(
            case
            for case in cases
            if case["gas_regime"] == "transition"
            and case["feedback_efficiency"] == 0.0
        ),
        next(
            case
            for case in cases
            if case["gas_regime"] == "transition"
            and case["feedback_efficiency"] == 3.0e-6
            and case["feedback_eta"] == 1.0
        ),
        next(
            case
            for case in cases
            if case["gas_regime"] == "transition"
            and case["feedback_efficiency"] == 3.0e-5
            and case["feedback_eta"] == 1.0
        ),
    ]
    for case in representative:
        with np.load(case["path"], allow_pickle=False) as data:
            axes[1, 1].plot(
                data["times_myr"],
                data["bondi_baryon_growth_limit_msun_myr"]
                / data["eddington_baryon_growth_limit_msun_myr"],
                label=f"eps_f={case['feedback_efficiency']:g}",
            )
    axes[1, 1].axhline(1.0, color="black", linestyle=":")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Time [Myr]")
    axes[1, 1].set_ylabel("Bondi / Eddington limit")
    axes[1, 1].set_title("Transition-regime limiter history, eta=1")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def crossing_threshold(cases: list[dict], regime: str, eta: float) -> float | None:
    selected = sorted(
        [
            case
            for case in cases
            if case["gas_regime"] == regime
            and case["feedback_efficiency"] > 0.0
            and case["feedback_eta"] == eta
        ],
        key=lambda case: case["feedback_efficiency"],
    )
    for lower, upper in zip(selected[:-1], selected[1:]):
        low_value = lower["extra_dark_matter_vs_no_eddington_msun"]
        high_value = upper["extra_dark_matter_vs_no_eddington_msun"]
        if low_value >= 0.0 and high_value < 0.0:
            fraction = low_value / (low_value - high_value)
            return float(
                10.0
                ** (
                    np.log10(lower["feedback_efficiency"])
                    + fraction
                    * np.log10(
                        upper["feedback_efficiency"]
                        / lower["feedback_efficiency"]
                    )
                )
            )
    return None


def main() -> None:
    cases, _ = load_cases()
    save_summary(cases)
    plot(cases)
    thresholds = {
        regime: {
            str(eta): crossing_threshold(cases, regime, eta)
            for eta in (0.5, 1.0)
        }
        for regime in ("transition", "eddington_saturated")
    }
    transition_prevention = {}
    for eta in (0.5, 1.0):
        prevented = sorted(
            [
                case
                for case in cases
                if case["gas_regime"] == "transition"
                and case["feedback_efficiency"] > 0.0
                and case["feedback_eta"] == eta
                and case["final_bondi_to_eddington_ratio"] < 1.0
            ],
            key=lambda case: case["feedback_efficiency"],
        )
        transition_prevention[str(eta)] = (
            prevented[0]["feedback_efficiency"] if prevented else None
        )
    statistics = {
        "all_cases_present": len(cases) == 18,
        "effective_to_self_binding_ratio": (
            cases[0]["feedback_binding_energy_erg"] / 1.187847719850361e50
        ),
        "feedback_reversal_log_interpolated_thresholds": thresholds,
        "sampled_transition_prevention_thresholds": transition_prevention,
        "maximum_mass_budget_residual_code": max(
            abs(case["mass_budget_residual_code"]) for case in cases
        ),
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
