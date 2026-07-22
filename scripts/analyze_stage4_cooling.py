"""Analyze finite cooling of the stage-4 Bondi feedback reservoir."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_cooling.tsv"
RESULTS = ROOT / "results" / "stage4" / "cooling"
SUMMARY = ROOT / "results" / "stage4" / "cooling_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "cooling_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_cooling.png"


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
            times = data["times_myr"]
            assembled = times >= metadata["assembly_time_myr"]
            bondi_limited = data["bondi_limited"]
            eddington_limited = assembled & ~bondi_limited
            ever_eddington = bool(np.any(eddington_limited))
            reclosure_time = float("nan")
            if ever_eddington:
                first_eddington = int(np.flatnonzero(eddington_limited)[0])
                for index in range(first_eddington + 1, len(times)):
                    if np.all(bondi_limited[index:]):
                        reclosure_time = float(times[index])
                        break
            if not ever_eddington:
                limiter_history = "bondi_throughout"
            elif np.isfinite(reclosure_time):
                limiter_history = "transient_eddington"
            else:
                limiter_history = "eddington_at_end"
            final_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
            peak_sound_speed = float(
                np.max(data["ambient_gas_sound_speed_km_s"])
            )
        cooling_label = row["cooling_time_myr"]
        cooling_time = (
            float("inf") if cooling_label == "inf" else float(cooling_label)
        )
        cases.append(
            {
                "task_id": task_id,
                "configuration": row["configuration"],
                "feedback_heating_fraction": float(
                    row["feedback_heating_fraction"]
                ),
                "feedback_efficiency": float(row["feedback_efficiency"]),
                "cooling_time_myr": cooling_time,
                "ever_eddington_after_assembly": ever_eddington,
                "limiter_history": limiter_history,
                "bondi_reclosure_time_myr": reclosure_time,
                "transition_time_myr": float(
                    metadata["bondi_to_eddington_transition_myr"]
                ),
                "final_bondi_to_eddington_ratio": final_ratio,
                "final_sound_speed_km_s": float(
                    metadata["final_ambient_gas_sound_speed_km_s"]
                ),
                "peak_sound_speed_km_s": peak_sound_speed,
                "final_density_fraction": float(
                    metadata["final_ambient_gas_density_msun_pc3"] / 300.0
                ),
                "retained_thermal_energy_fraction": float(
                    metadata["retained_thermal_energy_fraction"]
                ),
                "final_thermal_energy_erg": float(
                    metadata["final_baryon_thermal_energy_erg"]
                ),
                "cooling_loss_energy_erg": float(
                    metadata["cumulative_cooling_loss_energy_erg"]
                ),
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_onto_bh_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
            }
        )
    return cases


def cooling_bracket(cases: list[dict], configuration: str) -> list[float | None]:
    selected = sorted(
        (case for case in cases if case["configuration"] == configuration),
        key=lambda case: case["cooling_time_myr"],
    )
    prevented = [case for case in selected if not case["ever_eddington_after_assembly"]]
    if not prevented:
        return [selected[-1]["cooling_time_myr"], None]
    upper = prevented[0]["cooling_time_myr"]
    lower = max(
        (
            case["cooling_time_myr"]
            for case in selected
            if case["cooling_time_myr"] < upper
            and case["ever_eddington_after_assembly"]
        ),
        default=None,
    )
    return [lower, upper]


def save_summary(cases: list[dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "mixed_threshold": "#1b9e77",
        "mixed_strong": "#006d5b",
        "heating_threshold": "#d73027",
        "heating_strong": "#8b1a1a",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for configuration, color in colors.items():
        selected = sorted(
            (case for case in cases if case["configuration"] == configuration),
            key=lambda case: case["cooling_time_myr"],
        )
        x = [
            300.0 if np.isinf(case["cooling_time_myr"]) else case["cooling_time_myr"]
            for case in selected
        ]
        axes[0, 0].plot(
            x,
            [case["final_bondi_to_eddington_ratio"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
        axes[0, 1].plot(
            x,
            [case["final_sound_speed_km_s"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
        axes[1, 0].plot(
            x,
            [case["retained_thermal_energy_fraction"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
        axes[1, 1].plot(
            x,
            [case["dark_matter_accreted_msun"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
    axes[0, 0].axhline(1.0, color="black", linestyle=":")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel("Final Bondi / Eddington limit")
    axes[0, 0].set_title("Cooling restores the limiter transition")
    axes[0, 1].set_ylabel("Final gas sound speed [km/s]")
    axes[0, 1].set_title("Thermal support")
    axes[1, 0].set_ylabel("Retained / injected thermal energy")
    axes[1, 0].set_title("Leaky thermal reservoir")
    axes[1, 1].set_ylabel("DM accreted in 2 Myr [M_sun]")
    axes[1, 1].set_title("Dark-channel response")
    for axis in axes.flat:
        axis.set_xscale("log")
        axis.set_xlabel("Cooling time [Myr]; 300 = no cooling")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save_summary(cases)
    plot(cases)
    configurations = sorted({case["configuration"] for case in cases})
    threshold_diagnostics = {}
    for configuration in configurations:
        lower, upper = cooling_bracket(cases, configuration)
        selected = [
            case for case in cases if case["configuration"] == configuration
        ]
        lower_case = next(
            case for case in selected if case["cooling_time_myr"] == lower
        )
        upper_case = next(
            case for case in selected if case["cooling_time_myr"] == upper
        )
        threshold_diagnostics[configuration] = {
            "last_transitioning_cooling_time_myr": lower,
            "first_preventing_cooling_time_myr": upper,
            "last_transition_time_myr": lower_case["transition_time_myr"],
            "last_reclosure_time_myr": (
                lower_case["bondi_reclosure_time_myr"]
                if np.isfinite(lower_case["bondi_reclosure_time_myr"])
                else None
            ),
            "first_preventing_final_sound_speed_km_s": upper_case[
                "final_sound_speed_km_s"
            ],
            "first_preventing_retained_thermal_energy_fraction": upper_case[
                "retained_thermal_energy_fraction"
            ],
            "first_preventing_dark_matter_accreted_msun": upper_case[
                "dark_matter_accreted_msun"
            ],
        }
    statistics = {
        "case_count": len(cases),
        "transition_prevention_cooling_time_bracket_myr": {
            configuration: cooling_bracket(cases, configuration)
            for configuration in configurations
        },
        "threshold_diagnostics": threshold_diagnostics,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
        "retained_thermal_energy_fraction_range": [
            min(case["retained_thermal_energy_fraction"] for case in cases),
            max(case["retained_thermal_energy_fraction"] for case in cases),
        ],
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(statistics, sort_keys=True))


if __name__ == "__main__":
    main()
