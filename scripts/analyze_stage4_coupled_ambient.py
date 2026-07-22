"""Analyze feedback-coupled Hernquist density and gas sound speed."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_coupled_ambient.tsv"
RESULTS = ROOT / "results" / "stage4" / "coupled_ambient"
SUMMARY = ROOT / "results" / "stage4" / "coupled_ambient_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "coupled_ambient_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_coupled_ambient.png"


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
            bondi_limited = data["bondi_limited"]
            assembled = times >= metadata["assembly_time_myr"]
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
            indicator = eddington_limited.astype(float)
            eddington_duration = float(
                np.sum(
                    0.5
                    * (indicator[1:] + indicator[:-1])
                    * np.diff(times)
                )
            )
            final_limiter_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
        cases.append(
            {
                "task_id": task_id,
                "ambient_model": row["ambient_model"],
                "feedback_efficiency": float(row["feedback_efficiency"]),
                "feedback_heating_fraction": float(
                    row["feedback_heating_fraction"]
                ),
                "feedback_expansion_fraction": float(
                    metadata["feedback_expansion_fraction"]
                ),
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
                    metadata["final_baryon_scale_radius_pc"]
                    / metadata["initial_baryon_scale_radius_pc"]
                ),
                "final_ambient_density_msun_pc3": float(
                    metadata["final_ambient_gas_density_msun_pc3"]
                ),
                "final_ambient_sound_speed_km_s": float(
                    metadata["final_ambient_gas_sound_speed_km_s"]
                ),
                "bondi_limited_sample_fraction": float(
                    metadata["post_assembly_bondi_limited_sample_fraction"]
                ),
                "transition_time_myr": float(
                    metadata["bondi_to_eddington_transition_myr"]
                ),
                "ever_eddington_after_assembly": ever_eddington,
                "limiter_history": limiter_history,
                "bondi_reclosure_time_myr": reclosure_time,
                "eddington_limited_duration_myr": eddington_duration,
                "final_bondi_to_eddington_ratio": final_limiter_ratio,
                "feedback_to_binding_ratio": float(
                    metadata["final_feedback_to_binding_ratio"]
                ),
                "total_feedback_to_binding_ratio": float(
                    metadata["final_total_feedback_to_binding_ratio"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "path": path,
            }
        )
    return cases


def save_summary(cases: list[dict]) -> None:
    rows = [{key: value for key, value in case.items() if key != "path"} for case in cases]
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prevention_threshold(cases: list[dict], heating_fraction: float) -> float | None:
    selected = sorted(
        (
            case
            for case in cases
            if case["ambient_model"] == "evolving"
            and case["feedback_heating_fraction"] == heating_fraction
            and case["feedback_efficiency"] > 0.0
        ),
        key=lambda case: case["feedback_efficiency"],
    )
    prevented = [case for case in selected if not case["ever_eddington_after_assembly"]]
    return prevented[0]["feedback_efficiency"] if prevented else None


def prevention_bracket(cases: list[dict], heating_fraction: float) -> list[float | None]:
    selected = sorted(
        (
            case
            for case in cases
            if case["ambient_model"] == "evolving"
            and case["feedback_heating_fraction"] == heating_fraction
            and case["feedback_efficiency"] > 0.0
        ),
        key=lambda case: case["feedback_efficiency"],
    )
    threshold = prevention_threshold(cases, heating_fraction)
    if threshold is None:
        return [selected[-1]["feedback_efficiency"], None]
    lower = max(
        (
            case["feedback_efficiency"]
            for case in selected
            if case["feedback_efficiency"] < threshold
            and case["ever_eddington_after_assembly"]
        ),
        default=None,
    )
    return [lower, threshold]


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {0.0: "#2166ac", 0.5: "#1b9e77", 1.0: "#d73027"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for heating_fraction in (0.0, 0.5, 1.0):
        for ambient_model, linestyle in (("fixed", "--"), ("evolving", "-")):
            selected = sorted(
                (
                    case
                    for case in cases
                    if case["ambient_model"] == ambient_model
                    and case["feedback_heating_fraction"] == heating_fraction
                    and case["feedback_efficiency"] > 0.0
                ),
                key=lambda case: case["feedback_efficiency"],
            )
            label = f"heat={heating_fraction:g}, {ambient_model}"
            axes[0, 0].plot(
                [case["feedback_efficiency"] for case in selected],
                [case["final_bondi_to_eddington_ratio"] for case in selected],
                color=colors[heating_fraction],
                linestyle=linestyle,
                marker="o",
                label=label,
            )
            axes[1, 1].plot(
                [case["feedback_efficiency"] for case in selected],
                [case["dark_matter_accreted_msun"] for case in selected],
                color=colors[heating_fraction],
                linestyle=linestyle,
                marker="o",
                label=label,
            )
        evolving = sorted(
            (
                case
                for case in cases
                if case["ambient_model"] == "evolving"
                and case["feedback_heating_fraction"] == heating_fraction
                and case["feedback_efficiency"] > 0.0
            ),
            key=lambda case: case["feedback_efficiency"],
        )
        axes[0, 1].plot(
            [case["feedback_efficiency"] for case in evolving],
            [case["final_ambient_density_msun_pc3"] / 300.0 for case in evolving],
            color=colors[heating_fraction],
            marker="o",
            label=f"heat={heating_fraction:g}",
        )
        axes[1, 0].plot(
            [case["feedback_efficiency"] for case in evolving],
            [case["final_ambient_sound_speed_km_s"] for case in evolving],
            color=colors[heating_fraction],
            marker="o",
            label=f"heat={heating_fraction:g}",
        )
    axes[0, 0].axhline(1.0, color="black", linestyle=":")
    axes[0, 0].set_ylabel("Final Bondi / Eddington limit")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_title("Does the limiter transition survive?")
    axes[0, 1].set_ylabel("Final ambient density / initial")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("Homologous dilution")
    axes[1, 0].set_ylabel("Final gas sound speed [km/s]")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Feedback heating")
    axes[1, 1].set_ylabel("DM accreted in 2 Myr [M_sun]")
    axes[1, 1].set_title("Impact on dark accretion")
    for axis in axes.flat:
        axis.set_xscale("log")
        axis.set_xlabel("Feedback efficiency")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save_summary(cases)
    plot(cases)
    dynamic = [case for case in cases if case["ambient_model"] == "evolving"]
    dynamic_control = next(
        case for case in dynamic if case["feedback_efficiency"] == 0.0
    )
    threshold_cases = {}
    for fraction in (0.0, 0.5, 1.0):
        threshold = prevention_threshold(cases, fraction)
        threshold_case = next(
            case
            for case in dynamic
            if case["feedback_heating_fraction"] == fraction
            and case["feedback_efficiency"] == threshold
        )
        threshold_cases[str(fraction)] = {
            "feedback_efficiency": threshold,
            "final_density_fraction": (
                threshold_case["final_ambient_density_msun_pc3"] / 300.0
            ),
            "final_sound_speed_km_s": threshold_case[
                "final_ambient_sound_speed_km_s"
            ],
            "dark_matter_suppression_fraction": 1.0
            - threshold_case["dark_matter_accreted_msun"]
            / dynamic_control["dark_matter_accreted_msun"],
            "baryon_accreted_onto_bh_msun": threshold_case[
                "baryon_accreted_onto_bh_msun"
            ],
        }
    statistics = {
        "case_count": len(cases),
        "transition_prevention_threshold_by_heating_fraction": {
            str(fraction): prevention_threshold(cases, fraction)
            for fraction in (0.0, 0.5, 1.0)
        },
        "transition_prevention_bracket_by_heating_fraction": {
            str(fraction): prevention_bracket(cases, fraction)
            for fraction in (0.0, 0.5, 1.0)
        },
        "threshold_case_diagnostics": threshold_cases,
        "transient_eddington_cases": [
            {
                "task_id": case["task_id"],
                "feedback_efficiency": case["feedback_efficiency"],
                "feedback_heating_fraction": case[
                    "feedback_heating_fraction"
                ],
                "transition_time_myr": case["transition_time_myr"],
                "reclosure_time_myr": case["bondi_reclosure_time_myr"],
            }
            for case in dynamic
            if case["limiter_history"] == "transient_eddington"
        ],
        "minimum_final_density_fraction": min(
            case["final_ambient_density_msun_pc3"] / 300.0 for case in dynamic
        ),
        "maximum_final_sound_speed_km_s": max(
            case["final_ambient_sound_speed_km_s"] for case in dynamic
        ),
        "dynamic_dark_matter_accreted_range_msun": [
            min(case["dark_matter_accreted_msun"] for case in dynamic),
            max(case["dark_matter_accreted_msun"] for case in dynamic),
        ],
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(statistics, sort_keys=True))


if __name__ == "__main__":
    main()
