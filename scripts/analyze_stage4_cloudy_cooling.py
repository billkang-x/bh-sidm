"""Analyze the Cloudy cooling and cooling-suppression stage-4 matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_cloudy_cooling.tsv"
RESULTS = ROOT / "results" / "stage4" / "cloudy_cooling"
SUMMARY = ROOT / "results" / "stage4" / "cloudy_cooling_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "cloudy_cooling_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_cloudy_cooling.png"


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
            ever_eddington = bool(np.any(assembled & ~bondi_limited))
            reclosure_time = float("nan")
            if ever_eddington:
                first = int(np.flatnonzero(assembled & ~bondi_limited)[0])
                for index in range(first + 1, len(times)):
                    if np.all(bondi_limited[index:]):
                        reclosure_time = float(times[index])
                        break
            if not ever_eddington:
                history = "bondi_throughout"
            elif np.isfinite(reclosure_time):
                history = "transient_eddington"
            else:
                history = "eddington_at_end"
            final_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
        cases.append(
            {
                "task_id": task_id,
                "configuration": row["configuration"],
                "feedback_heating_fraction": float(
                    row["feedback_heating_fraction"]
                ),
                "feedback_efficiency": float(row["feedback_efficiency"]),
                "metallicity_solar": float(row["metallicity_solar"]),
                "cooling_rate_multiplier": float(
                    row["cooling_rate_multiplier"]
                ),
                "limiter_history": history,
                "transition_time_myr": float(
                    metadata["bondi_to_eddington_transition_myr"]
                ),
                "reclosure_time_myr": reclosure_time,
                "final_bondi_to_eddington_ratio": final_ratio,
                "final_sound_speed_km_s": float(
                    metadata["final_ambient_gas_sound_speed_km_s"]
                ),
                "final_gas_density_msun_pc3": float(
                    metadata["final_ambient_gas_density_msun_pc3"]
                ),
                "final_temperature_k": float(
                    metadata["final_ambient_gas_temperature_k"]
                ),
                "final_cooling_time_myr": float(
                    metadata["final_physical_cooling_time_myr"]
                ),
                "minimum_cooling_time_myr": float(
                    metadata["minimum_physical_cooling_time_myr"]
                ),
                "retained_thermal_energy_fraction": float(
                    metadata["retained_thermal_energy_fraction"]
                ),
                "scale_radius_expansion_factor": float(
                    metadata["final_baryon_scale_radius_pc"]
                    / metadata["initial_baryon_scale_radius_pc"]
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


def save_summary(cases: list[dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)


def suppression_bracket(
    cases: list[dict], configuration: str, metallicity: float
) -> list[float | None]:
    selected = sorted(
        (
            case
            for case in cases
            if case["configuration"] == configuration
            and case["metallicity_solar"] == metallicity
        ),
        key=lambda case: case["cooling_rate_multiplier"],
        reverse=True,
    )
    prevented = [case for case in selected if case["limiter_history"] == "bondi_throughout"]
    if not prevented:
        return [selected[-1]["cooling_rate_multiplier"], None]
    first = prevented[0]["cooling_rate_multiplier"]
    transitioning = [
        case["cooling_rate_multiplier"]
        for case in selected
        if case["cooling_rate_multiplier"] > first
        and case["limiter_history"] != "bondi_throughout"
    ]
    return [first, min(transitioning) if transitioning else None]


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "heating_strong": "#d73027",
        "heating_extreme": "#8b1a1a",
        "mixed_strong": "#1b9e77",
        "mixed_extreme": "#006d5b",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for configuration, color in colors.items():
        selected = sorted(
            (
                case
                for case in cases
                if case["configuration"] == configuration
                and case["cooling_rate_multiplier"] == 1.0
            ),
            key=lambda case: case["metallicity_solar"],
        )
        axes[0, 0].plot(
            [case["metallicity_solar"] for case in selected],
            [case["final_bondi_to_eddington_ratio"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
        axes[0, 1].plot(
            [case["metallicity_solar"] for case in selected],
            [case["final_sound_speed_km_s"] for case in selected],
            color=color,
            marker="o",
            label=configuration,
        )
    for configuration in ("heating_extreme", "mixed_extreme"):
        for metallicity, linestyle in ((0.0, "-"), (1.0, "--")):
            selected = sorted(
                (
                    case
                    for case in cases
                    if case["configuration"] == configuration
                    and case["metallicity_solar"] == metallicity
                ),
                key=lambda case: case["cooling_rate_multiplier"],
            )
            label = f"{configuration}, Z={metallicity:g}"
            axes[1, 0].plot(
                [case["cooling_rate_multiplier"] for case in selected],
                [case["final_bondi_to_eddington_ratio"] for case in selected],
                color=colors[configuration],
                linestyle=linestyle,
                marker="o",
                label=label,
            )
            axes[1, 1].plot(
                [case["cooling_rate_multiplier"] for case in selected],
                [case["dark_matter_accreted_msun"] for case in selected],
                color=colors[configuration],
                linestyle=linestyle,
                marker="o",
                label=label,
            )
    axes[0, 0].axhline(1.0, color="black", linestyle=":")
    axes[0, 0].set_xscale("symlog", linthresh=1.0e-4)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel("Final Bondi / Eddington limit")
    axes[0, 0].set_title("Optically thin metallicity response")
    axes[0, 1].set_xscale("symlog", linthresh=1.0e-4)
    axes[0, 1].set_ylabel("Final sound speed [km/s]")
    axes[0, 1].set_title("Cloudy thermal balance")
    axes[1, 0].axhline(1.0, color="black", linestyle=":")
    axes[1, 0].set_xscale("log")
    axes[1, 0].invert_xaxis()
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylabel("Final Bondi / Eddington limit")
    axes[1, 0].set_title("Required cooling suppression")
    axes[1, 1].set_xscale("log")
    axes[1, 1].invert_xaxis()
    axes[1, 1].set_ylabel("DM accreted in 2 Myr [M_sun]")
    axes[1, 1].set_title("Dark-channel response")
    for row in range(2):
        for column in range(2):
            axis = axes[row, column]
            axis.set_xlabel(
                "Z / Z_sun" if row == 0 else "Cooling-rate multiplier"
            )
            axis.grid(alpha=0.25)
            axis.legend(frameon=False, fontsize=7)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save_summary(cases)
    plot(cases)
    optically_thin = [
        case for case in cases if case["cooling_rate_multiplier"] == 1.0
    ]
    statistics = {
        "case_count": len(cases),
        "optically_thin_case_count": len(optically_thin),
        "optically_thin_prevention_count": sum(
            case["limiter_history"] == "bondi_throughout"
            for case in optically_thin
        ),
        "optically_thin_prevention_by_configuration": {
            configuration: sum(
                case["configuration"] == configuration
                and case["limiter_history"] == "bondi_throughout"
                for case in optically_thin
            )
            for configuration in (
                "no_feedback_control",
                "heating_strong",
                "heating_extreme",
                "mixed_strong",
                "mixed_extreme",
            )
        },
        "cooling_suppression_bracket": {
            f"{configuration}_Z{metallicity:g}": suppression_bracket(
                cases, configuration, metallicity
            )
            for configuration in ("heating_extreme", "mixed_extreme")
            for metallicity in (0.0, 1.0)
        },
        "minimum_physical_cooling_time_myr": min(
            case["minimum_cooling_time_myr"] for case in cases
        ),
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
