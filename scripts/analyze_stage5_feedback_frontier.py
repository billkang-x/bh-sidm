"""Analyze physical cooling and mixed feedback on the stage-5 frontier."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_feedback_frontier.tsv"
RESULTS = ROOT / "results" / "stage5" / "feedback_frontier"
CEILING_PATHS = {
    1.0: ROOT / "results" / "stage5" / "baryon_timing_refinement" / "task_005.npz",
    1.25: ROOT / "results" / "stage5" / "baryon_timing_refinement" / "task_006.npz",
    1.5: ROOT / "results" / "stage5" / "baryon_timing_refinement" / "task_007.npz",
}
SUMMARY = ROOT / "results" / "stage5" / "feedback_frontier_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "feedback_frontier_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_feedback_frontier.png"


def limiter_history(data: np.lib.npyio.NpzFile, assembly_time_myr: float) -> str:
    times = data["times_myr"]
    bondi_limited = data["bondi_limited"]
    assembled = times >= assembly_time_myr
    eddington = np.flatnonzero(assembled & ~bondi_limited)
    if not len(eddington):
        return "bondi_throughout"
    first = int(eddington[0])
    for index in range(first + 1, len(times)):
        if np.all(bondi_limited[index:]):
            return "transient_eddington"
    return "eddington_at_end"


def load_cases() -> tuple[list[dict], dict[float, dict]]:
    ceilings = {}
    for assembly_time, path in CEILING_PATHS.items():
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        ceilings[assembly_time] = {
            "final_mass_msun": float(metadata["final_black_hole_mass_msun"]),
            "dark_matter_msun": float(metadata["dark_matter_accreted_msun"]),
            "baryon_msun": float(metadata["baryon_accreted_onto_bh_msun"]),
        }
    optimal_ceiling = ceilings[1.25]
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    for manifest_row in manifest:
        task_id = int(manifest_row["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            final_mass = float(metadata["final_black_hole_mass_msun"])
            dark_mass = float(metadata["dark_matter_accreted_msun"])
            eddington_samples = np.flatnonzero(~data["bondi_limited"])
            first_eddington_time = (
                float(data["times_myr"][eddington_samples[0]])
                if len(eddington_samples)
                else float("nan")
            )
            final_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
            history = limiter_history(data, float(manifest_row["assembly_time_myr"]))
        assembly_time = float(manifest_row["assembly_time_myr"])
        matched_ceiling = ceilings[assembly_time]
        cases.append(
            {
                "task_id": task_id,
                "configuration": manifest_row["configuration"],
                "gas_regime": manifest_row["gas_regime"],
                "gas_density_msun_pc3": float(manifest_row["gas_density_msun_pc3"]),
                "assembly_time_myr": assembly_time,
                "feedback_efficiency": float(manifest_row["feedback_efficiency"]),
                "feedback_heating_fraction": float(manifest_row["feedback_heating_fraction"]),
                "metallicity_solar": float(manifest_row["metallicity_solar"]),
                "final_black_hole_mass_msun": final_mass,
                "ceiling_mass_retention_fraction": final_mass / optimal_ceiling["final_mass_msun"],
                "matched_mass_retention_fraction": final_mass / matched_ceiling["final_mass_msun"],
                "dark_matter_accreted_msun": dark_mass,
                "dark_accretion_retention_fraction": dark_mass / optimal_ceiling["dark_matter_msun"],
                "matched_dark_retention_fraction": dark_mass / matched_ceiling["dark_matter_msun"],
                "baryon_accreted_onto_bh_msun": float(metadata["baryon_accreted_onto_bh_msun"]),
                "limiter_history": history,
                "first_eddington_time_myr": first_eddington_time,
                "final_bondi_to_eddington_ratio": final_ratio,
                "final_density_fraction": float(
                    metadata["final_ambient_gas_density_msun_pc3"]
                    / metadata["gas_density_msun_pc3"]
                ),
                "final_sound_speed_km_s": float(metadata["final_ambient_gas_sound_speed_km_s"]),
                "scale_radius_expansion_factor": float(
                    metadata["final_baryon_scale_radius_pc"]
                    / metadata["initial_baryon_scale_radius_pc"]
                ),
                "feedback_binding_energy_erg": float(
                    metadata["feedback_binding_energy_erg"]
                ),
                "expansion_feedback_to_binding_ratio": float(
                    metadata["final_feedback_to_binding_ratio"]
                ),
                "total_feedback_to_binding_ratio": float(
                    metadata["final_total_feedback_to_binding_ratio"]
                ),
                "retained_thermal_energy_fraction": float(metadata["retained_thermal_energy_fraction"]),
                "minimum_cooling_time_myr": float(metadata["minimum_physical_cooling_time_myr"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    return cases, ceilings


def select(cases: list[dict], **criteria: object) -> dict:
    matches = [
        case
        for case in cases
        if all(case[key] == value for key, value in criteria.items())
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one case for {criteria}, found {len(matches)}")
    return matches[0]


def save(cases: list[dict], ceilings: dict[float, dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    central = [
        case
        for case in cases
        if case["configuration"] == "mixed_sweep"
    ]
    mechanism_controls = [
        case
        for case in cases
        if case["configuration"] in ("pure_expansion", "pure_heating")
    ]
    metallicity_controls = [
        case for case in cases if case["configuration"] == "solar_mixed"
    ]
    statistics = {
        "case_count": len(cases),
        "eddington_ceilings_by_assembly_time": {
            f"{assembly_time:g}": ceiling
            for assembly_time, ceiling in ceilings.items()
        },
        "central_mixed_response": {
            regime: [
                case
                for case in central
                if case["gas_regime"] == regime
            ]
            for regime in ("transition", "dense")
        },
        "epsilon_1e-3": {
            regime: select(
                cases,
                configuration="mixed_sweep",
                gas_regime=regime,
                feedback_efficiency=1.0e-3,
            )
            for regime in ("transition", "dense")
        },
        "minimum_retention_case": min(
            cases, key=lambda case: case["ceiling_mass_retention_fraction"]
        ),
        "minimum_matched_retention_case": min(
            cases, key=lambda case: case["matched_mass_retention_fraction"]
        ),
        "mechanism_controls": mechanism_controls,
        "metallicity_controls": metallicity_controls,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
        "minimum_physical_cooling_time_myr": min(
            case["minimum_cooling_time_myr"] for case in cases
        ),
        "limiter_history_counts": {
            history: sum(case["limiter_history"] == history for case in cases)
            for history in ("bondi_throughout", "transient_eddington", "eddington_at_end")
        },
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii"
    )


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    colors = {"transition": "#1b9e77", "dense": "#b2182b"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for regime, color in colors.items():
        selected = sorted(
            (
                case
                for case in cases
                if case["configuration"] == "mixed_sweep"
                and case["gas_regime"] == regime
            ),
            key=lambda case: case["feedback_efficiency"],
        )
        efficiencies = [case["feedback_efficiency"] for case in selected]
        axes[0, 0].plot(
            efficiencies,
            [case["ceiling_mass_retention_fraction"] for case in selected],
            color=color, marker="o", label=regime,
        )
        axes[0, 1].plot(
            efficiencies,
            [case["dark_accretion_retention_fraction"] for case in selected],
            color=color, marker="o", label=regime,
        )
        axes[1, 0].plot(
            efficiencies,
            [case["final_bondi_to_eddington_ratio"] for case in selected],
            color=color, marker="o", label=regime,
        )
    for regime, color in colors.items():
        for efficiency, linestyle in ((0.0, ":"), (1.0e-3, "--"), (1.0e-2, "-")):
            selected = sorted(
                (
                    case
                    for case in cases
                    if case["gas_regime"] == regime
                    and case["feedback_heating_fraction"] == 0.5
                    and case["metallicity_solar"] == 0.0
                    and case["feedback_efficiency"] == efficiency
                    and case["assembly_time_myr"] in (1.0, 1.25, 1.5)
                ),
                key=lambda case: case["assembly_time_myr"],
            )
            axes[1, 1].plot(
                [case["assembly_time_myr"] for case in selected],
                [case["ceiling_mass_retention_fraction"] for case in selected],
                color=color, linestyle=linestyle, marker="o",
                label=f"{regime}, eps={efficiency:g}",
            )
    for axis in axes.flat:
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=7)
    for axis in axes[0, :]:
        axis.set_xscale("symlog", linthresh=1.0e-5)
        axis.set_xlabel("Feedback efficiency")
    axes[1, 0].set_xscale("symlog", linthresh=1.0e-5)
    axes[1, 0].set_yscale("log")
    axes[1, 0].axhline(1.0, color="black", linestyle=":")
    axes[1, 0].set_xlabel("Feedback efficiency")
    axes[1, 1].set_xlabel("Assembly time [Myr]")
    axes[0, 0].set_ylabel("Final mass / no-feedback ceiling")
    axes[0, 1].set_ylabel("Dark accretion / no-feedback ceiling")
    axes[1, 0].set_ylabel("Final Bondi / Eddington")
    axes[1, 1].set_ylabel("Final mass / no-feedback ceiling")
    axes[0, 0].set_title("Physical-cooling mass retention")
    axes[0, 1].set_title("Dark-channel retention")
    axes[1, 0].set_title("Gas limiter response")
    axes[1, 1].set_title("Assembly-time robustness")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases, ceilings = load_cases()
    save(cases, ceilings)
    plot(cases)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
