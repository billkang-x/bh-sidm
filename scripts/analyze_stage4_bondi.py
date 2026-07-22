"""Analyze Bondi-to-Eddington transitions in the stage-4 gas matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage4_bondi.tsv"
RESULTS = ROOT / "results" / "stage4" / "bondi"
SUMMARY = ROOT / "results" / "stage4" / "bondi_summary.csv"
STATISTICS = ROOT / "results" / "stage4" / "bondi_statistics.json"
FIGURE = ROOT / "results" / "stage4" / "figures" / "stage4_bondi.png"


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
            initial_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][0]
                / data["eddington_baryon_growth_limit_msun_myr"][0]
            )
            final_ratio = float(
                data["bondi_baryon_growth_limit_msun_myr"][-1]
                / data["eddington_baryon_growth_limit_msun_myr"][-1]
            )
            transition_time = float(
                metadata["bondi_to_eddington_transition_myr"]
            )
            transition_mass = float("nan")
            if np.isfinite(transition_time):
                transition_index = int(
                    np.argmin(np.abs(data["times_myr"] - transition_time))
                )
                transition_mass = float(
                    data["black_hole_mass_msun"][transition_index]
                )
            critical_mass = 100.0 / initial_ratio
        cases.append(
            {
                "task_id": task_id,
                "gas_density_msun_pc3": float(row["gas_density_msun_pc3"]),
                "gas_sound_speed_km_s": float(row["gas_sound_speed_km_s"]),
                "initial_bondi_to_eddington_ratio": initial_ratio,
                "final_bondi_to_eddington_ratio": final_ratio,
                "analytic_critical_black_hole_mass_msun": critical_mass,
                "bondi_limited_sample_fraction": float(
                    metadata["post_assembly_bondi_limited_sample_fraction"]
                ),
                "transition_time_myr": transition_time,
                "transition_black_hole_mass_msun": transition_mass,
                "transition_mass_relative_error": (
                    (transition_mass - critical_mass) / critical_mass
                    if np.isfinite(transition_mass)
                    else float("nan")
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
                "dark_fraction_of_black_hole_growth": float(
                    metadata["dark_fraction_of_black_hole_growth"]
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


def matrix(cases: list[dict], field: str) -> np.ndarray:
    densities = sorted({case["gas_density_msun_pc3"] for case in cases})
    sounds = sorted({case["gas_sound_speed_km_s"] for case in cases})
    return np.array(
        [
            [
                next(
                    case[field]
                    for case in cases
                    if case["gas_density_msun_pc3"] == density
                    and case["gas_sound_speed_km_s"] == sound
                )
                for density in densities
            ]
            for sound in sounds
        ]
    )


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    densities = sorted({case["gas_density_msun_pc3"] for case in cases})
    sounds = sorted({case["gas_sound_speed_km_s"] for case in cases})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)
    panels = [
        ("initial_bondi_to_eddington_ratio", "Initial Bondi / Eddington", "magma", True),
        ("bondi_limited_sample_fraction", "Post-assembly Bondi-limited fraction", "viridis", False),
        ("baryon_accreted_onto_bh_msun", "Retained baryon mass [M_sun]", "plasma", True),
    ]
    for axis, (field, title, color_map, logarithmic) in zip(axes.flat[:3], panels):
        values = matrix(cases, field)
        shown = np.log10(values) if logarithmic else values
        image = axis.imshow(shown, origin="lower", cmap=color_map, aspect="auto")
        axis.set_xticks(range(len(densities)), [f"{value:g}" for value in densities])
        axis.set_yticks(range(len(sounds)), [f"{value:g}" for value in sounds])
        axis.set_xlabel("Gas density [M_sun/pc^3]")
        axis.set_ylabel("Sound speed [km/s]")
        axis.set_title(title)
        for i in range(len(sounds)):
            for j in range(len(densities)):
                axis.text(j, i, f"{values[i, j]:.3g}", ha="center", va="center", color="white")
        fig.colorbar(image, ax=axis, label="log10" if logarithmic else None)

    representative = sorted(
        cases,
        key=lambda case: case["initial_bondi_to_eddington_ratio"],
    )
    selected = [representative[0], representative[len(representative) // 2], representative[-1]]
    for case in selected:
        with np.load(case["path"], allow_pickle=False) as data:
            axes[1, 1].plot(
                data["times_myr"],
                data["bondi_baryon_growth_limit_msun_myr"]
                / data["eddington_baryon_growth_limit_msun_myr"],
                label=(
                    f"rho={case['gas_density_msun_pc3']:g}, "
                    f"c_s={case['gas_sound_speed_km_s']:g}"
                ),
            )
    axes[1, 1].axhline(1.0, color="black", linestyle=":")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Time [Myr]")
    axes[1, 1].set_ylabel("Bondi / Eddington limit")
    axes[1, 1].set_title("SIDM-driven limiter evolution")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save_summary(cases)
    plot(cases)
    transition_cases = [
        case
        for case in cases
        if case["initial_bondi_to_eddington_ratio"] < 1.0
        and np.isfinite(case["transition_time_myr"])
    ]
    statistics = {
        "all_cases_present": len(cases) == 9,
        "maximum_mass_budget_residual_code": max(
            abs(case["mass_budget_residual_code"]) for case in cases
        ),
        "initially_eddington_limited_cases": sum(
            case["initial_bondi_to_eddington_ratio"] >= 1.0 for case in cases
        ),
        "sidm_driven_transition_cases": len(transition_cases),
        "transition_time_range_myr": (
            [
                min(case["transition_time_myr"] for case in transition_cases),
                max(case["transition_time_myr"] for case in transition_cases),
            ]
            if transition_cases
            else None
        ),
        "maximum_absolute_transition_mass_relative_error": (
            max(
                abs(case["transition_mass_relative_error"])
                for case in transition_cases
            )
            if transition_cases
            else None
        ),
        "always_bondi_limited_cases": sum(
            case["final_bondi_to_eddington_ratio"] < 1.0 for case in cases
        ),
        "baryon_accreted_range_msun": [
            min(case["baryon_accreted_onto_bh_msun"] for case in cases),
            max(case["baryon_accreted_onto_bh_msun"] for case in cases),
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
