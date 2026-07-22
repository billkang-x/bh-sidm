"""Analyze the resolved stage-5 SIDM transport screening matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
MANIFEST = ROOT / "hpc" / "stage5_transport.tsv"
RESULTS = ROOT / "results" / "stage5" / "transport"
SUMMARY = ROOT / "results" / "stage5" / "transport_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "transport_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_transport.png"

from sidm_bh.constants import M_SUN_CGS, PC_CGS
from sidm_bh.sidm import conductivity_cgs, knudsen_number


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
        seed = float(row["black_hole_seed_msun"])
        final_mass = float(metadata["final_black_hole_mass_msun"])
        cases.append(
            {
                "task_id": task_id,
                "halo_mass_msun": float(row["halo_mass_msun"]),
                "halo_redshift": float(row["halo_redshift"]),
                "halo_concentration": float(row["halo_concentration"]),
                "sigma_over_m_cm2_g": float(row["sigma_over_m_cm2_g"]),
                "black_hole_seed_msun": seed,
                "final_black_hole_mass_msun": final_mass,
                "black_hole_growth_factor": final_mass / seed,
                "dark_matter_accreted_msun": float(
                    metadata["dark_matter_accreted_msun"]
                ),
                "baryon_accreted_onto_bh_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "dark_fraction_of_growth": float(
                    metadata["dark_fraction_of_black_hole_growth"]
                ),
                "r_min_over_influence_radius": float(
                    metadata["r_min_over_black_hole_influence_radius"]
                ),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "steps": int(metadata["steps"]),
            }
        )
    return cases


def save(cases: list[dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    best_by_halo_redshift = {}
    for halo_mass in sorted({case["halo_mass_msun"] for case in cases}):
        for redshift in sorted({case["halo_redshift"] for case in cases}):
            selected = [
                case
                for case in cases
                if case["halo_mass_msun"] == halo_mass
                and case["halo_redshift"] == redshift
            ]
            best_by_halo_redshift[
                f"M{halo_mass:.0e}_z{redshift:g}"
            ] = max(selected, key=lambda case: case["black_hole_growth_factor"])

    best_sigma_counts = {
        f"{sigma:g}": sum(
            case["sigma_over_m_cm2_g"] == sigma
            for case in best_by_halo_redshift.values()
        )
        for sigma in sorted({case["sigma_over_m_cm2_g"] for case in cases})
    }
    best_concentration_counts = {
        f"{concentration:g}": sum(
            case["halo_concentration"] == concentration
            for case in best_by_halo_redshift.values()
        )
        for concentration in sorted(
            {case["halo_concentration"] for case in cases}
        )
    }
    transport_regime = {}
    representative = [
        case
        for case in cases
        if case["halo_mass_msun"] == 1.0e9
        and case["halo_redshift"] == 30.0
        and case["halo_concentration"] == 8.0
    ]
    for case in representative:
        path = RESULTS / f"task_{case['task_id']:03d}.npz"
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            radius_index = int(
                np.argmin(
                    np.abs(
                        data["radii_pc"]
                        - metadata["initial_baryon_scale_radius_pc"]
                    )
                )
            )
            assembly_index = int(
                np.argmin(
                    np.abs(data["times_myr"] - metadata["assembly_time_myr"])
                )
            )
            states = {}
            for label, time_index in (("initial", 0), ("assembled", assembly_index)):
                density_cgs = (
                    data["density_msun_pc3"][time_index, radius_index]
                    * M_SUN_CGS
                    / PC_CGS**3
                )
                dispersion_cms = (
                    data["velocity_dispersion_km_s"][time_index, radius_index]
                    * 1.0e5
                )
                states[label] = {
                    "time_myr": float(data["times_myr"][time_index]),
                    "radius_pc": float(data["radii_pc"][radius_index]),
                    "knudsen_number": float(
                        knudsen_number(
                            density_cgs,
                            dispersion_cms,
                            case["sigma_over_m_cm2_g"],
                        )
                    ),
                    "conductivity_cgs": float(
                        conductivity_cgs(
                            density_cgs,
                            dispersion_cms,
                            case["sigma_over_m_cm2_g"],
                        )
                    ),
                }
        transport_regime[f"sigma_{case['sigma_over_m_cm2_g']:g}"] = states
    statistics = {
        "case_count": len(cases),
        "best_global_case": max(
            cases, key=lambda case: case["black_hole_growth_factor"]
        ),
        "best_by_halo_mass_and_redshift": best_by_halo_redshift,
        "best_sigma_counts": best_sigma_counts,
        "best_concentration_counts": best_concentration_counts,
        "representative_transport_regime_M1e9_z30_c8": transport_regime,
        "target_reached_counts": {
            f"{target:.0e}": sum(
                case["final_black_hole_mass_msun"] >= target
                for case in cases
            )
            for target in (1.0e5, 1.0e6, 1.0e7)
        },
        "r_min_over_influence_radius_range": [
            min(case["r_min_over_influence_radius"] for case in cases),
            max(case["r_min_over_influence_radius"] for case in cases),
        ],
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
    masses = sorted({case["halo_mass_msun"] for case in cases})
    redshifts = sorted({case["halo_redshift"] for case in cases})
    colors = {3.0: "#2166ac", 5.0: "#1b9e77", 8.0: "#b2182b"}
    fig, axes = plt.subplots(
        len(masses),
        len(redshifts),
        figsize=(12, 12),
        sharex=True,
        constrained_layout=True,
    )
    for row_index, halo_mass in enumerate(masses):
        for column_index, redshift in enumerate(redshifts):
            axis = axes[row_index, column_index]
            for concentration in sorted(
                {case["halo_concentration"] for case in cases}
            ):
                selected = sorted(
                    (
                        case
                        for case in cases
                        if case["halo_mass_msun"] == halo_mass
                        and case["halo_redshift"] == redshift
                        and case["halo_concentration"] == concentration
                    ),
                    key=lambda case: case["sigma_over_m_cm2_g"],
                )
                axis.plot(
                    [case["sigma_over_m_cm2_g"] for case in selected],
                    [case["black_hole_growth_factor"] for case in selected],
                    color=colors[concentration],
                    marker="o",
                    label=f"c={concentration:g}",
                )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.grid(alpha=0.25)
            axis.set_title(f"M={halo_mass:.0e}, z={redshift:g}")
            if row_index == len(masses) - 1:
                axis.set_xlabel("sigma/m [cm2/g]")
            if column_index == 0:
                axis.set_ylabel("2 Myr growth factor")
            if row_index == 0 and column_index == 0:
                axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    cases = load_cases()
    save(cases)
    plot(cases)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
