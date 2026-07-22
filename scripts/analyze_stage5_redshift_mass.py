"""Analyze the first stage-5 cosmological redshift-mass anchor matrix."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
MANIFEST = ROOT / "hpc" / "stage5_redshift_mass.tsv"
RESULT_SETS = {
    "scaled_rmin": ROOT / "results" / "stage5" / "redshift_mass",
    "fixed_rmin": ROOT / "results" / "stage5" / "redshift_mass_fixed_rmin",
}
SUMMARY = ROOT / "results" / "stage5" / "redshift_mass_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "redshift_mass_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_redshift_mass.png"


def load_cases() -> list[dict]:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    cases = []
    for boundary_protocol, results in RESULT_SETS.items():
        for row in manifest:
            task_id = int(row["task_id"])
            path = results / f"task_{task_id:03d}.npz"
            if not path.exists():
                raise FileNotFoundError(path)
            with np.load(path, allow_pickle=False) as data:
                metadata = json.loads(str(data["metadata_json"]))
            dark_mass = float(metadata["dark_matter_accreted_msun"])
            influence_radius_pc = (
                float(row["black_hole_seed_msun"])
                / float(row["halo_mass_msun"])
                * float(metadata["halo_virial_radius_pc"])
            )
            cases.append(
                {
                "task_id": task_id,
                "boundary_protocol": boundary_protocol,
                "halo_mass_msun": float(row["halo_mass_msun"]),
                "halo_redshift": float(row["halo_redshift"]),
                "halo_concentration": float(row["halo_concentration"]),
                "halo_scale_radius_pc": float(metadata["halo_scale_radius_pc"]),
                "halo_scale_density_msun_pc3": float(
                    metadata["halo_scale_density_msun_pc3"]
                ),
                "final_black_hole_mass_msun": float(
                    metadata["final_black_hole_mass_msun"]
                ),
                "dark_matter_accreted_msun": dark_mass,
                "baryon_accreted_onto_bh_msun": float(
                    metadata["baryon_accreted_onto_bh_msun"]
                ),
                "dark_fraction_of_growth": float(
                    metadata["dark_fraction_of_black_hole_growth"]
                ),
                "dark_accreted_halo_fraction": dark_mass
                / float(row["halo_mass_msun"]),
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
                "r_min_pc": float(metadata["r_min_pc"]),
                "r_min_over_rs": float(metadata["r_min_over_rs"]),
                "r_min_over_influence_radius": float(metadata["r_min_pc"])
                / influence_radius_pc,
                "cells": int(metadata["cells"]),
                "steps": int(metadata["steps"]),
                }
            )
    return cases


def grid(
    cases: list[dict],
    field: str,
    boundary_protocol: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cases = [
        case for case in cases if case["boundary_protocol"] == boundary_protocol
    ]
    masses = np.array(sorted({case["halo_mass_msun"] for case in cases}))
    redshifts = np.array(sorted({case["halo_redshift"] for case in cases}))
    values = np.empty((len(masses), len(redshifts)))
    for i, mass in enumerate(masses):
        for j, redshift in enumerate(redshifts):
            selected = [
                case
                for case in cases
                if case["halo_mass_msun"] == mass
                and case["halo_redshift"] == redshift
            ]
            if len(selected) != 1:
                raise RuntimeError("stage-5 anchor grid is incomplete")
            values[i, j] = selected[0][field]
    return masses, redshifts, values


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        ("scaled_rmin", "dark_matter_accreted_msun", "Scaled r_min: log10 DM accreted"),
        ("fixed_rmin", "dark_matter_accreted_msun", "Fixed r_min: log10 DM accreted"),
        ("boundary_ratio", "dark_matter_accreted_msun", "log10 scaled / fixed DM"),
        ("fixed_rmin", "dark_fraction_of_growth", "Fixed r_min: dark growth fraction"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for axis, (protocol, field, title) in zip(axes.flat, panels):
        if protocol == "boundary_ratio":
            masses, redshifts, scaled = grid(cases, field, "scaled_rmin")
            _, _, fixed = grid(cases, field, "fixed_rmin")
            shown = np.log10(scaled / fixed)
        else:
            masses, redshifts, values = grid(cases, field, protocol)
            shown = (
                values
                if field == "dark_fraction_of_growth"
                else np.log10(values)
            )
        image = axis.imshow(shown, origin="lower", aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(redshifts)), [f"{value:g}" for value in redshifts])
        axis.set_yticks(
            range(len(masses)),
            [f"{value:.0e}" for value in masses],
        )
        axis.set_xlabel("Halo redshift")
        axis.set_ylabel("M200c [M_sun]")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.85)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def save(cases: list[dict]) -> None:
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    best = {
        protocol: max(
            (
                case
                for case in cases
                if case["boundary_protocol"] == protocol
            ),
            key=lambda case: case["final_black_hole_mass_msun"],
        )
        for protocol in RESULT_SETS
    }
    concentration = 3.9200662204075956
    shape = math.log1p(concentration) - concentration / (1.0 + concentration)
    characteristic_over_critical = (
        200.0 * concentration**3 / (3.0 * shape)
    )
    from sidm_bh.cosmology import FlatLambdaCDM

    cosmology = FlatLambdaCDM()
    target_critical_density = 3.7 / characteristic_over_critical
    expansion_squared = (
        target_critical_density / cosmology.critical_density_msun_pc3(0.0)
    )
    equivalent_anchor_redshift = (
        (expansion_squared - cosmology.omega_lambda) / cosmology.omega_matter
    ) ** (1.0 / 3.0) - 1.0
    scaled = {
        (case["halo_mass_msun"], case["halo_redshift"]): case
        for case in cases
        if case["boundary_protocol"] == "scaled_rmin"
    }
    fixed = {
        (case["halo_mass_msun"], case["halo_redshift"]): case
        for case in cases
        if case["boundary_protocol"] == "fixed_rmin"
    }
    boundary_ratios = {
        key: scaled[key]["dark_matter_accreted_msun"]
        / fixed[key]["dark_matter_accreted_msun"]
        for key in scaled
    }

    def slopes(protocol: str, independent: str) -> dict[str, float]:
        selected = [
            case for case in cases if case["boundary_protocol"] == protocol
        ]
        if independent == "mass":
            groups = sorted({case["halo_redshift"] for case in selected})
            return {
                f"z{group:g}": float(
                    np.polyfit(
                        np.log10(
                            [
                                case["halo_mass_msun"]
                                for case in selected
                                if case["halo_redshift"] == group
                            ]
                        ),
                        np.log10(
                            [
                                case["dark_matter_accreted_msun"]
                                for case in selected
                                if case["halo_redshift"] == group
                            ]
                        ),
                        1,
                    )[0]
                )
                for group in groups
            }
        groups = sorted({case["halo_mass_msun"] for case in selected})
        return {
            f"M{group:.0e}": float(
                np.polyfit(
                    np.log10(
                        [
                            1.0 + case["halo_redshift"]
                            for case in selected
                            if case["halo_mass_msun"] == group
                        ]
                    ),
                    np.log10(
                        [
                            case["dark_matter_accreted_msun"]
                            for case in selected
                            if case["halo_mass_msun"] == group
                        ]
                    ),
                    1,
                )[0]
            )
            for group in groups
        }

    statistics = {
        "case_count": len(cases),
        "best_case": best,
        "dark_accretion_power_law_slopes": {
            protocol: {
                "with_halo_mass_at_fixed_redshift": slopes(protocol, "mass"),
                "with_one_plus_redshift_at_fixed_mass": slopes(protocol, "redshift"),
            }
            for protocol in RESULT_SETS
        },
        "scaled_to_fixed_rmin_dark_accretion_ratio": {
            "minimum": min(boundary_ratios.values()),
            "maximum": max(boundary_ratios.values()),
        },
        "resolved_case_count_rmin_le_influence_radius": {
            protocol: sum(
                case["boundary_protocol"] == protocol
                and case["r_min_over_influence_radius"] <= 1.0
                for case in cases
            )
            for protocol in RESULT_SETS
        },
        "legacy_anchor_equivalent_m200c_redshift": equivalent_anchor_redshift,
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )


def main() -> None:
    cases = load_cases()
    save(cases)
    plot(cases)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
