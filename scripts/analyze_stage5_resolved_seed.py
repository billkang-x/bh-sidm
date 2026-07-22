"""Analyze the stage-5 resolved seed-to-halo-ratio anchor matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage5_resolved_seed.tsv"
RESULTS = ROOT / "results" / "stage5" / "resolved_seed"
SUMMARY = ROOT / "results" / "stage5" / "resolved_seed_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "resolved_seed_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_resolved_seed.png"


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
        halo_mass = float(row["halo_mass_msun"])
        final_mass = float(metadata["final_black_hole_mass_msun"])
        dark_mass = float(metadata["dark_matter_accreted_msun"])
        cases.append(
            {
                "task_id": task_id,
                "halo_mass_msun": halo_mass,
                "halo_redshift": float(row["halo_redshift"]),
                "black_hole_seed_msun": seed,
                "seed_to_halo_ratio": seed / halo_mass,
                "halo_scale_density_msun_pc3": float(
                    metadata["halo_scale_density_msun_pc3"]
                ),
                "final_black_hole_mass_msun": final_mass,
                "black_hole_growth_factor": final_mass / seed,
                "dark_matter_accreted_msun": dark_mass,
                "dark_accreted_halo_fraction": dark_mass / halo_mass,
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


def grid(cases: list[dict], field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
                raise RuntimeError("resolved-seed grid is incomplete")
            values[i, j] = selected[0][field]
    return masses, redshifts, values


def plot(cases: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        ("final_black_hole_mass_msun", "log10 final BH mass", True),
        ("black_hole_growth_factor", "log10 2 Myr growth factor", True),
        ("dark_fraction_of_growth", "Dark fraction of growth", False),
        ("dark_accreted_halo_fraction", "log10 DM accreted / halo", True),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for axis, (field, title, logarithmic) in zip(axes.flat, panels):
        masses, redshifts, values = grid(cases, field)
        shown = np.log10(values) if logarithmic else values
        image = axis.imshow(shown, origin="lower", aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(redshifts)), [f"{value:g}" for value in redshifts])
        axis.set_yticks(range(len(masses)), [f"{value:.0e}" for value in masses])
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
    statistics = {
        "case_count": len(cases),
        "best_final_mass_case": max(
            cases, key=lambda case: case["final_black_hole_mass_msun"]
        ),
        "largest_growth_factor_case": max(
            cases, key=lambda case: case["black_hole_growth_factor"]
        ),
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


def main() -> None:
    cases = load_cases()
    save(cases)
    plot(cases)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
