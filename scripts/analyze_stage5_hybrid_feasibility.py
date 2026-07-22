"""Map a simulated 2 Myr dark burst onto conservative LRD time budgets."""

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

from sidm_bh.cosmology import FlatLambdaCDM
from sidm_bh.stage5 import LRD_TARGET_MASSES_MSUN, required_eddington_activity


ANCHORS = ROOT / "results" / "stage5" / "resolved_seed_summary.csv"
OUTPUT = ROOT / "results" / "stage5" / "hybrid_feasibility.csv"
STATISTICS = ROOT / "results" / "stage5" / "hybrid_feasibility_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_hybrid_feasibility.png"
OBSERVATION_REDSHIFTS = (11.0, 8.0, 6.0, 4.0)
DARK_BURST_DURATION_MYR = 2.0


def load_rows() -> list[dict]:
    with ANCHORS.open(newline="", encoding="ascii") as stream:
        anchors = list(csv.DictReader(stream))
    cosmology = FlatLambdaCDM()
    rows = []
    for anchor in anchors:
        formation_redshift = float(anchor["halo_redshift"])
        formation_age = cosmology.age_myr(formation_redshift)
        post_burst_age = formation_age + DARK_BURST_DURATION_MYR
        seed = float(anchor["black_hole_seed_msun"])
        post_burst_mass = float(anchor["final_black_hole_mass_msun"])
        for observation_redshift in OBSERVATION_REDSHIFTS:
            if observation_redshift >= formation_redshift:
                continue
            observation_age = cosmology.age_myr(observation_redshift)
            if observation_age <= post_burst_age:
                continue
            total_window = observation_age - formation_age
            post_burst_window = observation_age - post_burst_age
            for target in LRD_TARGET_MASSES_MSUN:
                seed_only_activity = required_eddington_activity(
                    seed,
                    target,
                    total_window,
                )
                post_burst_activity = required_eddington_activity(
                    post_burst_mass,
                    target,
                    post_burst_window,
                )
                rows.append(
                    {
                        "halo_mass_msun": float(anchor["halo_mass_msun"]),
                        "formation_redshift": formation_redshift,
                        "observation_redshift": observation_redshift,
                        "seed_mass_msun": seed,
                        "post_dark_burst_mass_msun": post_burst_mass,
                        "target_mass_msun": target,
                        "available_post_burst_time_myr": post_burst_window,
                        "seed_only_required_activity": seed_only_activity,
                        "post_burst_required_activity": post_burst_activity,
                        "activity_reduction_factor": (
                            post_burst_activity / seed_only_activity
                            if seed_only_activity > 0.0
                            else 0.0
                        ),
                        "seed_only_feasible": seed_only_activity <= 1.0,
                        "post_burst_feasible": post_burst_activity <= 1.0,
                        "dark_burst_changes_feasibility": (
                            seed_only_activity > 1.0
                            and post_burst_activity <= 1.0
                        ),
                    }
                )
    return rows


def save(rows: list[dict]) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    statistics = {
        "row_count": len(rows),
        "dark_burst_duration_myr": DARK_BURST_DURATION_MYR,
        "seed_only_feasible_count": sum(row["seed_only_feasible"] for row in rows),
        "post_burst_feasible_count": sum(row["post_burst_feasible"] for row in rows),
        "dark_burst_changes_feasibility_count": sum(
            row["dark_burst_changes_feasibility"] for row in rows
        ),
        "activity_reduction_factor_range": [
            min(row["activity_reduction_factor"] for row in rows),
            max(row["activity_reduction_factor"] for row in rows),
        ],
        "changed_cases": [
            row for row in rows if row["dark_burst_changes_feasibility"]
        ],
    }
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )


def grid(
    rows: list[dict],
    observation_redshift: float,
    target_mass: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    masses = np.array(sorted({row["halo_mass_msun"] for row in rows}))
    formation_redshifts = np.array(
        sorted({row["formation_redshift"] for row in rows})
    )
    values = np.full((len(masses), len(formation_redshifts)), np.nan)
    for i, mass in enumerate(masses):
        for j, formation_redshift in enumerate(formation_redshifts):
            selected = [
                row
                for row in rows
                if row["halo_mass_msun"] == mass
                and row["formation_redshift"] == formation_redshift
                and row["observation_redshift"] == observation_redshift
                and row["target_mass_msun"] == target_mass
            ]
            if selected:
                values[i, j] = selected[0]["post_burst_required_activity"]
    return masses, formation_redshifts, values


def plot(rows: list[dict]) -> None:
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    panels = (
        (11.0, 1.0e6),
        (11.0, 1.0e7),
        (8.0, 1.0e6),
        (8.0, 1.0e7),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for axis, (observation_redshift, target_mass) in zip(axes.flat, panels):
        masses, formation_redshifts, values = grid(
            rows,
            observation_redshift,
            target_mass,
        )
        image = axis.imshow(
            np.ma.masked_invalid(values),
            origin="lower",
            aspect="auto",
            cmap="RdYlGn_r",
            vmin=0.0,
            vmax=max(1.5, float(np.nanmax(values))),
        )
        axis.set_xticks(
            range(len(formation_redshifts)),
            [f"{value:g}" for value in formation_redshifts],
        )
        axis.set_yticks(range(len(masses)), [f"{value:.0e}" for value in masses])
        axis.set_xlabel("Formation redshift")
        axis.set_ylabel("M200c [M_sun]")
        axis.set_title(
            f"to {target_mass:.0e} M_sun by z={observation_redshift:g}"
        )
        fig.colorbar(image, ax=axis, shrink=0.85, label="Required f_Edd x duty")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    rows = load_rows()
    save(rows)
    plot(rows)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
