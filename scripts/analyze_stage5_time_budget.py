"""Build the stage-5 redshift and Eddington-activity feasibility budget."""

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


FORMATION_REDSHIFTS = (30.0, 25.0, 20.0, 15.0, 10.0)
OBSERVATION_REDSHIFTS = (11.0, 8.0, 6.0, 4.0)
SEED_MASSES_MSUN = (1.0, 10.0, 100.0, 1.0e3, 1.0e5)
CSV_PATH = ROOT / "results" / "stage5" / "time_budget.csv"
STATISTICS_PATH = ROOT / "results" / "stage5" / "time_budget_statistics.json"
FIGURE_PATH = ROOT / "results" / "stage5" / "figures" / "stage5_time_budget.png"


def build_rows() -> list[dict]:
    cosmology = FlatLambdaCDM()
    rows = []
    for formation_redshift in FORMATION_REDSHIFTS:
        for observation_redshift in OBSERVATION_REDSHIFTS:
            if observation_redshift >= formation_redshift:
                continue
            available = cosmology.elapsed_time_myr(
                formation_redshift,
                observation_redshift,
            )
            for seed_mass in SEED_MASSES_MSUN:
                for target_mass in LRD_TARGET_MASSES_MSUN:
                    activity = required_eddington_activity(
                        seed_mass,
                        target_mass,
                        available,
                    )
                    rows.append(
                        {
                            "formation_redshift": formation_redshift,
                            "observation_redshift": observation_redshift,
                            "available_time_myr": available,
                            "seed_mass_msun": seed_mass,
                            "target_mass_msun": target_mass,
                            "required_eddington_activity": activity,
                            "sub_eddington_feasible": activity <= 1.0,
                        }
                    )
    return rows


def save(rows: list[dict]) -> None:
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    cosmology = FlatLambdaCDM()
    statistics = {
        "row_count": len(rows),
        "cosmology": {
            "hubble_km_s_mpc": cosmology.hubble_km_s_mpc,
            "omega_matter": cosmology.omega_matter,
            "omega_lambda": cosmology.omega_lambda,
            "radiation_included": False,
        },
        "lrd_target_masses_msun": list(LRD_TARGET_MASSES_MSUN),
        "cosmic_ages_myr": {
            f"z{redshift:g}": cosmology.age_myr(redshift)
            for redshift in sorted(
                set(FORMATION_REDSHIFTS + OBSERVATION_REDSHIFTS),
                reverse=True,
            )
        },
        "sub_eddington_feasible_count": sum(
            row["sub_eddington_feasible"] for row in rows
        ),
    }
    STATISTICS_PATH.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="ascii",
    )


def plot(rows: list[dict]) -> None:
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    cosmology = FlatLambdaCDM()
    redshifts = np.linspace(4.0, 30.0, 300)
    ages = np.array([cosmology.age_myr(redshift) for redshift in redshifts])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axes[0, 0].plot(redshifts, ages, color="#222222")
    axes[0, 0].invert_xaxis()
    axes[0, 0].set_xlabel("Redshift")
    axes[0, 0].set_ylabel("Cosmic age [Myr]")
    axes[0, 0].set_title("Planck matter + Lambda time budget")

    panels = (
        (axes[0, 1], 1.0e5),
        (axes[1, 0], 1.0e6),
        (axes[1, 1], 1.0e7),
    )
    colors = {1.0: "#2166ac", 100.0: "#1b9e77", 1.0e5: "#b2182b"}
    for axis, target_mass in panels:
        for seed_mass in (1.0, 100.0, 1.0e5):
            selected = [
                row
                for row in rows
                if row["formation_redshift"] == 20.0
                and row["seed_mass_msun"] == seed_mass
                and row["target_mass_msun"] == target_mass
            ]
            selected.sort(key=lambda row: row["observation_redshift"])
            axis.plot(
                [row["observation_redshift"] for row in selected],
                [row["required_eddington_activity"] for row in selected],
                marker="o",
                color=colors[seed_mass],
                label=f"seed={seed_mass:g} M_sun",
            )
        axis.axhline(1.0, color="black", linestyle=":")
        axis.set_yscale("log")
        axis.set_xlabel("Observation redshift")
        axis.set_ylabel("Required f_Edd x duty")
        axis.set_title(f"z_form=20 to {target_mass:.0e} M_sun")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE_PATH, dpi=180)
    plt.close(fig)


def main() -> None:
    rows = build_rows()
    save(rows)
    plot(rows)
    print(STATISTICS_PATH.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
