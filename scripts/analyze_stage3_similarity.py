"""Compare fixed physical scales with two self-similar companions."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "hpc" / "stage3_similarity_matrix.tsv"
RESULTS = ROOT / "results" / "stage3" / "similarity"
SCALING_SUMMARY = ROOT / "results" / "stage3" / "scaling_summary.csv"
SUMMARY = ROOT / "results" / "stage3" / "similarity_summary.csv"
STATISTICS = ROOT / "results" / "stage3" / "similarity_statistics.json"
FIGURE = ROOT / "results" / "stage3" / "figures" / "stage3_similarity.png"


def slope(rows: list[dict], field: str) -> float:
    ordered = sorted(rows, key=lambda row: row["halo_mass_msun"])
    return float(
        np.polyfit(
            np.log10([row["halo_mass_msun"] for row in ordered]),
            np.log10([row[field] for row in ordered]),
            1,
        )[0]
    )


def load_existing_protocol() -> list[dict]:
    with SCALING_SUMMARY.open(newline="", encoding="utf-8") as stream:
        source = list(csv.DictReader(stream))
    selected = []
    for row in source:
        is_control = row["case_type"] == "control"
        is_baryon = (
            float(row["baryon_fraction"]) == 0.05
            and float(row["assembly_multiplier"]) == 1.25
        )
        if not (is_control or is_baryon):
            continue
        selected.append(
            {
                "protocol": "fixed_physical_seed_boundary_sigma",
                "case_type": "control" if is_control else "baryon",
                "halo_mass_msun": float(row["halo_mass_msun"]),
                "black_hole_mass_msun": 100.0,
                "sigma_over_m_cm2_g": 50.0,
                "r_min_pc": 0.005,
                "assembly_time_myr": float(row["assembly_time_myr"]),
                "accreted_dark_matter_msun": float(
                    row["accreted_dark_matter_msun"]
                ),
                "mass_budget_residual_code": float(
                    row["mass_budget_residual_code"]
                ),
            }
        )
    if len(selected) != 6:
        raise RuntimeError(f"expected 6 existing cases, found {len(selected)}")
    return selected


def load_companions() -> list[dict]:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    rows = []
    for source in manifest:
        task_id = int(source["task_id"])
        path = RESULTS / f"task_{task_id:03d}.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            growth = float(
                data["black_hole_mass_msun"][-1]
                - float(metadata["black_hole_seed_msun"])
            )
        rows.append(
            {
                "protocol": source["protocol"],
                "case_type": source["case_type"],
                "halo_mass_msun": float(source["halo_mass_msun"]),
                "black_hole_mass_msun": float(source["black_hole_mass_msun"]),
                "sigma_over_m_cm2_g": float(
                    source["sigma_over_m_cm2_g"]
                ),
                "r_min_pc": float(source["r_min_pc"]),
                "assembly_time_myr": float(source["assembly_time_myr"]),
                "accreted_dark_matter_msun": growth,
                "mass_budget_residual_code": float(
                    metadata["mass_budget_residual_code"]
                ),
            }
        )
    return rows


def main() -> None:
    rows = load_existing_protocol() + load_companions()
    protocols = sorted({row["protocol"] for row in rows})
    for protocol in protocols:
        controls = {
            row["halo_mass_msun"]: row["accreted_dark_matter_msun"]
            for row in rows
            if row["protocol"] == protocol and row["case_type"] == "control"
        }
        for row in rows:
            if row["protocol"] != protocol:
                continue
            row["accreted_fraction_of_halo"] = (
                row["accreted_dark_matter_msun"] / row["halo_mass_msun"]
            )
            row["enhancement_over_matched_control"] = (
                row["accreted_dark_matter_msun"]
                / controls[row["halo_mass_msun"]]
            )

    with SUMMARY.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    statistics = {
        "all_cases_present": len(rows) == 18,
        "maximum_mass_budget_residual_code": max(
            abs(row["mass_budget_residual_code"]) for row in rows
        ),
        "growth_mass_slopes": {},
        "fully_self_similar_fractional_spread": {},
    }
    for protocol in protocols:
        statistics["growth_mass_slopes"][protocol] = {}
        for case_type in ("control", "baryon"):
            selected = [
                row
                for row in rows
                if row["protocol"] == protocol
                and row["case_type"] == case_type
            ]
            statistics["growth_mass_slopes"][protocol][case_type] = slope(
                selected,
                "accreted_dark_matter_msun",
            )
            if protocol == "fully_dimensionless_self_similar":
                fractions = np.array(
                    [row["accreted_fraction_of_halo"] for row in selected]
                )
                statistics["fully_self_similar_fractional_spread"][case_type] = (
                    float((np.max(fractions) - np.min(fractions)) / np.mean(fractions))
                )
    STATISTICS.write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    labels = {
        "fixed_physical_seed_boundary_sigma": "Fixed physical scales",
        "scaled_seed_boundary_fixed_sigma": "Scaled seed + boundary",
        "fully_dimensionless_self_similar": "Fully dimensionless similar",
    }
    colors = {
        "fixed_physical_seed_boundary_sigma": "#d62828",
        "scaled_seed_boundary_fixed_sigma": "#3b4cc0",
        "fully_dimensionless_self_similar": "#2a9d8f",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    for protocol in protocols:
        baryon = sorted(
            [
                row
                for row in rows
                if row["protocol"] == protocol and row["case_type"] == "baryon"
            ],
            key=lambda row: row["halo_mass_msun"],
        )
        axes[0].plot(
            [row["halo_mass_msun"] for row in baryon],
            [row["accreted_fraction_of_halo"] for row in baryon],
            marker="o",
            color=colors[protocol],
            label=labels[protocol],
        )
        axes[1].plot(
            [row["halo_mass_msun"] for row in baryon],
            [row["enhancement_over_matched_control"] for row in baryon],
            marker="o",
            color=colors[protocol],
            label=labels[protocol],
        )
    for axis in axes:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Halo mass [M_sun]")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Accreted dark mass / halo mass")
    axes[0].set_title("Breaking and restoring self-similarity")
    axes[1].set_ylabel("Enhancement over matched control")
    axes[1].set_title("Baryon amplification")
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(SUMMARY)
    print(STATISTICS)
    print(FIGURE)
    print(json.dumps(statistics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
