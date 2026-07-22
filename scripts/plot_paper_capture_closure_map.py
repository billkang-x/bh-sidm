"""Plot the sparse Bondi-capture closure map from existing stage-5 summaries."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage5"
OUTPUT = ROOT / "paper" / "figures" / "capture_closure_map.png"


def read_rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def finite(value: str) -> bool:
    return math.isfinite(float(value))


def classify_boundary_group(rows: list[dict[str, str]]) -> str:
    crossings = [finite(row["time_to_1e7_msun_myr"]) for row in rows]
    if all(crossings):
        return "both_boundaries"
    if any(crossings):
        return "boundary_ambiguous"
    return "no_crossing"


def main() -> None:
    capture = read_rows("capture_sensitivity_summary.csv")
    threshold = read_rows("seed_threshold_summary.csv")
    threshold_final = read_rows("seed_threshold_final_summary.csv")
    light = read_rows("light_seed_boundary_summary.csv")

    classes: list[tuple[float, float, str]] = []
    for row in capture:
        status = "one_boundary" if finite(row["time_to_1e7_msun_myr"]) else "no_crossing"
        classes.append(
            (
                float(row["black_hole_seed_msun"]),
                float(row["dark_bondi_lambda"]),
                status,
            )
        )

    lambda_quarter: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in light:
        if row["axis"] == "inner_boundary":
            lambda_quarter[float(row["black_hole_seed_msun"])].append(row)
    for row in threshold:
        lambda_quarter[float(row["black_hole_seed_msun"])].append(row)
    for row in threshold_final:
        if row["axis"] == "inner_boundary":
            lambda_quarter[float(row["black_hole_seed_msun"])].append(row)
    for seed, rows in sorted(lambda_quarter.items()):
        classes.append((seed, 0.25, classify_boundary_group(rows)))

    styles = {
        "no_crossing": ("o", "#777777", "No crossing at tested boundary"),
        "boundary_ambiguous": ("D", "#e68613", "Boundary ambiguous"),
        "both_boundaries": ("s", "#2a9d5b", "Crossing at both boundaries"),
        "one_boundary": ("^", "#2878b5", "Crossing at one tested boundary"),
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for status, (marker, color, label) in styles.items():
        selected = [(seed, lam) for seed, lam, item in classes if item == status]
        if not selected:
            continue
        axes[0].scatter(
            [item[0] for item in selected],
            [item[1] for item in selected],
            marker=marker,
            color=color,
            s=62,
            label=label,
            zorder=3,
        )
    axes[0].set_xscale("log")
    axes[0].set_ylim(0.19, 0.31)
    axes[0].set_yticks((0.20, 0.25, 0.30))
    axes[0].set_xlabel(r"Seed mass [$M_\odot$]")
    axes[0].set_ylabel(r"Capture coefficient $\lambda_B$")
    axes[0].set_title("Sparse capture-closure map")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].grid(alpha=0.25)
    axes[0].text(
        0.98,
        0.04,
        "Discrete calculations; no interpolation",
        transform=axes[0].transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )

    for lam, color in ((0.20, "#777777"), (0.30, "#2878b5")):
        selected = sorted(
            (row for row in capture if float(row["dark_bondi_lambda"]) == lam),
            key=lambda row: float(row["black_hole_seed_msun"]),
        )
        axes[1].plot(
            [float(row["black_hole_seed_msun"]) for row in selected],
            [float(row["final_black_hole_mass_msun"]) for row in selected],
            marker="o",
            color=color,
            label=rf"$\lambda_B={lam:.2f}$",
        )
    axes[1].axhline(1.0e7, color="#d62728", linestyle="--", label=r"$10^7M_\odot$")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"Seed mass [$M_\odot$]")
    axes[1].set_ylabel(r"Final black-hole mass [$M_\odot$]")
    axes[1].set_title("Sensitivity to the capture coefficient")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
