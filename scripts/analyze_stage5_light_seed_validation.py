"""Validate the dark Bondi reservoir against resolved direct-flux runs."""

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

from sidm_bh.constants import G_CGS, M_SUN_CGS, MYR_CGS, PC_CGS
MANIFEST = ROOT / "hpc" / "stage5_light_seed_validation.tsv"
RESULTS = ROOT / "results" / "stage5" / "light_seed_validation"
SUMMARY = ROOT / "results" / "stage5" / "light_seed_validation_summary.csv"
STATISTICS = ROOT / "results" / "stage5" / "light_seed_validation_statistics.json"
FIGURE = ROOT / "results" / "stage5" / "figures" / "stage5_light_seed_validation.png"
DIRECT = {
    0.05208333333333334: ROOT / "results" / "stage5" / "trusted_peak_convergence" / "task_000.npz",
    0.10416666666666667: ROOT / "results" / "stage5" / "trusted_peak_convergence" / "task_001.npz",
    5.0 / 24.0: ROOT / "results" / "stage5" / "compactness_peak" / "task_002.npz",
}


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / (0.5 * (abs(first) + abs(second)))


def effective_lambda(data: np.lib.npyio.NpzFile) -> float:
    mass = data["black_hole_mass_msun"]
    rate = data["dark_matter_accretion_rate_msun_myr"]
    density = data["density_msun_pc3"][:, 0] * M_SUN_CGS / PC_CGS**3
    sound = np.sqrt(5.0 / 3.0) * data["velocity_dispersion_km_s"][:, 0] * 1.0e5
    physical_rate = (
        4.0
        * np.pi
        * G_CGS**2
        * (mass * M_SUN_CGS) ** 2
        * density
        / sound**3
        / (M_SUN_CGS / MYR_CGS)
    )
    ratio = np.divide(rate, physical_rate, out=np.full_like(rate, np.nan), where=physical_rate > 0.0)
    times = data["times_myr"]
    return float(np.nanmedian(ratio[times >= 1.25]))


def main() -> None:
    with MANIFEST.open(newline="", encoding="ascii") as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    direct = {}
    direct_lambda = {}
    for boundary, path in DIRECT.items():
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            direct[boundary] = metadata
            direct_lambda[boundary] = effective_lambda(data)

    cases = []
    for row in manifest:
        task_id = int(row["task_id"])
        boundary = float(row["reference_r_min_over_influence"])
        bondi_lambda = float(row["dark_bondi_lambda"])
        with np.load(RESULTS / f"task_{task_id:03d}.npz", allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        direct_metadata = direct[boundary]
        cases.append(
            {
                "task_id": task_id,
                "reference_r_min_over_influence": boundary,
                "dark_bondi_lambda": bondi_lambda,
                "direct_final_mass_msun": float(direct_metadata["final_black_hole_mass_msun"]),
                "reservoir_final_mass_msun": float(metadata["final_black_hole_mass_msun"]),
                "final_mass_relative_difference": relative_difference(
                    float(direct_metadata["final_black_hole_mass_msun"]),
                    float(metadata["final_black_hole_mass_msun"]),
                ),
                "dark_capture_fraction": float(metadata["dark_capture_fraction_of_available_supply"]),
                "final_dark_reservoir_msun": float(metadata["final_inner_dark_matter_reservoir_msun"]),
                "time_to_1e7_msun_myr": float(metadata["time_to_1e7_msun_myr"]),
                "mass_budget_residual_code": float(metadata["mass_budget_residual_code"]),
                "steps": int(metadata["steps"]),
            }
        )
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    calibrated = [case for case in cases if case["dark_bondi_lambda"] == 0.25]
    statistics = {
        "case_count": len(cases),
        "direct_effective_lambda_after_assembly": {
            f"{boundary:.12g}": value for boundary, value in sorted(direct_lambda.items())
        },
        "lambda_025_maximum_final_mass_relative_difference": max(
            case["final_mass_relative_difference"] for case in calibrated
        ),
        "lambda_025_passes_five_percent_validation": bool(
            all(case["final_mass_relative_difference"] < 0.05 for case in calibrated)
        ),
        "maximum_mass_budget_residual_code": max(
            case["mass_budget_residual_code"] for case in cases
        ),
    }
    STATISTICS.write_text(json.dumps(statistics, indent=2, sort_keys=True), encoding="ascii")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    calibrated.sort(key=lambda case: case["reference_r_min_over_influence"])
    axes[0].plot(
        [case["reference_r_min_over_influence"] for case in calibrated],
        [case["direct_final_mass_msun"] for case in calibrated],
        marker="o", label="resolved boundary flux",
    )
    axes[0].plot(
        [case["reference_r_min_over_influence"] for case in calibrated],
        [case["reservoir_final_mass_msun"] for case in calibrated],
        marker="s", label="Bondi reservoir",
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Feeding radius / reference influence radius")
    axes[0].set_ylabel("Final black-hole mass [M_sun]")
    axes[0].legend()

    sensitivity = sorted(
        (case for case in cases if np.isclose(case["reference_r_min_over_influence"], 0.10416666666666667)),
        key=lambda case: case["dark_bondi_lambda"],
    )
    axes[1].plot(
        [case["dark_bondi_lambda"] for case in sensitivity],
        [case["reservoir_final_mass_msun"] for case in sensitivity],
        marker="o",
    )
    axes[1].axhline(sensitivity[0]["direct_final_mass_msun"], color="black", linestyle="--")
    axes[1].set_xlabel("Dark Bondi lambda")
    axes[1].set_ylabel("Final black-hole mass [M_sun]")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(STATISTICS.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
