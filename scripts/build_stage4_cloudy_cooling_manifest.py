"""Build the stage-4 Cloudy cooling and trapping sensitivity matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_cloudy_cooling.tsv"
CONFIGURATIONS = [
    ("heating_strong", 1.0, 1.0e-4),
    ("heating_extreme", 1.0, 1.0e-3),
    ("mixed_strong", 0.5, 1.4e-4),
    ("mixed_extreme", 0.5, 1.0e-3),
]
METALLICITIES = [0.0, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
SUPPRESSION_MULTIPLIERS = [1.0e-2, 1.0e-4, 1.0e-6]
REFINEMENT_MULTIPLIERS = [3.0e-5, 1.0e-5, 3.0e-6]


def main() -> None:
    rows = []
    rows.append(
        {
            "task_id": 0,
            "configuration": "no_feedback_control",
            "feedback_heating_fraction": 1.0,
            "feedback_efficiency": 0.0,
            "metallicity_solar": 0.0,
            "cooling_rate_multiplier": 1.0,
        }
    )
    for name, heating_fraction, efficiency in CONFIGURATIONS:
        for metallicity in METALLICITIES:
            rows.append(
                {
                    "task_id": len(rows),
                    "configuration": name,
                    "feedback_heating_fraction": heating_fraction,
                    "feedback_efficiency": efficiency,
                    "metallicity_solar": metallicity,
                    "cooling_rate_multiplier": 1.0,
                }
            )
    for name, heating_fraction, efficiency in (
        CONFIGURATIONS[1],
        CONFIGURATIONS[3],
    ):
        for metallicity in (0.0, 1.0):
            for multiplier in SUPPRESSION_MULTIPLIERS:
                rows.append(
                    {
                        "task_id": len(rows),
                        "configuration": name,
                        "feedback_heating_fraction": heating_fraction,
                        "feedback_efficiency": efficiency,
                        "metallicity_solar": metallicity,
                        "cooling_rate_multiplier": multiplier,
                    }
                )
    name, heating_fraction, efficiency = CONFIGURATIONS[1]
    for metallicity in (0.0, 1.0):
        for multiplier in REFINEMENT_MULTIPLIERS:
            rows.append(
                {
                    "task_id": len(rows),
                    "configuration": name,
                    "feedback_heating_fraction": heating_fraction,
                    "feedback_efficiency": efficiency,
                    "metallicity_solar": metallicity,
                    "cooling_rate_multiplier": multiplier,
                }
            )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 39:
        raise RuntimeError(f"expected 39 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
