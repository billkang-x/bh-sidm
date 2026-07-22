"""Build finite-cooling scans for feedback-regulated Bondi accretion."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_cooling.tsv"
CONFIGURATIONS = [
    ("mixed_threshold", 0.5, 7.0e-5),
    ("mixed_strong", 0.5, 1.4e-4),
    ("heating_threshold", 1.0, 5.0e-5),
    ("heating_strong", 1.0, 1.0e-4),
]
COOLING_TIMES_MYR = [
    1.0e-3,
    3.0e-3,
    1.0e-2,
    3.0e-2,
    1.0e-1,
    3.0e-1,
    1.0,
    3.0,
    "inf",
]
REFINEMENT_CASES = [
    ("mixed_strong", 0.5, 1.4e-4, 0.15),
    ("mixed_strong", 0.5, 1.4e-4, 0.20),
    ("heating_strong", 1.0, 1.0e-4, 0.50),
    ("heating_strong", 1.0, 1.0e-4, 0.70),
    ("mixed_threshold", 0.5, 7.0e-5, 1.50),
    ("mixed_threshold", 0.5, 7.0e-5, 2.00),
    ("heating_threshold", 1.0, 5.0e-5, 4.00),
    ("heating_threshold", 1.0, 5.0e-5, 5.00),
    ("heating_threshold", 1.0, 5.0e-5, 7.00),
    ("heating_threshold", 1.0, 5.0e-5, 10.0),
    ("heating_threshold", 1.0, 5.0e-5, 30.0),
    ("heating_threshold", 1.0, 5.0e-5, 100.0),
]


def main() -> None:
    rows = []
    for name, heating_fraction, efficiency in CONFIGURATIONS:
        for cooling_time in COOLING_TIMES_MYR:
            rows.append(
                {
                    "task_id": len(rows),
                    "configuration": name,
                    "feedback_heating_fraction": heating_fraction,
                    "feedback_efficiency": efficiency,
                    "cooling_time_myr": cooling_time,
                }
            )
    for name, heating_fraction, efficiency, cooling_time in REFINEMENT_CASES:
        rows.append(
            {
                "task_id": len(rows),
                "configuration": name,
                "feedback_heating_fraction": heating_fraction,
                "feedback_efficiency": efficiency,
                "cooling_time_myr": cooling_time,
            }
        )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 48:
        raise RuntimeError(f"expected 48 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
