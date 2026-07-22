"""Build matched frozen/evolving Bondi-ambient feedback experiments."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_coupled_ambient.tsv"
EFFICIENCIES = [
    1.0e-7,
    3.0e-7,
    1.0e-6,
    3.0e-6,
    1.0e-5,
    3.0e-5,
    1.0e-4,
    3.0e-4,
    1.0e-3,
]
HEATING_FRACTIONS = [0.0, 0.5, 1.0]
AMBIENT_MODELS = ["fixed", "evolving"]
THRESHOLD_CASES = [
    (1.3e-4, 0.0),
    (1.7e-4, 0.0),
    (2.2e-4, 0.0),
    (4.0e-5, 0.5),
    (5.0e-5, 0.5),
    (7.0e-5, 0.5),
    (4.0e-5, 1.0),
    (5.0e-5, 1.0),
    (7.0e-5, 1.0),
]


def main() -> None:
    rows = []
    for ambient_model in AMBIENT_MODELS:
        rows.append(
            {
                "task_id": len(rows),
                "ambient_model": ambient_model,
                "feedback_efficiency": 0.0,
                "feedback_heating_fraction": 0.0,
                "feedback_eta": 0.5,
            }
        )
    for efficiency in EFFICIENCIES:
        for heating_fraction in HEATING_FRACTIONS:
            for ambient_model in AMBIENT_MODELS:
                rows.append(
                    {
                        "task_id": len(rows),
                        "ambient_model": ambient_model,
                        "feedback_efficiency": efficiency,
                        "feedback_heating_fraction": heating_fraction,
                        "feedback_eta": 0.5,
                    }
                )
    for efficiency, heating_fraction in THRESHOLD_CASES:
        rows.append(
            {
                "task_id": len(rows),
                "ambient_model": "evolving",
                "feedback_efficiency": efficiency,
                "feedback_heating_fraction": heating_fraction,
                "feedback_eta": 0.5,
            }
        )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 65:
        raise RuntimeError(f"expected 65 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
