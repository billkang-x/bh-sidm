"""Refine the feedback efficiency where baryonic catalysis reverses."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_feedback_refinement.tsv"
EFFICIENCIES = [1.0e-7, 3.0e-7, 1.0e-6, 3.0e-6]
EXPONENTS = [0.5, 1.0]


def main() -> None:
    rows = []
    for efficiency in EFFICIENCIES:
        for exponent in EXPONENTS:
            rows.append(
                {
                    "task_id": len(rows),
                    "feedback_efficiency": efficiency,
                    "feedback_eta": exponent,
                }
            )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 8:
        raise RuntimeError(f"expected 8 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
