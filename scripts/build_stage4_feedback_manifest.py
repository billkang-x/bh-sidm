"""Build the stage-4 Hernquist-expansion feedback matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_feedback.tsv"
FEEDBACK_EFFICIENCIES = [1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2]
EXPANSION_EXPONENTS = [0.5, 1.0]


def main() -> None:
    rows = [
        {
            "task_id": 0,
            "feedback_efficiency": 0.0,
            "feedback_eta": 0.5,
        }
    ]
    for feedback_efficiency in FEEDBACK_EFFICIENCIES:
        for feedback_eta in EXPANSION_EXPONENTS:
            rows.append(
                {
                    "task_id": len(rows),
                    "feedback_efficiency": feedback_efficiency,
                    "feedback_eta": feedback_eta,
                }
            )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 9:
        raise RuntimeError(f"expected 9 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
