"""Build the stage-4 no-feedback seed and Eddington-ratio matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_no_feedback.tsv"
SEEDS_MSUN = [10.0, 100.0, 1000.0]
EDDINGTON_RATIOS = [0.0, 0.1, 1.0]


def main() -> None:
    rows = []
    for seed_mass in SEEDS_MSUN:
        for eddington_ratio in EDDINGTON_RATIOS:
            rows.append(
                {
                    "task_id": len(rows),
                    "black_hole_mass_msun": seed_mass,
                    "eddington_ratio": eddington_ratio,
                    "duty_cycle": 1.0,
                    "radiative_efficiency": 0.1,
                    "feedback_efficiency": 0.0,
                    "feedback_eta": 0.5,
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
