"""Build the controlled stage-3 concentration/time/cross-section scan."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage3_matrix.tsv"
SCALE_RATIOS = [0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3]
ASSEMBLY_TIMES_MYR = [0.0, 0.05, 0.2, 0.5, 1.0]
CROSS_SECTIONS = [10.0, 30.0, 50.0, 100.0]


def main() -> None:
    rows = []
    for cross_section in CROSS_SECTIONS:
        rows.append(
            {
                "task_id": len(rows),
                "baryon_fraction": 0.0,
                "scale_radius_over_rs": 0.1,
                "assembly_time_myr": 0.0,
                "sigma_over_m_cm2_g": cross_section,
            }
        )
    for cross_section in CROSS_SECTIONS:
        for scale_ratio in SCALE_RATIOS:
            for assembly_time in ASSEMBLY_TIMES_MYR:
                rows.append(
                    {
                        "task_id": len(rows),
                        "baryon_fraction": 0.05,
                        "scale_radius_over_rs": scale_ratio,
                        "assembly_time_myr": assembly_time,
                        "sigma_over_m_cm2_g": cross_section,
                    }
                )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 164:
        raise RuntimeError(f"expected 164 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
