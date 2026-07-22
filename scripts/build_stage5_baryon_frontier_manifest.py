"""Build the stage-5 baryon screen on the resolved transport frontier."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_baryon_frontier.tsv"
REDSHIFTS = (10.0, 20.0, 30.0)
BARYON_FRACTIONS = (0.01, 0.05, 0.16)
SCALE_RADII_OVER_RS = (0.003, 0.01, 0.03)
ASSEMBLY_TIMES_MYR = (0.2, 0.65, 1.5)


def main() -> None:
    rows = []
    for redshift in REDSHIFTS:
        for baryon_fraction in BARYON_FRACTIONS:
            for scale_radius in SCALE_RADII_OVER_RS:
                for assembly_time in ASSEMBLY_TIMES_MYR:
                    rows.append(
                        {
                            "task_id": len(rows),
                            "halo_redshift": redshift,
                            "baryon_fraction": baryon_fraction,
                            "scale_radius_over_rs": scale_radius,
                            "assembly_time_myr": assembly_time,
                        }
                    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 81:
        raise RuntimeError(f"expected 81 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
