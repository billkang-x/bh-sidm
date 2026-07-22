"""Build the first stage-5 cosmological halo anchor matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_redshift_mass.tsv"
HALO_MASSES_MSUN = (1.0e6, 1.0e7, 1.0e8, 1.0e9)
REDSHIFTS = (10.0, 15.0, 20.0, 25.0, 30.0)


def main() -> None:
    rows = []
    for halo_mass in HALO_MASSES_MSUN:
        for redshift in REDSHIFTS:
            rows.append(
                {
                    "task_id": len(rows),
                    "halo_mass_msun": halo_mass,
                    "halo_redshift": redshift,
                    "halo_concentration": 4.0,
                    "black_hole_seed_msun": 100.0,
                    "baryon_fraction": 0.05,
                    "scale_radius_over_rs": 0.01,
                    "sigma_over_m_cm2_g": 50.0,
                }
            )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 20:
        raise RuntimeError(f"expected 20 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
