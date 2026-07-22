"""Build the resolved stage-5 SIDM transport screening matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_transport.tsv"
HALO_MASSES_MSUN = (1.0e6, 1.0e7, 1.0e8, 1.0e9)
REDSHIFTS = (10.0, 20.0, 30.0)
CONCENTRATIONS = (3.0, 5.0, 8.0)
CROSS_SECTIONS_CM2_G = (10.0, 30.0, 50.0, 100.0)
SEED_TO_HALO_RATIO = 1.0e-4


def main() -> None:
    rows = []
    for halo_mass in HALO_MASSES_MSUN:
        for redshift in REDSHIFTS:
            for concentration in CONCENTRATIONS:
                for cross_section in CROSS_SECTIONS_CM2_G:
                    rows.append(
                        {
                            "task_id": len(rows),
                            "halo_mass_msun": halo_mass,
                            "halo_redshift": redshift,
                            "halo_concentration": concentration,
                            "black_hole_seed_msun": (
                                SEED_TO_HALO_RATIO * halo_mass
                            ),
                            "sigma_over_m_cm2_g": cross_section,
                            "baryon_fraction": 0.05,
                            "scale_radius_over_rs": 0.01,
                        }
                    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 144:
        raise RuntimeError(f"expected 144 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
