"""Build the stage-4 ambient gas density and sound-speed matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_bondi.tsv"
GAS_DENSITIES_MSUN_PC3 = [30.0, 300.0, 3000.0]
SOUND_SPEEDS_KM_S = [5.0, 10.0, 30.0]


def main() -> None:
    rows = []
    for sound_speed in SOUND_SPEEDS_KM_S:
        for density in GAS_DENSITIES_MSUN_PC3:
            rows.append(
                {
                    "task_id": len(rows),
                    "gas_density_msun_pc3": density,
                    "gas_sound_speed_km_s": sound_speed,
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
