"""Build feedback tests with Bondi supply and effective binding energy."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage4_effective_feedback.tsv"
GAS_REGIMES = [
    ("transition", 300.0, 10.0),
    ("eddington_saturated", 3000.0, 10.0),
]
EFFICIENCIES = [1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5]
EXPONENTS = [0.5, 1.0]


def main() -> None:
    rows = []
    for regime, density, sound_speed in GAS_REGIMES:
        rows.append(
            {
                "task_id": len(rows),
                "gas_regime": regime,
                "gas_density_msun_pc3": density,
                "gas_sound_speed_km_s": sound_speed,
                "feedback_efficiency": 0.0,
                "feedback_eta": 0.5,
            }
        )
        for efficiency in EFFICIENCIES:
            for exponent in EXPONENTS:
                rows.append(
                    {
                        "task_id": len(rows),
                        "gas_regime": regime,
                        "gas_density_msun_pc3": density,
                        "gas_sound_speed_km_s": sound_speed,
                        "feedback_efficiency": efficiency,
                        "feedback_eta": exponent,
                    }
                )
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 18:
        raise RuntimeError(f"expected 18 tasks, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
