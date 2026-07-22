"""Build the stage-5 physical-cooling feedback-frontier matrix."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "hpc" / "stage5_feedback_frontier.tsv"
GAS_REGIMES = {
    "transition": 0.3,
    "dense": 300.0,
}
EFFICIENCIES = (0.0, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)


def main() -> None:
    rows: list[dict] = []

    def add(
        configuration: str,
        gas_regime: str,
        assembly_time_myr: float,
        feedback_efficiency: float,
        feedback_heating_fraction: float = 0.5,
        metallicity_solar: float = 0.0,
    ) -> None:
        rows.append(
            {
                "task_id": len(rows),
                "configuration": configuration,
                "gas_regime": gas_regime,
                "gas_density_msun_pc3": GAS_REGIMES[gas_regime],
                "assembly_time_myr": assembly_time_myr,
                "feedback_efficiency": feedback_efficiency,
                "feedback_heating_fraction": feedback_heating_fraction,
                "metallicity_solar": metallicity_solar,
            }
        )

    for gas_regime in GAS_REGIMES:
        for efficiency in EFFICIENCIES:
            add("mixed_sweep", gas_regime, 1.25, efficiency)
        for assembly_time in (1.0, 1.5):
            for efficiency in (0.0, 1.0e-3, 1.0e-2):
                add("timing_control", gas_regime, assembly_time, efficiency)

    for efficiency in (1.0e-3, 1.0e-2):
        add("pure_expansion", "dense", 1.25, efficiency, 0.0)
        add("pure_heating", "dense", 1.25, efficiency, 1.0)
        add("solar_mixed", "dense", 1.25, efficiency, 0.5, 1.0)

    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 34:
        raise RuntimeError(f"expected 34 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
