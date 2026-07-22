"""Build the controlled six-parameter SIDM applicability map."""

from __future__ import annotations

import csv
from math import log, sqrt
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sidm_bh.constants import G_CGS, M_SUN_CGS, MYR_CGS, PC_CGS
from sidm_bh.cosmology import FlatLambdaCDM
from sidm_bh.halos import NFWProfile


OUTPUT = ROOT / "hpc" / "stage5_applicability_map.tsv"
REFERENCE_SEED_MSUN = 1.0e5
SCREEN_BOUNDARY = 5.0 / 24.0
LATIN_HYPERCUBE_SIZE = 64
RANDOM_SEED = 20260720

MODELS = (
    ("constant_sigma1", "constant", 1.0, 0.0),
    ("vd_low_transport", "rutherford", 30.0, 10.0),
    ("vd_matched_transport", "rutherford", 30.0, 30.0),
    ("vd_high_transport", "rutherford", 30.0, 100.0),
)

BASE = {
    "halo_mass_msun": 1.0e9,
    "halo_redshift": 30.0,
    "halo_concentration": 8.0,
    "black_hole_seed_msun": 1.0e5,
    "scale_radius_over_rs": 0.02,
    "assembly_time_myr": 1.25,
}

AXIS_LEVELS = {
    "halo_mass_msun": (1.0e8, 3.0e8, 1.0e9, 3.0e9, 1.0e10),
    "halo_redshift": (15.0, 20.0, 25.0, 30.0),
    "halo_concentration": (4.0, 6.0, 8.0, 10.0, 12.0),
    "black_hole_seed_msun": (
        1.0e2,
        1.0e3,
        1.0e4,
        2.0e4,
        3.0e4,
        4.0e4,
        1.0e5,
    ),
    "scale_radius_over_rs": (0.005, 0.01, 0.02, 0.04, 0.08),
    "assembly_time_myr": (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75),
}


def latin_hypercube(size: int, dimensions: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    samples = np.empty((size, dimensions))
    for dimension in range(dimensions):
        samples[:, dimension] = (
            rng.permutation(size) + rng.random(size)
        ) / size
    return samples


def map_latin_hypercube(samples: np.ndarray) -> list[dict[str, float]]:
    rows = []
    for values in samples:
        rows.append(
            {
                "halo_mass_msun": 10.0 ** (8.0 + 2.0 * values[0]),
                "halo_redshift": 15.0 + 15.0 * values[1],
                "halo_concentration": 4.0 + 8.0 * values[2],
                "black_hole_seed_msun": 10.0 ** (2.0 + 3.0 * values[3]),
                "scale_radius_over_rs": 10.0 ** (
                    log(0.005, 10.0)
                    + values[4] * log(0.08 / 0.005, 10.0)
                ),
                "assembly_time_myr": 0.25 + 1.5 * values[5],
            }
        )
    return rows


def dynamical_time_myr(parameters: dict[str, float]) -> tuple[float, float]:
    cosmology = FlatLambdaCDM()
    virial_radius_pc = cosmology.spherical_overdensity_radius_pc(
        parameters["halo_mass_msun"],
        parameters["halo_redshift"],
    )
    profile = NFWProfile.from_mass_concentration(
        parameters["halo_mass_msun"],
        virial_radius_pc,
        parameters["halo_concentration"],
    )
    match_radius_over_rs = min(
        max(8.0 * parameters["scale_radius_over_rs"], 0.05),
        0.4,
    )
    match_radius_pc = match_radius_over_rs * profile.scale_radius_pc
    enclosed_mass_msun = (
        profile.enclosed_mass_msun(match_radius_pc)
        + parameters["black_hole_seed_msun"]
    )
    time_myr = sqrt(
        (match_radius_pc * PC_CGS) ** 3
        / (G_CGS * enclosed_mass_msun * M_SUN_CGS)
    ) / MYR_CGS
    return match_radius_over_rs, time_myr


def numerical_domain(parameters: dict[str, float]) -> tuple[float, float, int]:
    concentration = parameters["halo_concentration"]
    r_min = (
        SCREEN_BOUNDARY
        * REFERENCE_SEED_MSUN
        / parameters["halo_mass_msun"]
        * concentration
    )
    r_max = concentration * (125.0 / 6.0)
    cells = int(round(256 * log(r_max / r_min) / log(1.0e6)))
    return r_min, r_max, cells


def main_effect_rows() -> list[tuple[str, float, dict[str, float]]]:
    rows: list[tuple[str, float, dict[str, float]]] = [
        ("baseline", 0.0, BASE.copy())
    ]
    seen = {tuple(BASE.values())}
    for axis, levels in AXIS_LEVELS.items():
        for level in levels:
            parameters = {**BASE, axis: level}
            key = tuple(parameters.values())
            if key in seen:
                continue
            seen.add(key)
            rows.append((axis, level, parameters))
    if len(rows) != 28:
        raise RuntimeError(f"expected 28 main-effect points, generated {len(rows)}")
    return rows


def main() -> None:
    lhs = map_latin_hypercube(
        latin_hypercube(LATIN_HYPERCUBE_SIZE, 6, RANDOM_SEED)
    )
    rows = []
    for model_label, model, sigma0, velocity in MODELS:
        design_points = [
            ("main_effect", axis, axis_value, parameters)
            for axis, axis_value, parameters in main_effect_rows()
        ]
        design_points.extend(
            ("latin_hypercube", "joint", index, parameters)
            for index, parameters in enumerate(lhs)
        )
        for design, axis, axis_value, parameters in design_points:
            match_radius, dynamical_time = dynamical_time_myr(parameters)
            r_min, r_max, cells = numerical_domain(parameters)
            rows.append(
                {
                    "task_id": len(rows),
                    "design": design,
                    "axis": axis,
                    "axis_value": axis_value,
                    "model_label": model_label,
                    "cross_section_model": model,
                    "sigma0_over_m_cm2_g": sigma0,
                    "velocity_scale_km_s": velocity,
                    **parameters,
                    "match_radius_over_rs": match_radius,
                    "initial_dynamical_time_myr": dynamical_time,
                    "assembly_over_dynamical_time": (
                        parameters["assembly_time_myr"] / dynamical_time
                    ),
                    "r_min_over_reference_influence": SCREEN_BOUNDARY,
                    "r_min_over_rs": r_min,
                    "r_max_over_rs": r_max,
                    "cells": cells,
                }
            )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if len(rows) != 368:
        raise RuntimeError(f"expected 368 cases, generated {len(rows)}")
    print(OUTPUT.resolve())


if __name__ == "__main__":
    main()
