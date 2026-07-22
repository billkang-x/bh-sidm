"""Diagnose physical timescale matching in the stage-3 HPC matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
RESULTS = ROOT / "results" / "stage3"
MATRIX = RESULTS / "hpc_matrix"
FIGURES = RESULTS / "figures"

from sidm_bh.mesh import SphericalGrid
from sidm_bh.state import FluidState
from sidm_bh.timescales import (
    inward_flux_median_radius_code,
    local_timescale_profiles_code,
)
from sidm_bh.units import SimulationScales


SCALES = SimulationScales(30.0, 3.7)
TIME_COLUMNS = [
    "dynamical_time_myr",
    "collision_time_myr",
    "conduction_radius_time_myr",
    "conduction_gradient_time_myr",
    "inflow_time_myr",
]


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def evaluate_position(
    profile,
    index: int,
    position: str,
    metadata: dict,
    sample_time_myr: float,
    target_time_myr: float,
    sample_index: int,
    state: FluidState,
    radii_pc: np.ndarray,
    baryon_fraction: float,
    feed_radius_pc: float,
) -> dict:
    def physical_time(values: np.ndarray) -> float:
        value = float(values[index])
        return SCALES.time_from_code(value) if np.isfinite(value) else float("inf")

    return {
        "task_id": metadata["task_id"],
        "sigma_over_m_cm2_g": metadata["sigma_over_m_cm2_g"],
        "scale_radius_over_rs": metadata["scale_radius_over_rs"],
        "assembly_time_myr": metadata["assembly_time_myr"],
        "phase": metadata["phase"],
        "target_time_myr": target_time_myr,
        "sample_time_myr": sample_time_myr,
        "sample_index": sample_index,
        "position": position,
        "radius_pc": float(radii_pc[index]),
        "feed_radius_pc": feed_radius_pc,
        "baryon_mass_fraction": baryon_fraction,
        "density_msun_pc3": SCALES.density_from_code(state.density[index]),
        "radial_velocity_km_s": SCALES.velocity_from_code(
            state.radial_velocity[index]
        ),
        "velocity_dispersion_km_s": SCALES.velocity_from_code(
            state.velocity_dispersion[index]
        ),
        "dynamical_time_myr": physical_time(profile.dynamical_code),
        "collision_time_myr": physical_time(profile.collision_code),
        "conduction_radius_time_myr": physical_time(
            profile.conduction_radius_code
        ),
        "conduction_gradient_time_myr": physical_time(
            profile.conduction_gradient_code
        ),
        "inflow_time_myr": physical_time(profile.inflow_code),
        "thermal_length_pc": SCALES.radius_from_code(
            profile.thermal_length_code[index]
        ),
        "knudsen_number": float(profile.knudsen_number[index]),
    }


def build_rows() -> list[dict]:
    rows = []
    for path in sorted(MATRIX.glob("task_*.npz")):
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            if metadata["baryon_fraction"] == 0.0:
                continue
            metadata["task_id"] = int(path.stem.split("_")[1])
            times = data["times_myr"]
            assembly_time = float(metadata["assembly_time_myr"])
            early_time = float(times[1])
            targets = {
                "mid_assembly": assembly_time / 2.0 if assembly_time > 0.0 else early_time,
                "assembly_end": assembly_time if assembly_time > 0.0 else early_time,
                "final": float(times[-1]),
            }
            grid = SphericalGrid.from_log_spacing(
                SCALES.radius_to_code(metadata["r_min_pc"]),
                SCALES.radius_to_code(metadata["r_max_pc"]),
                int(metadata["cells"]),
            )
            radii_pc = data["radii_pc"]
            hernquist_radius_pc = (
                float(metadata["scale_radius_over_rs"])
                * SCALES.radius_scale_pc
            )
            hernquist_index = nearest_index(radii_pc, hernquist_radius_pc)
            sigma_code = SCALES.sigma_over_m_to_code(
                metadata["sigma_over_m_cm2_g"]
            )
            for phase, target_time in targets.items():
                sample_index = nearest_index(times, target_time)
                state = FluidState(
                    density=(
                        data["density_msun_pc3"][sample_index]
                        / SCALES.density_scale_msun_pc3
                    ),
                    radial_velocity=(
                        data["radial_velocity_km_s"][sample_index]
                        / SCALES.velocity_scale_km_s
                    ),
                    velocity_dispersion=(
                        data["velocity_dispersion_km_s"][sample_index]
                        / SCALES.velocity_scale_km_s
                    ),
                )
                mass_fraction = float(data["baryon_mass_fraction"][sample_index])
                baryon_mass_code = (
                    data["full_baryon_enclosed_mass_msun"]
                    * mass_fraction
                    / SCALES.mass_scale_msun
                )
                black_hole_mass_code = SCALES.mass_to_code(
                    float(data["black_hole_mass_msun"][sample_index])
                )
                profile = local_timescale_profiles_code(
                    state,
                    grid,
                    sigma_code,
                    black_hole_mass_code=black_hole_mass_code,
                    baryon_enclosed_mass_code=baryon_mass_code,
                )
                feed_radius_code = inward_flux_median_radius_code(
                    state,
                    grid,
                    maximum_radius_code=1.0,
                )
                if np.isfinite(feed_radius_code):
                    feed_radius_pc = SCALES.radius_from_code(feed_radius_code)
                    feed_index = nearest_index(
                        grid.centers_code,
                        feed_radius_code,
                    )
                else:
                    feed_radius_pc = float("nan")
                    feed_index = hernquist_index
                base_metadata = {
                    **metadata,
                    "phase": phase,
                }
                rows.append(
                    evaluate_position(
                        profile,
                        hernquist_index,
                        "hernquist_radius",
                        base_metadata,
                        float(times[sample_index]),
                        target_time,
                        sample_index,
                        state,
                        radii_pc,
                        mass_fraction,
                        feed_radius_pc,
                    )
                )
                rows.append(
                    evaluate_position(
                        profile,
                        feed_index,
                        "inward_flux_median",
                        base_metadata,
                        float(times[sample_index]),
                        target_time,
                        sample_index,
                        state,
                        radii_pc,
                        mass_fraction,
                        feed_radius_pc,
                    )
                )
    return rows


def initial_row_at_radius(
    task_id: int,
    radius_pc: float,
    position: str,
    feed_radius_pc: float,
) -> dict:
    path = MATRIX / f"task_{task_id:03d}.npz"
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        metadata.update(task_id=task_id, phase="initial")
        grid = SphericalGrid.from_log_spacing(
            SCALES.radius_to_code(metadata["r_min_pc"]),
            SCALES.radius_to_code(metadata["r_max_pc"]),
            int(metadata["cells"]),
        )
        state = FluidState(
            density=data["density_msun_pc3"][0] / SCALES.density_scale_msun_pc3,
            radial_velocity=(
                data["radial_velocity_km_s"][0] / SCALES.velocity_scale_km_s
            ),
            velocity_dispersion=(
                data["velocity_dispersion_km_s"][0]
                / SCALES.velocity_scale_km_s
            ),
        )
        profile = local_timescale_profiles_code(
            state,
            grid,
            SCALES.sigma_over_m_to_code(metadata["sigma_over_m_cm2_g"]),
            black_hole_mass_code=SCALES.mass_to_code(
                float(data["black_hole_mass_msun"][0])
            ),
            baryon_enclosed_mass_code=np.zeros(grid.num_cells),
        )
        index = nearest_index(data["radii_pc"], radius_pc)
        return evaluate_position(
            profile,
            index,
            position,
            metadata,
            0.0,
            0.0,
            0,
            state,
            data["radii_pc"],
            0.0,
            feed_radius_pc,
        )


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    with (RESULTS / "timescale_profiles_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with (RESULTS / "hpc_optima_summary.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        optima = list(csv.DictReader(stream))
    row_lookup = {
        (
            int(row["task_id"]),
            row["phase"],
            row["position"],
        ): row
        for row in rows
    }
    ratio_rows = []
    for optimum in optima:
        assembly_time = float(optimum["optimal_assembly_time_myr"])
        if assembly_time <= 0.0:
            continue
        task_id = int(optimum["task_id"])
        midpoint_rows = [
            row_lookup[(task_id, "mid_assembly", "hernquist_radius")],
            row_lookup[(task_id, "mid_assembly", "inward_flux_median")],
        ]
        feed_radius_pc = float(midpoint_rows[1]["feed_radius_pc"])
        initial_rows = [
            initial_row_at_radius(
                task_id,
                float(midpoint_rows[0]["radius_pc"]),
                "initial_hernquist_radius",
                feed_radius_pc,
            ),
            initial_row_at_radius(
                task_id,
                feed_radius_pc,
                "initial_future_feed_radius",
                feed_radius_pc,
            ),
        ]
        for row in midpoint_rows + initial_rows:
            position = row["position"]
            ratio_row = {
                "task_id": task_id,
                "sigma_over_m_cm2_g": float(optimum["sigma_over_m_cm2_g"]),
                "scale_radius_over_rs": float(optimum["scale_radius_over_rs"]),
                "optimal_assembly_time_myr": assembly_time,
                "position": position,
                "radius_pc": row["radius_pc"],
                "feed_radius_pc": row["feed_radius_pc"],
            }
            for column in TIME_COLUMNS:
                value = float(row[column])
                ratio_row[column] = value
                ratio_row[f"assembly_over_{column}"] = (
                    assembly_time / value if np.isfinite(value) and value > 0.0 else 0.0
                )
            ratio_rows.append(ratio_row)
    with (RESULTS / "timescale_optimum_ratios.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ratio_rows[0]))
        writer.writeheader()
        writer.writerows(ratio_rows)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    colors = {
        "dynamical_time_myr": "#3b4cc0",
        "collision_time_myr": "#2a9d8f",
        "conduction_radius_time_myr": "#d62828",
        "conduction_gradient_time_myr": "#f4a261",
        "inflow_time_myr": "#6c757d",
    }
    labels = {
        "dynamical_time_myr": "Dynamical",
        "collision_time_myr": "Collision",
        "conduction_radius_time_myr": "Conduction (L=r)",
        "conduction_gradient_time_myr": "Conduction (gradient L)",
        "inflow_time_myr": "Inflow",
    }
    positions = [
        "initial_hernquist_radius",
        "initial_future_feed_radius",
        "hernquist_radius",
        "inward_flux_median",
    ]
    titles = {
        "initial_hernquist_radius": "Initial state at Hernquist radius",
        "initial_future_feed_radius": "Initial state at future feed radius",
        "hernquist_radius": "Assembly midpoint at Hernquist radius",
        "inward_flux_median": "Assembly midpoint at feed radius",
    }
    for axis, position in zip(
        axes.flat,
        positions,
        strict=True,
    ):
        selected = [row for row in ratio_rows if row["position"] == position]
        assembly = np.array([row["optimal_assembly_time_myr"] for row in selected])
        for column in TIME_COLUMNS:
            local = np.array([row[column] for row in selected], dtype=float)
            finite = np.isfinite(local) & (local > 0.0)
            axis.scatter(
                local[finite],
                assembly[finite],
                color=colors[column],
                label=labels[column],
                s=38,
            )
        limits = [1.0e-4, 1.0e3]
        axis.plot(limits, limits, color="black", linestyle=":", linewidth=1.0)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(limits)
        axis.set_ylim(0.04, 2.0)
        axis.set_xlabel("Local timescale [Myr]")
        axis.set_ylabel("Optimal assembly time [Myr]")
        axis.grid(alpha=0.2)
        axis.set_title(titles[position])
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "stage3_timescale_matching.png", dpi=180)
    plt.close(fig)

    causal_rows = [
        row
        for row in ratio_rows
        if row["position"] == "initial_future_feed_radius"
    ]
    statistics = {}
    for column in TIME_COLUMNS[:-1]:
        local = np.array([row[column] for row in causal_rows], dtype=float)
        assembly = np.array(
            [row["optimal_assembly_time_myr"] for row in causal_rows],
            dtype=float,
        )
        finite = np.isfinite(local) & (local > 0.0)
        log_ratio = np.log10(assembly[finite] / local[finite])
        correlation = float(
            np.corrcoef(np.log10(local[finite]), np.log10(assembly[finite]))[0, 1]
        )
        statistics[column] = {
            "median_assembly_over_local": float(
                np.median(assembly[finite] / local[finite])
            ),
            "minimum_assembly_over_local": float(
                np.min(assembly[finite] / local[finite])
            ),
            "maximum_assembly_over_local": float(
                np.max(assembly[finite] / local[finite])
            ),
            "rms_log10_ratio_dex": float(np.sqrt(np.mean(log_ratio**2))),
            "maximum_absolute_log10_ratio_dex": float(np.max(np.abs(log_ratio))),
            "log10_pearson_correlation": correlation,
        }
    (RESULTS / "timescale_match_statistics.json").write_text(
        json.dumps(statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    representative = [
        (7, "Global maximum"),
        (128, "High-cross-section compact optimum"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for axis, (task_id, title) in zip(axes, representative, strict=True):
        path = MATRIX / f"task_{task_id:03d}.npz"
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            grid = SphericalGrid.from_log_spacing(
                SCALES.radius_to_code(metadata["r_min_pc"]),
                SCALES.radius_to_code(metadata["r_max_pc"]),
                int(metadata["cells"]),
            )
            state = FluidState(
                density=(
                    data["density_msun_pc3"][0]
                    / SCALES.density_scale_msun_pc3
                ),
                radial_velocity=(
                    data["radial_velocity_km_s"][0]
                    / SCALES.velocity_scale_km_s
                ),
                velocity_dispersion=(
                    data["velocity_dispersion_km_s"][0]
                    / SCALES.velocity_scale_km_s
                ),
            )
            profile = local_timescale_profiles_code(
                state,
                grid,
                SCALES.sigma_over_m_to_code(metadata["sigma_over_m_cm2_g"]),
                black_hole_mass_code=SCALES.mass_to_code(
                    float(data["black_hole_mass_msun"][0])
                ),
                baryon_enclosed_mass_code=np.zeros(grid.num_cells),
            )
            radius = data["radii_pc"]
            mask = radius <= 30.0
            axis.loglog(
                radius[mask],
                SCALES.time_from_code(profile.dynamical_code[mask]),
                color="#3b4cc0",
                label="Dynamical",
            )
            axis.loglog(
                radius[mask],
                SCALES.time_from_code(profile.collision_code[mask]),
                color="#2a9d8f",
                label="Collision",
            )
            axis.loglog(
                radius[mask],
                SCALES.time_from_code(profile.conduction_radius_code[mask]),
                color="#d62828",
                label="Conduction (L=r)",
            )
            matching = next(
                row
                for row in causal_rows
                if int(row["task_id"]) == task_id
            )
            assembly_time = float(matching["optimal_assembly_time_myr"])
            hernquist_radius = (
                float(metadata["scale_radius_over_rs"])
                * SCALES.radius_scale_pc
            )
            feed_radius = float(matching["feed_radius_pc"])
            axis.axhline(
                assembly_time,
                color="black",
                linestyle="--",
                label="Optimal assembly time",
            )
            axis.axvline(
                hernquist_radius,
                color="#6c757d",
                linestyle="--",
                label="Hernquist radius",
            )
            axis.axvline(
                feed_radius,
                color="#f4a261",
                linestyle=":",
                linewidth=2.0,
                label="Future feed radius",
            )
            axis.set_title(title)
            axis.set_xlabel("Radius [pc]")
            axis.set_ylabel("Initial local timescale [Myr]")
            axis.grid(alpha=0.2)
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "stage3_initial_timescale_profiles.png", dpi=180)
    plt.close(fig)

    for position in positions:
        selected = [row for row in ratio_rows if row["position"] == position]
        print(position)
        for column in TIME_COLUMNS:
            ratios = np.array(
                [row[f"assembly_over_{column}"] for row in selected],
                dtype=float,
            )
            ratios = ratios[ratios > 0.0]
            if len(ratios) == 0:
                print(column, "unavailable")
                continue
            print(
                column,
                "median=",
                float(np.median(ratios)),
                "range=",
                (float(ratios.min()), float(ratios.max())),
            )
    print(RESULTS / "timescale_profiles_summary.csv")
    print(RESULTS / "timescale_optimum_ratios.csv")
    print(RESULTS / "timescale_match_statistics.json")


if __name__ == "__main__":
    main()
