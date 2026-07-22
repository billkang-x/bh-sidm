"""Analyze extended assembly times and the heat-off mechanism check."""

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
HEAT_OFF = RESULTS / "heat_off"
EXTENDED = RESULTS / "extended_assembly"
FIGURES = RESULTS / "figures"

from sidm_bh.conduction import cell_conductivity_code
from sidm_bh.mesh import SphericalGrid
from sidm_bh.state import FluidState
from sidm_bh.timescales import inward_flux_median_radius_code
from sidm_bh.timescales import local_timescale_profiles_code
from sidm_bh.units import SimulationScales


SCALES = SimulationScales(30.0, 3.7)
BASE_ASSEMBLY_TIMES = [0.0, 0.05, 0.2, 0.5, 1.0]
MAIN_TASKS = {
    10.0: [4, 5, 6, 7, 8],
    30.0: [44, 45, 46, 47, 48],
    50.0: [84, 85, 86, 87, 88],
    100.0: [124, 125, 126, 127, 128],
}
CONTROL_TASKS = {10.0: 0, 30.0: 1, 50.0: 2, 100.0: 3}
OPTIMUM_PATHS = {
    10.0: RESULTS / "refined_assembly" / "sigma10_t0.6.npz",
    30.0: MATRIX / "task_047.npz",
    50.0: RESULTS / "refined_assembly" / "sigma50_t0.65.npz",
    100.0: RESULTS / "refined_assembly" / "sigma100_t0.75.npz",
}


def load_summary(path: Path) -> tuple[dict, float]:
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        growth = float(data["black_hole_mass_msun"][-1] - 100.0)
    return metadata, growth


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    heat_off_control_metadata, heat_off_control = load_summary(
        HEAT_OFF / "control.npz"
    )
    del heat_off_control_metadata
    controls = {0.0: heat_off_control}
    for sigma, task_id in CONTROL_TASKS.items():
        _, controls[sigma] = load_summary(MATRIX / f"task_{task_id:03d}.npz")

    rows = []
    heat_off_names = {
        0.0: "fb0.05_a0.01_t0.npz",
        0.05: "fb0.05_a0.01_t0.05.npz",
        0.2: "fb0.05_a0.01_t0.2.npz",
        0.5: "fb0.05_a0.01_t0.5.npz",
        1.0: "fb0.05_a0.01_t1.0.npz",
    }
    for sigma in [0.0, 10.0, 30.0, 50.0, 100.0]:
        for index, assembly_time in enumerate(BASE_ASSEMBLY_TIMES):
            if sigma == 0.0:
                path = HEAT_OFF / heat_off_names[assembly_time]
            else:
                path = MATRIX / f"task_{MAIN_TASKS[sigma][index]:03d}.npz"
            metadata, growth = load_summary(path)
            rows.append(
                {
                    "sigma_over_m_cm2_g": sigma,
                    "assembly_time_myr": assembly_time,
                    "conduction": metadata["conduction"],
                    "accreted_dark_matter_msun": growth,
                    "enhancement_over_matched_control": growth / controls[sigma],
                    "peak_accretion_rate_msun_myr": metadata[
                        "peak_accretion_rate_msun_myr"
                    ],
                    "final_accretion_rate_msun_myr": metadata[
                        "final_accretion_rate_msun_myr"
                    ],
                }
            )
    for path in sorted(EXTENDED.glob("sigma*_t*.npz")) + sorted(
        (RESULTS / "refined_assembly").glob("sigma*_t*.npz")
    ):
        metadata, growth = load_summary(path)
        sigma = float(metadata["sigma_over_m_cm2_g"])
        rows.append(
            {
                "sigma_over_m_cm2_g": sigma,
                "assembly_time_myr": float(metadata["assembly_time_myr"]),
                "conduction": metadata["conduction"],
                "accreted_dark_matter_msun": growth,
                "enhancement_over_matched_control": growth / controls[sigma],
                "peak_accretion_rate_msun_myr": metadata[
                    "peak_accretion_rate_msun_myr"
                ],
                "final_accretion_rate_msun_myr": metadata[
                    "final_accretion_rate_msun_myr"
                ],
            }
        )
    rows.sort(key=lambda row: (row["sigma_over_m_cm2_g"], row["assembly_time_myr"]))
    with (RESULTS / "mechanism_followup_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    colors = {0.0: "#6c757d", 10.0: "#d62828", 30.0: "#8a5cf5", 50.0: "#f4a261", 100.0: "#3b4cc0"}
    for sigma in [0.0, 10.0, 30.0, 50.0, 100.0]:
        selected = [row for row in rows if row["sigma_over_m_cm2_g"] == sigma]
        times = [row["assembly_time_myr"] for row in selected]
        label = "Heat off" if sigma == 0.0 else f"sigma/m={sigma:g}"
        axes[0].plot(
            times,
            [row["accreted_dark_matter_msun"] for row in selected],
            marker="o",
            color=colors[sigma],
            label=label,
        )
        axes[1].plot(
            times,
            [row["enhancement_over_matched_control"] for row in selected],
            marker="o",
            color=colors[sigma],
            label=label,
        )
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    axes[0].set_ylabel("Accreted dark mass at 2 Myr [M_sun]")
    axes[1].set_ylabel("Enhancement over matched control")
    for axis in axes:
        axis.set_xlabel("Assembly time [Myr]")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    axes[0].set_title("Compact-potential growth")
    axes[1].set_title("Transport changes both amplitude and optimum")
    fig.savefig(FIGURES / "stage3_mechanism_followup.png", dpi=180)
    plt.close(fig)

    flux_rows = []
    refined_timescale_rows = []
    for sigma, path in OPTIMUM_PATHS.items():
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
            times = data["times_myr"]
            sample_index = int(
                np.argmin(np.abs(times - metadata["assembly_time_myr"] / 2.0))
            )
            grid = SphericalGrid.from_log_spacing(
                SCALES.radius_to_code(metadata["r_min_pc"]),
                SCALES.radius_to_code(metadata["r_max_pc"]),
                int(metadata["cells"]),
            )
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
            conductivity = cell_conductivity_code(
                state,
                SCALES.sigma_over_m_to_code(sigma),
            )
            temperature_gradient = np.gradient(
                state.velocity_dispersion**2,
                grid.centers_code,
                edge_order=2,
            )
            outward_flux = -conductivity * temperature_gradient
            feed_radius_code = inward_flux_median_radius_code(
                state,
                grid,
                maximum_radius_code=1.0,
            )
            positions = {
                "hernquist_radius": SCALES.radius_to_code(
                    metadata["scale_radius_pc"]
                ),
                "inward_flux_median": feed_radius_code,
            }
            for position, radius_code in positions.items():
                index = int(np.argmin(np.abs(grid.centers_code - radius_code)))
                normalized_flux = outward_flux[index] / (
                    state.density[index] * state.velocity_dispersion[index] ** 3
                )
                flux_rows.append(
                    {
                        "case_id": path.stem,
                        "sigma_over_m_cm2_g": sigma,
                        "assembly_time_myr": metadata["assembly_time_myr"],
                        "sample_time_myr": float(times[sample_index]),
                        "position": position,
                        "radius_pc": SCALES.radius_from_code(
                            grid.centers_code[index]
                        ),
                        "temperature_gradient_code": temperature_gradient[index],
                        "outward_heat_flux_code": outward_flux[index],
                        "normalized_outward_heat_flux": normalized_flux,
                    }
                )
            profile = local_timescale_profiles_code(
                FluidState(
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
                ),
                grid,
                SCALES.sigma_over_m_to_code(sigma),
                black_hole_mass_code=SCALES.mass_to_code(
                    float(data["black_hole_mass_msun"][0])
                ),
                baryon_enclosed_mass_code=np.zeros(grid.num_cells),
            )
            feed_index = int(
                np.argmin(np.abs(grid.centers_code - feed_radius_code))
            )
            refined_timescale_rows.append(
                {
                    "case_id": path.stem,
                    "sigma_over_m_cm2_g": sigma,
                    "optimal_assembly_time_myr": metadata["assembly_time_myr"],
                    "future_feed_radius_pc": SCALES.radius_from_code(
                        grid.centers_code[feed_index]
                    ),
                    "initial_dynamical_time_myr": SCALES.time_from_code(
                        profile.dynamical_code[feed_index]
                    ),
                    "initial_collision_time_myr": SCALES.time_from_code(
                        profile.collision_code[feed_index]
                    ),
                    "initial_conduction_time_myr": SCALES.time_from_code(
                        profile.conduction_radius_code[feed_index]
                    ),
                }
            )
    with (RESULTS / "compact_heat_flux_diagnostics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flux_rows[0]))
        writer.writeheader()
        writer.writerows(flux_rows)

    with (RESULTS / "refined_optimum_timescales.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(refined_timescale_rows[0]))
        writer.writeheader()
        writer.writerows(refined_timescale_rows)

    refined_statistics = {}
    for column in (
        "initial_dynamical_time_myr",
        "initial_collision_time_myr",
        "initial_conduction_time_myr",
    ):
        local = np.array([row[column] for row in refined_timescale_rows])
        assembly = np.array(
            [row["optimal_assembly_time_myr"] for row in refined_timescale_rows]
        )
        ratio = assembly / local
        refined_statistics[column] = {
            "median_assembly_over_local": float(np.median(ratio)),
            "minimum_assembly_over_local": float(np.min(ratio)),
            "maximum_assembly_over_local": float(np.max(ratio)),
            "rms_log10_ratio_dex": float(
                np.sqrt(np.mean(np.log10(ratio) ** 2))
            ),
        }
    (RESULTS / "refined_optimum_match_statistics.json").write_text(
        json.dumps(refined_statistics, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    fig, axis = plt.subplots(figsize=(6.3, 4.8), constrained_layout=True)
    timescale_columns = [
        ("initial_dynamical_time_myr", "Dynamical", "#3b4cc0"),
        ("initial_collision_time_myr", "Collision", "#2a9d8f"),
        ("initial_conduction_time_myr", "Conduction (L=r)", "#d62828"),
    ]
    for column, label, color in timescale_columns:
        local = np.array([row[column] for row in refined_timescale_rows])
        assembly = np.array(
            [row["optimal_assembly_time_myr"] for row in refined_timescale_rows]
        )
        axis.scatter(local, assembly, color=color, s=55, label=label)
    limits = [0.07, 1.3]
    axis.plot(limits, limits, color="black", linestyle=":")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(limits)
    axis.set_ylim(0.4, 1.0)
    axis.set_xlabel("Initial timescale at future feed radius [Myr]")
    axis.set_ylabel("Refined optimal assembly time [Myr]")
    axis.set_title("Compact-potential optimum tracks dynamical time")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    fig.savefig(FIGURES / "stage3_refined_timescale_matching.png", dpi=180)
    plt.close(fig)

    optima = {}
    for sigma in [0.0, 10.0, 30.0, 50.0, 100.0]:
        selected = [row for row in rows if row["sigma_over_m_cm2_g"] == sigma]
        best = max(selected, key=lambda row: row["accreted_dark_matter_msun"])
        optima[str(sigma)] = best
    (RESULTS / "mechanism_followup_optima.json").write_text(
        json.dumps(optima, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(RESULTS / "mechanism_followup_summary.csv")
    print(RESULTS / "compact_heat_flux_diagnostics.csv")
    print(RESULTS / "refined_optimum_timescales.csv")
    print(RESULTS / "refined_optimum_match_statistics.json")
    print(json.dumps(optima, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
