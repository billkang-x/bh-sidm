"""Analyze stage-3 Hernquist assembly-history experiments."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "stage3"
FIGURES = RESULTS / "figures"
ASSEMBLY_TIMES = [0.0, 0.05, 0.2, 0.5, 1.0]
SCALE_RATIOS = [0.01, 0.1]
COLORS = ["#3b4cc0", "#2a9d8f", "#e9c46a", "#f4a261", "#d62828"]


def load(path: Path) -> dict:
    with np.load(path) as data:
        result = {name: data[name].copy() for name in data.files}
    result["metadata"] = json.loads(str(result.pop("metadata_json")))
    return result


def time_label(value: float) -> str:
    return "0" if value == 0.0 else str(value)


def assembly_path(scale_ratio: float, assembly_time: float, r_min: float = 0.005) -> Path:
    return RESULTS / (
        f"assembly_fb0.05_a{scale_ratio}_t{time_label(assembly_time)}_rmin{r_min}.npz"
    )


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    control = load(RESULTS / "control_fb0_rmin0.005.npz")
    control_growth = float(control["black_hole_mass_msun"][-1] - 100.0)
    static = {
        scale: load(RESULTS / f"fb0.05_a{scale}_rmin0.005.npz")
        for scale in SCALE_RATIOS
    }
    assembly = {
        (scale, time): load(assembly_path(scale, time))
        for scale in SCALE_RATIOS
        for time in ASSEMBLY_TIMES
    }

    rows = []
    for scale in SCALE_RATIOS:
        static_growth = float(static[scale]["black_hole_mass_msun"][-1] - 100.0)
        rows.append(
            {
                "scale_radius_over_rs": scale,
                "protocol": "static_equilibrium",
                "assembly_time_myr": "",
                "final_black_hole_mass_msun": static_growth + 100.0,
                "accreted_dark_matter_msun": static_growth,
                "enhancement_over_control": static_growth / control_growth,
                "ratio_to_static_equilibrium": 1.0,
                "peak_accretion_rate_msun_myr": static[scale]["metadata"][
                    "peak_accretion_rate_msun_myr"
                ],
                "final_accretion_rate_msun_myr": static[scale]["metadata"][
                    "final_accretion_rate_msun_myr"
                ],
            }
        )
        for time in ASSEMBLY_TIMES:
            data = assembly[(scale, time)]
            growth = float(data["black_hole_mass_msun"][-1] - 100.0)
            rows.append(
                {
                    "scale_radius_over_rs": scale,
                    "protocol": data["metadata"]["baryon_protocol"],
                    "assembly_time_myr": time,
                    "final_black_hole_mass_msun": growth + 100.0,
                    "accreted_dark_matter_msun": growth,
                    "enhancement_over_control": growth / control_growth,
                    "ratio_to_static_equilibrium": growth / static_growth,
                    "peak_accretion_rate_msun_myr": data["metadata"][
                        "peak_accretion_rate_msun_myr"
                    ],
                    "final_accretion_rate_msun_myr": data["metadata"][
                        "final_accretion_rate_msun_myr"
                    ],
                }
            )
    with (RESULTS / "assembly_history_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for axis, scale in zip(axes, SCALE_RATIOS, strict=True):
        growth = [
            assembly[(scale, time)]["black_hole_mass_msun"][-1] - 100.0
            for time in ASSEMBLY_TIMES
        ]
        static_growth = static[scale]["black_hole_mass_msun"][-1] - 100.0
        axis.plot(ASSEMBLY_TIMES, growth, color="#d62828", marker="o", label="Assembled")
        axis.axhline(static_growth, color="#264653", linestyle="--", label="Static equilibrium")
        axis.axhline(control_growth, color="#6c757d", linestyle=":", label="No baryons")
        axis.set_title(f"a_b/r_s = {scale}")
        axis.set_xlabel("Assembly time [Myr]")
        axis.set_ylabel("Accreted dark mass at 2 Myr [M_sun]")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    fig.savefig(FIGURES / "stage3_assembly_time_response.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    for axis, scale in zip(axes, SCALE_RATIOS, strict=True):
        axis.plot(
            control["times_myr"],
            control["black_hole_mass_msun"],
            color="#6c757d",
            linestyle=":",
            label="No baryons",
        )
        axis.plot(
            static[scale]["times_myr"],
            static[scale]["black_hole_mass_msun"],
            color="#264653",
            linestyle="--",
            label="Static equilibrium",
        )
        for time, color in zip(ASSEMBLY_TIMES, COLORS, strict=True):
            data = assembly[(scale, time)]
            label = "Instant turn-on" if time == 0.0 else f"T_asm={time} Myr"
            axis.semilogy(
                data["times_myr"],
                data["black_hole_mass_msun"],
                color=color,
                label=label,
            )
        axis.set_title(f"a_b/r_s = {scale}")
        axis.set_xlabel("Time [Myr]")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Black-hole mass [M_sun]")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "stage3_assembly_growth.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    radial_protocols = [
        ("Static equilibrium", None, "#264653", "--"),
        ("Instant turn-on", 0.0, "#3b4cc0", "-"),
        ("T_asm=0.5 Myr", 0.5, "#d62828", "-"),
        ("T_asm=1.0 Myr", 1.0, "#2a9d8f", "-"),
    ]
    for column, scale in enumerate(SCALE_RATIOS):
        for label, time, color, linestyle in radial_protocols:
            data = static[scale] if time is None else assembly[(scale, time)]
            radii = data["radii_pc"]
            mask = radii <= 30.0
            axes[0, column].loglog(
                radii[mask],
                data["density_msun_pc3"][-1, mask]
                / control["density_msun_pc3"][-1, mask],
                color=color,
                linestyle=linestyle,
                label=label,
            )
            axes[1, column].semilogx(
                radii[mask],
                data["radial_velocity_km_s"][-1, mask],
                color=color,
                linestyle=linestyle,
                label=label,
            )
        axes[0, column].axhline(1.0, color="#6c757d", linewidth=0.8, linestyle=":")
        axes[0, column].set_title(f"a_b/r_s = {scale}")
        axes[0, column].set_ylabel("Density / no-baryon density")
        axes[1, column].axhline(0.0, color="#6c757d", linewidth=0.8, linestyle=":")
        axes[1, column].set_xlabel("Radius [pc]")
        axes[1, column].set_ylabel("Radial velocity [km/s]")
        for row in range(2):
            axes[row, column].grid(alpha=0.2)
    axes[0, 1].legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "stage3_assembly_radial_response.png", dpi=180)
    plt.close(fig)

    boundary_radii = [0.0025, 0.005, 0.01]
    static_prefixes = {0.01: "enhance", 0.1: "suppress"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    boundary_rows = []
    for axis, scale in zip(axes, SCALE_RATIOS, strict=True):
        static_factors = []
        assembly_factors = []
        for radius in boundary_radii:
            suffix = str(radius)
            if radius == 0.005:
                control_case = control
                static_case = static[scale]
            else:
                control_case = load(RESULTS / f"control_fb0_a0.01_rmin{suffix}.npz")
                static_case = load(
                    RESULTS
                    / f"{static_prefixes[scale]}_fb0.05_a{scale}_rmin{suffix}.npz"
                )
            assembly_case = load(assembly_path(scale, 0.5, radius))
            control_delta = control_case["black_hole_mass_msun"][-1] - 100.0
            static_factor = (
                static_case["black_hole_mass_msun"][-1] - 100.0
            ) / control_delta
            assembly_factor = (
                assembly_case["black_hole_mass_msun"][-1] - 100.0
            ) / control_delta
            static_factors.append(static_factor)
            assembly_factors.append(assembly_factor)
            boundary_rows.append(
                {
                    "scale_radius_over_rs": scale,
                    "r_min_pc": radius,
                    "static_equilibrium_enhancement": static_factor,
                    "assembly_t0.5_enhancement": assembly_factor,
                }
            )
        axis.plot(boundary_radii, static_factors, marker="o", linestyle="--", color="#264653", label="Static equilibrium")
        axis.plot(boundary_radii, assembly_factors, marker="o", color="#d62828", label="T_asm=0.5 Myr")
        axis.axhline(1.0, color="#6c757d", linestyle=":")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(f"a_b/r_s = {scale}")
        axis.set_xlabel("Inner boundary [pc]")
        axis.set_ylabel("Dark-growth enhancement")
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    fig.savefig(FIGURES / "stage3_assembly_boundary.png", dpi=180)
    plt.close(fig)
    with (RESULTS / "assembly_boundary_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(boundary_rows[0]))
        writer.writeheader()
        writer.writerows(boundary_rows)

    print((RESULTS / "assembly_history_summary.csv").resolve())
    print((RESULTS / "assembly_boundary_summary.csv").resolve())


if __name__ == "__main__":
    main()
