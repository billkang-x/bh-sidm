"""Plot stage-2 heat-flow comparisons and radial histories."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nfw-on", type=Path, required=True)
    parser.add_argument("--nfw-off", type=Path, required=True)
    parser.add_argument("--sis-on", type=Path, required=True)
    parser.add_argument("--sis-off", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {name: data[name].copy() for name in data.files}


def plot_growth(cases: dict[str, tuple[dict, dict]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for column, (profile, (heat_on, heat_off)) in enumerate(cases.items()):
        for data, label, style in (
            (heat_on, "Heat flow", "-"),
            (heat_off, "No heat flow", "--"),
        ):
            axes[0, column].plot(
                data["times_myr"], data["black_hole_mass_msun"], style, label=label
            )
            rate = np.where(data["accretion_rate_msun_myr"] > 0.0,
                            data["accretion_rate_msun_myr"], np.nan)
            axes[1, column].semilogy(data["times_myr"], rate, style, label=label)
        axes[0, column].set_title(profile)
        axes[0, column].set_ylabel("Black-hole mass [M_sun]")
        axes[1, column].set_ylabel("Accretion rate [M_sun/Myr]")
        axes[1, column].set_xlabel("Time [Myr]")
        axes[0, column].grid(alpha=0.25)
        axes[1, column].grid(alpha=0.25)
        axes[0, column].legend(frameon=False)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_profiles(profile: str, heat_on: dict, heat_off: dict, output: Path) -> None:
    target_times = np.array([0.0, 0.4, 0.8, 1.2, 1.6, 2.0])
    indices = [int(np.argmin(np.abs(heat_on["times_myr"] - time))) for time in target_times]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(indices)))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    radius = heat_on["radii_pc"]
    for index, time, color in zip(indices, target_times, colors, strict=True):
        axes[0].loglog(radius, heat_on["density_msun_pc3"][index], color=color,
                       label=f"{time:.1f} Myr")
        axes[1].semilogx(radius, heat_on["radial_velocity_km_s"][index], color=color)
        axes[2].semilogx(radius, heat_on["velocity_dispersion_km_s"][index], color=color)
    axes[0].loglog(
        heat_off["radii_pc"], heat_off["density_msun_pc3"][-1], "k--",
        linewidth=1.2, label="2.0 Myr, no heat"
    )
    axes[1].semilogx(
        heat_off["radii_pc"], heat_off["radial_velocity_km_s"][-1], "k--",
        linewidth=1.2
    )
    axes[2].semilogx(
        heat_off["radii_pc"], heat_off["velocity_dispersion_km_s"][-1], "k--",
        linewidth=1.2
    )
    axes[0].set_ylabel("Density [M_sun/pc^3]")
    axes[1].set_ylabel("Radial velocity [km/s]")
    axes[2].set_ylabel("Velocity dispersion [km/s]")
    for axis in axes:
        axis.set_xlabel("Radius [pc]")
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{profile}: MC-Roe radial evolution with heat flow")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    nfw_on, nfw_off = load(args.nfw_on), load(args.nfw_off)
    sis_on, sis_off = load(args.sis_on), load(args.sis_off)
    plot_growth(
        {"NFW": (nfw_on, nfw_off), "SIS": (sis_on, sis_off)},
        args.output_dir / "stage2_growth_comparison.png",
    )
    plot_profiles("NFW", nfw_on, nfw_off, args.output_dir / "stage2_nfw_profiles.png")
    plot_profiles("SIS", sis_on, sis_off, args.output_dir / "stage2_sis_profiles.png")
    print(args.output_dir.resolve())


if __name__ == "__main__":
    main()
