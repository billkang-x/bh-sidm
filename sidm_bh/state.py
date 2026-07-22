"""Fluid state containers for the 1D SIDM prototype solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FluidState:
    """Primitive SIDM variables in dimensionless code units."""

    density: np.ndarray
    radial_velocity: np.ndarray
    velocity_dispersion: np.ndarray

    def __post_init__(self) -> None:
        density = np.asarray(self.density, dtype=float)
        radial_velocity = np.asarray(self.radial_velocity, dtype=float)
        velocity_dispersion = np.asarray(self.velocity_dispersion, dtype=float)
        if density.ndim != 1:
            raise ValueError("density must be one-dimensional")
        if radial_velocity.shape != density.shape:
            raise ValueError("radial_velocity shape must match density")
        if velocity_dispersion.shape != density.shape:
            raise ValueError("velocity_dispersion shape must match density")
        if np.any(density <= 0.0):
            raise ValueError("density values must be positive")
        if np.any(velocity_dispersion <= 0.0):
            raise ValueError("velocity_dispersion values must be positive")

        object.__setattr__(self, "density", density)
        object.__setattr__(self, "radial_velocity", radial_velocity)
        object.__setattr__(self, "velocity_dispersion", velocity_dispersion)

    @property
    def pressure(self) -> np.ndarray:
        return self.density * self.velocity_dispersion**2

    @property
    def specific_energy(self) -> np.ndarray:
        return 1.5 * self.velocity_dispersion**2 + 0.5 * self.radial_velocity**2

    @property
    def conservative(self) -> np.ndarray:
        return np.vstack(
            [
                self.density,
                self.density * self.radial_velocity,
                self.density * self.specific_energy,
            ]
        )

    @classmethod
    def from_pressure(
        cls,
        density: np.ndarray,
        radial_velocity: np.ndarray,
        pressure: np.ndarray,
    ) -> "FluidState":
        density_array = np.asarray(density, dtype=float)
        pressure_array = np.asarray(pressure, dtype=float)
        if np.any(pressure_array <= 0.0):
            raise ValueError("pressure values must be positive")
        velocity_dispersion = np.sqrt(pressure_array / density_array)
        return cls(density_array, radial_velocity, velocity_dispersion)

    @classmethod
    def from_conservative(cls, conservative: np.ndarray) -> "FluidState":
        conservative_array = np.asarray(conservative, dtype=float)
        if conservative_array.ndim != 2 or conservative_array.shape[0] != 3:
            raise ValueError("conservative must have shape (3, n)")

        density = conservative_array[0]
        momentum = conservative_array[1]
        energy_density = conservative_array[2]
        if np.any(density <= 0.0):
            raise ValueError("conservative density must be positive")

        radial_velocity = momentum / density
        specific_energy = energy_density / density
        thermal_energy = specific_energy - 0.5 * radial_velocity**2
        velocity_dispersion_squared = (2.0 / 3.0) * thermal_energy
        if np.any(velocity_dispersion_squared <= 0.0):
            raise ValueError("conservative energy implies non-positive pressure")

        return cls(density, radial_velocity, np.sqrt(velocity_dispersion_squared))
