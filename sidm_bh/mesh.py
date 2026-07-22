"""Spherical finite-volume grids for the 1D prototype solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SphericalGrid:
    """Dimensionless spherical finite-volume grid.

    Volumes and areas omit the common 4 pi factor, matching the baseline
    paper's discretization with V_i = (r_o^3 - r_i^3) / 3 and A = r^2.
    """

    interfaces_code: np.ndarray

    def __post_init__(self) -> None:
        interfaces = np.asarray(self.interfaces_code, dtype=float)
        if interfaces.ndim != 1:
            raise ValueError("interfaces_code must be one-dimensional")
        if len(interfaces) < 2:
            raise ValueError("at least two interfaces are required")
        if np.any(interfaces <= 0.0):
            raise ValueError("interfaces must be positive")
        if np.any(np.diff(interfaces) <= 0.0):
            raise ValueError("interfaces must be strictly increasing")
        object.__setattr__(self, "interfaces_code", interfaces)

    @classmethod
    def from_log_spacing(
        cls,
        r_min_code: float,
        r_max_code: float,
        num_cells: int,
    ) -> "SphericalGrid":
        if r_min_code <= 0.0:
            raise ValueError("r_min_code must be positive")
        if r_max_code <= r_min_code:
            raise ValueError("r_max_code must exceed r_min_code")
        if num_cells <= 0:
            raise ValueError("num_cells must be positive")
        return cls(np.geomspace(r_min_code, r_max_code, num_cells + 1))

    @property
    def num_cells(self) -> int:
        return len(self.interfaces_code) - 1

    @property
    def inner_edges_code(self) -> np.ndarray:
        return self.interfaces_code[:-1]

    @property
    def outer_edges_code(self) -> np.ndarray:
        return self.interfaces_code[1:]

    @property
    def centers_code(self) -> np.ndarray:
        return ((self.inner_edges_code**3 + self.outer_edges_code**3) / 2.0) ** (1.0 / 3.0)

    @property
    def widths_code(self) -> np.ndarray:
        return self.outer_edges_code - self.inner_edges_code

    @property
    def interface_areas_code(self) -> np.ndarray:
        return self.interfaces_code**2

    @property
    def cell_volumes_code(self) -> np.ndarray:
        return (self.outer_edges_code**3 - self.inner_edges_code**3) / 3.0

