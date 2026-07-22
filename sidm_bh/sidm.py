"""SIDM transport closures used by the accretion prototype."""

from __future__ import annotations

import math

import numpy as np

from .constants import G_CGS

CONDUCTIVITY_A = math.sqrt(16.0 / math.pi)
CONDUCTIVITY_B = 25.0 * math.sqrt(math.pi) / 32.0
CONDUCTIVITY_C = 0.75


def rutherford_momentum_transfer_ratio(
    relative_speed_over_transition: float | np.ndarray,
) -> float | np.ndarray:
    """Return the Rutherford transfer cross section divided by ``sigma_0``.

    The differential cross section is proportional to
    ``[1 + (v/w)^2 sin(theta/2)^2]^-2`` and the result is normalized to one
    in the low-relative-speed limit.
    """

    speed_ratio = np.asarray(relative_speed_over_transition, dtype=float)
    if np.any(speed_ratio < 0.0):
        raise ValueError("relative speed ratio cannot be negative")
    y = speed_ratio**2
    result = np.empty_like(y)
    small = y < 1.0e-4
    # Avoid cancellation between log1p(y) and y / (1 + y).
    result[small] = (
        1.0
        - (4.0 / 3.0) * y[small]
        + 1.5 * y[small] ** 2
        - 1.6 * y[small] ** 3
    )
    regular = ~small
    result[regular] = 2.0 * (
        np.log1p(y[regular]) - y[regular] / (1.0 + y[regular])
    ) / y[regular] ** 2
    if result.ndim == 0:
        return float(result)
    return result


def rutherford_viscosity_cross_section_ratio(
    relative_speed_over_transition: float | np.ndarray,
) -> float | np.ndarray:
    """Return the viscosity cross section normalized to its isotropic limit."""

    speed_ratio = np.asarray(relative_speed_over_transition, dtype=float)
    if np.any(speed_ratio < 0.0):
        raise ValueError("relative speed ratio cannot be negative")
    y = speed_ratio**2
    result = np.empty_like(y)
    small = y < 1.0e-3
    result[small] = (
        1.0
        - y[small]
        + 0.9 * y[small] ** 2
        - 0.8 * y[small] ** 3
        + (5.0 / 7.0) * y[small] ** 4
    )
    regular = ~small
    result[regular] = 6.0 * (
        (y[regular] + 2.0) * np.log1p(y[regular]) - 2.0 * y[regular]
    ) / y[regular] ** 3
    if result.ndim == 0:
        return float(result)
    return result


def maxwellian_viscosity_cross_section_ratio(
    dispersion_over_transition: float | np.ndarray,
    velocity_power: int,
    quadrature_order: int = 384,
) -> float | np.ndarray:
    """Return the normalized ``<sigma_visc v_rel^p>`` thermal average.

    The input is the local one-dimensional velocity dispersion divided by
    the Rutherford transition speed ``w``.  The logarithmic integration
    variable resolves the low-relative-speed tail even when the local
    dispersion is much larger than ``w``.  The conduction closure uses
    ``p=3`` in the LMFP term and ``p=5`` in the SMFP term.
    """

    ratio = np.asarray(dispersion_over_transition, dtype=float)
    if np.any(ratio < 0.0):
        raise ValueError("dispersion ratio cannot be negative")
    if velocity_power not in (3, 5):
        raise ValueError("velocity_power must be three or five")
    if quadrature_order < 8:
        raise ValueError("quadrature_order must be at least eight")
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    lower_log_t = -35.0
    upper_log_t = np.log(50.0)
    log_t = (
        0.5 * (upper_log_t - lower_log_t) * nodes
        + 0.5 * (upper_log_t + lower_log_t)
    )
    log_weights = 0.5 * (upper_log_t - lower_log_t) * weights
    t = np.exp(log_t)
    relative_speed_ratio = 2.0 * ratio[..., np.newaxis] * np.sqrt(t)
    viscosity_ratio = rutherford_viscosity_cross_section_ratio(
        relative_speed_ratio
    )
    t_power = 0.5 * (velocity_power + 3.0)
    normalization = math.gamma(0.5 * (velocity_power + 3.0))
    effective = np.sum(
        log_weights * t**t_power * np.exp(-t) * viscosity_ratio,
        axis=-1,
    ) / normalization
    if effective.ndim == 0:
        return float(effective)
    return effective


def effective_cross_section_ratio_tables(
    log10_dispersion_ratio_min: float = -6.0,
    log10_dispersion_ratio_max: float = 4.0,
    num_points: int = 4097,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tabulate the LMFP ``K3`` and SMFP ``K5`` Rutherford averages."""

    if log10_dispersion_ratio_max <= log10_dispersion_ratio_min:
        raise ValueError("table upper bound must exceed lower bound")
    if num_points < 2:
        raise ValueError("num_points must be at least two")
    log_ratio = np.linspace(
        log10_dispersion_ratio_min,
        log10_dispersion_ratio_max,
        num_points,
    )
    dispersion_ratio = 10.0**log_ratio
    k3 = maxwellian_viscosity_cross_section_ratio(
        dispersion_ratio, velocity_power=3
    )
    k5 = maxwellian_viscosity_cross_section_ratio(
        dispersion_ratio, velocity_power=5
    )
    return (
        log_ratio,
        np.asarray(k3, dtype=float),
        np.asarray(k5, dtype=float),
    )


def conductivity_cgs(
    density_cgs: float,
    velocity_dispersion_cms: float,
    sigma_over_m_cm2_g: float,
    calibration_c: float = CONDUCTIVITY_C,
) -> float:
    """Return SIDM effective conductivity in cgs units.

    This implements the SMFP/LMFP interpolation used by the baseline paper.
    The velocity dispersion is the one-dimensional dispersion v.
    """

    if density_cgs <= 0.0:
        raise ValueError("density_cgs must be positive")
    if velocity_dispersion_cms <= 0.0:
        raise ValueError("velocity_dispersion_cms must be positive")
    if sigma_over_m_cm2_g <= 0.0:
        raise ValueError("sigma_over_m_cm2_g must be positive")
    if calibration_c <= 0.0:
        raise ValueError("calibration_c must be positive")

    smfp_inverse = sigma_over_m_cm2_g / (CONDUCTIVITY_B * velocity_dispersion_cms)
    lmfp_inverse = (
        4.0
        * math.pi
        * G_CGS
        / (
            CONDUCTIVITY_A
            * calibration_c
            * density_cgs
            * velocity_dispersion_cms**3
            * sigma_over_m_cm2_g
        )
    )
    return 1.5 / (smfp_inverse + lmfp_inverse)


def conductivity_code(
    density_code: float,
    velocity_dispersion_code: float,
    sigma_over_m_code: float,
    calibration_c: float = CONDUCTIVITY_C,
) -> float:
    """Return dimensionless SIDM conductivity from the paper's Eq. (9)."""

    if density_code <= 0.0:
        raise ValueError("density_code must be positive")
    if velocity_dispersion_code <= 0.0:
        raise ValueError("velocity_dispersion_code must be positive")
    if sigma_over_m_code <= 0.0:
        raise ValueError("sigma_over_m_code must be positive")
    if calibration_c <= 0.0:
        raise ValueError("calibration_c must be positive")

    smfp_inverse = sigma_over_m_code / (
        CONDUCTIVITY_B * velocity_dispersion_code
    )
    lmfp_inverse = 1.0 / (
        CONDUCTIVITY_A
        * calibration_c
        * density_code
        * velocity_dispersion_code**3
        * sigma_over_m_code
    )
    return 1.5 / (smfp_inverse + lmfp_inverse)


def mean_free_path_cgs(density_cgs: float, sigma_over_m_cm2_g: float) -> float:
    if density_cgs <= 0.0:
        raise ValueError("density_cgs must be positive")
    if sigma_over_m_cm2_g <= 0.0:
        raise ValueError("sigma_over_m_cm2_g must be positive")
    return 1.0 / (density_cgs * sigma_over_m_cm2_g)


def gravitational_scale_height_cgs(
    density_cgs: float,
    velocity_dispersion_cms: float,
) -> float:
    if density_cgs <= 0.0:
        raise ValueError("density_cgs must be positive")
    if velocity_dispersion_cms <= 0.0:
        raise ValueError("velocity_dispersion_cms must be positive")
    return math.sqrt(velocity_dispersion_cms**2 / (4.0 * math.pi * G_CGS * density_cgs))


def knudsen_number(
    density_cgs: float,
    velocity_dispersion_cms: float,
    sigma_over_m_cm2_g: float,
) -> float:
    return mean_free_path_cgs(
        density_cgs,
        sigma_over_m_cm2_g,
    ) / gravitational_scale_height_cgs(density_cgs, velocity_dispersion_cms)
