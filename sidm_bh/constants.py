"""Physical constants and unit conversions used by the prototype models."""

from __future__ import annotations

G_CGS = 6.67430e-8
C_CGS = 2.99792458e10
M_SUN_CGS = 1.98847e33
PC_CGS = 3.0856775814913673e18
KPC_CGS = 1.0e3 * PC_CGS
MPC_CGS = 1.0e6 * PC_CGS
MYR_CGS = 1.0e6 * 365.25 * 24.0 * 3600.0
SIGMA_T_CGS = 6.6524587321e-25
M_PROTON_CGS = 1.67262192369e-24
K_BOLTZMANN_CGS = 1.380649e-16


def msun_to_g(mass_msun: float) -> float:
    return mass_msun * M_SUN_CGS


def g_to_msun(mass_g: float) -> float:
    return mass_g / M_SUN_CGS


def pc_to_cm(radius_pc: float) -> float:
    return radius_pc * PC_CGS


def cm_to_pc(radius_cm: float) -> float:
    return radius_cm / PC_CGS


def myr_to_s(time_myr: float) -> float:
    return time_myr * MYR_CGS


def s_to_myr(time_s: float) -> float:
    return time_s / MYR_CGS
