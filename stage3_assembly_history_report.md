# Stage 3 Hernquist Assembly-History Experiment

Date: 2026-07-11

Status: representative assembly-time scan complete

## Question and protocol

The static-equilibrium matrix showed that compact Hernquist potentials enhance SIDM accretion while an extended potential can suppress it. This experiment tests whether that conclusion survives a physical mass-assembly history.

All assembly runs start from the same NFW dark-matter density and the same hydrostatic pressure in the `DM + 100 M_sun BH` potential, with no baryonic force initially. The final Hernquist enclosed-mass profile is then multiplied by

```text
x = t / T_asm
f_asm = 3 x^2 - 2 x^3       for 0 < x < 1,
f_asm = 0                    at t <= 0,
f_asm = 1                    at t >= T_asm.
```

The smoothstep law has zero slope at both endpoints and reaches exactly the same final potential in every finite-assembly run. `T_asm = 0` is the instantaneous turn-on limit. All other numerical settings match the 2 Myr MC-Roe static experiment: 256 logarithmic cells, `r_min = 0.005 pc`, CFL `0.2`, entropy fix `0.1`, implicit conduction, and `sigma/m = 50 cm^2/g`.

## Main results

The matched no-baryon control accretes `22.12 M_sun` of dark matter in 2 Myr.

| a_b/r_s | Protocol | T_asm [Myr] | Accreted DM [M_sun] | Control enhancement | Ratio to static equilibrium |
|---:|:---|---:|---:|---:|---:|
| 0.01 | static equilibrium | - | 663.56 | 29.99 | 1.00 |
| 0.01 | instant turn-on | 0.00 | 981.00 | 44.34 | 1.48 |
| 0.01 | smooth assembly | 0.05 | 988.01 | 44.66 | 1.49 |
| 0.01 | smooth assembly | 0.20 | 1108.50 | 50.11 | 1.67 |
| 0.01 | smooth assembly | 0.50 | 1258.23 | 56.88 | 1.90 |
| 0.01 | smooth assembly | 1.00 | 1219.38 | 55.12 | 1.84 |
| 0.10 | static equilibrium | - | 15.57 | 0.70 | 1.00 |
| 0.10 | instant turn-on | 0.00 | 71.54 | 3.23 | 4.60 |
| 0.10 | smooth assembly | 0.05 | 69.30 | 3.13 | 4.45 |
| 0.10 | smooth assembly | 0.20 | 64.59 | 2.92 | 4.15 |
| 0.10 | smooth assembly | 0.50 | 57.25 | 2.59 | 3.68 |
| 0.10 | smooth assembly | 1.00 | 47.88 | 2.16 | 3.08 |

![Assembly-time response](results/stage3/figures/stage3_assembly_time_response.png)

## Findings

1. The compact case has a non-monotonic history response. Its largest 2 Myr growth occurs near `T_asm = 0.5 Myr`, not in the instantaneous limit. It accretes about `90%` more dark matter than the full-potential static-equilibrium run.
2. The extended case changes sign. The static-equilibrium protocol suppresses growth to `0.70` of control, whereas every assembled-potential run enhances it by `2.16-3.23`.
3. The sign reversal is caused by the initial thermodynamic protocol, not by a different final Hernquist mass. Static equilibrium preloads the SIDM with the pressure needed to support the external potential. Assembly starts from the lower-pressure no-baryon state and drives contraction and inward motion as the potential deepens.
4. At `0.1 pc` in the extended case, the final density is `6.30` times the no-baryon value for static equilibrium and `15.89` times for `T_asm = 0.5 Myr`. The corresponding radial velocities are `-0.212` and `-0.280 km/s`.

![Growth histories](results/stage3/figures/stage3_assembly_growth.png)

![Final radial response](results/stage3/figures/stage3_assembly_radial_response.png)

## Inner-boundary check

The representative `T_asm = 0.5 Myr` runs were repeated at all three established inner boundaries.

| a_b/r_s | r_min [pc] | Static-equilibrium enhancement | Assembly enhancement |
|---:|---:|---:|---:|
| 0.01 | 0.0025 | 22.27 | 42.58 |
| 0.01 | 0.0050 | 29.99 | 56.88 |
| 0.01 | 0.0100 | 43.54 | 77.66 |
| 0.10 | 0.0025 | 0.77 | 2.45 |
| 0.10 | 0.0050 | 0.70 | 2.59 |
| 0.10 | 0.0100 | 0.65 | 2.80 |

The absorbing radius remains the main normalization uncertainty, but neither the stronger compact response nor the extended-case sign reversal changes across the tested boundary range.

![Assembly boundary sensitivity](results/stage3/figures/stage3_assembly_boundary.png)

## Numerical validation and limitations

- The time-dependent fast solver is cross-validated against the reference MC-Roe implementation at the full conservative-state level.
- Completed runs have relative mass-budget residuals at the established `10^-15-10^-13` scale.
- The baryons remain an externally prescribed collisionless potential. There is no baryonic gas, cooling, feedback, star formation, or baryonic black-hole accretion.
- The smoothstep law is a controlled history experiment rather than a cosmological assembly model.
- Compact cases still have rising accretion rates at 2 Myr, so the reported maximum is a fixed-time response rather than an asymptotic optimum.

## Next decision

The concentration-assembly-time-cross-section surface has now been run on Paracloud and is documented in `stage3_hpc_parameter_scan_report.md`. The next calculation should diagnose local dynamical and conductive times to test the timescale-matching interpretation before expanding in baryon fraction and halo mass.
