# Stage 3 Static Hernquist Baryon Pilot

Date: 2026-07-11

Status: first matched experiment matrix complete

## Matched experiment protocol

- Dark halo: NFW with `rho_s = 3.7 M_sun/pc^3` and `r_s = 30 pc`.
- Nominal halo mass used to normalize baryons: `1e6 M_sun`.
- Black-hole seed: `100 M_sun`.
- SIDM cross section: `50 cm^2/g`.
- Domain: `0.005-5000 pc`, 256 logarithmic cells.
- Evolution: 2 Myr, MC-Roe, CFL `0.2`, entropy fix `0.1`, implicit conduction.
- Static baryons: Hernquist profile with `M_b = f_b M_halo`.

Each case begins with the same NFW dark-matter density. Its pressure is independently integrated to hydrostatic equilibrium in the combined `DM + BH + static baryon` potential. The no-baryon control uses the same protocol without the Hernquist term. Thus the comparison does not contain an impulsive `t=0` switch-on of the baryonic force.

The primary enhancement statistic is

```text
Enhancement_DM = [M_BH,b(2 Myr) - M_seed]
                 / [M_BH,0(2 Myr) - M_seed].
```

The matched no-baryon control reaches `122.12 M_sun`, corresponding to `22.12 M_sun` of accreted dark matter.

## Main 3x3 matrix

| f_b | a_b/r_s | M_BH(2 Myr) [M_sun] | Accreted DM [M_sun] | Enhancement_DM |
|---:|---:|---:|---:|---:|
| 0.01 | 0.001 | 759.67 | 659.67 | 29.82 |
| 0.01 | 0.010 | 371.67 | 271.67 | 12.28 |
| 0.01 | 0.100 | 119.54 | 19.54 | 0.88 |
| 0.05 | 0.001 | 1457.04 | 1357.04 | 61.34 |
| 0.05 | 0.010 | 763.56 | 663.56 | 29.99 |
| 0.05 | 0.100 | 115.57 | 15.57 | 0.70 |
| 0.16 | 0.001 | 2295.58 | 2195.58 | 99.25 |
| 0.16 | 0.010 | 1025.42 | 925.42 | 41.83 |
| 0.16 | 0.100 | 111.15 | 11.15 | 0.50 |

![Static-baryon enhancement matrix](results/stage3/figures/stage3_static_enhancement.png)

## Time dependence

![Static-baryon growth histories](results/stage3/figures/stage3_static_growth.png)

Compact baryons (`a_b/r_s <= 0.01`) strongly enhance dark accretion. The enhancement rises with baryon fraction and reaches approximately `100x` for the most compact `f_b = 0.16` case. In contrast, extended baryons (`a_b/r_s = 0.1`) suppress dark growth, with stronger suppression at larger baryon fraction.

The `a_b/r_s = 0.01` cases still have rising accretion rates at 2 Myr. Their endpoint is therefore a fixed-time comparison, not a saturated final mass.

## Inner-boundary robustness

Representative enhanced and suppressed cases were repeated with matched no-baryon controls at three absorbing radii.

| r_min [pc] | Enhanced case factor | Suppressed case factor |
|---:|---:|---:|
| 0.0025 | 22.27 | 0.77 |
| 0.0050 | 29.99 | 0.70 |
| 0.0100 | 43.54 | 0.65 |

Enhanced case: `f_b = 0.05`, `a_b/r_s = 0.01`. Suppressed case: `f_b = 0.05`, `a_b/r_s = 0.1`.

![Boundary robustness](results/stage3/figures/stage3_boundary_robustness.png)

The absorbing radius changes the numerical magnitude of the enhancement but not its sign. Any quoted enhancement constraint must therefore retain matched boundaries and report a boundary range.

## Interpretation

1. Baryonic concentration is at least as important as total baryonic mass.
2. A deep compact potential can increase gravitational inflow by one to two orders of magnitude.
3. A broad static potential can instead raise the equilibrium thermal support and conductive redistribution enough to reduce central dark supply.
4. The existence of both enhancement and suppression argues against using baryon fraction alone as the stage-3 control parameter.

## Limitations and next tests

- The baryonic potential is static and externally prescribed.
- No baryonic gas is accreted by the black hole.
- No baryonic cooling, feedback, contraction history, or mass assembly is included.
- The initial NFW density is held fixed while its equilibrium pressure changes with the external potential.
- The full 3x3 matrix has only one nominal inner boundary; representative points establish sign robustness but not a complete boundary-marginalized surface.

The equilibrium-static protocol has now been compared against a gradually assembled Hernquist potential. The representative assembled-potential scan and its sign reversal are documented in `stage3_assembly_history_report.md`. A full concentration-assembly-time-boundary matrix, followed by cross-section variation, is large enough to justify Paracloud Slurm execution once credentials are provided through the secure workflow.
