# Stage 2 Baseline Validation Report

Date: 2026-07-11

Status: complete

## Scope and acceptance criteria

Stage 2 required:

1. A time-dependent 1D spherical SIDM fluid solver with self-gravity, black-hole accretion, and thermal conduction.
2. SIS and NFW 2 Myr baseline reproduction.
3. Full heat-flow on/off comparisons.
4. Grid, inner-boundary, outer-boundary, timestep, and Riemann entropy-fix sensitivity checks.
5. Reproducible histories and radial profiles.

All five criteria are now satisfied. The solver uses nonuniform-grid MC reconstruction, Roe fluxes, explicit Euler source integration, and implicit frozen-coefficient conduction. Positivity fallback was disabled in the validation runs.

## Baseline endpoint reproduction

| Profile | M_BH(0) [M_sun] | M_BH(2 Myr) [M_sun] | Paper endpoint |
|---|---:|---:|---:|
| NFW, heat flow | 100 | 129.08 | approximately 130 |
| SIS, heat flow | 100 | 9646.8 | approximately 1e4 |

The mass histories, accretion histories, and 201 full radial snapshots are stored in `results/stage2/`.

## Heat-flow comparison

| Profile | Heat flow | M_BH(2 Myr) [M_sun] | Final rate [M_sun/Myr] | Inner density [M_sun/pc^3] |
|---|---|---:|---:|---:|
| NFW | On | 129.08 | 2.84 | 744.9 |
| NFW | Off | 234.36 | 57.00 | 19123.7 |
| SIS | On | 9646.8 | 155.7 | 49350 |
| SIS | Off | 32528.3 | 15886 | 5.15e6 |

Turning off conduction raises the final black-hole mass by factors of `1.82` for NFW and `3.37` for SIS. The no-conduction inner density is higher by factors of approximately `25.7` and `104`, respectively. This reproduces the physical conclusion of Figures 3 and 5: conductive heat transport depletes and heats the central supply, strongly suppressing late accretion.

![Growth and accretion comparison](results/stage2/figures/stage2_growth_comparison.png)

## Six-snapshot radial evolution

The profile plots show `0.0`, `0.4`, `0.8`, `1.2`, `1.6`, and `2.0 Myr` heat-flow states, plus the final no-heat state.

![NFW radial profiles](results/stage2/figures/stage2_nfw_profiles.png)

![SIS radial profiles](results/stage2/figures/stage2_sis_profiles.png)

## Outer-boundary sensitivity

The number of cells was adjusted to preserve approximately constant logarithmic resolution as `r_max` changed.

| Profile | Fixed r_min [pc] | r_max values [pc] | M_BH range [M_sun] | Relative span |
|---|---:|---|---:|---:|
| NFW | 0.005 | 2500, 5000, 10000 | 129.07551-129.07573 | 1.7e-6 |
| SIS | 0.005 | 500, 1000, 2000 | 11698.32-11699.98 | 1.4e-4 |

The outer boundary is negligible at the tested locations. This contrasts with the inner absorbing boundary, which changes the NFW endpoint by approximately `4%` for a factor-two shift and changes the SIS normalization more strongly.

## Other sensitivity results

- NFW 256-to-512 grid change: `0.05%` in final black-hole mass.
- SIS 128-to-256 grid change at fixed `r_min = 0.005 pc`: `1.1%`.
- NFW CFL `0.1-0.4`: below `0.001%` endpoint variation.
- SIS CFL `0.2-0.8`: approximately `0.012%` endpoint variation.
- Roe entropy fix `0-0.2`: below `0.001%` for NFW and approximately `0.02%` for SIS.
- Relative mass-budget residuals: approximately `1e-15` to `1e-13`.

Detailed sensitivity tables are in `baseline_2myr_results.md` and `convergence_results.md`.

## Remaining limitations

1. The inner absorbing boundary is the dominant numerical systematic, especially for SIS.
2. Primitive-variable MC reconstruction and the Harten entropy fix are deliberate robustness extensions beyond a literal reading of Appendix A.
3. Spherical symmetry, absent angular momentum, and the effective LMFP closure remain physical limitations rather than stage-2 implementation gaps.

## Stage decision

Stage 2 is complete. Stage 3 may begin with static Hernquist baryonic potentials, provided every baryonic run uses a matched no-baryon control with identical grid, boundaries, and numerical parameters. At least two inner-boundary values must be retained when quoting a baryonic enhancement factor.
