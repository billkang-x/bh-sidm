# Stage 3 HPC Parameter Scan

Date: 2026-07-12

Status: concentration-assembly-cross-section matrix and selected boundary checks complete

## Scope

The scan holds the NFW halo, `f_b = 0.05`, `M_BH,seed = 100 M_sun`, 2 Myr duration, 256 logarithmic cells, MC-Roe solver, implicit conduction, CFL `0.2`, and entropy fix `0.1` fixed. It varies

- `a_b/r_s = 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3`;
- `T_asm = 0, 0.05, 0.2, 0.5, 1.0 Myr`;
- `sigma/m = 10, 30, 50, 100 cm^2/g`.

Four no-baryon controls use the same cross sections. The main matrix therefore contains 160 assembled-baryon runs plus 4 matched controls.

## HPC execution and validation

- Cluster project: private remote workspace (path intentionally omitted from the public release).
- Main Slurm job: `40734527`, packaged as 8 workers to respect the account QOS submission limit.
- Boundary Slurm job: `40734548`, packaged as 5 workers.
- All workers completed with exit code `0:0`; peak memory was approximately `111 MB` per worker.
- The 164 main NPZ files cover all task IDs and match the manifest metadata exactly.
- All histories reach 2 Myr, contain finite states, and have monotonic black-hole mass.
- The largest mass-budget residual is `1.46e-12` in code units.
- The HPC and local `sigma/m = 50` no-baryon control masses agree exactly.

## Matched controls

| sigma/m [cm^2/g] | Accreted DM at 2 Myr [M_sun] |
|---:|---:|
| 10 | 33.4573 |
| 30 | 23.7867 |
| 50 | 22.1226 |
| 100 | 23.7514 |

The no-baryon response is non-monotonic in cross section, with minimum 2 Myr growth near `sigma/m = 50 cm^2/g` for this halo and numerical setup. Enhancement factors must therefore use a cross-section-matched control rather than a single common denominator.

## Three-dimensional response

![HPC enhancement surface](results/stage3/figures/stage3_hpc_enhancement_surface.png)

All 160 assembled-baryon cases exceed their matched no-baryon controls. The weakest case is still enhanced by `1.138`: `a_b/r_s = 0.3`, `T_asm = 1 Myr`, and `sigma/m = 100 cm^2/g`. Thus the static-equilibrium suppression found earlier does not reappear anywhere in the assembled-potential matrix tested here.

The global maximum is

```text
sigma/m = 10 cm^2/g
a_b/r_s = 0.01
T_asm = 0.5 Myr
accreted DM = 3422.10 M_sun
enhancement = 102.28
```

Selected maxima over assembly time are:

| a_b/r_s | sigma=10 | sigma=30 | sigma=50 | sigma=100 |
|---:|---:|---:|---:|---:|
| 0.01 | 102.28 @ 0.5 | 94.75 @ 0.5 | 56.88 @ 0.5 | 26.20 @ 1.0 |
| 0.02 | 54.20 @ 0.5 | 33.21 @ 0.5 | 20.57 @ 0.5 | 11.13 @ 0.5 |
| 0.04 | 14.23 @ 0.5 | 10.60 @ 0.0 | 7.96 @ 0.0 | 5.16 @ 0.0 |
| 0.10 | 3.27 @ 0.0 | 3.45 @ 0.0 | 3.23 @ 0.0 | 2.70 @ 0.0 |
| 0.30 | 1.31 @ 0.0 | 1.30 @ 0.0 | 1.31 @ 0.0 | 1.29 @ 0.0 |

Each entry is `maximum enhancement @ T_asm [Myr]`.

![HPC optimum map](results/stage3/figures/stage3_hpc_optimum_map.png)

## Main findings

1. Concentration remains the dominant control parameter. The response falls from order `100` at `a_b/r_s = 0.01` to order unity at `0.3`.
2. A delayed optimum exists only for compact potentials. At `a_b/r_s <= 0.02`, finite assembly near `0.5 Myr` generally outperforms instantaneous turn-on; for `sigma/m = 100` and `a_b/r_s = 0.01`, the sampled optimum shifts to `1 Myr`.
3. At `a_b/r_s >= 0.07`, slower assembly monotonically reduces 2 Myr growth because the full force acts later and the compact transient advantage is absent.
4. Cross-section dependence is not monotonic. Compact-potential enhancement decreases strongly between `sigma/m = 10` and `100`, while the extended `a_b/r_s = 0.3` maximum remains near `1.3` throughout.

![Cross-section response](results/stage3/figures/stage3_hpc_cross_section_response.png)

The movement of the compact-potential optimum with cross section is consistent with, but does not yet prove, a matching between the baryonic assembly time and an SIDM transport or contraction time.

## Boundary follow-up

| Case | r_min=0.0025 pc | r_min=0.005 pc | r_min=0.01 pc |
|:---|---:|---:|---:|
| Global maximum: `sigma=10, a/r_s=0.01, T=0.5` | 99.76 | 102.28 | 103.38 |
| High-sigma compact: `100, 0.01, 1.0` | 20.61 | 26.20 | 36.43 |
| Weak extended: `100, 0.3, 1.0` | 1.1388 | 1.1380 | 1.1413 |

![HPC boundary follow-up](results/stage3/figures/stage3_hpc_boundary_followup.png)

The global maximum and weak extended enhancement are highly stable under the tested boundary variation. The high-cross-section compact normalization remains boundary-sensitive, but its strong-enhancement sign is unambiguous.

## Limitations and next test

- Only `f_b = 0.05` and one NFW halo are covered by the three-dimensional scan.
- The smoothstep potential is controlled and reproducible but not a cosmological baryon assembly history.
- The 2 Myr endpoint is not asymptotic for the strongest compact cases.
- Baryons remain an external potential without gas dynamics, cooling, feedback, or black-hole baryonic accretion.

The local timescale and heat-off diagnostic is now complete and documented in `stage3_timescale_mechanism_report.md`. The compact optimum tracks the initial dynamical time at the future supply radius, while outward conduction removes compressional heat and amplifies the response. The next test should examine whether this scaling persists across baryon fraction and halo mass.
