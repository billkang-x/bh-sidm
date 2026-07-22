# MC-Roe 2 Myr Baseline and Sensitivity Results

Date: 2026-07-11

## Method

- Hydrodynamics: nonuniform-grid primitive MC reconstruction and Roe flux.
- Thermal transport: frozen-coefficient implicit SIDM conduction.
- Source integration: explicit Euler, matching the baseline paper.
- Positivity fallback: disabled.
- Black-hole growth: inner-boundary SIDM mass flux, updated every accepted step.
- Long runs: Numba kernel cross-validated against the reference implementation to relative field differences below `2e-14`.

The NFW natural scales are `r0 = 30 pc` and `rho0 = 3.7 M_sun/pc^3`. The SIS scales use `r0 = 1 pc` and `rho0 = rho_SIS(r0)`, giving `v0 = sqrt(2) c_s`; the initial SIS velocity dispersion is exactly constant at `c_s = 4.2 km/s`.

## Full baseline histories

### NFW

Parameters: `r_min = 0.005 pc`, `r_max = 5000 pc`, `M_BH(0) = 100 M_sun`, `sigma/m = 50 cm^2/g`, 256 cells, CFL `0.2`, entropy fix `0.1`.

| Time [Myr] | M_BH [M_sun] |
|---:|---:|
| 0.0 | 100.0000 |
| 0.4 | 119.4806 |
| 0.8 | 123.9100 |
| 1.2 | 126.2245 |
| 1.6 | 127.8232 |
| 2.0 | 129.0756 |

The peak and final accretion rates are `124.27` and `2.840 M_sun/Myr`. A 512-cell run gives `M_BH(2 Myr) = 129.1406 M_sun`, a `0.050%` increase. This reproduces the approximately `130 M_sun` endpoint in Figure 4.

### SIS

Parameters: `r_min = 0.001 pc`, `r_max = 1000 pc`, `M_BH(0) = 100 M_sun`, `sigma/m = 50 cm^2/g`, 128 cells, CFL `0.8`, entropy fix `0.1`.

| Time [Myr] | M_BH [M_sun] | Accretion rate [M_sun/Myr] |
|---:|---:|---:|
| 0.0 | 100.0000 | 0.000 |
| 0.4 | 5790.7034 | 10472.24 |
| 0.8 | 8904.6524 | 4895.79 |
| 1.2 | 9502.9103 | 222.14 |
| 1.6 | 9580.7448 | 176.13 |
| 2.0 | 9646.7746 | 155.69 |

The transient peak accretion rate is `8.47e4 M_sun/Myr`. The final mass reproduces the approximately `1e4 M_sun` endpoint in Figure 2. The nominal run required about `6.07e6` timesteps; its relative mass-budget residual was `5.3e-14`.

## NFW sensitivity

All entries are 256-cell, CFL `0.2`, entropy-fix `0.1` runs unless the varied column states otherwise.

| Variation | Value | M_BH(2 Myr) [M_sun] | Peak rate [M_sun/Myr] | Final rate [M_sun/Myr] |
|---|---:|---:|---:|---:|
| Baseline | - | 129.0756 | 124.272 | 2.840 |
| Grid cells | 512 | 129.1406 | 124.390 | 2.846 |
| Inner boundary [pc] | 0.0025 | 123.7876 | 119.852 | 2.539 |
| Inner boundary [pc] | 0.0100 | 134.3634 | 130.591 | 3.166 |
| CFL | 0.1 | 129.0758 | 124.269 | 2.840 |
| CFL | 0.4 | 129.0752 | 124.277 | 2.840 |
| Entropy fix | 0.0 | 129.0756 | 124.272 | 2.840 |
| Entropy fix | 0.2 | 129.0756 | 124.272 | 2.840 |

Changing the inner boundary by a factor of two changes the final NFW mass by approximately `-4.1%/+4.1%`. CFL and entropy-fix effects are below `0.001%` for the final mass.

## SIS sensitivity

### Inner boundary

These runs use 128 cells, CFL `0.8`, and entropy fix `0.1`.

| r_min [pc] | M_BH(2 Myr) [M_sun] | Peak rate [M_sun/Myr] | Final rate [M_sun/Myr] |
|---:|---:|---:|---:|
| 0.001 | 9646.8 | 84674 | 155.7 |
| 0.002 | 10304.9 | 53580 | 266.7 |
| 0.005 | 11699.1 | 30037 | 577.2 |
| 0.010 | 13292.3 | 18453 | 1103.9 |

The SIS endpoint remains of order `1e4 M_sun`, but its normalization is strongly controlled by the absorbing radius. The transient peak rate is even more boundary-sensitive and should not be treated as a robust physical observable.

### Grid, CFL, and entropy fix

For computational tractability these checks use `r_min = 0.005 pc`.

| Variation | Value | M_BH(2 Myr) [M_sun] |
|---|---:|---:|
| Grid cells | 128 | 11699.1 |
| Grid cells | 256 | 11827.2 |
| CFL | 0.2 | 11700.5 |
| CFL | 0.4 | 11700.0 |
| CFL | 0.8 | 11699.1 |
| Entropy fix | 0.0 | 11698.7 |
| Entropy fix | 0.1 | 11699.1 |
| Entropy fix | 0.2 | 11701.1 |

At fixed boundary, the 128-to-256 grid change is about `1.1%`; CFL changes the final mass by about `0.012%`, and the entropy-fix range changes it by about `0.02%`.

## Conclusions

1. MC-Roe reproduces both paper endpoints: approximately `130 M_sun` for NFW and approximately `1e4 M_sun` for SIS after 2 Myr.
2. CFL and Roe entropy-fix choices are negligible over the tested ranges.
3. Spatial convergence is adequate at 256 cells for NFW and at the percent level for the tested SIS boundary.
4. The inner absorbing boundary is the dominant systematic uncertainty, particularly for SIS. Baryonic enhancement claims must therefore be reported at multiple `r_min` values or normalized against a matched no-baryon control using the same boundary.
