# Stage 3 Baryon-Fraction and Halo-Mass Scaling

Date: 2026-07-12

Status: 60-case predictive scaling matrix and refinement complete

## Controlled protocol

The halo family is anchored to the paper baseline
`M_halo = 1e6 M_sun`, `rho_s = 3.7 M_sun/pc^3`, and `r_s = 30 pc`.
The anchor encloses `1e6 M_sun` at concentration `c = 3.9201`. Halo mass
is varied self-similarly at fixed scale density and concentration:

```text
rho_s = constant
r_s = 30 pc (M_halo / 1e6 M_sun)^(1/3)
r_vir = c r_s
```

This is a controlled fixed-density halo family, not a redshift-dependent
cosmological concentration relation. The black-hole seed remains `100 M_sun`
and the physical inner boundary remains `0.005 pc`. The outer boundary scales
with `r_s`, while the number of cells is increased from 256 to 270 and 284 to
preserve the baseline logarithmic resolution.

All baryon runs use `a_b/r_s = 0.01`, `sigma/m = 50 cm^2/g`, MC-Roe fluxes,
implicit conduction, and a 2 Myr endpoint. The matrix covers
`M_halo = 1e6, 1e7, 1e8 M_sun` and `f_b = 0.01, 0.05, 0.16`, with one matched
no-baryon control per halo mass.

## Pre-run assembly-time predictor

The earlier compact-potential results give future feed radii near
`3-4.5 a_b`. The new predictor fixes this information once, before the new
runs, as

```text
r_supply,pred = 3.5 a_b
T_pred = t_dyn(r_supply,pred, t=0)
```

No radius or timescale from the new evolved solutions is used to set
`T_pred`.

| M_halo [M_sun] | r_s [pc] | sampled r_supply [pc] | T_pred [Myr] | r_max [pc] | cells |
|---:|---:|---:|---:|---:|---:|
| 1e6 | 30.000 | 1.075 | 0.5641 | 5000 | 256 |
| 1e7 | 64.633 | 2.300 | 0.5936 | 10772 | 270 |
| 1e8 | 139.248 | 4.925 | 0.5952 | 23208 | 284 |

The near-constant high-mass value is the expected result for a fixed-density
self-similar halo. The small low-mass offset is caused by the fixed physical
black-hole seed.

The first matrix sampled assembly-time multipliers `0, 0.5, 1, 2`. A second
21-case matrix added `0.25/0.75` for `f_b=0.01`, `0.75/1.25` for `f_b=0.05`,
and `1.5/2.5/3.0` for `f_b=0.16`. All assembly times remain below the 2 Myr
endpoint. Except for the instantaneous low-mass, low-baryon case, every peak
is bracketed on both sides. Continuous peak locations below use a local
quadratic fit to log accreted mass over the three adjacent samples.

## Refined optima

| M_halo [M_sun] | f_b | discrete T/T_pred | fitted T/T_pred | fitted T [Myr] | accreted DM [M_sun] | enhancement |
|---:|---:|---:|---:|---:|---:|---:|
| 1e6 | 0.01 | 0.00 | 0.000 | 0.000 | 506.51 | 22.90 |
| 1e6 | 0.05 | 1.25 | 1.197 | 0.675 | 1276.11 | 57.68 |
| 1e6 | 0.16 | 2.00 | 1.984 | 1.119 | 2005.05 | 90.63 |
| 1e7 | 0.01 | 0.50 | 0.466 | 0.277 | 365.81 | 23.84 |
| 1e7 | 0.05 | 1.25 | 1.242 | 0.737 | 689.31 | 44.92 |
| 1e7 | 0.16 | 2.00 | 2.080 | 1.235 | 1247.21 | 81.27 |
| 1e8 | 0.01 | 0.50 | 0.380 | 0.226 | 168.89 | 9.38 |
| 1e8 | 0.05 | 0.75 | 0.869 | 0.517 | 412.28 | 22.90 |
| 1e8 | 0.16 | 1.50 | 1.586 | 0.944 | 956.43 | 53.13 |

![Scaling matrix](results/stage3/figures/stage3_scaling.png)

## Main finding: the clock depends on baryon strength

The original dynamical predictor succeeds to within 25 percent for all three
`f_b = 0.05` halos, but for none of the `f_b = 0.01` or `0.16` halos. Median
fitted peak multipliers are

```text
f_b = 0.01: T_opt/T_pred = 0.380
f_b = 0.05: T_opt/T_pred = 1.197
f_b = 0.16: T_opt/T_pred = 1.984
```

A least-squares power law through these three medians is

```text
T_opt = 1.06 t_dyn(3.5 a_b) (f_b / 0.05)^0.603.
```

Thus the dark-matter dynamical time sets the basic scale, but it is not a
complete assembly clock. A stronger baryon potential favors slower assembly,
consistent with a longer interval over which a larger supply region can
contract while conduction removes compression heat. This is an empirical
three-point scaling at fixed compactness and cross section, not yet a
universal law.

## Growth-amplitude scaling

At fixed halo mass, optimal accreted mass grows sublinearly with baryon
fraction:

```text
Delta M_BH,opt proportional to f_b^alpha
alpha = 0.501, 0.439, 0.621 for M_halo = 1e6, 1e7, 1e8 M_sun.
```

At fixed baryon fraction, optimal accreted mass decreases weakly with halo
mass, with slopes `-0.24`, `-0.25`, and `-0.16`. The matched controls accrete
`22.12`, `15.35`, and `18.00 M_sun`, giving an almost flat control slope of
`-0.045`.

This negative halo-mass slope is not a cosmological prediction. Holding the
seed mass and physical absorbing boundary fixed makes both quantities smaller
in halo units as mass increases, breaking exact self-similarity. A scaled-seed
control is required before interpreting this trend physically.

## Numerical audit

- Slurm jobs `40734715` and `40734740` completed all array elements with
  exit code `0:0`.
- All 60 expected cases are present; all remote error logs are empty.
- NFW concentration spread across saved metadata is exactly zero.
- The largest absolute mass-budget residual is `1.82e-12` in code units.
- The local test suite contains 63 passing tests.

## Self-similarity closure test

The apparent negative halo-mass trend was separated with two companion
protocols at `f_b=0.05` and `T_asm=1.25 T_pred`:

1. Scale the seed as `M_BH proportional to M_halo` and the inner boundary as
   `r_min proportional to r_s`, while keeping the physical cross section at
   `50 cm^2/g`.
2. Apply the same seed and boundary scaling, and additionally scale the
   physical cross section as `1/r_s` so its code-unit value is constant.

| Protocol | control growth slope | baryon growth slope | enhancement at 1e6/1e7/1e8 M_sun |
|:---|---:|---:|:---|
| Fixed physical seed, boundary, sigma | -0.045 | -0.255 | 57.68 / 44.92 / 21.87 |
| Scaled seed and boundary, fixed physical sigma | 1.068 | 0.744 | 57.68 / 24.70 / 13.02 |
| Fully dimensionless self-similar | 1.000 | 1.000 | 57.68 / 57.68 / 57.68 |

![Self-similarity closure](results/stage3/figures/stage3_similarity.png)

Scaling the seed and boundary reverses the negative absolute growth trend;
therefore that trend was dominated by fixed physical seed and boundary
scales. The remaining sublinear slope at fixed physical cross section is a
transport effect: the code-unit cross section changes with `r_s`. Once that
quantity is also held fixed, accreted mass scales exactly with halo mass. The
fractional spread in `Delta M/M_halo` is `2.46e-14` for baryon runs and
`1.53e-14` for controls, providing a stringent end-to-end similarity check.

Slurm job `40734760` completed all 12 companion cases with `0:0` exit codes
and empty error logs. Its largest mass-budget residual is `5.30e-13` code
units.

## Stage-3 conclusion

The supply-radius correlation has become a genuinely pre-run predictor, but
only after including baryon strength. The resulting picture is now:

1. Halo dynamics supplies the base clock.
2. Baryon strength shifts the preferred fraction of that clock.
3. Conductive heat removal sets the large amplification established by the
   previous heat-off experiment.

The requested Hernquist static/quasi-static baryon-potential scope is now
complete. It includes matched controls, concentration and assembly histories,
cross-section and boundary sensitivity, timescale and heat-flow mechanism
tests, baryon-fraction and halo-mass scaling, and the similarity audit above.
The results should not be generalized to Plummer or gas-core profiles without
new simulations.

The project can now enter stage 4: dynamic baryon accretion and feedback.
