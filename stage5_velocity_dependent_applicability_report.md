# Stage 5 Velocity-dependent SIDM and Applicability Map

## Scope and status

This extension tests whether the trusted constant-cross-section result
survives a physically motivated velocity-dependent SIDM transport closure,
then maps its domain of applicability across halo mass, redshift,
concentration, seed mass, baryon compactness, and assembly time.

The velocity-dependent implementation, 46-case calibration, 368-case
controlled applicability map, and 80-case targeted boundary/grid refinement
are complete. The final results below use threshold crossing as the robust
classifier; they do not interpret the sampled success fractions as a cosmic
population probability.

## Transport closure

The differential cross section is the Born/Rutherford form

```text
d sigma/d Omega = sigma0/(4 pi)
                  [1 + (v_rel/w)^2 sin^2(theta/2)]^-2.
```

The implementation follows the separate thermal averages advocated by
Outmezguine et al. (2022, arXiv:2204.06568):

```text
K_p = <sigma_visc v_rel^p> / <sigma_visc v_rel^p>_(w -> infinity)
kappa_SMFP = (3/2) B v / (sigma0 K5)
kappa_LMFP = (3/2) A C rho v^3 (sigma0 K3).
```

The full conductivity retains the baseline harmonic LMFP/SMFP interpolation.
Both `K3` and `K5` are integrated over a local Maxwell relative-speed
distribution. A logarithmic quadrature resolves the low-speed tail when
`v/w` is large, and a 4097-point logarithmic lookup table is interpolated
inside the Numba evolution kernel.

This is stronger than substituting a single cross section evaluated at the
virial velocity. It is still an effective fluid closure: the LMFP `K3`
weighting is an ansatz calibrated in the gravothermal literature and has not
been fully established by anisotropic N-body simulations.

## Verification

- The analytic Rutherford momentum-transfer and normalized viscosity cross
  sections recover unity in the low-speed limit and decrease monotonically.
- Both Maxwell averages recover `K3=K5=1` as `w -> infinity`.
- Lookup interpolation agrees with direct quadrature to `2e-4` relative
  tolerance over the tested range.
- A very large transition speed reproduces the constant-cross-section MC-Roe
  solution to `2e-10` relative tolerance.
- The full local suite passes 104/104 tests; the remote Python 3.12 environment
  passes all 12 targeted conduction and velocity-dependent tests.

## Pre-declared calibration

The physical anchor is the trusted stage-5 model:

```text
M200=1e9 M_sun, z=30, c=8, M_seed=1e5 M_sun
f_b=0.16, a_b/r_s=0.020, T_asm=1.25 Myr
physical cooling, evolving Bondi ambient, mixed feedback
influence-gated dark capture, t_end=2 Myr.
```

The matrix contains constant controls at `sigma/m=1,10,100 cm2/g` and a
Rutherford grid with `sigma0/m=1,3,10,30,100 cm2/g` and
`w=10,30,100,300 km/s`. Every model was run at feeding radii equal to
`0.208` and `0.104` times the reference seed influence radius. All 46 tasks
in job `40782036` completed with `0:0` exits and empty error logs; the maximum
mass-budget residual is `1.92e-11` in code units.

The three velocity-dependent models pre-declared for the applicability map
hold `sigma0/m=30 cm2/g` fixed and vary `w`:

| Model | `w` [km/s] | Virial LMFP `sigma0 K3/m` | Initial inner LMFP `sigma0 K3/m` | Smaller-boundary final mass [M_sun] | Crosses `1e7` at both boundaries |
|---|---:|---:|---:|---:|---|
| Low transport | 10 | 0.0265 | 0.000582 | `1.816e6` | no |
| Matched transport | 30 | 0.6998 | 0.0248 | `4.901e6` | no |
| High transport | 100 | 8.768 | 0.875 | `1.065e7` | yes, at 1.86 Myr on the smaller boundary |
| Constant control | n/a | 1.0 | 1.0 | `1.560e7` | yes, at 1.65 Myr on the smaller boundary |

The result is non-monotonic. Only the constant control and the
`sigma0/m=30 cm2/g, w=100 km/s` model cross `1e7 M_sun` at both completed
boundaries. The velocity-dependent model falls from an initial inner LMFP
effective cross section of `0.875` to `0.0689 cm2/g` by 2 Myr, so its success
cannot be represented by one fixed effective cross section.

Matching at the virial velocity is not sufficient. The `w=30 km/s` model has
a virial-scale LMFP value near unity but an initial inner value of only
`0.0248 cm2/g` and fails the target. The relevant predictor must include the
radial and temporal `K3/K5` history through baryonic assembly.

None of the 23 calibration models satisfies the pre-declared 5% terminal-mass
boundary criterion between the two tested radii. The constant control differs
by `7.71%`; the low, matched, and high velocity-dependent brackets differ by
`36.9%`, `16.8%`, and `22.2%`. Accordingly, target crossing is retained as the
scientific classifier, while exact terminal masses remain provisional.

## Controlled applicability design

Four microphysics models are evaluated independently. Each contains 28
one-axis points about the trusted anchor and 64 fixed-seed Latin-hypercube
points, for 368 screen calculations in total.

| Parameter | Controlled range |
|---|---|
| `M200` | `1e8-1e10 M_sun` |
| Redshift | `15-30` |
| Concentration | `4-12` |
| Seed mass | `1e2-1e5 M_sun` |
| `a_b/r_s` | `0.005-0.08` |
| `T_asm` | `0.25-1.75 Myr` |

The analysis also records `T_asm/t_dyn` at a consistently defined initial
matching radius. After the screen, the 16 points closest to the target in
each microphysics model are selected without manual intervention. All 64
receive a smaller-boundary run and the four closest per model also receive a
double-grid run, producing 80 refinement tasks. Final classes are
`robust_success`, `robust_failure`, and `boundary_ambiguous`.

The first applicability chunk is job `40782909`; coordinator job `40783158`
chains all remaining chunks, map analysis, refinement, and final numerical
classification while respecting the account's 50-element submission limit.

## Final six-parameter map

All 368 screen tasks and all 80 targeted refinement tasks completed. The
maximum mass-budget residual is `1.96e-10` in the screen and `4.54e-11` in the
refinement. Two screen points reached the original step cap; rerunning them
with `8e8` steps produced valid NPZ files before the final analysis.

The screen counts are:

| Microphysics model | Screen points | Screen crossings | Maximum final mass [M_sun] |
|---|---:|---:|---:|
| Constant `sigma/m=1` | 92 | 14 | `2.43e7` |
| Rutherford, low transport (`30,10`) | 92 | 0 | `5.54e6` |
| Rutherford, matched (`30,30`) | 92 | 2 | `1.88e7` |
| Rutherford, high transport (`30,100`) | 92 | 16 | `6.03e7` |

The screen crossings are not all physically accepted. After smaller-inner-
boundary and double-grid checks, the audited points divide as follows:

| Microphysics model | Audited points | Robust failures | Boundary ambiguous | Robust successes |
|---|---:|---:|---:|---:|
| Constant `sigma/m=1` | 16 | 8 | 2 | 6 |
| Rutherford, low transport | 16 | 16 | 0 | 0 |
| Rutherford, matched | 16 | 15 | 1 | 0 |
| Rutherford, high transport | 16 | 7 | 7 | 2 |

The robust success points are concentrated in the high-redshift,
intermediate-to-heavy-seed part of the map. In the constant control, accepted
examples include `z=25` at the anchor seed, `a_b/r_s=0.01`, and assembly times
`0.5-1.75 Myr`, plus one Latin-hypercube point with
`M_seed=7.62e4 M_sun`. In the high-transport velocity-dependent model, the
two accepted points are the anchor halo with `M_seed=1e5 M_sun`,
`a_b/r_s=0.02`, and assembly times `0.5` and `1.5 Myr`.

The most important negative result is stronger than the calibration alone:
the low-transport model has no screen crossing anywhere in the tested
six-dimensional domain, and the matched model has no boundary-robust crossing
in its audited neighborhood. Thus a velocity-dependent cross section can
remove the rapid-growth channel even when its nominal low-speed normalization
is large.

None of the 64 audited points passes the 5% terminal-mass criterion in the
full boundary-plus-grid sense. One low-transport failure happens to pass the
terminal-mass criterion, which does not change its threshold classification.
The scientifically stable observable from this map is therefore the
existence of a (10^7 M_sun) crossing, not the terminal mass to several
significant figures.

The descriptive logistic fits on the 64 Latin-hypercube points give
five-fold AUC `0.75` for the constant control and `0.885` for the high-
transport model. These values are useful for ranking sensitivities, but the
success class is sparse (2 and 3 positive LHS points respectively), so they
are not population-level predictive models. Across this controlled map the
seed mass and redshift are the clearest positive directions; assembly timing
and baryon compactness select the width of the successful window.

The machine-readable outputs are:

- `results/stage5/applicability_map_summary.csv`
- `results/stage5/applicability_map_statistics.json`
- `results/stage5/applicability_refinement_summary.csv`
- `results/stage5/applicability_refinement_statistics.json`

The final figures are [stage5_applicability_main_effects.png](results/stage5/figures/stage5_applicability_main_effects.png), [stage5_applicability_joint_map.png](results/stage5/figures/stage5_applicability_joint_map.png), and [stage5_applicability_refinement.png](results/stage5/figures/stage5_applicability_refinement.png).

## Updated conclusion

The velocity-dependent extension does not overturn the heavy-seed result, but
it narrows its physical domain. A high-transport Rutherford model can retain
a small, high-redshift, intermediate/heavy-seed success region, whereas
lower-transition-speed models lose it because the inner `K3/K5` transport
history is suppressed. The project can therefore claim a conditional,
seed-selective SIDM growth channel, not a generic velocity-dependent route
from light seeds to LRD masses.

![Velocity-dependent calibration](results/stage5/figures/stage5_velocity_calibration.png)
