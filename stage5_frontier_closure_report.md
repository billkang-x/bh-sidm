# Stage 5 Frontier Closure and Numerical Convergence

## Scope

This study closes the local stage-5 frontier around the physical-cooling,
evolving-Bondi, mixed-feedback model.  The common physical anchor is

```text
M200 = 1e9 M_sun, z = 30, M_seed = 1e5 M_sun
f_b = 0.16, T_asm = 1.25 Myr
rho_b,inf(0) = 0.3 M_sun pc^-3, c_s(0) = 10 km s^-1
epsilon_f = 1e-3, chi_heat = 0.5
Cloudy cooling, MC reconstruction, Roe flux, t_end = 2 Myr.
```

The calculation has four parts: 26 one-axis closure cases, 72 interaction
cases, 40 completed numerical-sensitivity cases around three candidate
models, and 15 targeted cases that close the moderate-concentration and
baryonic-compactness branches.  Physical inner and outer radii are held fixed
when concentration is varied.

## Parameter closure

The one-axis cross-section scan at `c=8` and `a_b/r_s=0.003` has an internal
maximum at `sigma/m=1 cm2 g^-1`: the final black-hole mass is
`8.2002e6 M_sun`, compared with `8.0875e6 M_sun` at `3 cm2 g^-1` and
`6.9958e6 M_sun` at `10 cm2 g^-1`.  The scan covers
`sigma/m=0.01-100 cm2 g^-1`, so this maximum is not a sampled edge.

Concentration and baryon fraction are monotonic over the tested intervals.
The apparent preference for the most compact baryon reservoirs in the first
scan is not reliable: the joint high-concentration, high-cross-section,
`a_b/r_s=0.0005` corner reaches a nominal `8.8890e7 M_sun`, but its initial
Hernquist radius is only twice the numerical inner radius.

Of the 72 interaction cases, 26 nominally exceed `1e7 M_sun` and nine exceed
`3e7 M_sun`; none exceeds `1e8 M_sun`.  These counts describe the raw
parameter surface, not a trusted scientific frontier.

## Rejection of the extreme branch

Grid, CFL, and Roe entropy-fix tests are individually well behaved even in
the extreme models.  The decisive diagnostic is the physical inner boundary:

| Model | Nominal mass [M_sun] | Fine-grid difference | CFL spread | Entropy-fix spread | Small-boundary difference | Status |
|---|---:|---:|---:|---:|---:|---|
| `c8, sigma1, a/r_s=0.003` | 8.200e6 | 0.104% | 0.00029% | 0.00017% | 3.24% | passes |
| `c12, sigma100, a/r_s=0.003` | 2.196e7 | 0.570% | 0.000006% | 0.047% | 63.6% | rejected |
| `c12, sigma100, a/r_s=0.0005` | 8.889e7 | 0.431% | 0.000026% | 0.036% | 56.5% | rejected |

The moderate `sigma/m=3 cm2 g^-1` branch also fails.  Its nominal masses are
`9.607e6 M_sun` at `c=10` and `1.150e7 M_sun` at `c=12`, but the two smallest
completed inner boundaries differ by 12.4% and 13.5%, respectively.  The
high-mass concentration branch is therefore an unresolved capture-boundary
effect, not evidence for robust growth toward `1e8 M_sun`.

Two optional `r_min/r_infl=0.026` calculations reached the pre-set
`1.6e8`-step limit.  They did not crash or enter the analysis.  The adjacent
completed boundary pair was already sufficient to accept the representative
model and reject both extreme models under the pre-declared 5% criterion.

## Trusted compactness peak

After rejecting the extreme branch, the resolved `c=8`,
`sigma/m=1 cm2 g^-1` model was rescanned in baryonic scale radius.  It has a
clear internal maximum:

| `a_b/r_s` | Final black-hole mass [M_sun] |
|---:|---:|
| 0.003 | 8.200e6 |
| 0.005 | 1.118e7 |
| 0.0075 | 1.359e7 |
| 0.010 | 1.514e7 |
| 0.015 | 1.665e7 |
| **0.020** | **1.685e7** |
| 0.030 | 1.473e7 |

For the best nominal case, `a_b=2.507 pc` initially.  It ends at
`1.6850e7 M_sun`, of which `1.6551e7 M_sun` is accreted dark matter and
`1.9883e5 M_sun` is retained baryonic mass.  Dark matter supplies 98.81% of
the total growth.  The Hernquist radius expands to `2.615 pc`, the one-zone
gas density falls from `0.300` to `0.264 M_sun pc^-3`, and the sound speed
rises from `10.0` to `11.99 km s^-1`.

The non-monotonic radius dependence matters physically: the optimum is not
the most compact reservoir.  A very compact potential is deep but localized;
an intermediate-scale potential couples to a larger part of the resolved
SIDM inflow region.

## Final convergence and threshold result

The `a_b/r_s=0.020` peak was tested at twice the grid resolution and two
smaller physical inner boundaries:

| Variant | Final mass [M_sun] | Time to `1e7 M_sun` [Myr] |
|---|---:|---:|
| baseline, 256 cells | 1.6850e7 | 1.53 |
| 512 cells | 1.6900e7 | 1.52 |
| `r_min/r_infl=0.104` | 1.5599e7 | 1.65 |
| `r_min/r_infl=0.052` | 1.4523e7 | 1.76 |

The 256-to-512-cell difference is 0.298%.  The two smallest inner-boundary
solutions differ by 7.14%, which fails the pre-declared 5% terminal-mass
criterion.  Therefore `1.685e7 M_sun` is a nominal peak, not a fully
converged terminal-mass prediction.

The threshold statement is stronger: every final peak variant crosses
`1e7 M_sun` within 1.52-1.76 Myr, and the most conservative completed
solution ends 45.2% above the threshold.  The robust stage-5 result is thus
the crossing of `1e7 M_sun` within 2 Myr, not the third significant figure of
the final mass.

## Light-seed qualification

The trusted peak above starts from a `1e5 M_sun` heavy seed. A subsequent
fixed-halo seed ladder separates dark-matter delivery across the numerical
feeding radius from Bondi capture by an unresolved black hole. The
mass-conserving influence-gated closure reproduces resolved heavy-seed runs
to better than `4e-6` in final mass.

Under the same optimized physical model, `1e2`, `1e3`, and `1e4 M_sun` seeds
end at only `103.3`, `1103`, and `1.33e4-1.72e4 M_sun`. The first two pass
the 5% boundary criterion; the `1e4 M_sun` terminal value is boundary
sensitive, but its failure to reach even `1e5 M_sun` is common to every
variant. At nominal capture coefficient `lambda=0.25`, robust `1e7 M_sun`
crossing begins near `4e4 M_sun`; the exact intermediate-seed threshold is
sensitive to feeding radius and `lambda`, whereas stellar-remnant seeds fail
throughout `lambda=0.20-0.30`.

The original threshold result therefore remains valid for a heavy seed but
must not be described as growth from a stellar-remnant light seed. Full
details are in `stage5_light_seed_boundary_report.md`.

## Scientific meaning of a longer run

Extending the present calculation merely to ask whether it reaches
`1e7 M_sun` has no additional scientific value: that event already occurs
inside the existing two-Myr window, including the most conservative boundary
variant.

A longer run becomes meaningful only as a post-crossing stability test.  A
matched 3-5 Myr pilot could determine whether the SIDM inflow turns over,
whether cooling-regulated Bondi feeding remains Eddington limited, and
whether feedback-driven Hernquist expansion eventually suppresses the dark
channel. A direct 10 Myr extrapolation with the present static halo, one-zone
gas closure, irreversible radius expansion, and 7.14% heavy-seed boundary
sensitivity would be a numerical continuation rather than a robust
prediction. The next stage-5 physics improvement should apply a
velocity-dependent SIDM cross section and the calibrated capture closure
across the halo-redshift parameter map before interpreting a long terminal
mass.

## Numerical audit

- Jobs `40761689`, `40761906`, and `40762091` completed all 98 closure and
  interaction tasks with `0:0` exits and empty error logs.
- Job `40762942` produced 40 of 42 planned convergence files.  The two missing
  ultra-small-boundary tasks stopped only at the configured maximum step
  count; all 40 accepted files are valid.
- Jobs `40764056`, `40764460`, `40765136`, and `40765277` completed all 15
  targeted cases with `0:0` exits and empty error logs.
- Jobs `40775904`, `40776966`, `40777024`, `40778274`, `40778815`, and
  `40779138` completed all 59 light-seed closure tasks with `0:0` exits and
  empty error logs.
- The maximum mass-budget residual in the final peak convergence set is
  `3.00e-11`; the nominal peak residual is `2.45e-12`.
- Machine-readable summaries are in `results/stage5/*_summary.csv` and
  `results/stage5/*_statistics.json`.

![Trusted compactness peak](results/stage5/figures/stage5_compactness_peak.png)

![Trusted peak convergence](results/stage5/figures/stage5_trusted_peak_convergence.png)
