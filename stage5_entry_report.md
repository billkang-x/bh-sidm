# Stage 5 Entry: Cosmological Time Budgets and Resolved Halo Anchors

Date: 2026-07-15

Status: observable target definition, cosmological halo construction,
boundary audit, resolved 2 Myr anchors, conservative feasibility map,
12 Myr no-feedback refinement, 144-case transport screen, and baryon-frontier
assembly-time, physical-cooling feedback, frontier closure, and fixed-halo
light-seed capture refinements complete

## Scope

Stage 5 upgrades the isolated fiducial calculation into a parameter-space
and observable-feasibility study. This entry calculation addresses four
prerequisites before a large scan:

1. replace the old single LRD mass target by an observation-aware mass
   bracket;
2. attach halo models to redshift through an explicit cosmology;
3. identify which cross-halo comparisons resolve the black-hole influence
   radius;
4. quantify how a directly simulated short dark-accretion episode changes
   the later Eddington duty-cycle requirement.

## LRD target bracket

The adopted targets are `1e5`, `1e6`, and `1e7 M_sun`. Recent
electron-scattering and bolometric-correction analyses place plausible LRD
black-hole masses near `1e5-1e7 M_sun`, while the direct dark-Bondi
comparison model targets `1e7 M_sun` in a `1e9 M_sun` halo within `500 Myr`.
The source list and rationale are frozen in `data/observations/README.md`.

These values are a model-dependent envelope, not three claims of measurement
precision.

## Cosmic-time budget

The new `FlatLambdaCDM` utility uses `H0=67.4 km/s/Mpc`,
`Omega_m=0.315`, and `Omega_Lambda=0.685`. Its matter-plus-Lambda analytic
ages are:

| Redshift | Cosmic age [Myr] |
|---:|---:|
| 30 | 99.84 |
| 25 | 129.98 |
| 20 | 179.06 |
| 15 | 269.23 |
| 11 | 414.45 |
| 10 | 472.21 |
| 8 | 637.91 |
| 6 | 929.47 |
| 4 | 1536.86 |

Radiation is omitted, so this is a percent-level time budget over
`4 <= z <= 30`, not a recombination-era cosmology.

For a `100 M_sun` seed formed at `z=20`, pure Eddington baryonic growth
requires `f_Edd x duty = 1.256` to reach `1e7 M_sun` by `z=8`, but only
`0.768` by `z=6`. A `1e5 M_sun` seed requires `0.979` to reach `1e7 M_sun`
by `z=11`. The high-redshift end of the LRD interval is therefore the
discriminating regime.

![Stage-5 time budget](results/stage5/figures/stage5_time_budget.png)

## Cosmological NFW construction

The stage-4 runner now supports an `M200c`, redshift, and concentration
triplet. It constructs

```text
r200 = [3 M200 / (4 pi 200 rho_critical(z))]^(1/3)
r_s = r200/c
rho_s = M200 / [4 pi r_s^3 f(c)].
```

The old fixed-density anchor, `rho_s=3.7 M_sun/pc^3` and `c=3.9201`, is
equivalent to an `M200c` halo at `z=25.44`. It cannot represent all redshifts
by changing halo mass alone. At `c=4`, `rho_s` changes from
`0.279 M_sun/pc^3` at `z=10` to `6.237 M_sun/pc^3` at `z=30`.

## Inner-boundary audit

Two matched 20-case matrices used a fixed `100 M_sun` seed over
`M200=1e6-1e9 M_sun` and `z=10-30`:

- scaled boundary: `r_min/r_s=1/6000`;
- fixed boundary: `r_min=0.005 pc`, with approximately fixed logarithmic
  resolution.

The inferred mass trend reverses:

| Boundary protocol | Power-law slope of 2 Myr dark mass with M200 |
|:---|:---|
| Fixed `r_min/r_s` | `+0.560` at z=10 to `+0.689` at z=30 |
| Fixed `r_min` | `-0.315` at z=10 to `-0.040` at z=30 |

The two dark-accretion masses differ by as much as a factor `621.9`. Only
`5/20` scaled-boundary and `6/20` fixed-boundary cases satisfy
`r_min <= r_influence`, where

```text
r_influence = G M_BH/V200^2 = (M_BH/M200) r200.
```

The apparent positive mass trend from the first matrix is therefore not a
physical stage-5 result. It is a boundary-resolution diagnostic.

![Stage-5 boundary audit](results/stage5/figures/stage5_redshift_mass.png)

## Resolved seed-fraction anchors

The trusted anchor matrix sets `M_seed/M200=1e-4`, giving seed masses
`1e2`, `1e3`, `1e4`, and `1e5 M_sun`, and retains `r_min/r_s=1/6000`.
All 20 cases then have

```text
r_min/r_influence = 0.4167.
```

Other parameters are `c=4`, `f_b=0.05`, `a_b/r_s=0.01`,
`sigma/m=50 cm^2/g`, `T_asm=0.65 Myr`, full Eddington baryon supply, and no
feedback. The 2 Myr results are:

- growth factors span `1.763-18.208`;
- dark matter supplies `92.80%-98.59%` of the added black-hole mass;
- five cases reach `1e5 M_sun`; none reaches `1e6 M_sun` in 2 Myr;
- the largest absolute mass is `5.642e5 M_sun` for
  `M200=1e9 M_sun`, `z=30`, and `M_seed=1e5 M_sun`;
- the largest relative growth is `18.21x` for
  `M200=1e6 M_sun`, `z=30`, and `M_seed=100 M_sun`.

Higher redshift robustly increases growth at fixed halo mass. At fixed seed
fraction, lower-mass halos have the larger fractional response, while the
largest halo still gives the largest absolute black-hole mass.

![Resolved stage-5 anchors](results/stage5/figures/stage5_resolved_seed.png)

## Conservative observable map

A deliberately conservative hybrid calculation allows the numerical dark
channel to operate for only the simulated first 2 Myr. It then switches dark
accretion off and asks what constant `f_Edd x duty` is required for baryons
to reach each LRD target by `z=11,8,6,4`.

Across 228 valid combinations:

- seed-only Eddington growth succeeds in 181;
- the 2 Myr dark head start followed by baryons succeeds in 195;
- the directly simulated dark burst changes 14 combinations from impossible
  (`f_Edd x duty > 1`) to possible.

For example, a `100 M_sun` seed in a `1e6 M_sun` halo formed at `z=20`
needs activity `1.005` to reach `1e6 M_sun` by `z=8` without the dark burst,
but only `0.751` after the simulated dark head start. This is a robust
observable benefit because it does not extrapolate dark accretion beyond the
computed 2 Myr.

![Conservative hybrid feasibility](results/stage5/figures/stage5_hybrid_feasibility.png)

## Ten-Myr no-feedback anchors

Three resolved `M200=1e9 M_sun`, `M_seed=1e5 M_sun` cases test whether the
early dark flow immediately saturates; the highest-redshift case is extended
to 12 Myr:

| Halo redshift | Duration [Myr] | Final mass [M_sun] | DM added [M_sun] | Baryons added [M_sun] | Dark fraction | Time to 1e6 [Myr] | Time to 1e7 [Myr] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 10 | 1.430e6 | 1.230e6 | 1.001e5 | 0.925 | 8.64 | none |
| 20 | 10 | 3.954e6 | 3.585e6 | 2.687e5 | 0.930 | 5.04 | none |
| 30 | 12 | 1.252e7 | 1.141e7 | 1.018e6 | 0.918 | 3.10 | 10.74 |

None reaches `1e7 M_sun` by exactly 10 Myr, but all terminal growth rates are
still increasing. The `z=30` refinement reaches `1e7 M_sun` at `10.74 Myr`
and ends at `1.252e7 M_sun` after 12 Myr. Dark matter contributes
`1.141e7 M_sun`, or `91.8%` of the added mass. The high-end LRD target is
therefore crossed directly by the time-dependent solver in the no-feedback
upper-envelope model; it is not an extrapolated crossing.

![Long-time anchors](results/stage5/figures/stage5_long_pilot.png)

## Resolved transport screen

The first main screen contains 144 resolved 2 Myr cases over

```text
M200 = 1e6, 1e7, 1e8, 1e9 M_sun
z = 10, 20, 30
c = 3, 5, 8
sigma/m = 10, 30, 50, 100 cm2/g.
```

All cases retain `M_seed/M200=1e-4`; their
`r_min/r_influence=0.208-0.556`. The optimum for every one of the 12
`M200-z` pairs is `c=8` and `sigma/m=10 cm2/g`. The global largest relative
growth is `78.42x` for `M200=1e6 M_sun`, `z=30`; the largest absolute mass is
`2.549e6 M_sun` for `M200=1e9 M_sun`, `z=30`. Across the complete screen,
46 cases reach `1e5 M_sun`, six reach `1e6 M_sun`, and none reaches
`1e7 M_sun` within 2 Myr.

The preference for the smallest scanned cross section has a transport
interpretation. At the representative `M200=1e9 M_sun`, `z=30`, `c=8`
baryon scale radius, the initial Knudsen numbers for cross sections
`10,30,50,100` are `0.0840,0.0280,0.0168,0.00840`; after assembly they are
all below `0.0035`. These central states are on the short-mean-free-path side
of the interpolation, where the effective conductivity decreases close to
`1/sigma`. Raising the cross section then traps compression heat and weakens
the dark inflow instead of accelerating it.

![Resolved transport screen](results/stage5/figures/stage5_transport.png)

## Baryon-frontier screen

The transport optimum was held fixed at `M200=1e9 M_sun`,
`M_seed=1e5 M_sun`, `c=8`, and `sigma/m=10 cm2/g`. An 81-case no-feedback
screen then varied

```text
z = 10, 20, 30
f_b = 0.01, 0.05, 0.16
a_b/r_s = 0.003, 0.01, 0.03
T_asm = 0.2, 0.65, 1.5 Myr.
```

The compact, high-baryon corner (`f_b=0.16`, `a_b/r_s=0.003`) maximizes
growth at all three redshifts. Its best 2 Myr masses are `1.832e6 M_sun` at
`z=10`, `3.674e6 M_sun` at `z=20`, and `6.957e6 M_sun` at `z=30`. Across
the full matrix, 34 cases reach `1e6 M_sun`, while none reaches
`1e7 M_sun`. The compact baryonic potential primarily catalyzes dark inflow:
dark matter supplies `98.25%-98.55%` of the added mass in the three
redshift optima.

![Baryon-frontier screen](results/stage5/figures/stage5_baryon_frontier.png)

## Assembly-time refinement

Because the coarse `z=20` and `z=30` optima occurred at the largest sampled
assembly time, ten additional cases resolved `T_asm=1.0-2.0 Myr` in
`0.25 Myr` intervals at `f_b=0.16` and `a_b/r_s=0.003`. Both redshifts have
a smooth internal maximum at `T_asm=1.25 Myr`:

| Redshift | Best final mass [M_sun] | DM added [M_sun] | Baryons added [M_sun] | Dark fraction |
|---:|---:|---:|---:|---:|
| 20 | 3.705e6 | 3.548e6 | 5.686e4 | 0.9842 |
| 30 | 7.008e6 | 6.802e6 | 1.060e5 | 0.9847 |

The `z=30` mass rises from `6.783e6 M_sun` at `1.0 Myr` to `7.008e6 M_sun`
at `1.25 Myr`, then falls monotonically to `6.273e6 M_sun` at `2.0 Myr`.
The optimum is therefore not a scan-edge artifact. It reveals a finite
assembly window in which the growing Hernquist potential most efficiently
drives dark inflow. Separating the roles of potential growth, compressional
heating, and later feedback requires the next controlled comparison; the
present no-feedback result alone does not identify the cause of the peak.

![Assembly-time refinement](results/stage5/figures/stage5_baryon_timing_refinement.png)

## Physical-cooling feedback frontier

Thirty-four cases restored the evolving Bondi ambient, optically thin Cloudy
cooling, and feedback split equally between heating and Hernquist expansion.
They compare the dense stage-4 gas anchor (`rho_inf=300 M_sun/pc3`) with a
seed-mass-scaled transition anchor (`0.3 M_sun/pc3`).

At the standard extreme test `epsilon_f=1e-3`, the dense and transition
cases retain `99.928%` and `99.825%` of the `7.008e6 M_sun` ceiling. Even at
`epsilon_f=1e-2`, they retain `99.332%` and `99.249%`. All cases enter the
Eddington branch; mixed feedback delays but never prevents the low-density
transition. The `T_asm=1.25 Myr` absolute-mass optimum remains unchanged.

The compact baryon reservoir has effective binding energy `1.056e57 erg`.
At `epsilon_f=1e-2`, mixed feedback expands it by only about `4.8%` and
suppresses dark accretion by `0.68%-0.72%`. Pure-heating controls retain
`99.997%` of the mass because Cloudy cooling reduces stored heat below
`0.021%`; pure expansion gives the strongest tested suppression but still
retains `98.712%`.

![Physical-cooling feedback frontier](results/stage5/figures/stage5_feedback_frontier.png)

The complete matrix and limitations are documented in
`stage5_feedback_frontier_report.md`.

## Frontier closure and trusted threshold crossing

A subsequent 26-case one-axis closure, 72-case interaction scan, and targeted
numerical-convergence campaign changed the location of the trusted frontier.
The raw high-concentration, high-cross-section corner reaches a nominal
`8.889e7 M_sun`, but changes by `56%-64%` under inner-boundary refinement and
is rejected.  Moderate `c=10-12`, `sigma/m=3 cm2 g^-1` candidates also fail
the 5% boundary criterion.

The surviving `c=8`, `sigma/m=1 cm2 g^-1` branch has an internal baryonic
scale-radius maximum at `a_b/r_s=0.020`.  Its nominal final mass is
`1.6850e7 M_sun`; a 512-cell run differs by only 0.298%.  The two smallest
inner-boundary solutions differ by 7.14%, so the exact terminal mass is not
yet converged to the pre-declared 5% standard.  Nevertheless, every tested
variant crosses `1e7 M_sun` at `1.52-1.76 Myr`, and the most conservative
completed solution ends at `1.4523e7 M_sun`.  Threshold crossing within 2 Myr
is therefore robust even though the terminal mass is not.

This result removes the motivation for a longer run whose only purpose is to
reach `1e7 M_sun`.  A 3-5 Myr extension is useful only as a post-crossing
stability and saturation test, preferably after improving the resolved
feeding/capture closure.  Full details are in
`stage5_frontier_closure_report.md`.

![Trusted peak convergence](results/stage5/figures/stage5_trusted_peak_convergence.png)

## Light-seed and capture closure

The trusted `1e5 M_sun` seed result was followed by a fixed-halo seed ladder
that separates SIDM supply across the numerical feeding radius from capture
by an unresolved black hole. A conservative central reservoir uses the
gamma=`5/3` dark Bondi coefficient `lambda=0.25` outside the influence region
and switches to direct flux capture only after
`r_feed/r_influence <= 0.104`. The gated closure reproduces resolved
heavy-seed runs to better than `4e-6` in final mass.

The result sharply limits the interpretation of the positive frontier.
Seeds of `1e2`, `1e3`, and `1e4 M_sun` end at approximately `103`, `1103`,
and `1.33e4-1.72e4 M_sun`, respectively, and never approach the LRD target.
The stellar-remnant cases pass the pre-declared 5% boundary and grid tests.
Their failure persists over `lambda=0.20-0.30` even though large amounts of
SIDM cross the feeding radius, because less than `0.2%` is captured.

At nominal `lambda=0.25`, both accepted feeding radii reach `1e7 M_sun` for
seeds of `4e4 M_sun` and above. Seeds `3e4-3.75e4 M_sun` are boundary
ambiguous, and the `4e4 M_sun` outcome itself changes across
`lambda=0.20-0.30`. The model therefore supports rapid amplification of an
intermediate/heavy seed, not direct formation of an LRD black hole from a
stellar-remnant seed. Full results are in
`stage5_light_seed_boundary_report.md`.

![Light-seed boundary closure](results/stage5/figures/stage5_light_seed_boundary.png)

## Interpretation and scan design

The stage-5 conclusion is now two-sided. A resolved early dark episode
materially relaxes the baryonic duty-cycle requirement, and the optimized
physical model rapidly amplifies intermediate/heavy seeds. The same model
does not rescue stellar-remnant seeds on a two-Myr timescale because capture,
not central delivery, is the bottleneck.

The remaining work should stay staged rather than becoming fully factorial:

1. replace the constant cross section by physically motivated
   velocity-dependent `sigma(v)/m` models and retest the seed threshold;
2. propagate the calibrated capture closure across `M_halo`, redshift, and
   concentration to obtain the final success/failure map;
3. replace irreversible Hernquist expansion by reversible radius evolution
   with the growing black hole included in the binding energy;
4. run a matched 3-5 Myr post-crossing pilot, and extend to 10 Myr only if the
   rate remains physically unsaturated.

This ordering prevents expensive feedback calculations in parameter regions
that already fail the redshift or resolution budget.

## Limits

- The 2-10 Myr anchors use static `M200`, concentration, and cosmological
  background; halo mergers and smooth cosmological growth are absent.
- The original no-feedback anchors use Eddington baryon supply. The feedback
  frontier adds a one-zone Bondi ambient but still lacks resolved gas
  hydrodynamics and angular momentum.
- The Hernquist assembly histories are prescribed and do not follow a
  cosmological galaxy model; the refined optimum is conditional on this
  smoothstep growth law.
- Fixed seeds now use an influence-gated Bondi reservoir. This is a
  calibrated one-zone capture closure, not resolved relativistic or
  phase-space transport inside the feeding radius.
- The SIDM cross section is constant. Velocity dependence remains a later
  stage-5 axis.
- The Planck time budget neglects radiation, while the fluid model neglects
  relativistic capture and resolved baryonic hydrodynamics.

## Numerical audit

- Jobs `40739143`, `40739246`, and `40739322` produced 60 valid 2 Myr files;
  all Slurm tasks exited `0:0` and all error logs were empty.
- Job `40758549` produced the three 10 Myr anchors with `0:0` exits and empty
  error logs.
- Job `40758751` completed the 12 Myr refinement with a `0:0` exit and an
  empty error log.
- Job `40758762` completed all 144 transport cases with `0:0` exits and empty
  error logs.
- Job `40759419` completed all 81 baryon-frontier cases with `0:0` exits and
  empty error logs.
- Job `40759476` completed all ten assembly-time refinements with `0:0` exits
  and empty error logs.
- Job `40761285` completed all 34 physical-cooling feedback cases with `0:0`
  exits and empty error logs.
- Jobs `40761689`, `40761906`, and `40762091` completed all 98 frontier
  closure and interaction cases with `0:0` exits and empty error logs.
- Job `40762942` produced 40 of 42 planned convergence cases; two optional
  ultra-small-boundary cases reached the configured step limit.  Jobs
  `40764056`, `40764460`, `40765136`, and `40765277` completed all 15 targeted
  closure and convergence cases.
- Jobs `40775904`, `40776966`, `40777024`, `40778274`, `40778815`, and
  `40779138` completed all 59 light-seed, gate-validation, threshold, and
  capture-sensitivity cases with `0:0` exits and empty error logs.
- Maximum mass-budget residuals are `2.51e-12` in the boundary audit,
  `5.67e-13` in the resolved matrix, `7.23e-13` in the long pilots, and
  `3.51e-12` in the transport screen, `6.90e-12` in the baryon-frontier
  screen, `1.28e-12` in the assembly-time refinement, and `2.42e-12` in the
  physical-cooling feedback matrix.
- The local test suite contains 98 passing tests after the light-seed capture
  additions.
