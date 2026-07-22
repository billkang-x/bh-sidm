# Stage 5 Physical-Cooling Feedback Frontier

Date: 2026-07-15

Status: 34-case Bondi, Cloudy-cooling, and feedback matrix complete

## Question

The no-feedback frontier reaches `7.0080835e6 M_sun` after 2 Myr for
`M200=1e9 M_sun`, `z=30`, `c=8`, `M_seed=1e5 M_sun`,
`sigma/m=10 cm2/g`, `f_b=0.16`, `a_b/r_s=0.003`, and
`T_asm=1.25 Myr`. This experiment restores the stage-4 dynamic Bondi
ambient, optically thin Cloudy cooling, and mixed heating/expansion feedback
to measure how much of that result survives.

## Matrix

The gas sound speed is initially `10 km/s`. Two densities separate the
relevant supply regimes:

- `rho_inf,0=300 M_sun/pc3`, the dense stage-4 anchor. The larger
  `1e5 M_sun` seed makes it strongly Eddington saturated.
- `rho_inf,0=0.3 M_sun/pc3`, scaled inversely with seed mass from the
  stage-4 transition anchor so that Bondi and Eddington supply are initially
  comparable.

The central 50/50 heating/expansion scan uses
`epsilon_f=0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2`. Timing controls use
`T_asm=1.0,1.25,1.5 Myr`. Pure-heating, pure-expansion, and solar-metallicity
controls isolate the mechanism at `epsilon_f=1e-3` and `1e-2`.

## Mass retention

The physical Bondi closure alone barely changes the ceiling. At zero
feedback the dense case ends at `7.007921e6 M_sun` (`99.9977%` retention),
and the transition-density case at `7.001893e6 M_sun` (`99.9117%`).

| Gas regime | epsilon_f | Final mass [M_sun] | Matched mass retained | DM retained | First Eddington time [Myr] |
|:---|---:|---:|---:|---:|---:|
| Dense | 1e-3 | 7.003061e6 | 99.928% | 99.928% | 0.05 |
| Transition | 1e-3 | 6.995839e6 | 99.825% | 99.862% | 0.65 |
| Dense | 1e-2 | 6.961247e6 | 99.332% | 99.320% | 0.06 |
| Transition | 1e-2 | 6.955483e6 | 99.249% | 99.280% | 0.69 |

All 34 cases enter the Eddington-limited branch and remain there. The low
density sequence delays the transition from `0.58 Myr` without feedback to
`0.69 Myr` at `epsilon_f=1e-2`, but cannot prevent it. By 2 Myr its final
Bondi/Eddington ratio is still `9.22` even in the strongest mixed case.

The original `T_asm=1.25 Myr` optimum also survives. At `epsilon_f=1e-2`,
the transition-density final masses for `T_asm=1.0,1.25,1.5 Myr` are
`6.737e6`, `6.955e6`, and `6.900e6 M_sun`. Relative to the no-feedback
case at the same assembly time, they retain `99.334%`, `99.249%`, and
`99.187%`.

![Physical-cooling feedback frontier](results/stage5/figures/stage5_feedback_frontier.png)

## Why feedback is weak

The compact `1.6e8 M_sun` Hernquist reservoir has an initial effective
Hernquist+NFW+black-hole binding energy of `1.056e57 erg`. At
`epsilon_f=1e-3`, the total deposited feedback energy reaches only
`1.94%-1.99%` of this binding energy; half is assigned to expansion. The
scale radius consequently grows by only `0.48%-0.50%`.

At `epsilon_f=1e-2`, the total feedback/binding ratio rises to
`19.2%-19.8%`. Mixed feedback expands the radius by `4.7%-4.8%`, lowers the
gas density by `12.9%-13.3%`, and suppresses the dark channel by only
`0.68%-0.72%`.

The mechanism controls confirm that expansion, not stored heat, causes the
remaining suppression:

| Channel at epsilon_f=1e-2 | Mass retained | Radius factor | Density factor | Final c_s [km/s] |
|:---|---:|---:|---:|---:|
| Pure heating | 99.997% | 1.000 | 0.999 | 10.74 |
| Mixed 50/50 | 99.332% | 1.048 | 0.867 | 10.65 |
| Pure expansion | 98.712% | 1.094 | 0.763 | 10.00 |

The minimum physical cooling time is `1.59e-4 Myr` (about 159 yr). Dense
pure heating retains only `0.021%` of injected thermal energy at
`epsilon_f=1e-2`. Changing from primordial to solar metallicity changes the
mixed final mass by less than `2 M_sun`; it cools the gas more efficiently
but does not alter the dark response.

## Interpretation

Within this closure, the `7e6 M_sun` result is robust: the standard
`epsilon_f=1e-3` stress test leaves more than `99.8%`, and even
`epsilon_f=1e-2` leaves more than `99.2%`, across both gas regimes. The
reason is physical within the model but conditional: the same compact,
high-mass baryon reservoir that catalyzes SIDM inflow is so deeply bound that
the accreted gas cannot expand it substantially in 2 Myr.

This does not establish that such a compact `f_b=0.16` reservoir is
cosmologically realizable. It establishes that ordinary optically thin
cooling and the stage-4 feedback prescription do not erase its SIDM
accretion consequence once that reservoir is assumed.

## Limits and next test

- The Bondi ambient is one-zone and its density normalization is not
  evaluated from a resolved gas profile at the Bondi radius.
- Cloudy cooling is optically thin, in equilibrium, and has no UV background,
  radiation trapping, molecular cooling, or angular momentum.
- Hernquist expansion is irreversible. The binding normalization also uses
  the initial black-hole mass, so both choices tend to overestimate late
  expansion as the black hole grows.
- The highest-growth baryon radius lies at the compact edge of the previous
  scan. A compactness/resolution consistency test is therefore more
  discriminating than another feedback-efficiency extension.

The next calculation should test the `a_b/r_s=0.003` frontier against a
resolved feeding radius and reversible radius evolution, then extend only
the surviving physical-cooling cases beyond 2 Myr.

## Numerical audit

- Slurm job `40761285` completed all 34 tasks with `0:0` exits and empty
  error logs.
- The canceled layout test `40761255` produced no scientific output and is
  excluded from the matrix.
- All 34 local NPZ files are present.
- The maximum mass-budget residual is `2.42e-12` code units.
- The analysis uses matched no-feedback ceilings at each assembly time.
