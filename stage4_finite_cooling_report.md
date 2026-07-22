# Stage 4 Finite-Cooling Feedback Reservoir

Date: 2026-07-13

Status: implementation, 2 Myr cooling scan, and threshold refinement complete

## Model

The accumulated feedback heat used by the dynamic Bondi sound speed is
replaced by a leaky one-zone reservoir:

```text
dE_th/dt = q_heat dotM_gas - E_th/t_cool

q_heat = chi_heat epsilon_f epsilon_r c^2

c_s^2 = c_s,0^2 + gamma (gamma - 1) E_th/M_b,rem.
```

For a constant gas rate over one MC-Roe step, the reservoir is updated
analytically:

```text
E_th,n+1 = E_th,n exp(-dt/t_cool)
           + q_heat dotM_gas t_cool [1-exp(-dt/t_cool)].
```

The implementation uses `expm1` for the small `dt/t_cool` limit. Omitting
`t_cool`, or explicitly supplying infinity, exactly recovers accumulated
no-cooling feedback. The output records injected heat, retained heat, and
cumulative cooling loss separately.

Cooling applies only to the thermal reservoir. The expansion-energy branch
remains cumulative and does not recontract after thermal losses.

## Experiment

All cases use the stage-4 transition baseline:
`M_BH,0=100 M_sun`, `M_halo=1e6 M_sun`, `M_b=5e4 M_sun`, `a_b=0.3 pc`,
`T_asm=0.65 Myr`, `rho_inf,0=300 M_sun/pc^3`, `c_s,0=10 km/s`,
`sigma/m=50 cm^2/g`, and a 2 Myr MC-Roe evolution.

Four feedback configurations were scanned over `t_cool=0.001-100 Myr` plus
the no-cooling limit:

- the first no-cooling prevention point for 50/50 heating and expansion,
  `epsilon_f=7e-5`;
- a stronger 50/50 case, `epsilon_f=1.4e-4`;
- the first no-cooling pure-heating prevention point, `epsilon_f=5e-5`;
- a stronger pure-heating case, `epsilon_f=1e-4`.

The broad scan and three refinements contain 48 cases.

## Critical cooling times

The bracket gives the last cooling time that still enters Eddington and the
first time that prevents the transition at every post-assembly sample.

| Configuration | chi_heat | epsilon_f | Critical t_cool [Myr] | Retained heat at first prevention | Final c_s [km/s] | DM accreted [M_sun] |
|:---|---:|---:|---:|---:|---:|---:|
| Mixed, threshold feedback | 0.5 | 7.0e-5 | 2.0-3.0 | 0.850 | 14.45 | 1205.27 |
| Mixed, strong feedback | 0.5 | 1.4e-4 | 0.15-0.20 | 0.237 | 12.58 | 1141.76 |
| Heating, threshold feedback | 1.0 | 5.0e-5 | 10-30 | 0.984 | 16.98 | 1291.30 |
| Heating, strong feedback | 1.0 | 1.0e-4 | 0.30-0.50 | 0.479 | 16.74 | 1290.72 |

![Finite cooling results](results/stage4/figures/stage4_cooling.png)

The no-cooling threshold points are not robust to ordinary thermal leakage.
The marginal pure-heating case requires a cooling time much longer than both
the `0.65 Myr` assembly time and the entire 2 Myr experiment. Its earlier
low feedback-efficiency threshold was therefore an almost adiabatic result.

Doubling the feedback changes the conclusion. Strong mixed feedback survives
when `t_cool` exceeds only `0.15-0.20 Myr`, while strong pure heating requires
`0.30-0.50 Myr`. Strong feedback can therefore prevent the transition with
finite cooling on a sub-Myr timescale.

## Transient Eddington windows

Several cases end Bondi limited but enter Eddington temporarily:

| Configuration | t_cool [Myr] | Bondi-to-Eddington [Myr] | Return to Bondi [Myr] |
|:---|---:|---:|---:|
| Mixed, strong | 0.15 | 1.36 | 1.76 |
| Mixed, threshold | 2.0 | 1.31 | 1.66 |
| Heating, threshold | 10 | 1.45 | 1.68 |

The strong pure-heating case at `t_cool=0.3 Myr` enters Eddington at
`1.28 Myr` and remains there through 2 Myr. A final Bondi/Eddington ratio
below unity is therefore insufficient to establish transition prevention;
the full limiter history is required.

## Dark-channel response

Pure heating remains cleanly separated from dark accretion. Across the
strong pure-heating cooling sequence, dark accretion changes only from
`1292.44` to `1288.25 M_sun`, while the gas limiter changes qualitatively.

Mixed feedback behaves differently. Increasing `t_cool` suppresses Bondi
fuel earlier, which reduces subsequent expansion-energy injection. In the
strong mixed sequence, dark accretion consequently rises from
`1133.09 M_sun` at `t_cool=0.001 Myr` to `1166.06 M_sun` without cooling.
This is the same feedback-fuel self-regulation found in the dynamic ambient
scan: more persistent heat can indirectly preserve dark accretion by
starving the expansion channel.

## Interpretation and limits

- Cooling restores Bondi-to-Eddington conversion unless the thermal
  retention time exceeds a feedback-dependent threshold.
- The relevant control is not `epsilon_f` or `t_cool` alone, but their
  combination and the heating/expansion partition.
- Pure heating changes the gas state with negligible dark-potential response;
  mixed feedback couples gas regulation to dark accretion.
- The cooling law is a prescribed exponential timescale, not a density- and
  temperature-dependent atomic cooling function.
- Thermal energy is instantaneously mixed through the remaining Hernquist
  reservoir; multiphase gas, radiative transfer, and angular momentum remain
  absent.
- The expansion branch is irreversible in the current closure. A model with
  recontraction could reduce the long-lived dark suppression of mixed cases.

The next useful extension is to replace prescribed `t_cool` with
`t_cool(rho,T,Z)` from an explicit cooling function, first in the same
one-zone reservoir and then, only if warranted, in resolved gas dynamics.

## Numerical audit

- Job `40738103` completed the 32 finite-cooling cases with successful `0:0`
  exits.
- Its four no-cooling rows initially received `inf` with a CRLF suffix from
  the Windows manifest, exposing an infinity-times-zero numerical form.
  No scientific result used those files. The parser, kernel, and manifest
  normalization were hardened, and job `40738125` replaced all four controls.
- Refinement jobs `40738132` and `40738143` completed all 12 additional cases
  with `0:0` exits.
- The three no-cooling configurations available from the previous scan agree
  to better than `2.5e-11 M_sun` in dark mass and `2.3e-13 M_sun` in baryon
  mass.
- The maximum SIDM mass-budget residual is `7.43e-13` code units.
- The local test suite contains 80 passing tests.
