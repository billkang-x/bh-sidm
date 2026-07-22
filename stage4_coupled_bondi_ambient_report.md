# Stage 4 Coupled Bondi Ambient and Hernquist Feedback

Date: 2026-07-13

Status: dynamic ambient implementation, matched controls, and threshold
refinement complete

## Closure implemented

The Bondi density and sound speed now follow the same finite Hernquist
reservoir and feedback energy budget that determine the baryonic potential:

```text
rho_inf(t) = rho_inf,0 [M_b,rem(t) / M_b,0] [a_b,0 / a_b(t)]^3

c_s^2(t) = c_s,0^2 + gamma (gamma - 1) E_heat(t) / M_b,rem(t)

E_heat = chi_heat epsilon_f epsilon_r Delta M_gas c^2
E_expand = (1 - chi_heat) epsilon_f epsilon_r Delta M_gas c^2

a_b(t) = a_b,0 [1 + E_expand(t) / E_bind,eff]^eta.
```

The heating and expansion fractions sum to one, so feedback energy is not
double counted. The Bondi rate is recomputed every solver step from the
current density and sound speed before applying the Eddington cap. This is a
minimal one-zone closure: it assumes homologous expansion, instantaneous
thermal mixing, and no radiative cooling. It is not a resolved gas
hydrodynamics model.

## Experiment

The fiducial transition case uses `M_BH,0=100 M_sun`, `M_halo=1e6 M_sun`,
`M_b=5e4 M_sun`, `a_b=0.3 pc`, `T_asm=0.65 Myr`,
`rho_inf,0=300 M_sun/pc^3`, `c_s,0=10 km/s`, `sigma/m=50 cm^2/g`, and a
2 Myr MC-Roe evolution. Feedback uses the effective initial Hernquist + NFW +
black-hole binding energy and `eta=0.5`.

There are 65 cases:

- matched frozen-ambient and evolving-ambient controls;
- `chi_heat=0`, `0.5`, and `1`;
- a broad feedback scan from `epsilon_f=1e-7` to `1e-3`;
- nine additional evolving-ambient cases near the transition-prevention
  thresholds.

The no-feedback evolving-ambient case still transitions at `0.89 Myr`,
accretes `1292.45 M_sun` of dark matter and `18.22 M_sun` of retained
baryons, and ends with Bondi/Eddington `=4.53`. Reservoir depletion alone is
therefore too weak to change the baseline result.

## Transition-prevention thresholds

"Prevented" means that no post-assembly output sample enters the
Eddington-limited branch, not merely that the final state is Bondi limited.

| Heating fraction | Last case that transitions | First case prevented | Final density / initial | Final sound speed [km/s] | DM suppression vs no feedback | Retained baryons [M_sun] |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 1.3e-4 | 1.7e-4 | 0.232 | 10.00 | 21.64% | 14.31 |
| 0.5 | 5.0e-5 | 7.0e-5 | 0.620 | 14.92 | 6.62% | 15.77 |
| 1.0 | 4.0e-5 | 5.0e-5 | 1.000 | 17.04 | 0.093% | 17.15 |

Thus potential expansion does answer the original question affirmatively:
it weakens dark accretion, dilutes the gas, and can prevent the
Bondi-to-Eddington transition. In this closure it requires
`epsilon_f` between `1.3e-4` and `1.7e-4`. Heating is more efficient at
blocking the gas transition because the Bondi rate scales as `c_s^-3`; pure
heating prevents it between `4e-5` and `5e-5` while leaving the dark channel
almost unchanged. The 50/50 model is an intermediate compromise.

![Coupled Bondi ambient results](results/stage4/figures/stage4_coupled_ambient.png)

## A transient Eddington window

The pure-expansion case at `epsilon_f=1.3e-4` crosses into Eddington at
`1.26 Myr`, remains there for approximately `0.71 Myr`, and returns to the
Bondi-limited branch at `1.97 Myr`. Its final Bondi/Eddington ratio is
`0.992`. This double crossing is absent from the fixed-ambient model and
shows that the limiter state is not necessarily monotonic once the gas
environment responds to feedback.

## Feedback self-regulates its own fuel

Dynamic gas suppression also reduces the subsequent feedback energy. At
`epsilon_f=1e-3`:

- pure expansion has `E_fb/E_bind=3.50` and `763.51 M_sun` of dark accretion
  with an evolving ambient, versus `8.33` and `585.26 M_sun` with a frozen
  ambient;
- the 50/50 case has `E_fb/E_bind=2.16` and `1011.08 M_sun` of dark accretion,
  versus `9.49` and `741.69 M_sun` with a frozen ambient.

The evolving model therefore weakens the extreme dark-matter suppression
relative to a frozen gas supply: dilution/heating shuts down Bondi inflow,
which removes the fuel that would have powered further expansion. This is a
negative feedback loop, not an enhancement relative to the no-feedback
case.

## Interpretation and limits

- Expansion couples the dark and gas channels: it lowers baryonic gravity
  and the Bondi density simultaneously.
- Heating mainly controls the gas limiter and can stop the Eddington switch
  without materially changing the dark potential.
- The inferred heating threshold is optimistic because the current closure
  has no radiative cooling or escape of thermal energy.
- The density law is homologous and is not evaluated from a resolved gas
  profile at the Bondi radius.
- The binding energy remains normalized to the initial black-hole mass; a
  growing black hole would increase confinement.
- Angular momentum and radiation transport remain absent.

The next high-value calculation is a cooling-time sensitivity test. It
should replace accumulated `E_heat` with a leaky thermal reservoir and test
whether the `4e-5` to `5e-5` heating threshold survives rapid cooling.

## Numerical audit

- Jobs `40738045` (evolving cases), `40738053` (matched frozen controls),
  `40738059` (high-efficiency extension), and `40738084` (threshold
  refinement) produced all 65 required results with successful `0:0` worker
  exits.
- The initial frozen-control launch in job `40738045` exposed a shell-only
  empty-array expansion error. Those tasks generated no results and were
  replaced by job `40738053` without changing the scientific configuration.
- The maximum SIDM mass-budget residual is `7.18e-13` code units.
- The local test suite contains 79 passing tests.
