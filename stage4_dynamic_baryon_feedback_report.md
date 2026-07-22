# Stage 4 Dynamic Baryon Accretion and Feedback

Date: 2026-07-12

Status: first no-feedback baseline and parametric expansion-feedback study complete

## Model implemented

The central mass is now evolved with two separately recorded channels:

```text
dM_BH/dt = dM_DM/dt + dM_b,BH/dt
```

The dark component is the measured MC-Roe inner-boundary flux. Baryonic gas
is drawn from a finite Hernquist reservoir at a prescribed fraction of the
Eddington inflow rate:

```text
dM_gas/dt  = f_Edd f_duty 4 pi G M_BH m_p
             / (epsilon_r sigma_T c)
dM_b,BH/dt = (1 - epsilon_r) dM_gas/dt.
```

The Eddington rate therefore responds to the full instantaneous black-hole
mass, including previously accreted dark matter. With `epsilon_r=0.1`, the
black-hole e-folding time at `f_Edd=f_duty=1` is `50.054 Myr`.

The reservoir is assembled with the stage-3 smoothstep history. Gas cannot be
accreted before it has assembled, and gas inflow reduces the remaining
Hernquist mass. The output records total black-hole mass, both instantaneous
rates, both retained mass contributions, total gas consumed, the remaining
reservoir, and the SIDM mass-budget residual.

## Feedback prescription

Feedback uses the minimal Hernquist expansion model selected in the project
blueprint:

```text
E_fb = epsilon_f epsilon_r Delta M_gas c^2
E_bind = G M_b^2 / (6 a_b,0)
a_b(t) = a_b,0 [1 + E_fb(t)/E_bind]^eta.
```

This changes the baryon enclosed-mass profile in the gravity source at every
step. It does not inject heat directly into the SIDM fluid.

## Controlled baseline

All runs use the stage-3 compact setup:

```text
M_halo = 1e6 M_sun
f_b = 0.05
a_b/r_s = 0.01
T_asm = 0.65 Myr
sigma/m = 50 cm^2/g
duration = 2 Myr
```

The no-feedback matrix varies seed mass over `10, 100, 1000 M_sun` and
Eddington ratio over `0, 0.1, 1`. The feedback matrix fixes the `100 M_sun`
seed and scans `epsilon_f=1e-7` through `1e-2` for `eta=0.5,1`.

## Regression to stage 3

With baryon accretion and feedback disabled, the new evolving-reservoir path
returns `1276.889397 M_sun` of dark accretion. Its final black-hole mass differs
from the stage-3 `T_asm=0.65 Myr` result by only `8.54e-7 M_sun`; the maximum
relative final-density difference is `4.87e-10`. The new path therefore
preserves the stage-3 solution.

## No-feedback results

The full-Eddington endpoints are:

| Seed [M_sun] | Final M_BH [M_sun] | Accreted DM [M_sun] | Retained baryons [M_sun] | Extra DM / baryon mass | Baryon growth / isolated Eddington | Dark fraction of growth |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 625.48 | 609.27 | 6.21 | 1.758 | 15.229 | 0.9899 |
| 100 | 1417.58 | 1297.78 | 19.80 | 1.055 | 4.857 | 0.9850 |
| 1000 | 3545.13 | 2466.75 | 78.38 | 0.371 | 1.923 | 0.9692 |

Here `extra DM` is measured against the matched `f_Edd=0` case for the same
seed. The `f_Edd=0` final masses are `608.35`, `1376.89`, and `3437.71 M_sun`.

![No-feedback channel separation](results/stage4/figures/stage4_no_feedback.png)

### Two-way catalysis

The coupling operates in both directions:

1. Retained baryon mass deepens the point-mass potential and produces
   `0.37-1.78 M_sun` of additional dark accretion per retained baryonic solar
   mass.
2. Dark accretion raises the Eddington limit. Consequently the retained
   baryon mass is `1.92-15.23` times larger than an isolated seed would gain
   at the same prescribed Eddington ratio.

The effect is strongest for the lightest seed because dark growth changes its
mass by the largest factor. Nevertheless, dark matter already supplies
`96.9-99.0%` of total growth in every full-Eddington run. There is no
post-assembly transition from baryon-dominated to dark-dominated growth in
this compact 2 Myr setup: it is dark-dominated from the start of the measured
post-assembly interval.

No case reaches `1e4 M_sun` within 2 Myr. The largest endpoint is
`3545.13 M_sun` for the `1000 M_sun` seed.

## Feedback results

The no-feedback `100 M_sun` run gains `1297.78 M_sun` of dark matter, compared
with `1276.89 M_sun` when Eddington baryon accretion is disabled. Expansion
feedback erases this positive difference at very low effective coupling:

| eta | Log-interpolated epsilon_f reversal | First sampled reversal | E_fb/E_bind there | a_b/a_b,0 there |
|---:|---:|---:|---:|---:|
| 0.5 | 2.04e-6 | 3e-6 | 0.0983 | 1.0480 |
| 1.0 | 1.08e-6 | 3e-6 | 0.0973 | 1.0973 |

At the first sampled reversal, feedback has supplied only about ten percent of
the adopted self-binding energy. A few-percent expansion of the compact
potential is already sufficient to remove the Eddington-induced dark
catalysis.

![Feedback and potential expansion](results/stage4/figures/stage4_feedback.png)

At the strongest sampled feedback, `epsilon_f=0.01` and `eta=1`, the scale
radius expands by a factor `85.0`; dark growth falls to `2.49%` and total
growth to `2.83%` of the no-feedback values. This extreme endpoint is a stress
test, not a calibrated physical model.

## Interpretation

The stage-4 result changes the simple sequential picture. In this compact
halo, Eddington baryon accretion does not first build most of the black hole
and then hand growth to SIDM. Instead, rapid SIDM inflow and Eddington growth
form a positive loop, but SIDM remains the dominant mass source. The same loop
is fragile because its dark component depends steeply on maintaining a very
compact baryon potential.

Thus the useful physical discriminator is not only Eddington ratio. It is the
competition between the mass-catalysis loop and the time-dependent central
compactness.

## Model limitations

- The gas supply is prescribed at an Eddington fraction; there is no Bondi
  supply calculation, gas pressure, cooling, angular momentum, or duty-cycle
  evolution.
- Feedback is converted instantly and globally into Hernquist expansion.
  There is no radiative transfer, delayed coupling, or hydrodynamic outflow.
- `E_bind` includes only Hernquist self-binding. Adding halo and central-mass
  binding would raise the effective reversal efficiency. The quoted
  `epsilon_f` thresholds are therefore model parameters, not observationally
  calibrated efficiencies.
- The 2 Myr endpoint is too short for ordinary isolated Eddington growth to
  reach the Little Red Dot mass range.

## Numerical audit and next test

- Slurm jobs `40734812`, `40734822`, and `40734833` completed with `0:0` exit
  codes and empty error logs.
- All `9 + 9 + 8` expected matrix files are present.
- The largest SIDM mass-budget residual is `6.67e-13` code units.
- The local suite contains 71 passing tests.

That refinement is now complete; see
`stage4_bondi_effective_binding_report.md`. The next step is to couple the
ambient Bondi density and sound speed to the evolving Hernquist reservoir and
feedback state before extending the model to `10-100 Myr` and comparing with
`1e4-1e6 M_sun` targets.
