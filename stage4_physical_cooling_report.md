# Stage 4 Density-Temperature-Metallicity Cooling

Date: 2026-07-13

Status: implementation, 2 Myr physical-cooling matrix, and trapping
sensitivity refinement complete

## Cooling closure

The prescribed cooling time is replaced by the no-UV-background Cloudy
equilibrium table distributed by Grackle. The pinned table covers
`10 <= T/K <= 1e9` and `-10 <= log10(n_H/cm^-3) <= 4`. Primordial and
solar-metal cooling coefficients are interpolated logarithmically and
combined as

```text
Lambda(T,n_H,Z) = Lambda_primordial + (Z/Z_sun) Lambda_metals,solar.
```

The local state and cooling time are

```text
n_H = X_H rho/m_H,
T = mu(T,n_H) m_H c_s^2/(gamma k_B),
t_cool = [rho c_s^2/(gamma(gamma-1))]/[n_H^2 Lambda].
```

Temperature is solved self-consistently with the tabulated equilibrium mean
molecular weight. The feedback reservoir remains

```text
dE_th/dt = q_heat dotM_gas - E_th/t_cool(rho,c_s,Z),
```

and its stiff, state-dependent update is solved implicitly by bounded
bisection in every MC-Roe step. The unheated `10 km/s` ambient state remains
a temperature floor: Cloudy cooling removes feedback heat, not the baseline
hydrostatic thermal support.

The imported table, source commit, checksum, and interpolation limits are
recorded in `data/cooling/README.md`.

## Experiment

All cases use the stage-4 transition baseline:
`M_BH,0=100 M_sun`, `M_halo=1e6 M_sun`, `M_b=5e4 M_sun`, `a_b=0.3 pc`,
`T_asm=0.65 Myr`, `rho_inf,0=300 M_sun/pc^3`, `c_s,0=10 km/s`,
`sigma/m=50 cm^2/g`, and a 2 Myr MC-Roe evolution.

The optically thin matrix contains a no-feedback control and four feedback
configurations over `Z/Z_sun=0, 0.001, 0.01, 0.1, 1`:

| Configuration | Heating fraction | Feedback efficiency |
|:---|---:|---:|
| Heating, strong | 1.0 | 1.0e-4 |
| Heating, extreme | 1.0 | 1.0e-3 |
| Mixed, strong | 0.5 | 1.4e-4 |
| Mixed, extreme | 0.5 | 1.0e-3 |

The extreme heating and mixed cases were additionally run with cooling-rate
multipliers `1e-2`, `1e-4`, and `1e-6` at `Z/Z_sun=0,1`. Pure heating was
refined at `3e-5`, `1e-5`, and `3e-6`, for 39 cases in total.

## Main result: expansion survives efficient cooling

Only the five optically thin mixed-extreme cases remain Bondi limited for the
full 2 Myr. The control and all 15 strong/extreme heating or mixed-strong
cases enter Eddington.

| Configuration | Z/Z_sun | Limiter history | Final Bondi/Eddington | Final c_s [km/s] | a_b/a_b,0 | Final rho/rho_0 |
|:---|---:|:---|---:|---:|---:|---:|
| No feedback | 0 | Eddington after 0.89 Myr | 4.531 | 10.000 | 1.000 | 1.000 |
| Heating, extreme | 0 | Eddington after 0.92 Myr | 4.127 | 10.316 | 1.000 | 1.000 |
| Heating, extreme | 1 | Eddington after 0.89 Myr | 4.441 | 10.067 | 1.000 | 1.000 |
| Mixed, extreme | 0 | Bondi throughout | 0.421 | 10.314 | 1.894 | 0.147 |
| Mixed, extreme | 1 | Bondi throughout | 0.439 | 10.057 | 1.913 | 0.143 |

Under standard cooling, extreme pure heating retains only
`0.033%-0.159%` of injected heat and barely changes the gas state. Extreme
mixed feedback also retains little thermal energy (`0.131%-0.754%`), but its
non-thermal expansion branch increases the Hernquist scale radius by about
`1.9x` and dilutes the ambient gas by `6.8-7.0x`. The dilution, rather than
persistent heating, prevents the Bondi-to-Eddington transition.

This same expansion reduces dark accretion from `1292.45 M_sun` in the
control to `864.29-868.80 M_sun`, a `32.8%-33.1%` reduction. The baryonic
mass retained by the black hole falls from `18.22 M_sun` to
`7.62-7.83 M_sun`.

## Pure-heating trapping requirement

Pure heating at `epsilon_f=1e-3` changes limiter class only between cooling
multipliers `1e-5` and `3e-6`, independently of whether `Z=0` or `Z=Z_sun`:

| Z/Z_sun | Cooling multiplier | Limiter history | Transition [Myr] | Final Bondi/Eddington | Final c_s [km/s] |
|---:|---:|:---|---:|---:|---:|
| 0 | 1e-5 | Eddington at end | 1.78 | 1.206 | 15.52 |
| 0 | 3e-6 | Bondi throughout | none | 0.707 | 18.51 |
| 1 | 1e-5 | Eddington at end | 1.75 | 1.241 | 15.37 |
| 1 | 3e-6 | Bondi throughout | none | 0.749 | 18.16 |

Thus a heating-only solution requires the optically thin cooling coefficient
to be reduced by roughly `1e5-3e5`. Ordinary metallicity variation does not
provide that suppression. This makes efficient trapping or a similarly
strong radiative-transfer effect a necessary condition for pure heating,
not a small correction.

![Physical cooling and trapping results](results/stage4/figures/stage4_cloudy_cooling.png)

## Interpretation

- The earlier no-cooling pure-heating threshold is not robust to physical
  atomic cooling. Even increasing the feedback efficiency to `1e-3` fails
  under the optically thin closure.
- Potential expansion is much more robust than thermal storage because its
  effect is cumulative and currently irreversible. It can prevent the gas
  transition while simultaneously suppressing the SIDM accretion channel.
- Metallicity has a secondary effect in this dense atomic regime. At solar
  metallicity it lowers the retained heat, but the extra gas supply can also
  inject slightly more expansion energy. It never changes limiter class in
  the scanned configurations.
- The relevant bifurcation is therefore feedback channel, not metallicity:
  heating is erased by cooling, whereas homologous expansion survives.

## Limits and next test

The Cloudy table assumes optically thin ionization equilibrium, no UV
background, and linear solar-pattern metal scaling. The baseline density
`n_H ~= 9.2e3 cm^-3` is close to the table's `1e4 cm^-3` upper edge. The
closure has no molecular cooling, self-shielding, optical depth, photon
diffusion, multiphase gas, angular momentum, or resolved hydrodynamic
outflow.

Most importantly, the expansion branch cannot recontract after feedback
power falls. The next discriminating experiment is therefore not a broader
metallicity scan. It is a reversible baryon-radius energy equation with
gravity-driven recontraction, run alongside an optical-depth estimate. This
will test whether the robust mixed-feedback branch is physical or partly an
artifact of irreversible homologous expansion.

## Numerical audit

- Base Slurm job `40738208` completed all 33 cases with `0:0` exits.
- Refinement job `40738246` completed all six cases with `0:0` exits and
  empty error logs.
- All 39 result files are present locally and the analysis consumes their
  full limiter histories, not only their final states.
- The maximum SIDM mass-budget residual is `7.34e-13` code units.
- The local `unittest` suite contains 86 passing tests.
