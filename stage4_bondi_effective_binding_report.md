# Stage 4 Bondi-Eddington Transition and Effective Binding

Date: 2026-07-13

Status: Bondi supply matrix and effective-binding feedback retest complete

## Refinement implemented

The prescribed Eddington inflow used in the first stage-4 matrix is replaced
by

```text
dotM_gas = min(dotM_Bondi, f_Edd dotM_Edd)

dotM_Bondi = 4 pi alpha G^2 M_BH^2 rho_b,inf
             / (c_s,b^2 + v_rel^2)^(3/2)

dotM_b,BH = (1 - epsilon_r) dotM_gas.
```

The code records the Bondi limit, Eddington limit, selected baryon rate, and
limiting branch at every output time. The ambient gas density and sound speed
are fixed parameters in this refinement; they are not yet evolved gas fields.

Feedback is also renormalized with an initial effective binding energy:

```text
E_bind,eff = E_self,Hernquist + E_NFW + E_BH.
```

For the fiducial `M_b=5e4 M_sun`, `a_b=0.3 pc`, and `M_BH=100 M_sun` case,

| Component | Binding energy [erg] |
|:---|---:|
| Hernquist self-gravity | 1.188e50 |
| NFW potential | 1.723e50 |
| Initial black hole | 1.425e48 |
| Total effective | 2.926e50 |

The effective normalization is `2.463` times the self-binding value used in
the first feedback scan.

## Bondi gas matrix

The no-feedback matrix uses `rho_b,inf=30,300,3000 M_sun/pc^3` and
`c_s,b=5,10,30 km/s`. All other parameters match the compact stage-4
`100 M_sun` seed baseline.

The matrix divides cleanly into three regimes:

- 3 cases are Eddington-limited from the initial seed mass.
- 3 cases remain Bondi-limited through 2 Myr.
- 3 cases start Bondi-limited and are pushed across the Eddington boundary by
  SIDM-driven black-hole growth.

![Bondi-Eddington matrix](results/stage4/figures/stage4_bondi.png)

## SIDM-driven limiter transitions

Because `dotM_Bondi/dotM_Edd` is proportional to
`M_BH rho_b,inf/c_s,b^3`, each ambient gas state has an analytic critical
black-hole mass. The three simulated transitions are:

| rho_b [M_sun/pc^3] | c_s [km/s] | Initial Bondi/Edd | Analytic M_crit [M_sun] | Simulated M at transition [M_sun] | Transition time [Myr] |
|---:|---:|---:|---:|---:|---:|
| 30 | 5 | 0.257 | 389.04 | 393.10 | 1.05 |
| 300 | 10 | 0.321 | 311.23 | 313.67 | 0.89 |
| 3000 | 30 | 0.119 | 840.32 | 850.88 | 1.61 |

The transition masses agree with the analytic thresholds to `0.8-1.3%`,
consistent with the 0.01 Myr output cadence. The transition is therefore
caused by the growing black-hole mass, rather than by assembly timing or an
output-classification artifact.

This is the clearest evidence so far for a genuinely sequential coupling:
SIDM accretion can move a seed from a gas-supply-limited state into an
Eddington-limited baryon-accretion state.

## Endpoint dependence on ambient gas

The retained baryon mass spans `0.173-19.800 M_sun` across the matrix.
Representative endpoints are:

| rho_b [M_sun/pc^3] | c_s [km/s] | Limiter history | Retained baryons [M_sun] | Accreted DM [M_sun] | Final M_BH [M_sun] |
|---:|---:|:---|---:|---:|---:|
| 30 | 30 | Bondi throughout | 0.173 | 1276.99 | 1377.17 |
| 300 | 10 | Bondi to Eddington at 0.89 Myr | 18.573 | 1293.73 | 1412.31 |
| 3000 | 10 | Eddington throughout | 19.800 | 1297.78 | 1417.58 |

The direct baryon channel is strongly gas-state dependent, while the dark
channel remains in the narrower `1276.99-1297.78 M_sun` range because the
same assembled Hernquist potential is present in every case.

## Effective-binding feedback retest

Feedback was rerun in two gas regimes:

- `transition`: `rho_b=300 M_sun/pc^3`, `c_s=10 km/s`;
- `eddington_saturated`: `rho_b=3000 M_sun/pc^3`, `c_s=10 km/s`.

The reversal thresholds, defined where dark growth falls below the matched
no-Eddington result, are:

| Gas regime | eta=0.5 | eta=1.0 |
|:---|---:|---:|
| Bondi-to-Eddington transition | 4.54e-6 | 2.49e-6 |
| Eddington saturated | 4.77e-6 | 2.67e-6 |

![Effective-binding feedback](results/stage4/figures/stage4_effective_feedback.png)

The thresholds rise by approximately the same factor as the binding-energy
normalization. Including NFW and initial-black-hole binding therefore weakens
feedback quantitatively but does not remove the compactness sensitivity.

In the transition gas case, feedback up to `epsilon_f=3e-5` does not prevent
the Bondi-to-Eddington transition. The strongest sampled case delays it only
from `0.89` to `0.90 Myr`. SIDM has already raised the black-hole mass above
the gas critical mass before the modest expansion substantially changes the
dark supply.

## Updated interpretation

The earlier statement that the compact setup is dark-dominated remains true
for the mass budget, but the Bondi refinement reveals a more useful sequence
in the accretion state:

1. The assembling compact baryon potential drives rapid SIDM inflow.
2. SIDM growth raises `M_BH` and hence `dotM_Bondi/dotM_Edd`.
3. In intermediate gas conditions, the baryon channel switches from Bondi to
   Eddington limited.
4. Baryon mass then deepens the point potential and modestly increases dark
   inflow, completing the positive loop.

Thus baryons need not dominate the black-hole mass to matter dynamically.
Their limiting regime can be changed by dark growth.

## Limitations and next decision

- `rho_b,inf` and `c_s,b` are fixed. Expansion feedback changes the Hernquist
  gravity but does not yet dilute or heat the ambient gas used in the Bondi
  rate.
- The effective binding energy uses the initial black-hole mass. The growing
  black hole would add further confinement.
- Angular momentum, cooling, radiative transfer, and gas depletion near the
  Bondi radius remain absent.
- The shared Hernquist potential means even the warmest, most diffuse ambient
  gas case retains the same global baryon gravity.

The next scientifically useful extension is to couple the Bondi ambient state
to the evolving Hernquist reservoir: evaluate density at a defined feeding
radius and let feedback raise the sound speed or lower that density. This will
test whether the SIDM-triggered transition survives self-consistent gas
dilution.

## Numerical audit

- Corrected Bondi job `40738000` completed all 9 cases with `0:0` exit codes.
- Effective-feedback job `40738013` completed all 18 cases with `0:0` exit
  codes.
- Remote error logs for both completed jobs are empty.
- The maximum SIDM mass-budget residual is `5.87e-13` code units.
- The local suite contains 75 passing tests.

An earlier submission (`40737992`) failed before integration because the
remote NumPy lacked a newer trapezoidal-integration API. It produced no result
files; the integration was replaced with an explicitly version-compatible
trapezoid sum before the successful rerun.
