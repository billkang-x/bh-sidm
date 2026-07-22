# Stage 3 Timescale and Transport Mechanism

Date: 2026-07-12

Status: local timescale diagnostic, refined assembly-time scan, and heat-off discriminant complete

## Diagnostic definitions

The diagnostic uses the same dimensionless mass and conductivity closure as the evolution code:

```text
t_dyn  = sqrt(r^3 / M_total(<r))
t_coll = 1 / (a rho (sigma/m) v)
t_cond = 3 rho L^2 / (2 kappa)
t_in   = r / |u_r|
Kn     = lambda / H
```

Here `a = sqrt(16/pi)`, `kappa` is the simulation's SMFP/LMFP-interpolated conductivity, and the primary conductive estimate uses `L = r`. A second estimate uses the local `v^2` gradient length clipped between one cell width and `r`.

The supply location is defined as the log-radius median of inward mass flux inside the NFW scale radius. For causal checks, this radius is measured at the assembly midpoint, but all compared timescales are then evaluated in the original `t = 0` no-baryon hydrostatic state.

## Coarse-matrix causal check

Nine `(sigma/m, a_b/r_s)` cells in the original matrix prefer a nonzero sampled assembly time. At their future supply radii, the ratios of the sampled optimum to initial local times are:

| Initial timescale | Median T_asm/t | Range | RMS log ratio [dex] | log correlation |
|:---|---:|---:|---:|---:|
| Dynamical | 0.985 | 0.775-1.408 | 0.065 | 0.767 |
| Conduction, L=r | 0.984 | 0.348-1.083 | 0.202 | 0.312 |
| Collision | 2.01 | 0.495-8.64 | 0.510 | -0.354 |

![Timescale matching](results/stage3/figures/stage3_timescale_matching.png)

The result disfavors direct matching to the particle collision time. Both dynamical and conductive times are relevant, but the initial dynamical time has substantially smaller scatter.

## Refined compact-potential optima

The strongest `a_b/r_s = 0.01` cases were extended to `T_asm = 2 Myr` and locally refined around each maximum.

| sigma/m [cm^2/g] | Refined T_opt [Myr] | Feed radius [pc] | Initial t_dyn [Myr] | Initial t_cond [Myr] | Initial t_coll [Myr] | Accreted DM [M_sun] | Enhancement |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0.60 | 1.019 | 0.545 | 1.046 | 0.812 | 3426.34 | 102.41 |
| 30 | 0.50 | 0.914 | 0.508 | 0.462 | 0.248 | 2253.71 | 94.75 |
| 50 | 0.65 | 1.135 | 0.584 | 0.558 | 0.177 | 1276.89 | 57.72 |
| 100 | 0.75 | 1.334 | 0.645 | 0.880 | 0.101 | 634.37 | 26.71 |

Across these refined maxima,

- `T_opt/t_dyn = 0.985-1.163`, with RMS log scatter `0.046 dex`;
- `T_opt/t_cond = 0.573-1.164`, with RMS log scatter `0.131 dex`;
- `T_opt/t_coll = 0.739-7.422`, with RMS log scatter `0.545 dex`.

![Refined timescale matching](results/stage3/figures/stage3_refined_timescale_matching.png)

The `sigma/m = 10` maximum is a broad plateau: results between `0.5` and `0.7 Myr` differ by approximately `0.12%`. This is therefore a dynamical matching window, not evidence for a narrow resonance.

![Initial radial timescales](results/stage3/figures/stage3_initial_timescale_profiles.png)

## Heat-off discriminant

The same compact potential was evolved with the conductive substep disabled. This is an adiabatic fluid diagnostic, not a collisionless-dark-matter model.

| Transport | T_opt [Myr] | Accreted DM [M_sun] | Matched-control enhancement |
|:---|---:|---:|---:|
| Heat off | 1.00 | 275.08 | 2.94 |
| sigma/m = 10 | 0.60 | 3426.34 | 102.41 |
| sigma/m = 30 | 0.50 | 2253.71 | 94.75 |
| sigma/m = 50 | 0.65 | 1276.89 | 57.72 |
| sigma/m = 100 | 0.75 | 634.37 | 26.71 |

![Mechanism follow-up](results/stage3/figures/stage3_mechanism_followup.png)

The finite-time optimum survives without conduction, but its enhancement is modest and its optimum shifts to a longer time. Dynamics therefore creates the assembly-time window, while conductive transport is essential for the one-to-two-order-of-magnitude amplification.

## Heat-flow direction

At every refined conducting optimum, the midpoint `v^2` gradient is negative at both the Hernquist radius and the inward-flux median radius. The corresponding `q_r = -kappa d(v^2)/dr` is positive: heat flows outward.

At the feed radius, the normalized outward heat flux decreases from approximately `1.96` at `sigma/m = 10` to `0.59` at `100`. The conductive term is therefore removing compressional heat from the inflowing region, allowing sustained contraction; it is not enhancing accretion by heating the center from outside.

## Interpretation

1. The assembly duration controls how far the halo can respond coherently. The optimal compact cases reach a supply region near `0.9-1.3 pc`, whose initial dynamical time tracks `T_opt`.
2. The cross section changes heat-removal efficiency and the resulting supply radius. This moves the broad optimal window and strongly changes its amplitude.
3. The local collision time does not predict the optimum and should not be used as the assembly clock.
4. The static-equilibrium suppression is absent because that protocol preloads thermal support; the assembled potential instead creates compressional heating followed by outward conductive cooling and contraction.

## Limitations and next test

- The future supply radius is measured from each evolved solution. Evaluating its timescale at `t = 0` avoids thermodynamic circularity, but the radius itself is not yet an independent pre-run prediction.
- The refined mechanism scan covers only `a_b/r_s = 0.01`, `f_b = 0.05`, and one NFW halo.
- The endpoint remains fixed at 2 Myr, and several strong cases still have large final accretion rates.

The next stage-3 test should vary baryon fraction and halo mass while choosing `T_asm` from the initial dynamical time at a consistently defined early supply radius. Success would turn the present correlation into a predictive scaling relation.
