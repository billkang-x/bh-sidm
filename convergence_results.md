# MC and Roe Grid-Convergence Results

Date: 2026-07-11

## Numerical setup

- Spatial schemes: piecewise-constant Rusanov, MC-Rusanov, and MC-Roe.
- MC reconstruction: nonuniform-grid primitive variables `rho`, `u`, and `p`.
- Roe solver: Appendix A Roe averages with an acoustic Harten entropy fix of `0.1 a`.
- Time integration: first-order explicit Euler plus implicit conduction splitting.
- CFL number: `0.15`.
- Roe positivity fallback: disabled for all convergence runs.
- Diagnostic radius: fixed physical radius `r = 200 pc`.
- SMFP case: `rho_s = 0.8 M_sun/pc^3`, `r_s = 10 kpc`, `sigma/m = 5 cm^2/g`, `t = 1.51 Myr`.
- LMFP case: `rho_s = 0.0194 M_sun/pc^3`, `r_s = 2.586 kpc`, `sigma/m = 5 cm^2/g`, `t = 22.65 Myr`.

`rho_cond/rho_init` measures the absolute evolution. `rho_cond/rho_ad` divides out the corresponding no-conduction run at the same resolution and scheme.

## Figure 1 early-snapshot convergence

| Scheme | Cells | SMFP rho_cond/rho_init | SMFP rho_cond/rho_ad | LMFP rho_cond/rho_init | LMFP rho_cond/rho_ad |
|---|---:|---:|---:|---:|---:|
| PC-Rusanov | 128 | 0.553557 | 0.716717 | 0.766244 | 0.944463 |
| PC-Rusanov | 256 | 0.630405 | 0.717739 | 0.845392 | 0.953853 |
| PC-Rusanov | 512 | 0.675237 | 0.725356 | 0.888488 | 0.964262 |
| MC-Rusanov | 128 | 0.720029 | 0.739967 | 0.923778 | 0.975237 |
| MC-Rusanov | 256 | 0.722832 | 0.739842 | 0.929303 | 0.976593 |
| MC-Rusanov | 512 | 0.724062 | 0.740184 | 0.932013 | 0.977806 |
| MC-Roe | 128 | 0.723045 | 0.742296 | 0.927365 | 0.978318 |
| MC-Roe | 256 | 0.724882 | 0.741972 | 0.931306 | 0.978628 |
| MC-Roe | 512 | 0.725654 | 0.741800 | 0.932853 | 0.978660 |

For MC-Roe, the 256-to-512 change in `rho_cond/rho_init` is about `0.11%` in the SMFP case and `0.17%` in the LMFP case. MC reconstruction removes most of the first-order Rusanov resolution drift; changing the Riemann solver from Rusanov to Roe produces a smaller correction.

## NFW black-hole accretion convergence

Setup: `rho_s = 3.7 M_sun/pc^3`, `r_s = 30 pc`, `M_BH(0) = 100 M_sun`, `sigma/m = 50 cm^2/g`, and `t = 0.02 Myr`.

| Scheme | Cells | M_BH [M_sun] | Peak accretion rate [M_sun/Myr] | Inner density [code] |
|---|---:|---:|---:|---:|
| PC-Rusanov | 128 | 101.626329 | 93.5994 | 8120.24 |
| PC-Rusanov | 256 | 101.776139 | 108.2134 | 8729.94 |
| PC-Rusanov | 512 | 101.835171 | 114.7090 | 8854.58 |
| MC-Roe | 128 | 101.870773 | 119.0750 | 8859.90 |
| MC-Roe | 256 | 101.880515 | 119.6400 | 8902.74 |
| MC-Roe | 512 | 101.882826 | 119.7681 | 8914.70 |

For MC-Roe, the 256-to-512 change is about `0.0023%` in black-hole mass and `0.11%` in peak accretion rate. All MC-Roe runs completed without positivity fallback, and their relative mass-budget residuals remained at approximately `10^-15`.

## Interpretation

1. The dominant error in the earlier prototype was piecewise-constant reconstruction, not the choice between Rusanov and Roe alone.
2. MC-Roe provides stable, tightly converged early-time NFW accretion histories at 256 cells.
3. Quantitative 2 Myr SIS/NFW reproduction should use MC-Roe with at least 256 cells and retain a 512-cell sensitivity run.
4. Primitive-variable MC reconstruction and the entropy fix are deliberate robustness choices; they differ slightly from a literal conservative-variable, entropy-uncorrected reading of Appendix A.
