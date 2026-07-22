# Stage 5 Light-seed and Inner-boundary Closure

## Scientific question

The previous trusted frontier used a `1e5 M_sun` seed in a
`1e9 M_sun`, `z=30` halo. This establishes rapid amplification of a heavy
seed, but does not test whether a stellar-remnant seed can enter the same
runaway. Directly preserving `r_min/r_influence` for a `1e2 M_sun` seed
would shrink the physical inner boundary by a factor of 1000. The global
MC-Roe CFL step scales approximately with that radius and makes a direct
two-Myr calculation impractical.

This study separates delivery across a resolved feeding radius from capture
by an unresolved black hole, then tests fixed seed masses in one fixed halo.
All other parameters retain the physical-cooling trusted frontier:

```text
M200 = 1e9 M_sun, z = 30, c = 8
f_b = 0.16, a_b/r_s = 0.020, T_asm = 1.25 Myr
sigma/m = 1 cm2/g
rho_b,inf(0) = 0.3 M_sun/pc3, c_s(0) = 10 km/s
epsilon_f = 1e-3, chi_heat = 0.5
Cloudy cooling, MC reconstruction, Roe flux, t_end = 2 Myr.
```

## Conservative capture closure

The inner fluid flux is now recorded as a supply rate rather than being
unconditionally added to the black hole:

```text
dotM_supply = max[-r_feed^2 rho(r_feed) u(r_feed), 0]
dotM_capture = lambda M_BH^2 rho(r_feed) / c_s(r_feed)^3
```

The second expression is the gamma=`5/3` Bondi rate in the paper's code
units. Supplied but uncaptured dark matter enters a conservative central
reservoir and contributes to the gravitational field. The adopted
`lambda=0.25` is the analytic gamma=`5/3` value. The three directly resolved
heavy-seed runs independently imply post-assembly effective values
`0.261`, `0.265`, and `0.271`.

A local Bondi formula is not applicable after the feeding radius lies inside
the black-hole influence region. The final closure is therefore influence
gated: Bondi capture operates while `r_feed/r_influence > 0.104`; after the
black hole grows past that threshold, current supply is captured directly
and the central reservoir drains on its local free-fall time.

A deliberately attempted pure-Bondi validation failed by construction: it
left `96.6%-98.8%` of the supplied dark matter in the reservoir even for the
already resolved `1e5 M_sun` seed. Those five diagnostic cases are retained
but are not used as the accepted model. The influence-gated closure matches
the three direct heavy-seed runs to at most `3.97e-6` in final mass; the two
smallest resolved boundaries agree to `2.90e-7`. Moving the gate from
`r_feed/r_influence=0.052` to `0.208` changes the calibrated final mass by
only `2.42e-6` fraction.

## Fixed-halo seed ladder

The main matrix varies only seed mass and feeding radius. The two smallest
feeding radii are `0.052` and `0.104` times the initial influence radius of
the reference `1e5 M_sun` seed. Logarithmic cell counts preserve the
baseline radial resolution, and every seed also has a double-resolution
test at the second radius.

| Seed [M_sun] | Conservative final mass [M_sun] | Growth factor | Boundary difference | Grid difference | `1e7 M_sun` in 2 Myr |
|---:|---:|---:|---:|---:|---|
| `1e2` | 103.254 | 1.033 | 1.24% | 0.0048% | no |
| `1e3` | 1103.17 | 1.103 | 1.97% | 0.0053% | no |
| `1e4` | 17187.6 | 1.719 | 14.4% | 0.175% | no at every boundary |
| `1e5` | `1.45235e7` | 145.2 | 7.14% | 0.335% | yes, at 1.76 Myr |

The `1e2` and `1e3 M_sun` results pass the pre-declared 5% boundary and grid
criteria. The `1e4 M_sun` terminal mass is boundary sensitive, but all
three boundaries end in the narrow scientific class `1.33e4-1.72e4 M_sun`;
none reaches even `1e5 M_sun`. The `1e5 M_sun` terminal mass retains the
known 7.14% boundary sensitivity, while all variants still robustly cross
`1e7 M_sun`.

## Delivery is not capture

The light-seed failure is not caused by a lack of central dark-matter
delivery. At the smallest feeding radius, the fluid supplies approximately

```text
seed 1e2:  1.42e4 M_sun supplied, 3.25 M_sun captured
seed 1e3:  5.28e4 M_sun supplied, 1.03e2 M_sun captured
seed 1e4:  6.74e6 M_sun supplied, 7.17e3 M_sun captured
```

The capture fractions are only `0.023%`, `0.195%`, and `0.106%`. Their gas
channels also remain Bondi limited for the full assembled interval, retaining
less than `1`, `0.11`, and `23 M_sun` of baryons, respectively. A central
baryonic potential can transport SIDM inward without guaranteeing that a
small black hole captures that supply.

## Nonlinear seed threshold

A 14-case refinement over `2e4-8e4 M_sun`, followed by nine boundary-plus-grid
cases over `3.25e4-3.75e4 M_sun`, finds a sharp transition:

- `2e4 M_sun` fails at both feeding radii, ending at
  `3.87e4-4.79e4 M_sun`.
- `3e4-3.75e4 M_sun` is boundary ambiguous. The smaller feeding radius
  enters direct capture and reaches `1e7 M_sun` at `1.85-1.88 Myr`; the
  larger radius remains below `9.59e4 M_sun`.
- At nominal `lambda=0.25`, seeds of `4e4 M_sun` and above cross
  `1e7 M_sun` at both radii. The `4e4 M_sun` crossing occurs at
  `1.71-1.85 Myr`.
- Double-resolution changes the non-runaway `3.25e4-3.75e4 M_sun` results by
  only `0.53%-0.55%`; the branch split is not a grid artifact.

The transition should not be quoted as one universal critical mass. A final
`lambda=0.20-0.30` sensitivity matrix gives:

| Seed [M_sun] | Final mass at `lambda=0.20` | Final mass at `lambda=0.30` | Classification |
|---:|---:|---:|---|
| `1e2` | 102.59 | 103.93 | robust failure |
| `1e3` | 1089.65 | 1113.26 | robust failure |
| `1e4` | 15587.6 | 19052.3 | robust failure |
| `3e4` | 58383.9 | 83770.3 | failure at the conservative radius |
| `4e4` | 85470.0 | `1.5423e7` | capture-model sensitive |

Thus the stellar-remnant light-seed conclusion is insensitive to the
calibrated `lambda` interval, while the intermediate-seed runaway threshold
is model dependent at a few times `1e4 M_sun`.

## Scientific conclusion

For this high-redshift halo and two-Myr physical-cooling baseline, baryon-aided
SIDM inflow does **not** turn `1e2-1e3 M_sun` stellar-remnant seeds into LRD
black holes. It also fails for a `1e4 M_sun` seed across all tested feeding
radii and capture coefficients. The rapid `1e7 M_sun` channel requires a
pre-existing intermediate or heavy seed near a few times `1e4 M_sun`, and
the exact threshold depends on unresolved capture physics.

This narrows the claim of the project. The robust positive result is rapid,
dark-dominated amplification of an already intermediate/heavy seed. The
robust negative result is that central SIDM delivery alone does not solve the
stellar-remnant light-seed problem on a two-Myr timescale.

## Limits and audit

- The reservoir is a one-zone subgrid capture model. It is mass conserving
  and calibrated to the resolved branch, but it does not resolve phase-space
  transport, angular momentum, or relativistic capture inside `r_feed`.
- The result applies to one optimized `M200=1e9 M_sun`, `z=30` halo and a
  constant SIDM cross section. It is not yet a population-level seed bound.
- Jobs `40775904`, `40776966`, `40777024`, `40778274`, `40778815`, and
  `40779138` completed all 59 tasks with `0:0` exits and empty error logs.
- The maximum mass-budget residual among accepted light-seed matrices is
  `5.49e-11`; all values remain far below the scientific mass differences.

![Light-seed boundary closure](results/stage5/figures/stage5_light_seed_boundary.png)

![Refined seed threshold](results/stage5/figures/stage5_seed_threshold.png)

![Capture sensitivity](results/stage5/figures/stage5_capture_sensitivity.png)
