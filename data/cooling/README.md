# Grackle Cloudy Cooling Data

`CloudyData_noUVB.h5` is copied unchanged from the official Grackle data
repository at submodule commit
`928696482fbe15d9bac4382de6134d95568f099c`:

https://github.com/grackle-project/grackle_data_files/blob/928696482fbe15d9bac4382de6134d95568f099c/input/CloudyData_noUVB.h5

SHA-256:

```text
0abe25cceeb5c0825381c5f17059982a9a2cdd27ce369a475c559fba6a8fa106
```

The table contains no-UV-background Cloudy equilibrium cooling coefficients
for primordial gas and the solar-metal increment, plus the equilibrium mean
molecular weight. Its axes cover `10 <= T/K <= 1e9` and
`-10 <= log10(n_H/cm^-3) <= 4`.

The project interpolates each cooling coefficient logarithmically, adds the
metal contribution linearly in `Z/Z_sun`, and clamps queries to the tabulated
density and temperature bounds. This is an optically thin equilibrium
closure, not a molecular network, self-shielding model, or radiative-transfer
calculation.
