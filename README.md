# SIDM black-hole growth with baryonic assembly

This repository contains the public, reproducible part of a one-dimensional
spherical study of black-hole growth in self-interacting dark-matter (SIDM)
haloes with an evolving Hernquist baryonic reservoir.  The calculation
combines conducting-fluid SIDM evolution, MC reconstruction, Roe fluxes,
Eddington-limited baryonic accretion, evolving Bondi gas conditions, optically
thin cooling, and a conservative small-seed capture closure.

The accompanying manuscript is:

> **Seed-mass dependence of baryon-assisted black-hole growth in
> self-interacting dark matter haloes**

The current study is a controlled local numerical model, not a cosmological
population calculation.  Its conclusions are restricted to the explored
parameter domain and to the stated inner-boundary and irreversible-feedback
closures.

## Repository contents

- `sidm_bh/`: simulation and physics modules.
- `scripts/`: analysis, validation, and figure-generation scripts.
- `tests/`: focused regression and numerical checks.
- `configs/`: compact configuration files used by the public workflows.
- `data/`: small tabulated inputs used by the model.
- `results/`: lightweight CSV/JSON summaries used by the manuscript.
- `paper/`: LaTeX source, bibliography, final figures, and the reviewed PDF.
- `*.md`: research notes and stage reports documenting the numerical study.

## Reproduction

The code targets Python 3.10+ with NumPy, SciPy, Matplotlib, and PyTest.  A
typical local setup is:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install numpy scipy matplotlib pytest
pytest -q
```

The manuscript-level `plot_paper_*.py` scripts read the compact summaries in
`results/` and write figures under `paper/figures/`.  Lower-level analysis
scripts document the full raw-data workflow and require the excluded snapshot
archive.  The manuscript can be built with a local MiKTeX or TeX Live
installation from `paper/` using `pdflatex`, `bibtex`, and two final
`pdflatex` passes.

## Data policy

The repository intentionally excludes raw time-series snapshots (`.npz` and
other data-heavy binary formats), local cluster submission files and package
wheels, generated caches, LaTeX intermediates, private correspondence, and
the source PDF used as the project prompt.  These products are not required to
inspect the code or the reported compact summaries.  No passwords, tokens,
SSH material, or cluster credentials are included.

The public summaries are sufficient to regenerate the manuscript figures, but
they are not a replacement for the excluded raw simulation snapshots.  The
paper does not interpret sampled success fractions as population probabilities.

## Citation and license

Please cite the accompanying manuscript and the references listed in
`paper/references.bib`.  No software license has been asserted yet; until the
authors specify one, the repository should be treated as publicly viewable but
not automatically reusable.
