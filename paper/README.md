# Manuscript build

The submission draft is `main.tex`; the compiled article is
`kang_sidm_baryon_accretion.pdf`.  The bibliography database contains 36
entries, all of which are cited in the current manuscript.  DOI and arXiv
metadata checks are recorded in `reference_audit.md`.

The manuscript was built on Windows with MiKTeX using:

```powershell
$env:Path = 'D:\CTEX\MiKTeX\miktex\bin\x64;' + $env:Path
pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex
bibtex.exe main
pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex
pdflatex.exe -interaction=nonstopmode -halt-on-error main.tex
```

All displayed figures are copied into `figures/` so that the paper builds
independently of the wider project directory.
