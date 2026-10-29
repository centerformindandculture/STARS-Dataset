# Computing environment

The analysis runs with R 4.3.1 or later and the recommended `cluster` package.
Regression fitting and replicate-weight covariance calculations use base R.
Session files in the output directories record the installed R and package versions.

Python 3.10 or later is required. Workbook preparation uses openpyxl 3.1.5;
PNG map rendering uses Pillow 12.3.0. CSV, SVG and HTML profile generation uses
Python's standard library. Install dependencies with:

```sh
Rscript scripts/setup.R
python -m pip install -r requirements.txt -r requirements-figures.txt
```

R uses UTC, C-locale collation, and Mersenne-Twister/Inversion/Rejection random
number generators. Clustering seeds are explicit in `R/contextual_clustering.R`.
Full-sample logistic fits use tolerance 1e-12 and at most 100 iterations.
Replicate score fits use convergence and score-residual checks.
