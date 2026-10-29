# Montana suicide analyses, 2012–2021

Code, aggregate data, tables, and figures for analyses of Montana suicide mortality,
individual characteristics, and county and temporal context.

Start with the [analysis index](docs/ANALYSIS_INDEX.md) for methods and links to results.
The [output inventory](outputs/INDEX.md) lists every generated file by analysis.

## Run the analyses

Run commands from the repository root. R 4.3.1 or later, the recommended `cluster`
package, and Python 3.10 or later are required. Python figure rendering uses Pillow.

The complete command is `python scripts/run_analyses.py aggregate` or
`python scripts/run_analyses.py full`. Supply `--rscript /path/to/Rscript` if R is
not on PATH. Add `--prepare` to rebuild clustering inputs from primary records and
source workbooks before running. The individual commands are listed below.

```sh
Rscript scripts/setup.R
python -m pip install -r requirements.txt -r requirements-figures.txt
Rscript scripts/reproduce.R aggregate
python scripts/build_county_profiles.py --png
python scripts/build_county_map_legend_only.py
Rscript scripts/validate.R aggregate
```

Aggregate mode fits county and temporal clustering and mortality comparisons from
bundled summary data and copies the descriptive summaries. It does not fit
individual regressions or recalculate descriptive summaries from primary records.

With the [primary inputs](docs/INPUTS.md) available:

```sh
Rscript scripts/reproduce.R full
python scripts/audit_demographic_coding.py
python scripts/build_county_profiles.py --png
python scripts/build_county_map_legend_only.py
Rscript scripts/validate.R full
python scripts/build_output_index.py
python scripts/check_integrity.py
```

`STARS_SOURCE_DIR` identifies the directory containing the primary-data folders;
its default is the repository's parent directory. `STARS_OUTPUT_DIR` changes the
output root; its default is `outputs`. Each command writes descriptive filenames
within analysis folders. CSV files contain either formatted presentation tables or
numeric estimates; `_long.csv` files contain numeric coefficients and exact p-values.

## Repository layout

The [three-category marital-status sensitivity](docs/MARITAL_THREE_CATEGORIES.md)
separates other, divorced and separated in parallel four-model sequences. It is
included in the full and regression workflows. After principal models are available,
run it separately with `Rscript scripts/run_marital_three_categories.R`.

| Directory | Contents |
| --- | --- |
| `R/` | Analysis modules, demographic recodes, covariance estimation, output paths |
| `scripts/` | Analysis commands, input preparation, figures, validation, indexes |
| `data/processed/clustering/` | Contextual features, county and temporal death counts and denominators |
| `data/processed/descriptive/` | Aggregate descriptive summaries |
| `data/external/census/` | Census source files and county boundary geometry |
| `outputs/descriptive/` | Source-cohort characteristics and missingness |
| `outputs/regressions/` | Principal models, sensitivity analyses, education comparisons, diagnostics |
| `outputs/clustering/` | Geographic, monthly and weekday clustering and mortality comparisons |
| `outputs/data_quality/` | Coding audit and 2020 comparisons |
| `docs/` | Analysis index, methods, input definitions and reproducibility instructions |

## Data access and reuse

STARS individual records and IPUMS microdata are not distributed. Obtain STARS
records through the study authors or data provider and ACS extracts through
[IPUMS USA](https://usa.ipums.org/). Aggregate data and generated results are
included for inspection and for analyses that do not require individual records.
Provider access and reuse terms apply; this repository does not grant rights to
the underlying datasets. Code and analysis attribution: Gabe Maayan.
