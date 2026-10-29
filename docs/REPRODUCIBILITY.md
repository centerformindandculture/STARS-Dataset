# Running and checking the analyses

Run commands from the repository root. Set `STARS_SOURCE_DIR` to the directory
containing the [primary-data folders](INPUTS.md) if they are not in the parent.
Set `STARS_OUTPUT_DIR` for a separate output destination.

| Command | Analyses |
| --- | --- |
| `Rscript scripts/reproduce.R aggregate` | Clustering and mortality comparisons from aggregate data; descriptive summary copies |
| `Rscript scripts/reproduce.R full` | Aggregate analyses plus source-cohort descriptives, individual regressions and 2020 comparisons |
| `Rscript scripts/run_clustering.R` | All three contextual clustering domains and mortality comparisons |
| `Rscript scripts/run_descriptives.R` | Clustering assignments and source-cohort descriptives |
| `Rscript scripts/run_regressions.R` | Principal, sensitivity, education and alternate-reference regressions |
| `Rscript scripts/compare_2020.R` | 2020 distribution and missingness comparisons |
| `python scripts/audit_demographic_coding.py` | Dataset-key and demographic-label audit |
| `python scripts/build_county_profiles.py --png` | County comparison CSV/HTML and map SVG/PNG |
| `python scripts/build_county_map_legend_only.py` | County map with enlarged legend and no other text |

The partial-education, joint-education and alternate-age-reference commands also
have individual runners under `scripts/`; they require the principal regression
outputs for sample/covariance consistency checks. The regression runner executes
them in dependency order.

`Rscript scripts/run_marital_three_categories.R` fits parallel four-model sequences
with separate divorced and separated categories after principal regression outputs
are available. The regression and full workflows include this analysis. See
[category definitions and contrast tests](MARITAL_THREE_CATEGORIES.md).

Validate with `Rscript scripts/validate.R aggregate` or `full`. Checks cover cohort
accounting, county coverage, rate denominators, fit convergence, resident and adult
eligibility, supported years, fixed samples, all 80 annual ACS refits, signed weights,
joint-education equivalence, Holm adjustments, and alternate-reference transformations.

Run `python scripts/build_output_index.py` to generate the output index and file
inventory. `python scripts/check_integrity.py` verifies the release checksums.
Use its `--write` option after preparing a release. Checksums cover code,
documentation, bundled aggregate inputs, and generated outputs, excluding the
checksum manifest itself and author-maintained document files.

Aggregate mode does not re-estimate individual regressions or validate them.
Full individual-level reproduction requires restricted primary records and the
specified ACS extract. Census/ACS input estimates, predictor measurement, and
cluster construction have uncertainties beyond those included in model intervals.
