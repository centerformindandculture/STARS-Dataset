# Three-category marital status

Four-model sequences use other recorded status (reference), divorced and separated
in Model 4. STARS DIVORCED and MARRIED (SEPARATED) define the two exposure groups;
IPUMS MARST 4 and 3 define the corresponding groups. Other contains married but
not separated, widowed and never-married/single records. Missing status stays missing.

The sequences parallel the principal regressions, with and without race/ethnicity.
Adults must have recorded Montana residence and meet the principal complete-case
criteria. All columns use identical records within each sequence. Non-race models
cover 2013–2019 and 2021; race-inclusive models cover 2017–2019 and 2021. Models 1–3
are checked against their principal counterparts; Model 4 replaces the combined
divorced/separated indicator with separate divorced and separated indicators.

Every model adjusts for categorical year. Model 4 also includes sex, age, veteran
status, farmer/rancher occupation and associate-degree attainment. Race/ethnicity
is included in every model of the race-inclusive sequence. Uncertainty combines
Poisson mortality-event and annual ACS SDR covariance with 80 signed person
replicate weights per annual sample, following [the regression inference method](REPLICATE_INFERENCE.md).

The divorced-versus-separated contrast tests equality of the two Model 4 log-odds
coefficients with a 1-df Wald test. The contrast odds ratio divides the divorced
odds ratio by the separated odds ratio. Raw and Holm-adjusted p-values are reported;
the two race specifications form one two-test family. Coefficient significance
stars and confidence intervals remain pointwise. Separate recorded categories
do not establish that either status was classified consistently across years.

Run `Rscript scripts/run_marital_three_categories.R` after the principal models,
or use `Rscript scripts/run_regressions.R` or the full analysis workflow.
Outputs are in `outputs/regressions/sensitivity/marital_three_categories/`:

- `sequential_marital_three_categories.csv`: four models without race.
- `sequential_marital_three_categories_with_race.csv`: four models with race.
- Corresponding `_long.csv` and `_notes.txt`: numeric estimates and table definitions.
- `divorced_separated_contrast.csv`: direct contrasts and equality tests.
- `marital_category_counts_by_year.csv`: source/year/category counts and weights.
- `principal_sample_comparison_checks.csv`: coefficient and uncertainty agreement for Models 1–3.
- `analysis_report.txt`: specification and contrast summary.
- Model, replicate and sample-accounting files, including covariance components.
