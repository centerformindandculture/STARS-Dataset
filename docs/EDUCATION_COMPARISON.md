# Joint education-indicator comparison

Fully adjusted models include both partial-attainment and associate-degree-or-higher
indicators, which form three mutually exclusive education groups. The reference
is other recorded below-associate attainment: no high school/below grade 9 or
high-school/GED. The same samples support both binary models and the joint model.

Two 1-df coefficient Wald tests evaluate whether partial attainment adds information
beyond the degree split and whether degree attainment adds information beyond the
partial split. Tests use combined mortality-event and annual ACS SDR covariance.
Holm adjustment covers the two incremental tests within each race specification;
race-inclusive tests form a separate sensitivity family. The partial-versus-degree
contrast has a pointwise confidence interval.

The indicators share a reference in the joint model, so their odds ratios differ
from the binary models' contrasts. The tests assess statistical information in
recorded categories; they do not establish measurement validity or causal effects.
An equivalent three-level factor verifies the fitted linear predictors.

Run `Rscript scripts/run_education_comparison.R` after principal and partial fits,
or use the complete regression runner. Results are in
`outputs/regressions/education_comparison`. Start with
`education_incremental_tests_formatted.csv`; exact p-values are in
`education_incremental_tests.csv` and group contrasts in `education_group_contrasts.csv`.
