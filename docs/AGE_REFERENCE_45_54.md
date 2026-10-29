# Age reference: 45–54 years

The principal and partial-attainment sequences, with and without race, are also
expressed relative to ages 45–54. This is an exact linear reparameterization:
the intercept, age coefficients, complete covariance components and annual ACS
replicate vectors are transformed together. Records, fitted probabilities and
non-age coefficients do not change. Each age test compares that group with 45–54.

Run `Rscript scripts/run_age_reference_45_54.R` after the principal and partial fits,
or use the complete regression runner. Results are in
`outputs/regressions/age_reference_45_54`. The four formatted sequential tables have
`_age45_54.csv` suffixes; `_long.csv` files retain numeric estimates. The
`reference_change_checks.csv` file verifies all 16 model transformations and
recomputes the SDR covariance from transformed replicate coefficients.
