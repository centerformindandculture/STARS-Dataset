# Partial-attainment models

These four-model sequences replace the Model 4 degree indicator with partial high
school or college attainment. STARS SOME HIGH SCHOOL and SOME COLLEGE and IPUMS
EDUCD 30/40/50/61/65/71 define partial attainment. Other recorded attainment,
including below-grade-9 schooling, high-school/GED and degrees, defines the reference.
NIU/unknown education stays missing. The indicator does not identify dropout or
distinguish current enrollment from discontinued education.

Principal and partial-attainment sequences use identical eligible samples;
Models 1–3 are checked for coefficient and standard-error equality. Race-inclusive
sequences require recorded race/ethnicity. Each model uses categorical year and
combined mortality-event/ACS covariance. Source classification uncertainty remains.

Run `Rscript scripts/run_regressions.R` or, after principal fits,
`Rscript scripts/run_partial_education.R`. Results are in
`outputs/regressions/partial_education`. Numeric `_long.csv` files accompany
formatted sequential tables, sample checks, crosswalks and covariance components.
