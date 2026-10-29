# Analysis variables and eligibility

| Variable | STARS | IPUMS ACS |
| --- | --- | --- |
| Sex | `DEM_GEN2GenderNum`: 0 male, 1 female | `SEX`: 1 male, 2 female |
| Age | `DEM_AGE1Age` | `AGE` |
| Veteran | `DEM_OCC5ArmedForcesEver`: 0 no, 1 yes | `VETSTAT`: 1 no, 2 yes; other codes missing |
| Farmer/rancher | `DEM_OCC3RancherFarmer`: 0 other, 1 farmer | `OCC2010`: 205 farmers/ranchers/other agricultural managers, other codes comparison |
| Degree attainment | ASSOCIATES/BACHELORS/MASTERS/DOCTORATE | `EDUCD`: 81/101/114/115/116 |
| Partial attainment | SOME HIGH SCHOOL/SOME COLLEGE | `EDUCD`: 30/40/50/61/65/71 |
| Divorced/separated | DIVORCED/MARRIED (SEPARATED) | `MARST`: 3/4 |
| Race/ethnicity | WHITE/AMERICAN INDIAN/OTHER; Hispanic flag enters Other | `RACE`: 1 White, 3 American Indian, other codes Other; `HISPAN` 1–4 enters Other |
| Year | `DEATH_2Year` | `YEAR` |
| Weight | 1 | `PERWT` |

Race/ethnicity combines heterogeneous categories and is not a race-only measure.
Education and marital-status labels show substantial recording discontinuities.
Missing or unsupported codes remain missing. `EDUCD` 0/1 are not treated as attainment.
The partial-attainment reference contains all other recorded education levels;
the degree reference contains all below-associate levels.

Regressions require recorded Montana residence (`DEM_RES09State = MONTANA`) and
age 18+. Age categories are 18–24, 25–34, 35–44, 45–54, 55–64, 65–74, 75–84 and 85+.
The reference is 18–24 except for age-25+ and age-45–54-reference analyses.
All sequential models use the Model 4 complete-case sample within each sequence.
Years must contain complete-case deaths and ACS comparisons; unsupported years
are excluded from both sources. Principal years are 2013–2019 and 2021;
race-inclusive years are 2017–2019 and 2021.

Descriptive summaries cover all 2,764 source records, including nonresident and
unknown-residence records. County mortality uses 2,650 county-matched Montana
resident deaths. Temporal mortality uses the same county-resident definition
within years with complete date recording.

County clustering feature definitions are in
`data/processed/clustering/county_feature_definitions.csv`. Death counts and
suicide rates do not enter county, month or weekday clustering.
