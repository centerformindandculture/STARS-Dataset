# Inputs and provenance

Primary paths are relative to `STARS_SOURCE_DIR` (default: the repository's parent).

| File | Use |
| --- | --- |
| `R - Data Exploration/20230831_Montana_Suicide_Database.csv` | 2,764 source records, 2012–2021; descriptives, mortality numerators, regressions and 2020 comparisons |
| `IPUMSData/ipums_2012-2019_2021_with_replicate_weights.csv` | 91,968 Montana ACS person records, 2012–2019 and 2021 |
| `IPUMSData/ipums_2012-2019_2021_with_replicate_weights.cbk.txt` | Extract sample and variable definitions |
| `20230528 MT Dataset Key.docx` | Demographic coding audit |

The IPUMS extract contains `YEAR`, `SAMPLE`, `STATEFIP`, `SERIAL`, `PERNUM`, `SEX`,
`AGE`, `MARST`, `RACE`, `HISPAN`, `EDUCD`, `VETSTAT`, `PERWT`, `OCC2010`, the
replicate-weight presence flag `REPWTP`, and weights `REPWTP1`–`REPWTP80`.
Each annual sample has `SAMPLE = YEAR*100+1` and `STATEFIP = 30`. The 2020 ACS
sample is excluded. Input checksums accompany regression and descriptive outputs.

STARS records were supplied by the Montana Department of Public Health and Human
Services. The pipeline begins with the supplied CSV; upstream manual cleaning
cannot be reconstructed from the available source files. No individual records,
identifiers or narrative fields are included in the repository.

## Context and population denominators

County features use 2007–2011 ACS five-year tables B01002, B02001, B19013, B17001,
B11001, B12001 and B15002. Uninsured percentage uses
`CensusCountyData/ACSST5Y2012_HealthInsuranceCoverage_County.xlsx`, sheet `DataValues`,
columns `County` and `PercentUninsured` (2008–2012 ACS). Estimates and margins of
error are retained in `data/processed/clustering/county_acs_source_cells.csv`.

Annual county populations and age/sex denominators use Census revised intercensal
estimates for 2012–2019 and vintage-2021 estimates for 2020–2021. URLs and checksums
are recorded in the clustering source manifests. County boundaries are 2020
U.S. Census Bureau TIGERweb geometry in NAD83 / Conus Albers (EPSG:5070).

Temporal context comes from the `Data` sheets in these workbooks:

| File | Measures |
| --- | --- |
| `MonthData/weatherresults.xlsx` | Freezing/hot days, precipitation, average/maximum/minimum temperature, humidity |
| `MonthData/CrimeByMonth2007.xlsx` | Monthly 2007 crime index |
| `DayData/carcrashes.xlsx` | Weekday crash proportions |
| `DayData/substancedeaths.xlsx` | Weekday substance-related mortality measure |
| `DayData/happiness.xlsx` | Digitized weekday happiness measure |

Workbook provenance identifies NOAA observations for 2012–2021,
[monthly Montana humidity](https://www.currentresults.com/Weather/Montana/humidity-january.php),
[Montana crime statistics](https://mbcc.mt.gov/_docs/Data/CIM/CIM-2006-07.pdf),
[NSC crash proportions](https://injuryfacts.nsc.org/motor-vehicle/overview/crashes-by-time-of-day-and-day-of-week/),
[CDC WONDER](https://wonder.cdc.gov/) Montana deaths in 2012–2020 with ICD-10 F10–F19,
and [Dodds et al.'s weekday happiness figure](https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0026752.g005).
These links document provenance; the supplied workbook values determine the analyses.
Weather measures retain their stored units. The precipitation unit and physical
temperature thresholds require caution because source notes are not fully consistent.
A constant unit conversion would not change z-standardized clustering.

## Prepare aggregate inputs

```sh
python scripts/prepare_clustering_inputs.py
python scripts/prepare_county_adjusted.py
Rscript scripts/run_descriptives.R
Rscript scripts/export_descriptive_aggregates.R
```

Preparation requires the primary CSV and source workbooks. Census downloads are
cached under `data/external/census`; a network connection is needed only for
uncached downloads. `prepare_temporal_context.py` exports external measures
without using suicide outcomes. Source files are read without modification.
