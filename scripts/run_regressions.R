#!/usr/bin/env Rscript
# Run year-matched regressions only, from the supplement root.
source("R/common.R")
output_dir <- results_root()
source("R/individual_regressions.R")
  local(source("R/marital_three_categories.R",local=TRUE))
  local(source("R/partial_education.R",local=TRUE))
  local(source("R/education_comparison.R",local=TRUE))
  local(source("R/age_reference_45_54.R",local=TRUE))
