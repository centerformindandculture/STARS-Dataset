#!/usr/bin/env Rscript
# Aggregate mode uses bundled summary data; full mode requires primary records.
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args))args[1] else "aggregate"
if(!mode %in% c("aggregate","full"))stop("Mode must be aggregate or full")
source("R/common.R")
output_dir <- results_root()
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))
source("R/contextual_clustering.R")
if(mode=="aggregate")source("R/aggregate_descriptive_statistics.R")
if(mode=="full") {
  source("R/descriptive_statistics.R")
  source("R/individual_regressions.R")
  local(source("R/marital_three_categories.R",local=TRUE))
  local(source("R/partial_education.R",local=TRUE))
  local(source("R/education_comparison.R",local=TRUE))
  local(source("R/age_reference_45_54.R",local=TRUE))
  source("R/compare_2020.R")
}
message("Completed ",mode," analyses in ",output_dir)
