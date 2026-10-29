source("R/common.R")
origin <- file.path(results_root(),"descriptive")
destination <- "data/processed/descriptive"
checks <- read_table(file.path(origin,"validation_checks.csv"))
stopifnot(all(checks$passed))
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
files <- list.files(origin,pattern="\\.(csv|txt)$",full.names=TRUE)
stopifnot(length(files)>20,all(file.copy(files,destination,overwrite=TRUE)))
cat("Exported",length(files),"source-cohort aggregate descriptive files.\n")
