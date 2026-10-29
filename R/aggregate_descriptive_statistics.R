# Public snapshots of final CSV-source descriptives; raw reconstruction is full mode.
cd_from <- "data/processed/descriptive"
cd_to <- result_path(output_dir,"descriptive")
dir.create(cd_to,recursive=TRUE,showWarnings=FALSE)
cd_files <- list.files(cd_from,pattern="\\.(csv|txt)$",full.names=TRUE)
stopifnot(length(cd_files)>20,all(file.copy(cd_files,cd_to,overwrite=TRUE)))
cd_verify <- read_table(file.path(cd_to,"validation_checks.csv"))
stopifnot(all(cd_verify$passed))
