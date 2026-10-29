# Output locations shared by analysis writers and cross-analysis readers.
results_root <- function() Sys.getenv("STARS_OUTPUT_DIR", unset="outputs")

result_path <- function(directory, ...) {
  path <- do.call(file.path, c(list(directory), list(...)))
  path <- gsub("\\\\", "/", path)
  regression_root <- paste0(gsub("\\\\", "/", results_root()), "/regressions")
  clustering_root <- paste0(gsub("\\\\", "/", results_root()), "/clustering")
  if (startsWith(path, paste0(regression_root, "/"))) {
    relative <- substring(path, nchar(regression_root)+2)
    # Explicit subdirectories (education analyses and alternate age reference)
    # already identify their analysis. Route shared regression writers by stem.
    parts <- strsplit(relative, "/", fixed=TRUE)[[1]]
    if (length(parts)==1 || parts[1]=="replicate_details") {
      name <- tail(parts,1)
      if (grepl("\\.[[:alnum:]]+$", name)) {
        folder <- if (grepl("^(age_sex_|age_effects_by_sex)",name)) "sensitivity/age_sex_interaction" else
          if (grepl("^(without_marital_|marital_adjustment_)",name)) "sensitivity/marital_adjustment" else
          if (grepl("^late_years",name)) "sensitivity/late_years" else
          if (grepl("^age_25_plus",name)) "sensitivity/age_25_plus" else
          if (grepl("^exclude_2021",name)) "sensitivity/exclude_2021" else
          if (grepl("^(minimal_available_|sample_selection_|retained_excluded_|selection_)",name)) "sensitivity/sample_selection" else
          if (grepl("^race_complete_",name)) "sensitivity/race_adjustment" else
          if (grepl("^(principal_model|principal_race_model|sequential_models|principal_sample_flow|principal_shared_sample|race_shared_sample)",name)) "principal" else
          "diagnostics"
        path <- file.path(regression_root, folder, relative)
      }
    }
  }
  if (startsWith(path, paste0(clustering_root, "/"))) {
    relative <- substring(path,nchar(clustering_root)+2)
    if (!grepl("/",relative,fixed=TRUE) && grepl("\\.[[:alnum:]]+$",relative)) {
      folder <- if(startsWith(relative,"county_")) "geographic" else
        if(startsWith(relative,"month_")) "monthly" else
        if(startsWith(relative,"weekday_")) "weekday" else "overview"
      path <- file.path(clustering_root,folder,relative)
    }
  }
  if (grepl("\\.[[:alnum:]]+$",basename(path)))
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  path
}
