# Supplementary adult regressions replacing broad attainment with partial schooling.
# Run through scripts/run_partial_education.R; no principal outputs are overwritten.
source("R/common.R",local=TRUE)
output_dir <- file.path(results_root(),"regressions","partial_education")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source("R/regression_education_inputs.R",local=TRUE)
source("R/replicate_inference.R",local=TRUE)
sample_counts <- missingness <- model_checks <- model_coverage <- model_warnings <- list()
source("R/regression_helpers.R",local=TRUE)
labels <- c(labels,EDUC_PARTIALPARTIAL="Partial high school or college (reference: other attainment)")
base_variables <- c("SUIC","SEX","AGE_ADULT","YEAR","VETSTAT","FARM","EDUC_PARTIAL","DIVSEP","weight")
residence_name <- "DEM_RES09State"
require_columns(raw,residence_name,"STARS source")
# Residence uses the same source field audited in the descriptive workflow.
stopifnot(length(residence_name)==1)
residence <- raw[[residence_name]]
for (race in c(FALSE,TRUE)) {
  variables <- c(base_variables,if(race) "RACE")
  complete <- complete.cases(ipums_stars[variables])
  available <- table(ipums_stars$YEAR[complete],factor(ipums_stars$SUIC[complete],levels=0:1))
  supported <- rownames(available)[rowSums(available>0)==2]
  common <- droplevels(ipums_stars[complete & as.character(ipums_stars$YEAR) %in% supported,])
  stem <- if(race) "sequential_partial_education_with_race" else "sequential_partial_education"
  terms <- c("SEX","AGE_ADULT",if(race) "RACE")
  formulas <- list(reformulate(terms,"SUIC"),reformulate(c(terms,"VETSTAT"),"SUIC"),
    reformulate(c(terms,"VETSTAT","FARM"),"SUIC"),
    reformulate(c(terms,"VETSTAT","FARM","EDUC_PARTIAL","DIVSEP"),"SUIC"))
  models <- setNames(lapply(seq_along(formulas),function(i)
    fit_logistic(formulas[[i]],common,paste0(stem,"_model",i))),paste0("Model_",1:4))
  publish_models(models,stem,race=race,specification=c(
    if(race) "Model 1: sex, age, race/ethnicity, categorical year." else "Model 1: sex, age, categorical year.",
    "Model 2 adds veteran status; Model 3 adds farmer/rancher occupation;",
    "Model 4 adds partial high school/college and combined divorced/separated."))
  notes_path <- result_path(output_dir,paste0(stem,"_notes.txt"))
  notes <- readLines(notes_path,warn=FALSE)
  notes <- notes[!grepl("^Attainment:|^Some college without",notes)]
  writeLines(c(notes,
    "Supplementary models use partial high school or college as the education indicator.",
    "Partial: STARS SOME HIGH SCHOOL or SOME COLLEGE; IPUMS EDUCD 30/40/50/61/65/71.",
    "Reference: all other recorded education, including no high school/below grade 9, high school/GED and degrees.",
    "NIU/unknown education remains missing; partial does not establish dropout, enrollment or causation.",
    "These tables share the confirmed-Montana-resident principal cohorts to isolate the education definition.",
    "Known non-Montana and unknown-residence deaths are excluded before complete-case selection."),notes_path)
  counts <- table(factor(ifelse(is.na(residence[common$SOURCE_ROW[common$SUIC==1]]),"UNKNOWN",
    ifelse(residence[common$SOURCE_ROW[common$SUIC==1]]=="MONTANA","MONTANA","NON_MONTANA")),
    levels=c("MONTANA","NON_MONTANA","UNKNOWN")))
  stopifnot(counts["NON_MONTANA"]==0,counts["UNKNOWN"]==0)
  write_table(data.frame(residence=names(counts),cases=as.integer(counts)),paste0(stem,"_residence_counts.csv"))
  checks <- lapply(1:3,function(i) {
    old_name <- if(race) paste0("principal_race_model",i) else paste0("principal_model",i)
    path <- result_path(dirname(output_dir),paste0(old_name,".csv"))
    if(!file.exists(path))return(data.frame(model=i,verified=FALSE,reason="Principal coefficient file absent"))
    old <- read_table(path);b <- coef(models[[i]])
    err <- max(abs(old$coefficient-b[match(old$term,names(b))]))
    se_err <- max(abs(old$std_error-sqrt(diag(vcov(models[[i]])))[match(old$term,names(b))]))
    stopifnot(is.finite(err),err<1e-8,is.finite(se_err),se_err<1e-8)
    data.frame(model=i,verified=TRUE,reason=paste("Maximum coefficient difference",format(err,scientific=TRUE),"; SE difference",format(se_err,scientific=TRUE)))
  })
  write_table(do.call(rbind,checks),paste0(stem,"_comparison_checks.csv"))
}
write_table(as.data.frame(table(year=raw$DEATH_2Year,original=raw$DEM_EDU1EducationAchieved,
  partial=all_stars$EDUC_PARTIAL)),"stars_partial_crosswalk_counts.csv")
write_table(as.data.frame(table(year=ipums_raw$YEAR,EDUCD=ipums_raw$EDUCD,
  partial=ipums$EDUC_PARTIAL,attainment=ipums$EDUC_ATTAIN)),"ipums_partial_crosswalk_counts.csv")
for (item in c("sample_counts","missingness","model_checks","model_coverage","model_warnings"))
  write_table(do.call(rbind,get(item)),paste0(item,".csv"))
write_table(do.call(rbind,replicate_checks),"replicate_model_checks.csv")
write_table(do.call(rbind,replicate_changes),"covariance_method_comparison.csv")
stopifnot(all(do.call(rbind,model_checks)$passed),all(do.call(rbind,replicate_checks)$passed))
paths <- c(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"),ipums_path,codebook_path)
write_table(data.frame(input=paths,md5=unname(tools::md5sum(paths))),"input_checksums.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
message("Partial education results: ",output_dir)
