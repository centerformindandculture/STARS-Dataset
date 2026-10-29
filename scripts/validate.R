#!/usr/bin/env Rscript
# Validate aggregate analyses or the complete primary-data workflow.
source("R/common.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "aggregate"
if (!mode %in% c("aggregate", "full")) stop("Mode must be aggregate or full")
output_dir <- results_root()
checks <- list()
check <- function(name, value) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(value))
}
outcome_checks <- read_table(result_path(file.path(output_dir,"clustering"),"validation_checks.csv"))
check("outcome-free analysis checks pass",nrow(outcome_checks)==28 && all(outcome_checks$passed))
outcome_selection <- read_table(result_path(file.path(output_dir,"clustering"),"cluster_selection.csv"))
check("three clustering domains",setequal(outcome_selection$domain,c("county","month","weekday")))
check("selected k in prespecified range",all(outcome_selection$selected_k %in% 3:4))
for (k in 3:4) {
  mortality <- read_table(result_path(file.path(output_dir,"clustering"),paste0("county_k",k,"_mortality.csv")))
  check(paste("county",k,"all counties counted"),sum(mortality$counties)==56)
  check(paste("county",k,"deaths reconcile"),sum(mortality$DEATHS)==2650)
  check(paste("county",k,"rates reconcile"),
    all(abs(mortality$rate_per_100000-mortality$DEATHS/mortality$person_years*1e5)<1e-8))
}
adjusted_checks <- read_table(result_path(file.path(output_dir,"clustering"),"county_adjusted_validation.csv"))
check("adjusted county analyses pass",nrow(adjusted_checks)==2 && all(adjusted_checks$passed) && all(adjusted_checks$counties==56))
current_checks <- read_table(result_path(output_dir,"descriptive/validation_checks.csv"))
check("source-cohort descriptives pass",all(current_checks$passed))
if (mode == "full") {
  comparison2020 <- read_table(result_path(output_dir,"data_quality","comparison_2020","validation_checks.csv"))
  check("2020 descriptive comparison checks",all(comparison2020$passed))

  matched <- result_path(output_dir,"regressions")
  models <- read_table(result_path(matched,"model_checks.csv"))
  reps <- read_table(result_path(matched,"replicate_model_checks.csv"))
  check("all regression replicate checks pass",all(reps$passed) && all(reps$replicates_per_year==80) &&
    all(reps$refits==80*reps$years) && all(reps$maximum_score_residual<1e-5))
  check("case weights fixed during ACS replication",all(reps$case_weights_fixed))
  check("negative replicate weights retained",all(reps$negative_replicates_preserved))
  check("regression covariance is two-source",all(grepl("ACS SDR",models$variance_method,fixed=TRUE)))

  expected <- c("principal_model1","principal_model2","principal_model3",
    "principal_model4","race_complete_model_without_race","principal_race_model4",
    "exclude_2021_model")
  check("all IPUMS models present",all(c(expected,
    "principal_race_model1","principal_race_model2","principal_race_model3",
    "without_marital_model","without_marital_model_with_race",
    "age_sex_interaction_model","age_sex_interaction_model_with_race",
    paste0("late_years_model",1:4),paste0("age_25_plus_model",1:4),
    paste0("age_25_plus_with_race_model",1:4),"minimal_available_model") %in% models$model))
  check("all regression records are adults",all(models$minimum_age>=18))
  check("all regression deaths are confirmed Montana residents",all(models$residence_verified) &&
    all(models$confirmed_montana_cases==models$cases) && all(models$non_montana_cases==0) &&
    all(models$unknown_residence_cases==0))
  for (i in seq_len(nrow(models))) check(paste(models$model[i],"fit checks"),models$passed[i])
  shared <- read_table(result_path(matched,"principal_shared_sample_checks.csv"))
  check("four primary models use identical records",nrow(shared)==4 &&
    all(shared$same_records) && length(unique(shared$n))==1)
  main_checks <- models[models$model %in%
    c("principal_model1","principal_model2","principal_model3","principal_model4"),]
  check("four primary model case counts identical",length(unique(main_checks$cases))==1)
  check("four primary model comparator counts identical",length(unique(main_checks$comparison_records))==1)
  estimates <- read_table(result_path(matched,"principal_model4.csv"))
  check("broader attainment replaces partial education",
    "EDUC_ATTAINASSOCIATE_OR_HIGHER" %in% estimates$term &&
    !any(grepl("^EDUCPARTIAL",estimates$term)))
  check("divorced and separated grouped in main model",
    "DIVSEPDIVORCED/SEPARATED" %in% estimates$term && !any(grepl("^MARST2",estimates$term)))
  marital_dir <- file.path(matched,"sensitivity","marital_three_categories")
  marital_models <- read_table(file.path(marital_dir,"model_checks.csv"))
  marital_reps <- read_table(file.path(marital_dir,"replicate_model_checks.csv"))
  marital_agreement <- read_table(file.path(marital_dir,"principal_sample_comparison_checks.csv"))
  marital_contrast <- read_table(file.path(marital_dir,"divorced_separated_contrast.csv"))
  check("eight three-category marital models pass with resident adults",nrow(marital_models)==8 &&
    all(marital_models$passed) && all(marital_models$residence_verified) && all(marital_models$minimum_age>=18))
  check("three-category marital ACS replicates pass",nrow(marital_reps)==8 && all(marital_reps$passed) &&
    all(marital_reps$refits==80*marital_reps$years) && all(marital_reps$negative_replicates_preserved))
  check("marital models 1-3 agree with principal sequences",nrow(marital_agreement)==6 &&
    all(marital_agreement$passed) && all(marital_agreement$coefficient_error<1e-8) &&
    all(marital_agreement$standard_error_error<1e-8))
  check("two divorced-separated contrasts with Holm adjustment",nrow(marital_contrast)==2 &&
    all(marital_contrast$df==1) && all(is.finite(marital_contrast$p)) &&
    max(abs(marital_contrast$p_holm-p.adjust(marital_contrast$p,"holm")))<1e-12)
  for(race in c(FALSE,TRUE)) {
    stem <- paste0("sequential_marital_three_categories",if(race)"_with_race")
    coefficients <- read_table(file.path(marital_dir,paste0(stem,"_long.csv")))
    category_terms <- coefficients[grepl("^MARST2",coefficients$term),]
    check(paste(stem,"separate categories enter only model4"),nrow(category_terms)==2 &&
      setequal(category_terms$term,c("MARST2DIVORCED","MARST2SEPARATED")) &&
      all(category_terms$model=="Model_4") && !any(grepl("^DIVSEP",coefficients$term)))
    model_names <- paste0(stem,"_model",1:4)
    sample <- marital_models[match(model_names,marital_models$model),]
    reference <- models[models$model==if(race)"principal_race_model4" else "principal_model4",]
    check(paste(stem,"uses principal cohort counts"),nrow(reference)==1 &&
      all(sample$cases==reference$cases) && all(sample$comparison_records==reference$comparison_records))
  }
  partial_dir <- result_path(matched,"partial_education")
  pm <- read_table(file.path(partial_dir,"model_checks.csv"))
  pr <- read_table(file.path(partial_dir,"replicate_model_checks.csv"))
  check("eight supplementary partial-education models pass",nrow(pm)==8 && all(pm$passed))
  check("partial-education deaths are confirmed Montana residents",all(pm$residence_verified) &&
    all(pm$confirmed_montana_cases==pm$cases) && all(pm$non_montana_cases==0) && all(pm$unknown_residence_cases==0))
  check("partial-education SDR checks pass",nrow(pr)==8 && all(pr$passed) &&
    all(pr$refits==80*pr$years) && all(pr$negative_replicates_preserved))
  for (race in c(FALSE,TRUE)) {
    stem <- paste0("sequential_partial_education",if(race) "_with_race" else "")
    pc <- read_table(file.path(partial_dir,paste0(stem,"_comparison_checks.csv")))
    pl <- read_table(file.path(partial_dir,paste0(stem,"_long.csv")))
    pw <- read_table(file.path(partial_dir,paste0(stem,".csv")))
    check(paste(stem,"models 1-3 reproduce principal estimates and SEs"),nrow(pc)==3 && all(pc$verified))
    check(paste(stem,"education enters model 4 only"),
      identical(as.character(pl$model[pl$term=="EDUC_PARTIALPARTIAL"]),"Model_4") &&
      !any(grepl("^EDUC_ATTAIN",pl$term)))
    check(paste(stem,"finite estimates and valid intervals"),all(is.finite(pl$p)) &&
      all(pl$p>=0 & pl$p<=1) && all(pl$ci_lower>0 & pl$ci_upper>=pl$ci_lower))
    check(paste(stem,"same counts across four columns"),
      length(unique(unlist(pw[pw$Predictor=="STARS cases",-1])))==1 &&
      length(unique(unlist(pw[pw$Predictor=="IPUMS comparisons",-1])))==1)
    notes <- readLines(file.path(partial_dir,paste0(stem,"_notes.txt")),warn=FALSE)
    check(paste(stem,"resident eligibility documented"),any(grepl("Only deaths with recorded MONTANA",notes,fixed=TRUE)))
    rc <- read_table(file.path(partial_dir,paste0(stem,"_residence_counts.csv")))
    check(paste(stem,"no nonresident or unknown-residence cases"),
      all(rc$cases[rc$residence!="MONTANA"]==0))
  }
  age45_dir <- result_path(matched,"age_reference_45_54")
  age45 <- read_table(file.path(age45_dir,"reference_change_checks.csv"))
  check("16 age45-reference transformations preserve fits",nrow(age45)==16 && all(age45$passed) &&
    all(age45$age_reference=="45-54") && all(age45$linear_predictor_error<1e-10) &&
    all(age45$sdr_covariance_error<1e-10) && all(age45$non_age_coefficient_error<1e-10))
  for(stem in c("sequential_models","sequential_models_with_race",
      "sequential_partial_education","sequential_partial_education_with_race")) {
    tab <- read_table(file.path(age45_dir,paste0(stem,"_age45_54.csv")))
    num <- read_table(file.path(age45_dir,paste0(stem,"_age45_54_long.csv")))
    check(paste(stem,"age45 reference and seven contrasts"),
      all(unlist(tab[tab$Predictor=="Age reference",-1])=="45-54") &&
      !any(num$term=="AGE_ADULT45-54") && sum(grepl("^AGE_ADULT",num$term))==28 &&
      sum(num$term=="AGE_ADULT18-24")==4)
  }
  education_dir <- result_path(matched,"education_comparison")
  em <- read_table(file.path(education_dir,"model_checks.csv"))
  er <- read_table(file.path(education_dir,"replicate_model_checks.csv"))
  ec <- read_table(file.path(education_dir,"comparison_checks.csv"))
  et <- read_table(file.path(education_dir,"education_incremental_tests.csv"))
  check("joint education models pass with confirmed residents",nrow(em)==2 && all(em$passed) &&
    all(em$residence_verified) && all(em$non_montana_cases==0) && all(em$unknown_residence_cases==0))
  check("joint education replicate inference passes",nrow(er)==2 && all(er$passed) &&
    all(er$refits==80*er$years) && all(er$negative_replicates_preserved))
  check("joint education samples and equivalent factor agree",nrow(ec)==2 && all(ec$passed) &&
    all(ec$same_samples_as_both_binary_models) && all(ec$current_input_checksums_match) &&
    all(ec$three_group_linear_predictor_error<1e-8))
  check("four incremental education tests present",nrow(et)==4 && all(et$df==1) &&
    all(table(et$analysis)==2) && all(is.finite(et$chi_square)) && all(et$chi_square>=0))
  for(scope in unique(et$analysis)) {
    ed <- et[et$analysis==scope,]
    check(paste(scope,"education Holm adjustment"),max(abs(ed$p_holm-p.adjust(ed$p,"holm")))<1e-12)
  }
  input <- read_table(result_path(matched,"input_checksums.csv"))
  check("main model uses designated August31 CSV",
    any(grepl("20230831_Montana_Suicide_Database.csv",input$input,fixed=TRUE)) &&
    !any(grepl(".sav",input$input,fixed=TRUE)))
  counts <- read_table(result_path(matched,"model_sample_counts_by_year.csv"))
  years <- c(2012:2019,2021)
  main <- counts[counts$model=="principal_model4",]
  check("main model has both sources for supported common-sample years",
    setequal(main$year,c(2013:2019,2021)) &&
    all(table(main$year,main$source)==1))
  check("2020 excluded from all IPUMS regression samples",!any(counts$year==2020))
  check("case weights equal one",all(counts$retained_weight_sum[counts$source=="STARS"]==
    counts$retained_records[counts$source=="STARS"]))
  check("sample accounting balances",all(counts$eligible_records==
    counts$retained_records+counts$excluded_records))
  sensitivity <- counts[counts$model=="exclude_2021_model",]
  check("pre2020 sensitivity uses seven supported years",setequal(sensitivity$year,2013:2019))
  coverage <- read_table(result_path(matched,"model_year_coverage.csv"))
  check("only years with both outcomes fitted",
    all(coverage$included == (coverage$complete_STARS > 0 & coverage$complete_IPUMS > 0)))
  for (i in seq_len(nrow(models))) {
    selected <- coverage[coverage$model==models$model[i],]
    check(paste(models$model[i],"year effects"),
      models$year_terms[i] == sum(selected$included)-1)
  }
  check("unsupported years have no retained records",
    all(counts$retained_records[!counts$year_supported] == 0))
}

result <- do.call(rbind, checks)
write_table(result, paste0(mode,"_checks.csv"), file.path(output_dir,"validation"))
print(result[!result$passed, , drop = FALSE])
if (!all(result$passed)) stop("Validation failed; see validation_checks.csv")
message(nrow(result), " checks passed.")
