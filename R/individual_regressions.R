# Publication CSV ranges use ASCII hyphens for reliable direct opening in Excel.
# IPUMS regressions: matching annual populations and deaths.
# source is the August 31 CSV. Education measures associate degree or higher;
# divorced and separated are combined. Two-source replicate-aware uncertainty is used.
parent_output_dir <- output_dir
output_dir <- result_path(output_dir, "regressions")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
analysis_years <- c(2012:2019, 2021L)
ipums_path <- primary_path("IPUMSData/ipums_2012-2019_2021_with_replicate_weights.csv")
codebook_path <- primary_path("IPUMSData/ipums_2012-2019_2021_with_replicate_weights.cbk.txt")
ipums_raw <- read_table(ipums_path)
replicate_names <- paste0("REPWTP",1:80)
require_columns(ipums_raw,c("REPWTP",replicate_names),"IPUMS replicate extract")
stopifnot(all(ipums_raw$REPWTP==1),all(is.finite(as.matrix(ipums_raw[replicate_names]))),
  grepl("REPWTP80",paste(readLines(codebook_path,warn=FALSE),collapse=" "),fixed=TRUE))
source("R/replicate_inference.R")
replicate_weight_summary <- do.call(rbind,lapply(sort(unique(ipums_raw$YEAR)),function(y) {
 d <- ipums_raw[ipums_raw$YEAR==y,]
 do.call(rbind,lapply(replicate_names,function(v)data.frame(year=y,replicate=v,
  weight_sum=sum(d[[v]]),minimum=min(d[[v]]),maximum=max(d[[v]]),negative=sum(d[[v]]<0),zero=sum(d[[v]]==0))))
}))
write_table(replicate_weight_summary,"replicate_weight_summary.csv")

require_columns(ipums_raw, c("YEAR", "SAMPLE", "STATEFIP", "SERIAL", "PERNUM"),
                "Year-matched IPUMS extract")
stopifnot(!anyNA(ipums_raw$YEAR), setequal(unique(ipums_raw$YEAR), analysis_years),
  all(ipums_raw$STATEFIP == 30), all(ipums_raw$SAMPLE == ipums_raw$YEAR * 100 + 1),
  all(is.finite(ipums_raw$PERWT) & ipums_raw$PERWT > 0),
  !anyDuplicated(ipums_raw[c("YEAR", "SAMPLE", "SERIAL", "PERNUM")]))
# Verify observed codes against the supplied extract's documented categories.
allowed <- list(SEX = 1:2, MARST = 1:6, RACE = 1:9, HISPAN = c(0:4,9),
                VETSTAT = c(0:2,9))
for (v in names(allowed)) {
  if (anyNA(ipums_raw[[v]]) || !all(ipums_raw[[v]] %in% allowed[[v]]))
    stop("Unexpected or missing IPUMS codes in ", v)
}
ipums <- harmonize_ipums(ipums_raw)
ipums$YEAR <- factor(ipums_raw$YEAR, levels = analysis_years)
ipums$EDUC_ATTAIN <- attainment_ipums(ipums_raw)
stars_csv_raw <- read_stars_csv()
stars_csv_all <- harmonize_stars(stars_csv_raw)
stars_csv_all$EDUC_ATTAIN <- attainment_stars(stars_csv_raw)
stars_residence <- stars_csv_raw$DEM_RES09State
resident <- montana_resident(stars_csv_raw)
write_table(regression_residence_audit(stars_csv_raw,analysis_years),"regression_residence_eligibility.csv")
stars_rows <- which(stars_csv_all$DYEAR %in% analysis_years & resident)
stars <- stars_csv_all[stars_rows, ]
stars$weight <- 1
stars$YEAR <- factor(stars$DYEAR, levels = analysis_years)
ipums$SOURCE_ROW <- seq_len(nrow(ipums))
stars$SOURCE_ROW <- stars_rows
shared <- intersect(names(ipums), names(stars))
ipums_stars <- rbind(ipums[shared], stars[shared])
# Apply the same adult eligibility and age cutpoints to both sources.
ipums_stars$AGE_ADULT <- cut(ipums_stars$AGE,
  breaks=c(17,24,34,44,54,64,74,84,Inf),
  labels=c("18-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"))
ipums_stars_all <- ipums_stars
adult_eligibility <- do.call(rbind,lapply(0:1,function(k) {
  d <- ipums_stars_all[ipums_stars_all$SUIC==k,]
  data.frame(source=c("IPUMS","STARS")[k+1L],input_records=nrow(d),
    age_missing=sum(is.na(d$AGE)),under_18=sum(d$AGE<18,na.rm=TRUE),
    adult_records=sum(d$AGE>=18,na.rm=TRUE))
}))
write_table(adult_eligibility,"adult_eligibility.csv")
ipums_stars <- ipums_stars[!is.na(ipums_stars$AGE) & ipums_stars$AGE>=18,]
stopifnot(all(ipums_stars$AGE>=18),!anyNA(ipums_stars$AGE_ADULT),
  levels(ipums_stars$AGE_ADULT)[1]=="18-24")
stopifnot(setequal(as.integer(as.character(stars$YEAR)), analysis_years))
write_table(data.frame(year = 2012:2021,
  deaths = as.integer(table(factor(stars_csv_all$DYEAR, levels = 2012:2021))),
  confirmed_resident_deaths = as.integer(table(factor(stars_csv_all$DYEAR[resident],levels=2012:2021))),
  included_in_ipums_models = 2012:2021 %in% analysis_years), "stars_year_eligibility.csv")
sample_counts <- list()
missingness <- list()
model_checks <- list()
model_coverage <- list()
model_warnings <- list()
source("R/regression_helpers.R")
main_variables <- c("SUIC","SEX","AGE_ADULT","YEAR","VETSTAT","FARM","EDUC_ATTAIN","DIVSEP","weight")
common_complete <- complete.cases(ipums_stars[main_variables])
by_year <- table(ipums_stars$YEAR[common_complete],
  factor(ipums_stars$SUIC[common_complete],levels=0:1))
common_years <- rownames(by_year)[rowSums(by_year>0)==2]
common <- droplevels(ipums_stars[common_complete &
  as.character(ipums_stars$YEAR) %in% common_years,])
flow <- list()
for (y in levels(ipums_stars$YEAR)) for (k in 0:1) {
  before <- ipums_stars[ipums_stars$YEAR==y & ipums_stars$SUIC==k,]
  retained <- common[common$YEAR==y & common$SUIC==k,]
  flow[[length(flow)+1]] <- data.frame(year=y,source=c("IPUMS","STARS")[k+1L],
    input_records=nrow(before),complete_on_model4=sum(complete.cases(before[main_variables])),
    year_supported=y %in% common_years,retained_records=nrow(retained))
}
write_table(do.call(rbind,flow),"principal_sample_flow.csv")
model1 <- fit_logistic(SUIC ~ SEX + AGE_ADULT,common,"principal_model1")
model2 <- fit_logistic(SUIC ~ SEX + AGE_ADULT + VETSTAT,common,"principal_model2")
model3 <- fit_logistic(SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM,common,"principal_model3")
ipums_norace <- fit_logistic(SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP,
  common,"principal_model4")
primary_models <- list(Model_1=model1,Model_2=model2,Model_3=model3,Model_4=ipums_norace)
same_rows <- vapply(primary_models,function(m)
  identical(rownames(model.frame(m)),rownames(model.frame(model1))),logical(1))
stopifnot(all(same_rows))
write_table(data.frame(model=names(primary_models),n=vapply(primary_models,nobs,integer(1)),
  same_records=unname(same_rows)),"principal_shared_sample_checks.csv")
# CSV-ready publication tables share one formatter and unadjusted Wald p-values.
publish_models(primary_models,"sequential_models")
# Compare race models on the identical race-complete sample.
race_variables <- c(main_variables,"RACE")
race_common <- ipums_stars[complete.cases(ipums_stars[race_variables]),]
ipums_race_comparison <- fit_logistic(
  SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP,
  race_common,"race_complete_model_without_race")
ipums_race <- fit_logistic(
  SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP + RACE,
  race_common,"principal_race_model4")
stopifnot(identical(rownames(model.frame(ipums_race_comparison)),rownames(model.frame(ipums_race))))
race_model1 <- fit_logistic(SUIC ~ SEX + AGE_ADULT + RACE,race_common,"principal_race_model1")
race_model2 <- fit_logistic(SUIC ~ SEX + AGE_ADULT + RACE + VETSTAT,race_common,"principal_race_model2")
race_model3 <- fit_logistic(SUIC ~ SEX + AGE_ADULT + RACE + VETSTAT + FARM,race_common,"principal_race_model3")
race_models <- list(Model_1=race_model1,Model_2=race_model2,Model_3=race_model3,Model_4=ipums_race)
publish_models(race_models,"sequential_models_with_race",race=TRUE)
write_table(data.frame(model=names(race_models),n=vapply(race_models,nobs,integer(1)),
  same_records=vapply(race_models,function(m)identical(rownames(model.frame(m)),
    rownames(model.frame(race_model1))),logical(1))),"race_shared_sample_checks.csv")

pre2020 <- droplevels(common[common$YEAR!="2021",])
pre2020_model <- fit_logistic(
  SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP,
  pre2020,"exclude_2021_model")
source("R/regression_sensitivities.R")
source("R/regression_selection.R")
source("R/regression_age_diagnostics.R")
write_table(as.data.frame(table(year=stars_csv_raw$DEATH_2Year,
  original=stars_csv_raw$DEM_EDU1EducationAchieved,attainment=stars_csv_all$EDUC_ATTAIN)),
  "stars_attainment_crosswalk_counts.csv")
write_table(as.data.frame(table(year=ipums_raw$YEAR,EDUCD=ipums_raw$EDUCD,
  attainment=ipums$EDUC_ATTAIN)),"ipums_attainment_crosswalk_counts.csv")
write_table(as.data.frame(table(year=stars_csv_raw$DEATH_2Year,
  original=stars_csv_raw$DEM_MAR1MaritalStatus,DIVSEP=stars_csv_all$DIVSEP)),
  "stars_marital_crosswalk_counts.csv")
write_table(do.call(rbind,sample_counts), "model_sample_counts_by_year.csv")
write_table(do.call(rbind,missingness), "model_missingness_by_year.csv")
write_table(do.call(rbind,replicate_checks),"replicate_model_checks.csv")
write_table(do.call(rbind,replicate_changes),"covariance_method_comparison.csv")
stopifnot(all(vapply(replicate_checks,function(x)x$passed,logical(1))))
checks <- do.call(rbind,model_checks)
write_table(checks, "model_checks.csv")
write_table(do.call(rbind,model_coverage),"model_year_coverage.csv")
write_table(do.call(rbind,model_warnings), "model_warnings.csv")
extract_counts <- as.data.frame(table(year=ipums_raw$YEAR))
names(extract_counts) <- c("year","records")
write_table(extract_counts,"extract_record_counts.csv")
# Aggregate coding audits contain no person-level identifiers.
for (v in c("EDUCD","OCC2010","VETSTAT","MARST")) {
  write_table(as.data.frame(table(year=ipums_raw$YEAR, code=ipums_raw[[v]])),
    paste0("ipums_codes_",tolower(v),".csv"))
}
input_paths <- c(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"),
                 ipums_path, codebook_path)
write_table(data.frame(input=input_paths,md5=unname(tools::md5sum(input_paths))),"input_checksums.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
main <- checks[checks$model=="principal_model4",]
writeLines(c("Year-matched IPUMS regression results",
  paste("Run time (UTC):",format(Sys.time(),tz="UTC",usetz=TRUE)),
  "IPUMS extract: ipums_2012-2019_2021_with_replicate_weights.csv; codebook checked alongside extract.",
  "Years: 2012-2019 and 2021. STARS 2020 deaths are excluded from these regressions.",
  "STARS source: 20230831_Montana_Suicide_Database.csv.",
  "All regressions are restricted to adults age 18 or older in both sources.",
  "STARS deaths require recorded MONTANA residence; non-Montana and unknown residence are excluded before sample selection.",
  "Age categories split older adults into 65-74, 75-84 and 85+; reference 18-24.",
  "Four-model table uses one complete-case sample for all columns.",
  "See ../sensitivity/age_sex_interaction/sensitivity_report.txt for marital and age-by-sex comparisons.",
  "See ../sensitivity/sample_selection/selection_sensitivity_report.txt for eligibility and sample-selection comparisons.",
  "Model 1: sex + age + year; Model 2 adds veteran; Model 3 adds farmer/rancher;",
  "Model 4 adds associate degree or higher and combined divorced/separated.",
  "EDUC_ATTAIN reference: less than associate degree (some college without degree included).",
  "STARS degree labels: ASSOCIATES, BACHELORS, MASTERS, DOCTORATE.",
  "IPUMS degree EDUCD codes: 81, 101, 114, 115, 116; NIU 0/1 stays missing.",
  "Reference categories: male; age 18-24; other marital status; less than associate degree;",
  "not in armed forces; other occupation. Year reference is the first supported year.",
  "STARS case weight = 1; IPUMS weight = PERWT. No division by nine or ten.",
  "Each year's cases are compared with its annual population, with year intercepts.",
  paste("Eligible STARS records:",nrow(stars)),
  paste("IPUMS input records:",nrow(ipums)),
  paste("Main model complete cases:",main$cases,"STARS and",main$comparison_records,"IPUMS"),
  paste("Main model year reference:",main$year_reference),
  "Some candidate years have no complete-case deaths for particular models.",
  "Such years are excluded from BOTH sources for that model; see model_year_coverage.csv.",
  "Predictor missingness and year exclusions are separate fields in the audit outputs.",
  paste("Model checks passed:",sum(checks$passed),"of",nrow(checks)),
  "All models use categorical year adjustment; sensitivity model omits 2021.",
  "Principal regressions use broad attainment; a separately labeled supplementary partial-education sequence is available.",
  "Divorced and separated are combined in the principal and supporting marital models.",
  "Combining categories does not resolve the documented source coding uncertainties.",
  "Uncertainty: Poisson case-process covariance plus year-blocked ACS SDR covariance; normal Wald intervals and tests.",
  "80 person replicate weights per annual ACS sample are supplied; signed values are preserved. No weights are divided by the number of years.",
  "Temporal alignment does not resolve measurement differences or missing-data bias."),
  result_path(output_dir,"analysis_report.txt"))
if (!all(checks$passed)) stop("Year-matched model checks failed; inspect model_checks.csv")
message("Year-matched IPUMS results: ",output_dir)
output_dir <- parent_output_dir
