# Sequential adult regressions with separate divorced and separated categories.
source("R/common.R",local=TRUE)
output_dir <- file.path(results_root(),"regressions","sensitivity","marital_three_categories")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source("R/regression_education_inputs.R",local=TRUE)
source("R/replicate_inference.R",local=TRUE)
sample_counts <- missingness <- model_checks <- model_coverage <- model_warnings <- list()
source("R/regression_helpers.R",local=TRUE)
labels <- c(labels,MARST2DIVORCED="Divorced (reference: other recorded status)",
  MARST2SEPARATED="Separated (reference: other recorded status)")
stopifnot(identical(is.na(ipums_stars$MARST2),is.na(ipums_stars$DIVSEP)),
  all(as.character(ipums_stars$MARST2[ipums_stars$DIVSEP=="OTHER" & !is.na(ipums_stars$DIVSEP)])=="OTHER"))
variables <- c("SUIC","SEX","AGE_ADULT","YEAR","VETSTAT","FARM","EDUC_ATTAIN","MARST2","weight")
tests <- comparisons <- groups <- list()
principal_dir <- file.path(results_root(),"regressions")
paths <- c(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"),ipums_path,codebook_path)
principal_inputs <- read_table(result_path(principal_dir,"input_checksums.csv"))
stopifnot(identical(unname(tools::md5sum(paths)),principal_inputs$md5))
for(race in c(FALSE,TRUE)) {
  scope <- if(race)"with_race" else "without_race"
  complete <- complete.cases(ipums_stars[c(variables,if(race)"RACE")])
  available <- table(ipums_stars$YEAR[complete],factor(ipums_stars$SUIC[complete],levels=0:1))
  supported <- rownames(available)[rowSums(available>0)==2]
  common <- droplevels(ipums_stars[complete & as.character(ipums_stars$YEAR) %in% supported,])
  stem <- paste0("sequential_marital_three_categories",if(race)"_with_race")
  baseline <- c("SEX","AGE_ADULT",if(race)"RACE")
  formulas <- list(reformulate(baseline,"SUIC"),reformulate(c(baseline,"VETSTAT"),"SUIC"),
    reformulate(c(baseline,"VETSTAT","FARM"),"SUIC"),
    reformulate(c(baseline,"VETSTAT","FARM","EDUC_ATTAIN","MARST2"),"SUIC"))
  models <- setNames(lapply(1:4,function(i)fit_logistic(formulas[[i]],common,
    paste0(stem,"_model",i))),paste0("Model_",1:4))
  publish_models(models,stem,race=race,specification=c(
    if(race)"Model 1: sex, age, race/ethnicity and categorical year." else "Model 1: sex, age and categorical year.",
    "Model 2 adds veteran status; Model 3 adds farmer/rancher occupation.",
    "Model 4 adds associate degree or higher and three-category marital status: other, divorced, separated."))
  notes_path <- file.path(output_dir,paste0(stem,"_notes.txt"))
  notes <- readLines(notes_path,warn=FALSE)
  notes <- notes[!grepl("^Divorced and separated are pooled",notes)]
  writeLines(c(notes,
    "Marital reference: OTHER, including never married/single, widowed, and married but not separated.",
    "STARS DIVORCED maps to divorced; MARRIED (SEPARATED) maps to separated.",
    "IPUMS MARST 4 maps to divorced; 3 maps to separated; 1/2/5/6 map to other.",
    "All four models share the corresponding principal complete-case cohort; only Model 4 changes the marital grouping.",
    "Recorded-category discontinuities and missingness are not corrected by separating these labels.",
    "The divorced-versus-separated contrast and its 1-df equality test are in divorced_separated_contrast.csv."),notes_path)
  for(i in 1:3) {
    principal_name <- paste0(if(race)"principal_race_model" else "principal_model",i)
    previous <- read_table(result_path(principal_dir,paste0(principal_name,".csv")))
    b <- coef(models[[i]]);se <- sqrt(diag(vcov(models[[i]])))
    error <- max(abs(previous$coefficient-b[match(previous$term,names(b))]))
    se_error <- max(abs(previous$std_error-se[match(previous$term,names(se))]))
    stopifnot(is.finite(error),error<1e-8,is.finite(se_error),se_error<1e-8)
    comparisons[[length(comparisons)+1]] <- data.frame(analysis=scope,model=i,
      coefficient_error=error,standard_error_error=se_error,passed=TRUE)
  }
  m <- models[[4]];b <- coef(m);V <- vcov(m)
  L <- setNames(rep(0,length(b)),names(b));L["MARST2DIVORCED"] <- 1;L["MARST2SEPARATED"] <- -1
  estimate <- sum(L*b);se <- sqrt(as.numeric(t(L)%*%V%*%L))
  tests[[scope]] <- data.frame(analysis=scope,comparison="Divorced versus separated",
    odds_ratio=exp(estimate),ci_lower=exp(estimate-qnorm(.975)*se),
    ci_upper=exp(estimate+qnorm(.975)*se),chi_square=(estimate/se)^2,df=1,
    p=pchisq((estimate/se)^2,1,lower.tail=FALSE),
    cases=sum(common$SUIC==1),comparison_records=sum(common$SUIC==0),years=paste(supported,collapse=";"))
  for(k in 0:1)for(y in levels(common$YEAR))for(g in levels(common$MARST2)) {
    d <- common[common$SUIC==k & common$YEAR==y & common$MARST2==g,]
    groups[[length(groups)+1]] <- data.frame(analysis=scope,source=c("IPUMS","STARS")[k+1],
      year=y,marital_status=g,records=nrow(d),weight_sum=sum(d$weight))
  }
}
contrast <- do.call(rbind,tests);contrast$p_holm <- p.adjust(contrast$p,"holm")
write_table(contrast,"divorced_separated_contrast.csv")
write_table(do.call(rbind,groups),"marital_category_counts_by_year.csv")
write_table(do.call(rbind,comparisons),"principal_sample_comparison_checks.csv")
for(item in c("sample_counts","missingness","model_checks","model_coverage","model_warnings"))
  write_table(do.call(rbind,get(item)),paste0(item,".csv"))
write_table(do.call(rbind,replicate_checks),"replicate_model_checks.csv")
write_table(do.call(rbind,replicate_changes),"covariance_method_comparison.csv")
stopifnot(all(do.call(rbind,model_checks)$passed),all(do.call(rbind,replicate_checks)$passed))
write_table(data.frame(input=paths,md5=unname(tools::md5sum(paths))),"input_checksums.csv")
writeLines(capture.output(sessionInfo()),file.path(output_dir,"sessionInfo.txt"))
writeLines(c("Three-category marital-status sensitivity analysis",
  "Categories: other (reference), divorced, separated. Missing status remains missing.",
  "Four sequential models share the principal sample within each race specification.",
  "Non-race years: 2013-2019 and 2021. Race-inclusive years: 2017-2019 and 2021.",
  "Model 4 adjusts for sex, age, veteran, farmer/rancher, degree attainment and categorical year; the race sequence also includes race/ethnicity.",
  "Uncertainty combines Poisson mortality-event covariance and annual ACS SDR covariance from 80 signed person replicate weights per year.",
  "Contrast OR >1 indicates higher adjusted odds for divorced than separated records.",
  "The equality test is a 1-df coefficient Wald test; the two race specifications form one two-test Holm family.",
  "Coefficient stars and contrast confidence intervals are pointwise and not multiplicity-adjusted.",
  "These comparisons do not resolve the uncertainty in recorded marital-status classifications.",
  "",capture.output(print(contrast,row.names=FALSE))),file.path(output_dir,"analysis_report.txt"))
print(contrast,row.names=FALSE)
message("Three-category marital results: ",output_dir)
