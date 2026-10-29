# Common regression fitting and publication-table functions.
fit_logistic <- function(formula, data, name) {
  formula <- update(formula, . ~ . + YEAR)
  variables <- unique(c(all.vars(formula), "weight"))
  complete <- complete.cases(data[variables])
  available <- table(data$YEAR[complete], factor(data$SUIC[complete], levels=0:1))
  supported <- rownames(available)[rowSums(available > 0) == 2]
  if (length(supported) < 2) stop("Fewer than two supported years for ", name)
  # Structural missingness can leave a year with population observations but
  # no usable deaths. Exclude that year from BOTH sources for this model.
  for (y in levels(data$YEAR)) model_coverage[[length(model_coverage)+1L]] <<- data.frame(
    model=name, year=y, complete_IPUMS=available[y,"0"],
    complete_STARS=available[y,"1"], included=y %in% supported)
  model_data <- droplevels(data[as.character(data$YEAR) %in% supported, ])
  used <- model_data[complete.cases(model_data[variables]), ]
  case_states <- stars_residence[used$SOURCE_ROW[used$SUIC==1]]
  residence_verified <- !anyNA(case_states) && all(case_states=="MONTANA")
  stopifnot(residence_verified,all(ipums_raw$STATEFIP[used$SOURCE_ROW[used$SUIC==0]]==30))
  observed <- character()
  model <- withCallingHandlers(
    glm(formula, data = model_data, weights = weight, family = binomial(), control = glm.control(epsilon=1e-12,maxit=100), na.action = na.omit),
    warning = function(w) { observed <<- c(observed, conditionMessage(w)); invokeRestart("muffleWarning") })
  model <- replicate_vcov(model,used,name)
  estimates <- coefficient_table(model, name)
  estimates$ci_lower <- exp(estimates$coefficient - qnorm(.975) * estimates$std_error)
  estimates$ci_upper <- exp(estimates$coefficient + qnorm(.975) * estimates$std_error)
  write_table(estimates, paste0(name, ".csv"))
  for (k in 0:1) for (y in levels(data$YEAR)) {
    before <- data[data$SUIC == k & data$YEAR == y, ]
    after <- used[used$SUIC == k & used$YEAR == y, ]
    sample_counts[[length(sample_counts)+1L]] <<- data.frame(
      model = name, source = c("IPUMS","STARS")[k+1L], year = y,
      eligible_records = nrow(before), retained_records = nrow(after),
      excluded_records = nrow(before)-nrow(after), year_supported = y %in% supported,
      eligible_weight_sum = sum(before$weight), retained_weight_sum = sum(after$weight))
    for (v in variables) missingness[[length(missingness)+1L]] <<- data.frame(
      model=name, source=c("IPUMS","STARS")[k+1L], year=y,
      variable=v, records=nrow(before), missing=sum(is.na(before[[v]])))
  }
  stats <- summary(model)$coefficients
  valid <- residence_verified && all(used$AGE>=18) && model$converged && model$rank == length(coef(model)) &&
    all(is.finite(stats)) && all(is.finite(estimates$ci_lower)) &&
    all(is.finite(estimates$ci_upper)) && nobs(model) == nrow(used) &&
    all(table(used$YEAR, used$SUIC) > 0) &&
    length(grep("^YEAR", names(coef(model)))) == nlevels(model_data$YEAR)-1L
  model_checks[[length(model_checks)+1L]] <<- data.frame(model=name,
    model_rows=nobs(model), cases=sum(used$SUIC == 1), comparison_records=sum(used$SUIC == 0),
    confirmed_montana_cases=sum(case_states=="MONTANA"), non_montana_cases=sum(!is.na(case_states) & case_states!="MONTANA"), unknown_residence_cases=sum(is.na(case_states)), residence_verified=residence_verified,
    weight_total=sum(model$prior.weights), converged=model$converged, rank=model$rank,
    parameters=length(coef(model)), variance_method=model$covariance_method, minimum_age=min(used$AGE), maximum_age=max(used$AGE), year_reference=levels(model_data$YEAR)[1], year_terms=length(grep("^YEAR",names(coef(model)))),
    passed=valid)
  model_warnings[[length(model_warnings)+1L]] <<- data.frame(model=name,
    warning=if(length(observed)) unique(observed) else "")
  writeLines(c(model$covariance_method,capture.output(summary(model))), result_path(output_dir,paste0(name,"_summary.txt")))
  model
}
significance_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
}
labels <- setNames(c("Female (reference: male)",
  "Age 25-34 (reference: 18-24)","Age 35-44","Age 45-54","Age 55-64",
  "Age 65-74","Age 75-84","Age 85 or older",
  "Ever in armed forces (reference: no)","Farmer/rancher (reference: other)",
  "Associate degree or higher (reference: less)",
  "Divorced/separated (reference: other recorded status)",
  "American Indian (reference: White)","Other race/ethnicity (reference: White)"),
  c("SEXFEMALE","AGE_ADULT25-34","AGE_ADULT35-44","AGE_ADULT45-54","AGE_ADULT55-64",
    "AGE_ADULT65-74","AGE_ADULT75-84","AGE_ADULT85+",
    "VETSTATAF","FARMFARMER","EDUC_ATTAINASSOCIATE_OR_HIGHER","DIVSEPDIVORCED/SEPARATED",
    "RACEAMERICAN INDIAN","RACEOTHER"))
publish_models <- function(models, stem, race = FALSE, specification = NULL, same_sample = TRUE) {
  final <- models[[length(models)]]
  terms <- setdiff(names(coef(final)),c("(Intercept)",grep("^YEAR",names(coef(final)),value=TRUE)))
  stopifnot(all(terms %in% names(labels)))
  frames <- lapply(models,model.frame)
  same <- vapply(frames,function(f)identical(rownames(f),rownames(frames[[1]])),logical(1))
  if(same_sample)stopifnot(all(same))
  wide <- data.frame(term=terms);long <- list()
  for (name in names(models)) {
    m <- models[[name]];b <- coef(m);se <- sqrt(diag(vcov(m)));st <- summary(m)$coefficients
    wide[[name]] <- vapply(terms,function(t)if(!t %in% names(b)) "" else
      sprintf("%.3f (%.3f-%.3f)%s",exp(b[t]),exp(b[t]-qnorm(.975)*se[t]),
        exp(b[t]+qnorm(.975)*se[t]),significance_stars(st[t,4])),character(1))
    long[[name]] <- data.frame(model=name,term=rownames(st),coefficient=st[,1],
      std_error=st[,2],p=st[,4],odds_ratio=exp(st[,1]),
      ci_lower=exp(st[,1]-qnorm(.975)*st[,2]),ci_upper=exp(st[,1]+qnorm(.975)*st[,2]))
  }
  footer <- data.frame(term=c("Year adjusted","STARS cases","IPUMS comparisons","Total observations"))
  for (name in names(models)) {
    f <- frames[[name]]
    footer[[name]] <- c("Yes",sum(f$SUIC==1),sum(f$SUIC==0),nrow(f))
  }
  wide$term <- unname(labels[wide$term]);publication <- rbind(wide,footer)
  names(publication)[1] <- "Predictor"
  write_table(publication,paste0(stem,".csv"))
  write_table(do.call(rbind,long),paste0(stem,"_long.csv"))
  notes <- c("Entries are odds ratios (two-source Wald 95% confidence intervals).",
    "* p < 0.05; ** p < 0.01; *** p < 0.001. Stars use unrounded, two-sided coefficient Wald p-values; no multiplicity adjustment.",
    if(same_sample) "All columns use the identical adult complete-case sample within this table." else
      "Columns use different adult samples; sample sizes and source counts are reported below.",
    paste("Eligible age groups in this table:",paste(levels(frames[[1]]$AGE_ADULT),collapse=", ")),
    if(!is.null(specification)) specification else c(
    if(race) "Race is included in every column; Model 1 also includes sex, age and categorical year." else
      "Model 1: sex, age, categorical year.",
    "Model 2 adds veteran status; Model 3 adds farmer/rancher occupation;",
    "Model 4 adds associate degree or higher and combined divorced/separated."),
    paste("Age reference:",levels(frames[[1]]$AGE_ADULT)[1],"for every age row. Year effects are included but not displayed."),
    paste("Study years in the common sample:",paste(levels(frames[[1]]$YEAR),collapse=", ")),
    "STARS case weights are 1; IPUMS weights are PERWT.",
    "Only deaths with recorded MONTANA residence are eligible; known non-Montana and unknown-residence records are excluded.",
    "Attainment: associate/bachelor/master/doctorate or IPUMS EDUCD 81/101/114/115/116.",
    "Some college without a degree is classified as less than associate.",
    "Divorced and separated are pooled; unresolved source coding limitations remain.",
    "Covariance combines a Poisson case-process sandwich with ACS sampling covariance from 80 SDR replicates per annual sample (4/80 scaling).",
    "Annual ACS samples and mortality counts are treated as independent; see docs/REPLICATE_INFERENCE.md.",
    "Complete-case selection, including age-specific IPUMS variable universes, can distort age comparisons; see age_diagnostics_report.txt.")
  if(race) notes <- c(notes,
    "Race reference: White. Recorded Hispanic ethnicity is pooled into Other for both sources, following the existing harmonization.",
    "Other combines several race/ethnicity categories. This is not a clean race-only classification.",
    "The race table additionally requires recorded race; cross-table differences can reflect this selection as well as race adjustment.")
  writeLines(notes,result_path(output_dir,paste0(stem,"_notes.txt")))
  invisible(publication)
}
