# Adult age-by-sex and same-sample marital-status sensitivities.
# Called by R/individual_regressions.R; shares the weights, age groups and year controls.
marital_no_race <- fit_logistic(
  SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN,
  common,"without_marital_model")
marital_race <- fit_logistic(
  SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + RACE,
  race_common,"without_marital_model_with_race")
marital_spec <- c("Both columns adjust for sex, adult age categories, veteran status, farming, attainment and categorical year.",
  "The With_marital column additionally includes divorced/separated status.",
  "Both columns retain the main Model 4 sample, including its requirement for recorded marital status, to isolate coefficient changes from adjustment.",
  "Removing marital adjustment in this comparison does not recover people excluded for missing marital data.")
publish_models(list(Without_marital=marital_no_race,With_marital=ipums_norace),
  "marital_adjustment_comparison",specification=marital_spec)
publish_models(list(Without_marital=marital_race,With_marital=ipums_race),
  "marital_adjustment_comparison_with_race",race=TRUE,
  specification=c(marital_spec,"Both columns additionally adjust for the harmonized race/ethnicity variable."))
interaction_no_race <- fit_logistic(
  SUIC ~ SEX * AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP,
  common,"age_sex_interaction_model")
interaction_race <- fit_logistic(
  SUIC ~ SEX * AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP + RACE,
  race_common,"age_sex_interaction_model_with_race")
interaction_tests <- list();contrast_outputs <- list();sex_cells <- list()
for (kind in c("without_race","with_race")) {
  fit <- if(kind=="without_race")interaction_no_race else interaction_race
  additive <- if(kind=="without_race")ipums_norace else ipums_race
  stopifnot(identical(rownames(model.frame(fit)),rownames(model.frame(additive))))
  b <- coef(fit);V <- vcov(fit);nm <- names(b)
  idx <- which(grepl(":",nm) & grepl("SEX",nm) & grepl("AGE_ADULT",nm))
  statistic <- as.numeric(t(b[idx])%*%solve(V[idx,idx,drop=FALSE],b[idx]))
  interaction_tests[[kind]] <- data.frame(analysis=kind,chi_square=statistic,df=length(idx),
    p=pchisq(statistic,length(idx),lower.tail=FALSE),n=nobs(fit),
    test="Joint coefficient Wald test; two-source case-process plus ACS SDR covariance")
  frame <- model.frame(fit);ages <- levels(frame$AGE_ADULT)
  numeric <- list();wide <- data.frame(Age_group=ages)
  for (sx in levels(frame$SEX)) {
    new <- frame[rep(1,length(ages)),,drop=FALSE]
    new$SEX <- factor(rep(sx,length(ages)),levels=levels(frame$SEX))
    new$AGE_ADULT <- factor(ages,levels=ages)
    X <- model.matrix(delete.response(terms(fit)),new,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
    stopifnot(identical(colnames(X),names(b)))
    displayed <- character(length(ages))
    for (j in seq_along(ages)) {
      if(j==1) {
        delta <- 0;se <- 0;pv <- NA_real_;displayed[j] <- "1.000 (reference)"
      } else {
        L <- X[j,]-X[1,];delta <- sum(L*b);se <- sqrt(as.numeric(t(L)%*%V%*%L))
        pv <- 2*pnorm(abs(delta/se),lower.tail=FALSE)
        displayed[j] <- sprintf("%.3f (%.3f-%.3f)%s",exp(delta),exp(delta-qnorm(.975)*se),
          exp(delta+qnorm(.975)*se),significance_stars(pv))
      }
      numeric[[length(numeric)+1L]] <- data.frame(analysis=kind,sex=sx,age_group=ages[j],
        reference_age=ages[1],log_odds_ratio=delta,std_error=se,odds_ratio=exp(delta),
        ci_lower=exp(delta-qnorm(.975)*se),ci_upper=exp(delta+qnorm(.975)*se),p=pv)
      cell <- frame[frame$SEX==sx & frame$AGE_ADULT==ages[j],]
      sex_cells[[length(sex_cells)+1L]] <- data.frame(analysis=kind,sex=sx,age_group=ages[j],
        STARS_cases=sum(cell$SUIC==1),IPUMS_records=sum(cell$SUIC==0),
        IPUMS_population_weight=sum(cell[["(weights)"]][cell$SUIC==0]))
    }
    wide[[sx]] <- displayed
  }
  stem <- if(kind=="without_race")"age_effects_by_sex" else "age_effects_by_sex_with_race"
  write_table(wide,paste0(stem,".csv"))
  write_table(do.call(rbind,numeric),paste0(stem,"_long.csv"))
  writeLines(c("Entries are age odds ratios within each sex, relative to age 18-24 of the same sex.",
    "Adults only; fitted as one Model 4 specification with all age-by-sex interaction terms.",
    "Covariates: veteran status, farming, attainment, divorced/separated and categorical year.",
    if(kind=="with_race")"Also adjusts for the existing harmonized race/ethnicity variable." else "Race is not included.",
    "Same sample as the corresponding additive Model 4; these are adjusted contrasts, not crude rates.",
    "* p<0.05; ** p<0.01; *** p<0.001; unadjusted two-sided contrast Wald tests.",
    "Two-source Wald intervals include ACS sampling and Poisson case-process covariance; coding uncertainty is not propagated.",
    "Sparse sex/age cells may have wide intervals; see age_sex_cell_counts.csv.",
    "The joint interaction test is in age_sex_interaction_tests.csv."),
    result_path(output_dir,paste0(stem,"_notes.txt")))
}
tests <- do.call(rbind,interaction_tests)
write_table(tests,"age_sex_interaction_tests.csv")
write_table(do.call(rbind,sex_cells),"age_sex_cell_counts.csv")
report <- c("Adult regression sensitivity analyses",
  "All regressions use adults age 18+ in both sources; reference age is 18-24.",
  "Older groups: 65-74, 75-84 and 85+.",
  "marital_adjustment_comparison.csv and its with_race companion compare omission/inclusion of marital status on identical records.",
  "Marital coding discontinuities remain unresolved; comparisons isolate adjustment, not data repair.",
  "The sex-specific tables derive within-sex age contrasts from Model 4 plus all age-by-sex interactions.",
  "Each interaction model shares records with the corresponding additive Model 4.",
  "Joint interaction Wald tests (two-source case-process plus ACS SDR covariance):")
for(i in seq_len(nrow(tests)))report <- c(report,sprintf("%s: chi-square(%d)=%.4f, p=%.6g; n=%d",
  tests$analysis[i],tests$df[i],tests$chi_square[i],tests$p[i],tests$n[i]))
report <- c(report,"","Older-age odds ratios with and without marital adjustment:")
for(kind in c("without_race","with_race")) {
  reduced <- if(kind=="without_race")marital_no_race else marital_race
  full <- if(kind=="without_race")ipums_norace else ipums_race
  for(t in c("AGE_ADULT65-74","AGE_ADULT75-84","AGE_ADULT85+"))report <- c(report,
    sprintf("%s %s: without marital %.3f; with marital %.3f",kind,t,exp(coef(reduced)[t]),exp(coef(full)[t])))
}
writeLines(report,result_path(output_dir,"sensitivity_report.txt"))
