# Aggregate diagnostics for age patterns; called within R/individual_regressions.R.
# No person-level output; preserve main table samples and specifications.
stopifnot(all(is.na(ipums$AGE) | ipums$AGE>=0 & ipums$AGE<=120),
  all(is.na(stars$AGE) | stars$AGE>=0 & stars$AGE<=120))
rows <- list()
for (source in c("IPUMS","STARS")) for (ag in levels(ipums_stars_all$AGE7)) {
  k <- if(source=="IPUMS")0 else 1
  base <- ipums_stars_all[ipums_stars_all$SUIC==k & !is.na(ipums_stars_all$AGE7) & ipums_stars_all$AGE7==ag,]
  kept <- common[common$SUIC==k & common$AGE7==ag,]
  rr <- data.frame(source=source,age_group=ag,eligible_records=nrow(base),retained_records=nrow(kept),
    eligible_weight=sum(base$weight),retained_weight=sum(kept$weight),
    retained_percent=100*nrow(kept)/nrow(base),
    retained_weight_percent=100*sum(kept$weight)/sum(base$weight))
  for(v in setdiff(main_variables,c("SUIC","weight","YEAR")))rr[[paste0("missing_",v)]] <- sum(is.na(base[[v]]))
  rows[[length(rows)+1]] <- rr
}
write_table(do.call(rbind,rows),"age_sample_retention.csv")
niu <- do.call(rbind,lapply(sort(unique(ipums_raw$AGE)),function(a) {
  x <- ipums_raw[ipums_raw$AGE==a,]
  data.frame(age=a,records=nrow(x),veteran_NIU=sum(x$VETSTAT==0),
    veteran_NIU_weight=sum(x$PERWT[x$VETSTAT==0]))
}))
write_table(niu,"age_ipums_veteran_universe.csv")
# Crude all-record source ratios, unaffected by covariate missingness.
# Comparisons are annual population-weight sums, not weighted survey-record counts.
breaks <- c(-Inf,17,24,34,44,54,64,74,84,Inf)
age_labels <- c("<18","18-24","25-34","35-44","45-54","55-64","65-74","75-84","85+")
coverage <- do.call(rbind,lapply(levels(ipums_stars_all$YEAR),function(y) {
 d <- ipums_stars_all[ipums_stars_all$YEAR==y & ipums_stars_all$SUIC==1,]
 data.frame(year=y,STARS_cases=nrow(d),known_age=sum(!is.na(d$AGE)),
  known_sex_age=sum(complete.cases(d[c("AGE","SEX")])),
  known_race_sex_age=sum(complete.cases(d[c("AGE","SEX","RACE")])))
}))
write_table(coverage,"age_source_year_coverage.csv")
crude <- list()
for (stratum in c("All","Male","White male")) {
  needed <- c("AGE",if(stratum!="All")"SEX",if(stratum=="White male")"RACE")
  available <- ipums_stars_all[ipums_stars_all$SUIC==1 & complete.cases(ipums_stars_all[needed]),]
  supported <- names(which(table(available$YEAR)>0))
  d <- ipums_stars_all[as.character(ipums_stars_all$YEAR) %in% supported,]
  if(stratum!="All")d <- d[!is.na(d$SEX) & d$SEX=="MALE",]
  if(stratum=="White male")d <- d[!is.na(d$RACE) & d$RACE=="WHITE",]
  group <- cut(d$AGE,breaks=breaks,labels=age_labels)
  for(ag in age_labels) {
    keep <- !is.na(group) & group==ag
    deaths <- sum(d$SUIC[keep]==1);pop <- sum(d$weight[keep & d$SUIC==0])
    crude[[length(crude)+1]] <- data.frame(stratum=stratum,age_group=ag,
      STARS_deaths=deaths,IPUMS_person_years=pop,crude_per_100000=deaths/pop*1e5,
      study_years=paste(supported,collapse=";"))
  }
}
crude <- do.call(rbind,crude);write_table(crude,"age_crude_source_rates.csv")
# Isolate covariate adjustment on the fixed main-table sample.
age_fits <- list()
age_fits[["Full available sex_age_year"]] <- glm(SUIC~SEX+AGE_ADULT+YEAR,
  data=ipums_stars[as.character(ipums_stars$YEAR) %in% coverage$year[coverage$known_sex_age>0],],
  weights=weight,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100),na.action=na.omit)
age_fits[["Common sample Model 1"]] <- model1
age_fits[["Common sample Model 3"]] <- model3
age_fits[["Common sample Model 3 plus education"]] <- glm(
  SUIC~SEX+AGE_ADULT+VETSTAT+FARM+EDUC_ATTAIN+YEAR,data=common,weights=weight,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100))
age_fits[["Common sample Model 3 plus marital"]] <- glm(
  SUIC~SEX+AGE_ADULT+VETSTAT+FARM+DIVSEP+YEAR,data=common,weights=weight,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100))
age_fits[["Common sample Model 4"]] <- ipums_norace
age_fits[["Race complete Model 4 without race"]] <- ipums_race_comparison
age_fits[["Race complete Model 4 with race"]] <- ipums_race
results <- list()
for (nm in names(age_fits)) {
  fit <- age_fits[[nm]]
  if(!inherits(fit,"stars_replicate_glm")) {
    matched <- match(rownames(model.frame(fit)),rownames(ipums_stars))
    stopifnot(!anyNA(matched))
    fit <- replicate_vcov(fit,ipums_stars[matched,,drop=FALSE],paste0("age_diagnostic_",match(nm,names(age_fits))))
  }
  st <- summary(fit)$coefficients
  stopifnot(fit$converged,fit$rank==length(coef(fit)),all(is.finite(st)))
  for(t in grep("^AGE_ADULT",rownames(st),value=TRUE))results[[length(results)+1]] <- data.frame(
    model=nm,term=t,odds_ratio=exp(st[t,1]),ci_lower=exp(st[t,1]-qnorm(.975)*st[t,2]),
    ci_upper=exp(st[t,1]+qnorm(.975)*st[t,2]),p=st[t,4],n=nobs(fit),
    cases=sum(model.frame(fit)$SUIC==1))
}
results <- do.call(rbind,results);write_table(results,"age_adjustment_diagnostics.csv")
report <- c("Age pattern diagnostics for year-matched STARS/IPUMS regressions",
  "All source years are 2012-2019 and 2021; 2020 is excluded.",
  "The main table is held on the Model 4 complete-case sample even in Model 1.",
  "Inspect age_sample_retention.csv for source-specific age selection and missingness.",
  "IPUMS veteran not-in-universe observations are not coded as nonveterans in these models.",
  "Age groups use the same numeric cutpoints for both sources; the 18-24 reference is explicit.",
  "Age coefficients are conditional associations with case versus population-record status, not predictions of future death.",
  "No sex-by-age or race-by-age interactions appear in the four-column tables.",
  "Principal regressions include adults age 18+ with separate 65-74, 75-84 and 85+ categories.",
  "The retention audit preserves the raw child counts to document their exclusion; no child enters a regression.",
  "Crude audit values use STARS cases and IPUMS population weights without Model 4 complete-case filtering.",
  "Each crude stratum excludes years with no cases having the necessary recorded age/sex/race in BOTH sources.",
  "Study years are listed in the crude CSV. Remaining within-year missing demographics can underestimate source rates.",
  "These are source-based descriptive rates, not independently validated official resident mortality rates.",
  "Marital status has documented source-coding discontinuities. Adjustment can transmit these errors into age coefficients.",
  "The race table retains the existing Hispanic-to-Other classification; White male refers to that harmonized category.",
  "Standard errors combine Poisson case-process and year-blocked ACS SDR covariance; coding uncertainty is not propagated.","",
  "Older-adult odds ratios relative to age 18-24:")
for(i in which(results$term %in% c("AGE_ADULT65-74","AGE_ADULT75-84","AGE_ADULT85+")))report <- c(report,sprintf("%s %s: %.3f (%.3f-%.3f), n=%d",
  results$model[i],results$term[i],results$odds_ratio[i],results$ci_lower[i],results$ci_upper[i],results$n[i]))
report <- c(report,"",paste("Crude source-based White male rates per 100,000 person-years; years:",crude$study_years[crude$stratum=="White male"][1]))
for(i in which(crude$stratum=="White male"))report <- c(report,sprintf("%s: deaths=%d, person-years=%.0f, rate=%.2f",
  crude$age_group[i],crude$STARS_deaths[i],crude$IPUMS_person_years[i],crude$crude_per_100000[i]))
writeLines(report,result_path(output_dir,"age_diagnostics_report.txt"))
