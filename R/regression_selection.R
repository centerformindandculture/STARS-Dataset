# Missingness, calendar-period and record-selection sensitivity analyses.
# Called by R/individual_regressions.R after main and interaction models; exports aggregates only.
selection_profiles <- list();selection_age <- list()
for (scope in c("all_eligible_years","supported_main_years")) for (source in c("STARS","IPUMS")) {
  z <- ipums_stars[ipums_stars$SUIC==if(source=="STARS")1 else 0,,drop=FALSE]
  if(scope=="supported_main_years")z <- z[as.character(z$YEAR) %in% common_years,,drop=FALSE]
  z$RETAINED <- rownames(z) %in% rownames(common)
  variables <- c("SEX","AGE_ADULT","YEAR","RACE","VETSTAT","FARM","EDUC_ATTAIN","DIVSEP")
  if(source=="STARS") {
    z$COUNTY <- clean_text(stars_csv_raw$DEM_RES03CountyCorrected[z$SOURCE_ROW])
    z$METHOD <- clean_text(stars_csv_raw$COD_Means[z$SOURCE_ROW])
    variables <- c(variables,"COUNTY","METHOD")
  }
  for(status in c(TRUE,FALSE)) {
    d <- z[z$RETAINED==status,,drop=FALSE]
    selection_age[[length(selection_age)+1]] <- data.frame(scope=scope,source=source,
      status=if(status)"Retained" else "Excluded",n=nrow(d),weight_total=sum(d$weight),
      mean_age=mean(d$AGE),sd_age=sd(d$AGE),median_age=median(d$AGE))
    for(v in variables) {
      values <- as.character(d[[v]])
      values[is.na(values) | values=="" | values=="NA"] <- "Missing"
      for(cat in sort(unique(values))) {
        n <- sum(values==cat);wt <- sum(d$weight[values==cat])
        selection_profiles[[length(selection_profiles)+1]] <- data.frame(scope=scope,source=source,
          status=if(status)"Retained" else "Excluded",variable=v,category=cat,
          n=n,records=nrow(d),percent_records=100*n/nrow(d),weight_sum=wt,
          percent_weight=100*wt/sum(d$weight))
      }
    }
  }
}
profile <- do.call(rbind,selection_profiles)
write_table(profile,"retained_excluded_profiles_long.csv")
write_table(do.call(rbind,selection_age),"retained_excluded_age_summary.csv")
case <- profile[profile$source=="STARS",]
keys <- unique(case[c("scope","variable","category")]);wide <- keys
for(status in c("Retained","Excluded")) {
  q <- case[case$status==status,]
  index <- match(paste(keys$scope,keys$variable,keys$category,sep="|"),paste(q$scope,q$variable,q$category,sep="|"))
  wide[[status]] <- ifelse(is.na(index),"0 (0.0%)",sprintf("%d (%.1f%%)",q$n[index],q$percent_records[index]))
}
write_table(wide,"retained_excluded_STARS_profiles.csv")
# Same four specifications, fixed complete-case sample within each sensitivity table.
fit_sequence <- function(data,prefix,race=FALSE) {
  block <- if(race)" + RACE" else ""
  fs <- c("SUIC ~ SEX + AGE_ADULT", "SUIC ~ SEX + AGE_ADULT + VETSTAT",
    "SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM",
    "SUIC ~ SEX + AGE_ADULT + VETSTAT + FARM + EDUC_ATTAIN + DIVSEP")
  result <- lapply(seq_along(fs),function(i)fit_logistic(as.formula(paste0(fs[i],block)),data,
    paste0(prefix,"_model",i)))
  names(result) <- paste0("Model_",1:4)
  publish_models(result,prefix,race=race)
  result
}
selection_late <- droplevels(common[as.character(common$YEAR) %in% c("2017","2018","2019","2021"),])
selection_late_models <- fit_sequence(selection_late,"late_years")
selection_age25 <- droplevels(common[common$AGE>=25,])
selection_age25_models <- fit_sequence(selection_age25,"age_25_plus")
selection_race25 <- droplevels(race_common[race_common$AGE>=25,])
selection_race25_models <- fit_sequence(selection_race25,"age_25_plus_with_race",race=TRUE)
minimal_available <- fit_logistic(SUIC~SEX+AGE_ADULT,ipums_stars,"minimal_available_model")
publish_models(list(Largest_available=minimal_available,Model4_common_sample=model1),
  "sample_selection_comparison",same_sample=FALSE,
  specification=c("Both columns adjust only for sex, adult age categories and categorical year.",
    "Largest_available uses all records complete on those variables in supported years.",
    "Model4_common_sample uses the identical records selected for principal Model 4.",
    "Coefficient differences therefore reflect record selection, including potential changes in year composition."))
selection_fits <- list(Main_Model4=ipums_norace,Late_years_Model4=selection_late_models[[4]],
  Age25plus_Model4=selection_age25_models[[4]],Age25plus_race_Model4=selection_race25_models[[4]],
  Minimal_largest_sample=minimal_available,Minimal_fixed_sample=model1)
selection_summary <- do.call(rbind,lapply(names(selection_fits),function(nm) {
  m <- selection_fits[[nm]];f <- model.frame(m)
  data.frame(model=nm,n=nrow(f),cases=sum(f$SUIC==1),comparison_records=sum(f$SUIC==0),
    age_reference=levels(f$AGE_ADULT)[1],minimum_age=min(ipums_stars$AGE[match(rownames(f),rownames(ipums_stars))]),
    years=paste(levels(f$YEAR),collapse=";"))
}))
write_table(selection_summary,"selection_sensitivity_samples.csv")
report <- c("Missingness and analytical-sample sensitivities",
 "Retained/excluded profiles describe eligible adults in 2012-2019 and 2021.",
 "Two scopes separate exclusions across all candidate years from exclusions within supported principal years.",
 "STARS percentages use record counts; IPUMS exports also include PERWT-weighted percentages.",
 "Missing categories are reported explicitly. These are descriptive selection diagnostics, not tests of missing-at-random assumptions.",
 "Late-year models use 2017, 2018, 2019 and 2021 without requiring race to be observed.",
 "Age25plus tables use age 25-34 as reference; their age odds ratios cannot be directly compared with an age-18-24 reference.",
 "All four models share records within each restricted table.",
 "The minimal-model comparison changes the sample while keeping sex, age and categorical-year specifications fixed.",
 "Intervals and p-values use combined case-process and ACS SDR covariance; no imputation is performed.","")
for(i in seq_len(nrow(selection_summary)))report <- c(report,sprintf("%s: cases=%d; comparisons=%d; age reference=%s; years=%s",
 selection_summary$model[i],selection_summary$cases[i],selection_summary$comparison_records[i],
 selection_summary$age_reference[i],selection_summary$years[i]))
writeLines(report,result_path(output_dir,"selection_sensitivity_report.txt"))
