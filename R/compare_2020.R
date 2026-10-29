# Exploratory comparison of 2020 STARS decedents with other study years.
# No population denominator or IPUMS comparator is used.
parent_2020_output_dir <- output_dir
output_dir <- result_path(output_dir,"data_quality","comparison_2020")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
raw2020 <- read_stars_csv()
d2020 <- harmonize_stars(raw2020)
d2020$EDUC_DETAILS <- clean_text(raw2020$DEM_EDU1EducationAchieved)
d2020$COUNTY <- county_name(raw2020$DEM_RES03CountyCorrected)
d2020$MONTH <- raw2020$DEATH_4MonthNum
d2020$WEEKDAY <- clean_text(raw2020$DEATH_5Day)
fields <- c("SEX","AGE7","EDUC_DETAILS","MARST2","VETSTAT","FARM","RACE","COD","COUNTY")
missing_fields <- c(fields,"AGE","MONTH","WEEKDAY")
main_fields <- c("VETSTAT","SEX","AGE7","DIVSEP","EDUC","FARM")
d2020$MAIN_COMPLETE <- complete.cases(d2020[main_fields])
d2020$RACE_COMPLETE <- complete.cases(d2020[c(main_fields,"RACE")])
yearly <- do.call(rbind,lapply(2012:2021,function(y) {
  z <- d2020[d2020$DYEAR==y,]
  data.frame(year=y,n=nrow(z),age_known=sum(!is.na(z$AGE)),
    age_mean=mean(z$AGE,na.rm=TRUE),age_sd=sd(z$AGE,na.rm=TRUE),
    main_complete=sum(z$MAIN_COMPLETE),race_complete=sum(z$RACE_COMPLETE))
}))
write_table(yearly,"annual_summary.csv")
annual_categories <- list()
for (y in 2012:2021) for (v in fields) {
  values <- clean_text(d2020[[v]][d2020$DYEAR==y])
  levels_all <- sort(unique(clean_text(d2020[[v]])))
  levels_all <- levels_all[!is.na(levels_all)]
  freq <- table(factor(values,levels=levels_all))
  annual_categories[[length(annual_categories)+1]] <- data.frame(
    year=y, variable=v, category=names(freq), n=as.integer(freq),
    known=sum(!is.na(values)), percent=if(sum(!is.na(values)))100*as.integer(freq)/sum(!is.na(values)) else NA_real_)
}
annual_categories <- do.call(rbind,annual_categories)
write_table(annual_categories,"category_counts_by_year.csv")

year_missing <- do.call(rbind,lapply(2012:2021,function(y) {
  z <- d2020[d2020$DYEAR==y,]
  do.call(rbind,lapply(missing_fields,function(v) data.frame(
    year=y,variable=v,n=nrow(z),missing=sum(is.na(z[[v]])),
    percent_missing=100*mean(is.na(z[[v]])))))
}))
write_table(year_missing,"missingness_by_year.csv")
contrasts <- list(all_other_years=c(2012:2019,2021L),
  comparable_coverage=c(2017:2019,2021L),recent_pre2020=2017:2019, adjacent_2021=2021L)
tests <- list(); categories <- list(); miss <- list(); ages <- list()
set.seed(2020)
for (contrast in names(contrasts)) {
  z <- d2020[d2020$DYEAR %in% c(2020,contrasts[[contrast]]),]
  group <- factor(ifelse(z$DYEAR==2020,"2020","comparison"),levels=c("2020","comparison"))
  a <- z$AGE[group=="2020"]; b <- z$AGE[group=="comparison"]
  a <- a[!is.na(a)];b <- b[!is.na(b)]
  age_test <- t.test(a,b)
  pooled <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
  ages[[length(ages)+1]] <- data.frame(contrast=contrast,n_2020=length(a),n_other=length(b),
    mean_2020=mean(a),sd_2020=sd(a),mean_other=mean(b),sd_other=sd(b),
    difference=mean(a)-mean(b),ci_lower=age_test$conf.int[1],ci_upper=age_test$conf.int[2],
    standardized_difference=(mean(a)-mean(b))/pooled,p=age_test$p.value)
  tests[[length(tests)+1]] <- data.frame(contrast=contrast,variable="AGE",
    n_2020=length(a),n_other=length(b),statistic=unname(age_test$statistic),
    df=unname(age_test$parameter),p=age_test$p.value,cramers_v=NA_real_,
    method="Welch t-test",minimum_expected=NA_real_,simulation_B=0)
  for (v in fields) {
    tab <- table(group,clean_text(z[[v]]))
    tab <- tab[,colSums(tab)>0,drop=FALSE]
    base <- suppressWarnings(chisq.test(tab,correct=FALSE))
    sparse <- any(base$expected<1) || mean(base$expected<5)>.2
    fit <- if(sparse) chisq.test(tab,simulate.p.value=TRUE,B=49999) else base
    tests[[length(tests)+1]] <- data.frame(contrast=contrast,variable=v,
      n_2020=sum(tab[1,]),n_other=sum(tab[2,]),statistic=unname(base$statistic),
      df=unname(base$parameter),p=fit$p.value,
      cramers_v=sqrt(unname(base$statistic)/sum(tab)),method=fit$method,
      minimum_expected=min(base$expected),simulation_B=if(sparse)49999 else 0)
    for (category in colnames(tab)) {
      n1 <- tab[1,category]; n0 <- tab[2,category]
      N1 <- sum(tab[1,]); N0 <- sum(tab[2,])
      categories[[length(categories)+1]] <- data.frame(contrast=contrast,variable=v,
        category=category,n_2020=n1,known_2020=N1,percent_2020=100*n1/N1,
        n_other=n0,known_other=N0,percent_other=100*n0/N0,
        percentage_point_difference=100*(n1/N1-n0/N0))
    }
  }
  for (v in missing_fields) {
    tab <- table(group,factor(is.na(z[[v]]),levels=c(FALSE,TRUE)))
    pv <- if(any(colSums(tab)==0)) 1 else fisher.test(tab)$p.value
    miss[[length(miss)+1]] <- data.frame(contrast=contrast,variable=v,
      n_2020=sum(group=="2020"),missing_2020=tab[1,2],
      percent_missing_2020=100*tab[1,2]/sum(tab[1,]),
      n_other=sum(group=="comparison"),missing_other=tab[2,2],
      percent_missing_other=100*tab[2,2]/sum(tab[2,]),p=pv)
  }
}
tests <- do.call(rbind,tests);categories <- do.call(rbind,categories)
miss <- do.call(rbind,miss); ages <- do.call(rbind,ages)
tests$p_holm <- ave(tests$p,tests$contrast,FUN=function(x)p.adjust(x,"holm"))
miss$p_holm <- ave(miss$p,miss$contrast,FUN=function(x)p.adjust(x,"holm"))
write_table(tests,"distribution_tests.csv")
write_table(categories,"category_comparisons.csv")
write_table(miss,"missingness_comparisons.csv")
write_table(ages,"age_comparisons.csv")
checks <- data.frame(check=c("all 2764 source records accounted for","298 records in 2020",
  "four contrasts with ten tests each","all p-values finite and within bounds",
  "Holm correction does not reduce p-values","category totals reconcile",
  "no 2020 records in comparison year sets"),
  passed=c(sum(yearly$n)==2764,yearly$n[yearly$year==2020]==298,
    nrow(tests)==40 && all(table(tests$contrast)==10),
    all(is.finite(tests$p) & tests$p>=0 & tests$p<=1),
    all(tests$p_holm>=tests$p),
    all(aggregate(n_2020~contrast+variable,categories,sum)$n_2020 ==
      tests$n_2020[match(paste(aggregate(n_2020~contrast+variable,categories,sum)$contrast,
        aggregate(n_2020~contrast+variable,categories,sum)$variable),
        paste(tests$contrast,tests$variable))]),
    all(!vapply(contrasts,function(x)2020 %in% x,logical(1)))))
write_table(checks,"validation_checks.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
source2020 <- primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv")
write_table(data.frame(input=source2020,md5=unname(tools::md5sum(source2020))),"input_checksum.csv")
lines <- c("Exploratory 2020 STARS comparison",
  "Source: 20230831_Montana_Suicide_Database.csv.",
  paste("Run time (UTC):",format(Sys.time(),tz="UTC",usetz=TRUE)),
  "298 records in 2020; 2466 records in the other nine years (mean 274 per year).",
  "This compares recorded decedent composition and missingness, not mortality rates.",
  "All-other-years comparison: 2012-2019 and 2021.",
  "Comparable-coverage comparison: 2017-2019 and 2021 (race is available from 2017 in this file).",
  "Recent pre2020 comparison: 2017-2019; avoids including 2021 as a comparator.",
  "Adjacent-year comparison: 2021 only; assesses whether recorded shifts persist after 2020.",
  "Ten distribution tests per contrast: age mean plus nine categorical distributions.",
  "Holm adjustment is applied across these ten tests separately within each contrast.",
  "Missingness tests form a separate Holm-adjusted family per contrast.",
  "Sparse category tables use fixed-margin Monte Carlo chi-square tests (49999 replicates;",
  "seed 2020). Other category tests use Pearson chi-square without continuity correction.",
  "Age uses Welch's t-test and a 95% confidence interval for the mean difference.",
  "Category percentages exclude missing values; missingness is reported separately.",
  "No differences being detected does not establish equivalence or absence of pandemic effects.",
  "These unadjusted exploratory comparisons cannot identify causal effects of COVID-19.",
  "Marital status uses divorced/separated/other recorded categories.",
  "County comparisons use recorded county names, not outcome-derived county clusters.",
  "The primary descriptive comparison is all other years; comparable-coverage and",
  "recent-pre2020 contrasts assess sensitivity to changes in data collection and period.",
  "", "Distribution test results (Holm p < .05):")
for (contrast in names(contrasts)) {
  hit <- tests$variable[tests$contrast==contrast & tests$p_holm<.05]
  lines <- c(lines,paste(contrast,":",if(length(hit))paste(hit,collapse=", ") else "none"))
}
lines <- c(lines,"","Missingness differences (Holm p < .05):")
for (contrast in names(contrasts)) {
  hit <- miss$variable[miss$contrast==contrast & miss$p_holm<.05]
  lines <- c(lines,paste(contrast,":",if(length(hit))paste(hit,collapse=", ") else "none"))
}
zero_hs <- annual_categories$year[
  annual_categories$variable=="EDUC_DETAILS" &
  annual_categories$category=="HIGH SCHOOL OR GED" &
  annual_categories$known>0 & annual_categories$n==0]
lines <- c(lines,"","Interpretation requiring coding verification:",
  paste("No HIGH SCHOOL OR GED records despite known education in years:",
    paste(zero_hs,collapse=", ")),
  "Recorded education and separation proportions change sharply around 2020 and",
  "similar patterns appear in 2021. These could reflect category harmonization",
  "or collection differences; pandemic explanations are not established.",
  "The comparison uses recorded labels and the documented demographic groupings.")
writeLines(lines,result_path(output_dir,"analysis_report.txt"))
note <- c("Suggested paper note",
  "",
  "The STARS dataset contained 298 deaths in 2020. In exploratory comparisons",
  "with the other study years, recorded education and marital-status distributions",
  "differed after Holm adjustment for ten distribution tests (both adjusted p<0.001).",
  "Among records with known values, high-school/GED completion was recorded for",
  "119/295 (40.3%) in 2020 versus 182/2205 (8.3%) in other years, and separation",
  "for 99/296 (33.4%) versus 146/2230 (6.5%). No other tested distribution",
  "differences remained significant after adjustment. Similar category patterns",
  "were present in 2021. Abrupt year-to-year changes in category use indicate",
  "that temporal consistency of education and marital-status coding requires",
  "verification before interpreting these differences as pandemic-related.",
  "These comparisons describe decedent composition rather than population",
  "mortality rates, and non-significance does not establish equivalence.")
writeLines(note,result_path(output_dir,"paper_note.txt"))
cat(paste(lines,collapse="\n"),"\n")
if(!all(checks$passed))stop("2020 comparison validation failed")
output_dir <- parent_2020_output_dir
