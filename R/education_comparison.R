# Joint-indicator comparisons of the two binary education groupings.
source("R/common.R",local=TRUE)
output_dir <- file.path(results_root(),"regressions","education_comparison")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source("R/regression_education_inputs.R",local=TRUE)
source("R/replicate_inference.R",local=TRUE)
sample_counts <- missingness <- model_checks <- model_coverage <- model_warnings <- list()
source("R/regression_helpers.R",local=TRUE)
stopifnot(!any(ipums_stars$EDUC_PARTIAL=="PARTIAL" & ipums_stars$EDUC_ATTAIN=="ASSOCIATE_OR_HIGHER",na.rm=TRUE))
base_variables <- c("SUIC","SEX","AGE_ADULT","YEAR","VETSTAT","FARM","DIVSEP","EDUC_ATTAIN","EDUC_PARTIAL","weight")
tests <- list();contrasts <- list();checks <- list();groups <- list()
for(race in c(FALSE,TRUE)) {
 scope <- if(race) "with_race" else "without_race"
 complete <- complete.cases(ipums_stars[c(base_variables,if(race)"RACE")])
 tab <- table(ipums_stars$YEAR[complete],factor(ipums_stars$SUIC[complete],levels=0:1))
 years <- rownames(tab)[rowSums(tab>0)==2]
 common <- droplevels(ipums_stars[complete & as.character(ipums_stars$YEAR) %in% years,])
 primary <- dirname(output_dir)
 current_paths <- c(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"),ipums_path,codebook_path)
 for(location in c(primary,result_path(primary,"partial_education"))) {
  old <- read_table(file.path(location,"input_checksums.csv"))
  stopifnot(identical(unname(tools::md5sum(current_paths)),old$md5))
 }
 principal_name <- if(race) "principal_race_model4" else "principal_model4"
 partial_name <- paste0("sequential_partial_education",if(race)"_with_race","_model4")
 bm <- read_table(result_path(primary,"model_checks.csv"));pm <- read_table(result_path(primary,"partial_education","model_checks.csv"))
 selected <- rbind(bm[bm$model==principal_name,],pm[pm$model==partial_name,])
 stopifnot(nrow(selected)==2,all(selected$cases==sum(common$SUIC==1)),all(selected$comparison_records==sum(common$SUIC==0)),all(selected$residence_verified))
 name <- paste0("joint_education_",scope)
 model <- fit_logistic(reformulate(c("SEX","AGE_ADULT","VETSTAT","FARM","DIVSEP","EDUC_ATTAIN","EDUC_PARTIAL",if(race)"RACE"),"SUIC"),common,name)
 b <- coef(model);V <- vcov(model)
 terms <- c("EDUC_PARTIALPARTIAL","EDUC_ATTAINASSOCIATE_OR_HIGHER")
 questions <- c("Partial attainment adds beyond the degree split","Degree attainment adds beyond the partial split")
 scope_tests <- do.call(rbind,lapply(1:2,function(i) {
  t <- terms[i];se <- sqrt(V[t,t]);chi <- b[t]^2/V[t,t]
  data.frame(analysis=scope,test=questions[i],term=t,chi_square=unname(chi),df=1,
   p=unname(pchisq(chi,1,lower.tail=FALSE)),odds_ratio=unname(exp(b[t])),
   ci_lower=unname(exp(b[t]-qnorm(.975)*se)),ci_upper=unname(exp(b[t]+qnorm(.975)*se)),
   cases=sum(common$SUIC==1),comparison_records=sum(common$SUIC==0),years=paste(years,collapse=";"))
 }))
 scope_tests$p_holm <- p.adjust(scope_tests$p,"holm");tests[[scope]] <- scope_tests
 for(comparison in c("Partial vs other below-associate","Degree vs other below-associate","Partial vs degree")) {
  L <- setNames(rep(0,length(b)),names(b))
  if(comparison!="Degree vs other below-associate")L[terms[1]] <- 1
  if(comparison=="Degree vs other below-associate")L[terms[2]] <- 1
  if(comparison=="Partial vs degree")L[terms[2]] <- -1
  estimate <- sum(L*b);se <- sqrt(as.numeric(t(L)%*%V%*%L))
  contrasts[[length(contrasts)+1]] <- data.frame(analysis=scope,comparison=comparison,odds_ratio=exp(estimate),ci_lower=exp(estimate-qnorm(.975)*se),ci_upper=exp(estimate+qnorm(.975)*se))
 }
 group <- ifelse(common$EDUC_PARTIAL=="PARTIAL","PARTIAL",ifelse(common$EDUC_ATTAIN=="ASSOCIATE_OR_HIGHER","DEGREE","OTHER_BELOW_ASSOCIATE"))
 for(source in 0:1)for(g in unique(group)) {
  d <- common[common$SUIC==source & group==g,]
  groups[[length(groups)+1]] <- data.frame(analysis=scope,source=c("IPUMS","STARS")[source+1],education_group=g,records=nrow(d),weight_sum=sum(d$weight))
 }
 common$EDUC_THREE <- factor(group,levels=c("OTHER_BELOW_ASSOCIATE","DEGREE","PARTIAL"))
 equivalent <- suppressWarnings(glm(reformulate(c("SEX","AGE_ADULT","VETSTAT","FARM","DIVSEP","EDUC_THREE",if(race)"RACE","YEAR"),"SUIC"),data=common,weights=weight,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100)))
 eta_error <- max(abs(predict(model,type="link")-predict(equivalent,type="link")))
 stopifnot(equivalent$converged,eta_error<1e-8)
 checks[[scope]] <- data.frame(analysis=scope,same_samples_as_both_binary_models=TRUE,current_input_checksums_match=TRUE,three_group_linear_predictor_error=eta_error,passed=TRUE)
}
result <- do.call(rbind,tests);group_results <- do.call(rbind,contrasts)
write_table(result,"education_incremental_tests.csv")
write_table(group_results,"education_group_contrasts.csv")
publication <- result[c("analysis","test")]
publication[["OR (95% CI)"]] <- sprintf("%.3f (%.3f-%.3f)",result$odds_ratio,result$ci_lower,result$ci_upper)
publication[["Wald chi-square (df=1)"]] <- round(result$chi_square,3)
publication[["Raw p"]] <- format.pval(result$p,digits=4,eps=1e-16)
publication[["Holm p"]] <- format.pval(result$p_holm,digits=4,eps=1e-16)
write_table(publication,"education_incremental_tests_formatted.csv")
write_table(do.call(rbind,groups),"education_group_sample_counts.csv")
write_table(do.call(rbind,checks),"comparison_checks.csv")
for(item in c("sample_counts","missingness","model_checks","model_coverage","model_warnings"))write_table(do.call(rbind,get(item)),paste0(item,".csv"))
write_table(do.call(rbind,replicate_checks),"replicate_model_checks.csv")
write_table(do.call(rbind,replicate_changes),"covariance_method_comparison.csv")
stopifnot(all(do.call(rbind,model_checks)$passed),all(do.call(rbind,replicate_checks)$passed))
paths <- c(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"),ipums_path,codebook_path)
write_table(data.frame(input=paths,md5=unname(tools::md5sum(paths))),"input_checksums.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
notes <- c("Joint education-indicator comparison, confirmed-Montana-resident adult samples.",
 "Both indicators enter the fully adjusted model: sex, adult age, veteran, farmer/rancher, divorced/separated and categorical year; race is added in the companion model.",
 "Age reference: 18-24. Education tests are invariant to changing the age reference.",
 "Education reference: other recorded below-associate attainment (STARS NO HIGH SCHOOL or HIGH SCHOOL OR GED; matching IPUMS below-grade-9 or high-school/GED codes).",
 "Partial: STARS SOME HIGH SCHOOL or SOME COLLEGE; IPUMS EDUCD 30/40/50/61/65/71. Degree: associate or higher.",
 "The mutually exclusive indicators create a three-group education classification.",
 "Testing the partial coefficient at zero tests whether partial adds beyond the degree-only model, which pools partial with other below-associate.",
 "Testing the degree coefficient at zero tests whether degree adds beyond the partial-only model, which pools degree with other below-associate.",
 "Each is a 1-df Wald test using combined Poisson mortality-event and year-blocked ACS SDR covariance; 80 signed replicates per annual sample, scaling 4/80.",
 "Raw and Holm-adjusted p-values are reported. The two incremental tests form one family within each sequence; race-inclusive tests are a separate sensitivity family.",
 "The two indicator ORs share the other-below-associate reference and differ from their binary-model contrasts.",
 "The partial-versus-degree contrast is descriptive with a pointwise 95% CI; contrast CIs are not multiplicity-adjusted.",
 "Sample equality, checksum consistency, resident eligibility and equivalence to a three-level factor are checked.",
 "These tests assess additional statistical information in recorded categories. They cannot establish causal effects or distinguish real patterns from recording/classification changes.")
writeLines(notes,result_path(output_dir,"analysis_notes.txt"))
writeLines(c(notes,"","Incremental tests:",capture.output(print(publication,row.names=FALSE)),"","Three-group contrasts:",capture.output(print(group_results,row.names=FALSE))),result_path(output_dir,"analysis_report.txt"))
print(publication,row.names=FALSE)
message("Education comparison results: ",output_dir)
