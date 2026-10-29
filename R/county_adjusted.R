# County mortality standardized to a common Montana age/year or age/sex/year mix.
# Uses public aggregate cells and the shared clustering output root.
ca_assign <- read_table(result_path(output_dir,"county_k3_assignments.csv"))
ca_rates <- list();ca_tests <- list();ca_check <- list()
for (ca_kind in c("age","age_sex")) {
  ca_d <- read_table(file.path("data/processed/clustering",
    paste0("county_year_",ca_kind,"_mortality.csv")))
  ca_d$CLUSTER <- factor(ca_assign$CLUSTER[match(ca_d$COUNTY,ca_assign$COUNTY)],levels=1:3)
  ca_d$YEAR <- factor(ca_d$YEAR)
  ca_ages <- c("<15","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+")
  ca_d$AGE_GROUP <- factor(ca_d$AGE_GROUP,levels=ca_ages)
  ca_d$STRATUM <- if(ca_kind=="age")ca_d$AGE_GROUP else interaction(ca_d$AGE_GROUP,ca_d$SEX,drop=TRUE)
  stopifnot(!anyNA(ca_d$CLUSTER),length(unique(ca_d$COUNTY))==56,
    all(ca_d$POPULATION>=0),all(ca_d$DEATHS[ca_d$POPULATION==0]==0))
  ca_d <- droplevels(ca_d[ca_d$POPULATION>0,])
  ca_fit <- glm(DEATHS~CLUSTER+STRATUM+YEAR+offset(log(POPULATION)),data=ca_d,family=poisson())
  ca_X <- model.matrix(ca_fit);ca_mu <- fitted(ca_fit);ca_b <- coef(ca_fit)
  ca_score <- rowsum(ca_X*as.numeric(ca_d$DEATHS-ca_mu),ca_d$COUNTY,reorder=FALSE)
  ca_G <- nrow(ca_score);ca_n <- nrow(ca_X);ca_p <- ncol(ca_X)
  ca_bread <- solve(crossprod(ca_X,ca_X*as.numeric(ca_mu)))
  ca_V <- ca_bread%*%crossprod(ca_score)%*%ca_bread * ca_G/(ca_G-1)*(ca_n-1)/(ca_n-ca_p)
  ca_idx <- grep("^CLUSTER",names(ca_b));ca_q <- length(ca_idx)
  ca_F <- as.numeric(t(ca_b[ca_idx])%*%solve(ca_V[ca_idx,ca_idx,drop=FALSE],ca_b[ca_idx])/ca_q)
  ca_pvalue <- pf(ca_F,ca_q,ca_G-1,lower.tail=FALSE)
  ca_tests[[ca_kind]] <- data.frame(adjustment=ca_kind,F=ca_F,df1=ca_q,df2=ca_G-1,p=ca_pvalue,
    first_year=min(as.integer(as.character(ca_d$YEAR))),last_year=max(as.integer(as.character(ca_d$YEAR))),
    n_cells=nrow(ca_d),deaths=sum(ca_d$DEATHS),person_years=sum(ca_d$POPULATION),
    se_method="County-clustered sandwich with finite-sample correction; t/F with 55 df")
  # Internal direct standard: pooled Montana person-years by stratum and year.
  ca_std <- aggregate(POPULATION~STRATUM+YEAR,ca_d,sum)
  ca_weights <- ca_std$POPULATION/sum(ca_std$POPULATION)
  for(ca_g in 1:3) {
    ca_new <- ca_std;ca_new$CLUSTER <- factor(ca_g,levels=1:3);ca_new$POPULATION <- 1
    ca_M <- model.matrix(delete.response(terms(ca_fit)),ca_new,
      contrasts.arg=ca_fit$contrasts,xlev=ca_fit$xlevels)
    ca_pred <- exp(as.numeric(ca_M%*%ca_b));ca_rate <- sum(ca_weights*ca_pred)
    ca_grad <- colSums(ca_M*as.numeric(ca_weights*ca_pred))
    ca_se_log <- sqrt(as.numeric(t(ca_grad)%*%ca_V%*%ca_grad))/ca_rate
    ca_lo <- exp(log(ca_rate)-qt(.975,ca_G-1)*ca_se_log)*1e5
    ca_hi <- exp(log(ca_rate)+qt(.975,ca_G-1)*ca_se_log)*1e5
    ca_part <- ca_d[ca_d$CLUSTER==ca_g,]
    ca_rates[[paste(ca_kind,ca_g)]] <- data.frame(adjustment=ca_kind,cluster=ca_g,
      standardized_rate_per_100000=ca_rate*1e5,ci_lower=ca_lo,ci_upper=ca_hi,
      matched_crude_rate_per_100000=sum(ca_part$DEATHS)/sum(ca_part$POPULATION)*1e5,
      deaths=sum(ca_part$DEATHS),person_years=sum(ca_part$POPULATION),
      first_year=min(as.integer(as.character(ca_d$YEAR))),last_year=max(as.integer(as.character(ca_d$YEAR))))
  }
  ca_se <- sqrt(diag(ca_V))
  write_table(data.frame(term=names(ca_b),coefficient=ca_b,std_error=ca_se,
    rate_ratio=exp(ca_b),ci_lower=exp(ca_b-qt(.975,55)*ca_se),
    ci_upper=exp(ca_b+qt(.975,55)*ca_se),p=2*pt(abs(ca_b/ca_se),55,lower.tail=FALSE)),
    paste0("county_k3_",ca_kind,"_adjusted_coefficients.csv"))
  ca_check[[ca_kind]] <- data.frame(analysis=ca_kind,passed=ca_fit$converged &&
    ca_fit$rank==length(ca_b) && all(is.finite(ca_b)) && all(is.finite(ca_V)) && all(diag(ca_V)>0),
    counties=ca_G,cells=nrow(ca_d),deaths=sum(ca_d$DEATHS))
}
ca_test <- do.call(rbind,ca_tests);ca_test$p_holm <- p.adjust(ca_test$p,"holm")
write_table(ca_test,"county_k3_adjusted_tests.csv")
write_table(do.call(rbind,ca_rates),"county_k3_adjusted_rates.csv")
write_table(do.call(rbind,ca_check),"county_adjusted_validation.csv")
stopifnot(all(vapply(ca_check,function(x)x$passed,logical(1))))
writeLines(c("Age-adjusted county rates: 2012-2021, known ages only.",
  "Age-and-sex-adjusted county rates: 2013-2021; 2012 has no recorded sex.",
  "Age bands: <15, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75-84, 85+.",
  "Population totals exactly reconcile to the denominators used for the crude county analysis.",
  "Each model includes county cluster, age (or joint age-sex stratum), categorical year and log population offset.",
  "Rates are model-standardized to the same pooled Montana stratum/year population distribution for all clusters.",
  "These are internally standardized rates, not rates standardized to the US 2000 population.",
  "Intervals use delta-method log-rate uncertainty, county-clustered sandwich with finite-sample correction covariance and t(55) limits.",
  "Omnibus comparisons use county-clustered F tests, with Holm adjustment across the two adjusted sensitivity analyses.",
  "Missing age/sex records are excluded from the relevant numerator; coverage accounting is retained separately.",
  "The adjusted models assume common cluster rate ratios across demographic strata and years.",
  "Population-estimate error and cluster-construction uncertainty are not propagated."),
  result_path(output_dir,"county_adjusted_notes.txt"))

cat("\nAdditional county demographic-standardized sensitivity tests:\n",
  paste(sprintf("%s: Holm p=%.6g; years %d-%d",ca_test$adjustment,ca_test$p_holm,
    ca_test$first_year,ca_test$last_year),collapse="\n"),"\nSee county_k3_adjusted_rates.csv and county_adjusted_notes.txt.\n",
  file=result_path(output_dir,"analysis_report.txt"),append=TRUE)
