# Exact reparameterization of resident-only four-model sequences.
# Coefficients AND covariance components/ACS replicate estimates are transformed.
source("R/common.R",local=TRUE)
base <- file.path(results_root(),"regressions")
output_dir <- result_path(base,"age_reference_45_54")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
checks <- list();inputs <- character()
ages <- c("18-24","25-34","35-44","45-54","55-64","65-74","75-84","85+")
star <- function(p)ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))
read_input <- function(path) {
  inputs <<- unique(c(inputs,path))
  read_table(path)
}
transform_model <- function(source_dir,old_name,new_name) {
  old <- read_input(result_path(source_dir,paste0(old_name,".csv")))
  bn <- old$term;b <- setNames(old$coefficient,bn)
  ref <- "AGE_ADULT45-54";newn <- sub(ref,"AGE_ADULT18-24",bn,fixed=TRUE)
  stopifnot(ref %in% bn,!"AGE_ADULT18-24" %in% bn,!anyDuplicated(newn))
  A <- diag(length(b));dimnames(A) <- list(newn,bn)
  A["(Intercept)",ref] <- 1
  for(t in grep("^AGE_ADULT",newn,value=TRUE)) {
    A[t,] <- 0
    if(t!="AGE_ADULT18-24")A[t,t] <- 1
    A[t,ref] <- -1
  }
  nb <- setNames(as.numeric(A%*%b),newn)
  V <- list()
  for(component in c("case_process","acs_sampling","combined")) {
    tab <- read_input(result_path(source_dir,"replicate_details",paste0(old_name,"_",component,"_covariance.csv")))
    M <- as.matrix(tab[match(bn,tab$term),bn,drop=FALSE]);rownames(M) <- bn
    stopifnot(all(is.finite(M)))
    if(component=="combined")stopifnot(max(abs(sqrt(diag(M))-old$std_error))<1e-10)
    V[[component]] <- A%*%M%*%t(A)
    write_table(data.frame(term=newn,V[[component]],check.names=FALSE),
      paste0(new_name,"_",component,"_covariance.csv"),result_path(output_dir,"replicate_details"))
  }
  stopifnot(max(abs(V$combined-V$case_process-V$acs_sampling))<1e-10)
  rep <- read_input(result_path(source_dir,"replicate_details",paste0(old_name,"_replicate_coefficients.csv")))
  rb <- as.matrix(rep[bn])%*%t(A)
  write_table(data.frame(rep[c("year","replicate")],rb,check.names=FALSE),
    paste0(new_name,"_replicate_coefficients.csv"),result_path(output_dir,"replicate_details"))
  # Independently recompute year-blocked SDR covariance after transformation.
  sdr <- matrix(0,length(b),length(b))
  for(y in unique(rep$year)) {
    selected <- rep$year==y;stopifnot(sum(selected)==80)
    diff <- sweep(rb[selected,,drop=FALSE],2,nb,"-")
    sdr <- sdr+(4/80)*crossprod(diff)
  }
  sdr_error <- max(abs(sdr-V$acs_sampling))
  # Linear predictors must agree for every age, with other coefficients held fixed.
  eta_error <- max(vapply(ages,function(age) {
    eta_old <- b["(Intercept)"]+if(age=="18-24")0 else b[paste0("AGE_ADULT",age)]
    eta_new <- nb["(Intercept)"]+if(age=="45-54")0 else nb[paste0("AGE_ADULT",age)]
    abs(eta_old-eta_new)
  },numeric(1)))
  other <- setdiff(bn,c("(Intercept)",grep("^AGE_ADULT",bn,value=TRUE)))
  other_error <- max(abs(nb[other]-b[other]))
  back_error <- max(abs(solve(A,nb)-b))
  passed <- sdr_error<1e-10 && eta_error<1e-10 && other_error<1e-10 && back_error<1e-10 &&
    all(diag(V$combined)>0) && all(is.finite(V$combined))
  stopifnot(passed)
  checks[[length(checks)+1]] <<- data.frame(model=new_name,source_model=old_name,
    age_reference="45-54",linear_predictor_error=eta_error,non_age_coefficient_error=other_error,
    back_transform_error=back_error,sdr_covariance_error=sdr_error,
    annual_replications=nrow(rep),passed=passed)
  se <- sqrt(diag(V$combined));p <- 2*pnorm(abs(nb/se),lower.tail=FALSE)
  result <- data.frame(term=newn,coefficient=nb,std_error=se,p=p,odds_ratio=exp(nb),
    ci_lower=exp(nb-qnorm(.975)*se),ci_upper=exp(nb+qnorm(.975)*se),row.names=NULL)
  write_table(result,paste0(new_name,".csv"))
  result
}
non_age_labels <- c(SEXFEMALE="Female (reference: male)",
  VETSTATAF="Ever in armed forces (reference: no)",FARMFARMER="Farmer/rancher (reference: other)",
  EDUC_ATTAINASSOCIATE_OR_HIGHER="Associate degree or higher (reference: less)",
  EDUC_PARTIALPARTIAL="Partial high school or college (reference: other attainment)",
  "DIVSEPDIVORCED/SEPARATED"="Divorced/separated (reference: other recorded status)",
  "RACEAMERICAN INDIAN"="American Indian (reference: White)",RACEOTHER="Other race/ethnicity (reference: White)")
for(partial in c(FALSE,TRUE))for(race in c(FALSE,TRUE)) {
  source_dir <- if(partial)result_path(base,"partial_education") else base
  oldstem <- paste0(if(partial)"sequential_partial_education" else "sequential_models",if(race)"_with_race")
  oldnames <- if(partial)paste0(oldstem,"_model",1:4) else if(race)
    c(paste0("principal_race_model",1:3),"principal_race_model4") else
    c(paste0("principal_model",1:3),"principal_model4")
  stem <- paste0(oldstem,"_age45_54")
  mc <- read_input(result_path(source_dir,"model_checks.csv"))
  mc <- mc[match(oldnames,mc$model),]
  stopifnot(nrow(mc)==4,!anyNA(mc$model),all(mc$residence_verified),
    all(mc$non_montana_cases==0),all(mc$unknown_residence_cases==0),
    length(unique(mc$cases))==1,length(unique(mc$comparison_records))==1)
  models <- lapply(1:4,function(i)transform_model(source_dir,oldnames[i],paste0(stem,"_model",i)))
  age_terms <- paste0("AGE_ADULT",setdiff(ages,"45-54"))
  final_terms <- models[[4]]$term
  terms <- c("SEXFEMALE",age_terms,setdiff(final_terms,c("(Intercept)","SEXFEMALE",
    grep("^AGE_ADULT|^YEAR",final_terms,value=TRUE))))
  labels <- c(setNames(paste("Age",setdiff(ages,"45-54"),"(reference: 45-54)"),age_terms),non_age_labels)
  stopifnot(all(terms %in% names(labels)),!"AGE_ADULT45-54" %in% final_terms)
  wide <- data.frame(Predictor=unname(labels[terms]));long <- list()
  for(i in 1:4) {
    d <- models[[i]];j <- match(terms,d$term)
    wide[[paste0("Model_",i)]] <- ifelse(is.na(j),"",sprintf("%.3f (%.3f-%.3f)%s",
      d$odds_ratio[j],d$ci_lower[j],d$ci_upper[j],star(d$p[j])))
    long[[i]] <- data.frame(model=paste0("Model_",i),d)
  }
  footer <- data.frame(Predictor=c("Age reference","Year adjusted","STARS cases","IPUMS comparisons","Total observations"))
  for(i in 1:4)footer[[paste0("Model_",i)]] <- c("45-54","Yes",mc$cases[i],mc$comparison_records[i],mc$model_rows[i])
  write_table(rbind(wide,footer),paste0(stem,".csv"))
  write_table(do.call(rbind,long),paste0(stem,"_long.csv"))
  notes_path <- result_path(source_dir,paste0(oldstem,"_notes.txt"))
  inputs <- unique(c(inputs,notes_path));notes <- readLines(notes_path,warn=FALSE)
  notes[grepl("^Age reference:",notes)] <- "Age reference: 45-54 for every age row. Year effects are included but not displayed."
  writeLines(c(notes,
    "Exact reparameterization of the confirmed-Montana-resident models; no records or covariates change.",
    "For age g, new log OR = original log OR(g vs 18-24) minus original log OR(45-54 vs 18-24).",
    "The intercept, complete covariance matrices, and every annual ACS replicate coefficient vector are transformed together.",
    "Fitted probabilities, likelihood, and non-age coefficients are unchanged; age p-values test contrasts with ages 45-54.",
    "See reference_change_checks.csv for exact-fit and replicated-covariance checks."),result_path(output_dir,paste0(stem,"_notes.txt")))
}
write_table(do.call(rbind,checks),"reference_change_checks.csv")
write_table(data.frame(input=inputs,md5=unname(tools::md5sum(inputs))),"input_checksums.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
message("Age-reference tables: ",output_dir,"; ",length(checks)," model transformations checked.")
