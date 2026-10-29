# Two-source covariance for mortality cases and annual ACS population comparisons.
# ACS: 80 SDR replicates per annual sample, varying one year at a time (independent
# annual-sample assumption). Cases are fixed during ACS replication.
# Death process: Poisson event-score sandwich. The two components are added.
# Signed replicate weights are preserved; no truncation or absolute values.
replicate_checks <- list();replicate_changes <- list()
replicate_vcov <- function(model,used,name) {
  stopifnot(identical(rownames(used),rownames(model.frame(model))))
  X <- model.matrix(model);b <- coef(model);y <- used$SUIC;w <- model$prior.weights
  stopifnot(all(w[y==1]==1),all(is.finite(X)),all(y %in% 0:1))
  key <- do.call(paste,c(as.data.frame(X),sep=":"));group <- match(key,unique(key))
  first <- !duplicated(group);A <- X[first,,drop=FALSE];G <- nrow(A);P <- ncol(A)
  case <- as.numeric(rowsum(as.numeric(y==1),group,reorder=FALSE))
  population <- as.numeric(rowsum(ifelse(y==0,w,0),group,reorder=FALSE))
  rr <- rowsum(as.matrix(ipums_raw[used$SOURCE_ROW[y==0],replicate_names]),group[y==0],reorder=FALSE)
  rep_population <- matrix(0,G,80);rep_population[as.integer(rownames(rr)),] <- rr
  year <- as.character(used$YEAR[first]);years <- sort(unique(year))
  fitted_p <- plogis(as.numeric(A%*%b))
  H <- crossprod(A,A*as.numeric((case+population)*fitted_p*(1-fitted_p)))
  bread <- chol2inv(chol(H))
  # Poisson event-process variation, conditional on the estimated comparison population.
  V_case <- bread%*%crossprod(A,A*as.numeric(case*(1-fitted_p)^2))%*%bread
  # Solve signed-weight score equations, initialized at the full-sample estimate.
  # Full-sample weights remain positive; signed totals occur only in variance replicates.
  score_fit <- function(total) {
    beta <- b
    for(iteration in 1:30) {
      prob <- plogis(as.numeric(A%*%beta))
      score <- as.numeric(crossprod(A,case-total*prob))
      info <- crossprod(A,A*as.numeric(total*prob*(1-prob)))
      R <- chol(info)
      step <- backsolve(R,forwardsolve(t(R),score))
      beta <- beta+step
      if(max(abs(step))<1e-9)break
    }
    residual <- max(abs(crossprod(A,case-total*plogis(as.numeric(A%*%beta)))))
    if(iteration==30 || residual>1e-5 || any(!is.finite(beta)))stop("Replicate score solve failed: ",name)
    list(beta=as.numeric(beta),residual=residual,iterations=iteration)
  }
  # Numerical reference check against ordinary binomial fitting where weights are positive.
  check_population <- population+.001*(rep_population[,1]-population)
  stopifnot(all(check_population>=0))
  reference <- suppressWarnings(glm.fit(A,cbind(case,check_population),family=binomial(),
    control=glm.control(epsilon=1e-12,maxit=100)))
  check_solution <- score_fit(case+check_population)$beta
  reference_error <- max(abs(reference$coefficients-check_solution))
  base_error <- max(abs(score_fit(case+population)$beta-b))
  stopifnot(reference$converged,reference_error<1e-6,base_error<1e-6)
  V_acs <- matrix(0,P,P);rep_records <- list();max_residual <- 0;max_iterations <- 0
  start <- proc.time()[3]
  for(yy in years) {
    selected <- year==yy;differences <- matrix(0,80,P)
    for(r in 1:80) {
      total <- case+population
      total[selected] <- case[selected]+rep_population[selected,r]
      result <- score_fit(total);differences[r,] <- result$beta-b
      max_residual <- max(max_residual,result$residual);max_iterations <- max(max_iterations,result$iterations)
      record <- data.frame(year=yy,replicate=r)
      for(j in seq_along(b))record[[names(b)[j]]] <- result$beta[j]
      rep_records[[length(rep_records)+1]] <- record
    }
    # 4/80 SDR multiplier with MSE centering at the full-sample estimate.
    V_acs <- V_acs+(4/80)*crossprod(differences)
  }
  dimnames(V_case) <- dimnames(V_acs) <- list(names(b),names(b))
  V <- (V_case+V_acs);V <- (V+t(V))/2
  old_V <- vcov(model);old_se <- sqrt(diag(old_V));se <- sqrt(diag(V))
  old_p <- 2*pnorm(abs(b/old_se),lower.tail=FALSE);new_p <- 2*pnorm(abs(b/se),lower.tail=FALSE)
  passed <- all(is.finite(V)) && all(diag(V)>0) && min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>-1e-10
  replicate_checks[[length(replicate_checks)+1]] <<- data.frame(model=name,years=length(years),
    replicates_per_year=80,refits=80*length(years),unique_covariate_cells=G,
    maximum_score_residual=max_residual,maximum_iterations=max_iterations,
    positive_weight_reference_error=reference_error,full_sample_score_error=base_error,
    case_weights_fixed=TRUE,negative_replicates_preserved=TRUE,passed=passed,
    elapsed_seconds=unname(proc.time()[3]-start))
  replicate_changes[[length(replicate_changes)+1]] <<- data.frame(model=name,term=names(b),
    coefficient=b,odds_ratio=exp(b),conventional_se=old_se,case_process_se=sqrt(diag(V_case)),
    acs_sampling_se=sqrt(diag(V_acs)),combined_se=se,se_ratio=se/old_se,
    conventional_p=old_p,combined_p=new_p,significance_05_changed=(old_p<.05)!=(new_p<.05))
  stopifnot(passed)
  write_table(do.call(rbind,rep_records),paste0(name,"_replicate_coefficients.csv"),
    directory=result_path(output_dir,"replicate_details"))
  for(component in c("case_process","acs_sampling","combined")) {
    value <- switch(component,case_process=V_case,acs_sampling=V_acs,combined=V)
    write_table(data.frame(term=rownames(value),value,check.names=FALSE),
      paste0(name,"_",component,"_covariance.csv"),directory=result_path(output_dir,"replicate_details"))
  }
  model$two_source_vcov <- V;model$case_process_vcov <- V_case;model$acs_sampling_vcov <- V_acs
  model$covariance_method <- "Poisson case-process sandwich + year-blocked ACS SDR (4/80); normal Wald inference"
  class(model) <- c("stars_replicate_glm",class(model))
  cat(sprintf("Replicate uncertainty: %s; %d year-block refits; %.1fs\n",name,80*length(years),proc.time()[3]-start))
  model
}
vcov.stars_replicate_glm <- function(object,...)object$two_source_vcov
summary.stars_replicate_glm <- function(object,...) {
  result <- summary.glm(object,...)
  result$coefficients[,2] <- sqrt(diag(object$two_source_vcov))
  result$coefficients[,3] <- coef(object)/result$coefficients[,2]
  result$coefficients[,4] <- 2*pnorm(abs(result$coefficients[,3]),lower.tail=FALSE)
  result$cov.unscaled <- result$cov.scaled <- object$two_source_vcov
  if(!is.null(result$correlation))result$correlation <- cov2cor(object$two_source_vcov)
  result$covariance_method <- object$covariance_method
  result
}
