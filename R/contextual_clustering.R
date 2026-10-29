# Contextual clustering and ecological mortality comparisons.
parent_of_output <- output_dir
output_dir <- result_path(output_dir,"clustering")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
inputs <- "data/processed/clustering"
check_results <- list()
of_check <- function(name,value) {
  check_results[[length(check_results)+1L]] <<- data.frame(check=name,passed=isTRUE(value))
}
county_context <- read_table(file.path(inputs,"county_context.csv"))
month_context <- read_table(file.path(inputs,"month_context.csv"))
weekday_context <- read_table(file.path(inputs,"weekday_context.csv"))
county_deaths <- read_table(file.path(inputs,"county_year_mortality.csv"))
month_deaths <- read_table(file.path(inputs,"month_year_mortality.csv"))
weekday_deaths <- read_table(file.path(inputs,"weekday_year_mortality.csv"))
of_check("all 56 counties unique",nrow(county_context)==56 && !anyDuplicated(county_context$COUNTY))
of_check("all county years represented",nrow(county_deaths)==560 && !anyDuplicated(county_deaths[c("COUNTY","YEAR")]))
of_check("2650 county-matched deaths",sum(county_deaths$DEATHS)==2650)
of_check("Petroleum included with zero deaths",sum(county_deaths$DEATHS[county_deaths$COUNTY=="Petroleum"])==0 &&
  sum(county_deaths$COUNTY=="Petroleum")==10)
of_check("positive population denominators",all(county_deaths$POPULATION>0))
of_check("complete temporal years only",setequal(month_deaths$YEAR,c(2018,2020,2021)) &&
  setequal(weekday_deaths$YEAR,c(2013:2016,2018,2020,2021)))
of_check("temporal outcome totals agree in common years",sum(month_deaths$DEATHS)==
  sum(weekday_deaths$DEATHS[weekday_deaths$YEAR %in% month_deaths$YEAR]))
of_check("temporal person years agree",abs(sum(month_deaths$PERSON_YEARS)-sum(weekday_deaths$PERSON_YEARS[weekday_deaths$YEAR %in% month_deaths$YEAR]))<1e-6)
ari <- function(a,b) {
  t <- table(a,b);choose2 <- function(x)x*(x-1)/2
  s <- sum(choose2(t));r <- sum(choose2(rowSums(t)));c <- sum(choose2(colSums(t)))
  e <- r*c/choose2(sum(t));den <- (r+c)/2-e
  if(abs(den)<1e-10)return(1)
  (s-e)/den
}
selections <- list();all_assignments <- list();stabilities <- list()
cluster_context <- function(context,key,features,domain,max_k,order_feature,seed) {
  stopifnot(all(vapply(context[features],is.numeric,logical(1))),!anyNA(context[features]))
  # Strict whitelists: no deaths, rates, or attributes of STARS decedents.
  of_check(paste(domain,"no outcome features"),
    !any(features %in% c("SUICOUNT","DEATHSPERCAP","deaths","deathsperday","medage","gunshot","vetstat")))
  x <- scale(context[features])
  of_check(paste(domain,"finite scaled features"),all(is.finite(x)))
  write_table(data.frame(variable=features,center=attr(x,"scaled:center"),
    scale=attr(x,"scaled:scale")),paste0(domain,"_scaling.csv"))
  write_table(cbind(context[key],as.data.frame(x)),paste0(domain,"_scaled_features.csv"))
  fits <- list();diagnostics <- list()
  for(k in 1:max_k) {
    set.seed(seed+k)
    fit <- kmeans(x,k,nstart=100,iter.max=200)
    fits[[as.character(k)]] <- fit
    sil <- if(k==1)NA_real_ else mean(cluster::silhouette(fit$cluster,dist(x))[,3])
    diagnostics[[k]] <- data.frame(k=k,wss=fit$tot.withinss,silhouette=sil,
      minimum_cluster_size=min(fit$size),maximum_cluster_size=max(fit$size))
  }
  set.seed(seed+100)
  gap <- cluster::clusGap(x,FUNcluster=function(x,k)list(cluster=kmeans(x,k,
    nstart=30,iter.max=200)$cluster),K.max=max_k,B=200,d.power=2,spaceH0="scaledPCA")
  diag <- do.call(rbind,diagnostics)
  diag$gap <- gap$Tab[,"gap"];diag$gap_se <- gap$Tab[,"SE.sim"]
  # Elbow: maximum perpendicular distance from line joining normalized endpoints.
  xx <- (diag$k-1)/(max_k-1)
  yy <- (diag$wss-min(diag$wss))/(max(diag$wss)-min(diag$wss))
  elbow <- diag$k[which.max((1-xx-yy)/sqrt(2))]
  gap_k <- cluster::maxSE(diag$gap,diag$gap_se,method="Tibs2001SEmax")
  silhouette_k <- diag$k[which.max(replace(diag$silhouette,is.na(diag$silhouette),-Inf))]
  candidates <- diag[diag$k %in% 3:4,]
  # smaller solution if its silhouette is within 0.02 of the best in this range.
  chosen <- min(candidates$k[candidates$silhouette>=max(candidates$silhouette)-.02])
  selected <- data.frame(domain=domain,n_units=nrow(x),n_features=ncol(x),
    selected_k=chosen,elbow_k=elbow,silhouette_k=silhouette_k,gap_1se_k=gap_k,
    selection_rule="Smallest k in 3:4 within 0.02 of best silhouette in 3:4")
  selections[[domain]] <<- selected
  write_table(diag,paste0(domain,"_k_diagnostics.csv"))
  profiles <- list();assignments <- list()
  for(k in 3:4) {
    original <- fits[[as.character(k)]]$cluster
    ordering <- order(tapply(context[[order_feature]],original,mean))
    relabel <- setNames(seq_len(k),as.character(ordering))
    assignment <- unname(relabel[as.character(original)])
    assignments[[as.character(k)]] <- assignment
    all_assignments[[paste(domain,k)]] <<- data.frame(domain=domain,k=k,
      unit=as.character(context[[key]]),cluster=assignment,preferred=k==chosen)
    write_table(data.frame(context[key],CLUSTER=assignment),
      paste0(domain,"_k",k,"_assignments.csv"))
    for(g in seq_len(k)) for(v in features) {
      values <- context[[v]][assignment==g]
      profiles[[length(profiles)+1L]] <- data.frame(k=k,cluster=g,variable=v,
        n=length(values),mean=mean(values),sd=if(length(values)>1)sd(values) else NA_real_)
    }
    # Perturbation diagnostic: remove one unit, refit with fixed full-data scaling.
    for(i in seq_len(nrow(x))) {
      set.seed(seed+k*1000+i)
      refit <- kmeans(x[-i,,drop=FALSE],k,nstart=100,iter.max=200)
      stabilities[[length(stabilities)+1L]] <<- data.frame(domain=domain,k=k,
        omitted_unit=as.character(context[[key]][i]),adjusted_rand=ari(original[-i],refit$cluster))
    }
  }
  write_table(do.call(rbind,profiles),paste0(domain,"_cluster_profiles.csv"))
  png(result_path(output_dir,paste0(domain,"_k_diagnostics.png")),width=1800,height=800,res=140)
  par(mfrow=c(1,3),mar=c(5,5,3,1))
  plot(diag$k,diag$wss,type="b",xlab="Number of clusters",ylab="Within-cluster sum of squares",
    main=paste(domain,"elbow"));abline(v=chosen,col="#a33",lty=2)
  plot(diag$k,diag$silhouette,type="b",xlab="Number of clusters",ylab="Mean silhouette",
    main="Silhouette");abline(v=chosen,col="#a33",lty=2)
  plot(diag$k,diag$gap,type="b",ylim=range(c(diag$gap-diag$gap_se,diag$gap+diag$gap_se)),
    xlab="Number of clusters",ylab="Gap statistic",main="Gap with simulation SE")
  arrows(diag$k,diag$gap-diag$gap_se,diag$k,diag$gap+diag$gap_se,angle=90,code=3,length=.03)
  abline(v=chosen,col="#a33",lty=2);invisible(dev.off())
  list(assignments=assignments,selected=chosen)
}
geo_features <- c("MEDAGE","AIAN_ALONE_PCT","MEDINCOME","POVERTY_PCT",
  "SINGLE_HH_PCT","DIVSEP_PCT","DEGREE_PCT","UNINSURED_PCT")
month_features <- c("freezedays","hotdays","prcp","tavg","tmax","tmin","humidity","crimeindex")
weekday_features <- c("crashperc","substancedeaths","happiness")
geo <- cluster_context(county_context,"COUNTY",geo_features,"county",8,"MEDINCOME",1700)
mon <- cluster_context(month_context,"DMONTH",month_features,"month",6,"tavg",2700)
day <- cluster_context(weekday_context,"DDAY",weekday_features,"weekday",5,"happiness",3700)
write_table(do.call(rbind,selections),"cluster_selection.csv")
write_table(do.call(rbind,all_assignments),"all_assignments.csv")
stability <- do.call(rbind,stabilities)
write_table(stability,"leave_one_out_stability.csv")
stab_summary <- aggregate(adjusted_rand~domain+k,stability,function(x)
  c(mean=mean(x),minimum=min(x),median=median(x)))
write_table(data.frame(domain=stab_summary$domain,k=stab_summary$k,stab_summary$adjusted_rand),
  "stability_summary.csv")
omnibus <- list();coef_outputs <- list();pairwise_outputs <- list();model_diagnostics <- list()
summaries <- list()
mortality_fit <- function(data,domain,k,metric,unit) {
  d <- data;d$CLUSTER <- factor(d$CLUSTER);d$YEAR <- factor(d$YEAR)
  d$EXPOSURE <- if(metric=="count")1 else if(domain=="county")d$POPULATION else d$PERSON_YEARS
  fit <- glm(DEATHS~CLUSTER+YEAR+offset(log(EXPOSURE)),data=d,family=poisson())
  X <- model.matrix(fit);mu <- fitted(fit)
  dispersion <- sum(residuals(fit,type="pearson")^2)/df.residual(fit)
  if(domain=="county") {
    score <- rowsum(X*as.numeric(d$DEATHS-mu),d[[unit]],reorder=FALSE)
    G <- nrow(score);n <- nrow(X);p <- ncol(X)
    bread <- solve(crossprod(X,X*as.numeric(mu)))
    V <- bread%*%crossprod(score)%*%bread * G/(G-1)*(n-1)/(n-p)
    df <- G-1;se_method <- "County-clustered sandwich with finite-sample correction; t/F with 55 df"
  } else {
    V <- vcov(fit)*max(1,dispersion)
    df <- df.residual(fit);se_method <- "Poisson covariance inflated by max(1, Pearson dispersion)"
  }
  b <- coef(fit);ses <- sqrt(diag(V));q <- length(levels(d$CLUSTER))-1L
  idx <- grep("^CLUSTER",names(b))
  Fstat <- as.numeric(t(b[idx])%*%solve(V[idx,idx,drop=FALSE],b[idx])/q)
  pv <- pf(Fstat,q,df,lower.tail=FALSE)
  omnibus[[length(omnibus)+1L]] <<- data.frame(domain=domain,k=k,metric=metric,
    preferred=k==selections[[domain]]$selected_k,F=Fstat,df1=q,df2=df,p=pv,se_method=se_method)
  estimates <- data.frame(domain=domain,k=k,metric=metric,term=names(b),
    coefficient=b,std_error=ses,ratio=exp(b),ci_lower=exp(b-qt(.975,df)*ses),
    ci_upper=exp(b+qt(.975,df)*ses),p=2*pt(abs(b/ses),df,lower.tail=FALSE))
  coef_outputs[[length(coef_outputs)+1L]] <<- estimates
  pair <- list()
  for(g1 in 1:(k-1))for(g2 in (g1+1):k) {
    L <- rep(0,length(b))
    if(g2!=1)L[match(paste0("CLUSTER",g2),names(b))] <- 1
    if(g1!=1)L[match(paste0("CLUSTER",g1),names(b))] <- -1
    delta <- sum(L*b);s <- sqrt(as.numeric(t(L)%*%V%*%L))
    pair[[length(pair)+1L]] <- data.frame(domain=domain,k=k,metric=metric,
      comparison=paste0("cluster",g2,"/cluster",g1),ratio=exp(delta),
      ci_lower=exp(delta-qt(.975,df)*s),ci_upper=exp(delta+qt(.975,df)*s),
      p=2*pt(abs(delta/s),df,lower.tail=FALSE))
  }
  pair <- do.call(rbind,pair);pair$p_holm <- p.adjust(pair$p,"holm")
  pairwise_outputs[[length(pairwise_outputs)+1L]] <<- pair
  model_diagnostics[[length(model_diagnostics)+1L]] <<- data.frame(domain=domain,k=k,
    metric=metric,n=nrow(d),deaths=sum(d$DEATHS),dispersion=dispersion,
    converged=fit$converged,rank=fit$rank,parameters=length(b),se_method=se_method)
  of_check(paste(domain,k,metric,"fit valid"),fit$converged &&
    fit$rank==length(b) && all(is.finite(b)) && all(is.finite(V)) && all(diag(V)>0))
  invisible(fit)
}
pool_rates <- function(d,domain,k,exposure) {
  r <- aggregate(d[c("DEATHS",exposure)],list(cluster=d$CLUSTER),sum)
  names(r)[3] <- "person_years"
  r$domain <- domain;r$k <- k;r$preferred <- k==selections[[domain]]$selected_k
  r$rate_per_100000 <- r$DEATHS/r$person_years*1e5
  r$rate_ci_lower <- ifelse(r$DEATHS==0,0,qchisq(.025,2*r$DEATHS)/2/r$person_years*1e5)
  r$rate_ci_upper <- qchisq(.975,2*(r$DEATHS+1))/2/r$person_years*1e5
  r
}
for(k in 3:4) {
  cd <- county_deaths
  cd$CLUSTER <- geo$assignments[[as.character(k)]][match(cd$COUNTY,county_context$COUNTY)]
  of_check(paste("county",k,"every county assigned"),!anyNA(cd$CLUSTER))
  mortality_fit(cd,"county",k,"count","COUNTY")
  mortality_fit(cd,"county",k,"rate","COUNTY")
  totals <- pool_rates(transform(cd,PERSON_YEARS=POPULATION),"county",k,"PERSON_YEARS")
  totals$counties <- as.integer(table(geo$assignments[[as.character(k)]]))
  totals$mean_annual_deaths_per_county <- totals$DEATHS/(totals$counties*10)
  totals$share_of_deaths <- totals$DEATHS/sum(totals$DEATHS)
  write_table(totals,paste0("county_k",k,"_mortality.csv"))
  write_table(cd,paste0("county_k",k,"_county_year_analysis.csv"))
  summaries[[paste0("county",k)]] <- totals
  for(domain in c("month","weekday")) {
    dd <- if(domain=="month")month_deaths else weekday_deaths
    cc <- if(domain=="month")mon else day
    keys <- if(domain=="month")month_context$DMONTH else weekday_context$DDAY
    dd$CLUSTER <- cc$assignments[[as.character(k)]][match(dd$UNIT,keys)]
    of_check(paste(domain,k,"all units assigned"),!anyNA(dd$CLUSTER))
    mortality_fit(dd,domain,k,"rate","UNIT")
    totals <- pool_rates(dd,domain,k,"PERSON_YEARS")
    totals$calendar_days <- aggregate(dd$DAYS,list(dd$CLUSTER),sum)$x
    totals$deaths_per_calendar_day <- totals$DEATHS/totals$calendar_days
    write_table(totals,paste0(domain,"_k",k,"_mortality.csv"))
    write_table(dd,paste0(domain,"_k",k,"_year_analysis.csv"))
    summaries[[paste0(domain,k)]] <- totals
  }
}
omni <- do.call(rbind,omnibus)
omni$p_holm_preferred <- NA_real_
omni$p_holm_preferred[omni$preferred] <- p.adjust(omni$p[omni$preferred],"holm")
write_table(omni,"mortality_omnibus_tests.csv")
write_table(do.call(rbind,coef_outputs),"mortality_model_coefficients.csv")
write_table(do.call(rbind,pairwise_outputs),"mortality_pairwise_comparisons.csv")
write_table(do.call(rbind,model_diagnostics),"mortality_model_diagnostics.csv")
# Descriptive crude rates with fixed-population Poisson intervals, not model uncertainty.
for(domain in c("county","month","weekday")) {
  k <- selections[[domain]]$selected_k
  totals <- summaries[[paste0(domain,k)]]
  png(result_path(output_dir,paste0(domain,"_preferred_mortality.png")),
    width=1000,height=750,res=140)
  par(mar=c(5,5,4,1))
  mids <- barplot(totals$rate_per_100000,names.arg=paste("Cluster",totals$cluster),
    ylim=c(0,max(totals$rate_ci_upper)*1.15),col="#497c98",
    ylab="Deaths per 100,000 person-years",main=paste(domain,"clusters: crude mortality"))
  arrows(mids,totals$rate_ci_lower,mids,totals$rate_ci_upper,angle=90,code=3,length=.06)
  invisible(dev.off())
}
write_table(do.call(rbind,check_results),"validation_checks.csv")
writeLines(capture.output(sessionInfo()),result_path(output_dir,"sessionInfo.txt"))
report <- c("Outcome-free clusters and ecological mortality comparisons",
  paste("Run time (UTC):",format(Sys.time(),tz="UTC",usetz=TRUE)),
  "Geographic coverage: all 56 Montana counties, including zero-death Petroleum.",
  "Outcome totals: 2650 matched county-resident deaths; 114 unmatched/out-of-state records excluded.",
  "County features: seven 2007-2011 ACS contextual variables plus 2008-2012 ACS uninsured percentage.",
  "No STARS-derived demographic summaries, death counts or mortality rates enter clustering.",
  "Monthly features: eight external weather/crime measures. Weekday features: external crash,",
  "substance-death and happiness measures; no suicide-derived age/veteran/outcome features.",
  "All features z-standardized; k-means uses 100 starts. Gap statistic: B=200, d.power=2.",
  "Elbow uses maximum distance from normalized WSS endpoint line. Silhouette and gap are also reported.",
  "Selection rule: smallest k in the prespecified 3-4 range within 0.02 of the best silhouette in that range.",
  "Unconstrained elbow, silhouette and gap optima are reported separately.",
  "Both k=3 and k=4 assignments, profiles and mortality comparisons are retained.",
  "No mortality result is used for k selection or labeling. Labels ordered by context only.",
  "Stability: leave-one-unit-out refits, using fixed full-data feature scaling; see stability_summary.csv.",
  "County mortality denominators: annual July population, summed as person-years.",
  "2012-2019: revised Census intercensal estimates; 2020-2021: vintage 2021 estimates.",
  "Pooled mortality rates are crude; demographic-standardized rates are in ../geographic/county_k3_adjusted_rates.csv.",
  "County rate model: Poisson cluster + year + log(population) offset, county-clustered covariance with finite-sample correction.",
  "County count model: same model without population offset; ratios compare annual deaths per county.",
  "Monthly mortality uses 2018, 2020, 2021 (complete month recording).",
  "Weekday mortality uses 2013-2016, 2018, 2020, 2021 (complete weekday recording).",
  "Temporal analyses use county-matched Montana residents and statewide person-time, adjusted",
  "for number of days represented and calendar year. These are not full-decade seasonal estimates.",
  "Temporal covariance is Poisson covariance multiplied by max(1, Pearson dispersion).",
  "Primary omnibus tests are Holm-adjusted across the four preferred comparisons (geographic counts,",
  "geographic rates, month rates, weekday rates). Pairwise p-values adjusted within each fit.",
  "Pooled rate intervals are descriptive Poisson intervals; Census/ACS uncertainty is not propagated.",
  "Baseline contextual profiles may change during the decade; small-county ACS estimates can be noisy.",
  "Month clustering has 12 units; weekday clustering has 7. Results are exploratory.",
  "Context profiles are descriptive; no post-clustering significance tests of defining features are used.",
  "")
for(domain in names(selections)) {
  s <- selections[[domain]]
  report <- c(report,sprintf("%s: preferred k=%d; elbow=%d; silhouette optimum=%d; gap 1SE=%d",
    domain,s$selected_k,s$elbow_k,s$silhouette_k,s$gap_1se_k))
}
report <- c(report,"","Preferred mortality omnibus tests:")
for(i in which(omni$preferred))report <- c(report,sprintf("%s %s k=%d: p=%.6g, Holm p=%.6g",
  omni$domain[i],omni$metric[i],omni$k[i],omni$p[i],omni$p_holm_preferred[i]))
report <- c(report,paste("Validation checks passed:",sum(vapply(check_results,function(x)x$passed,logical(1))),
  "of",length(check_results)))
writeLines(report,result_path(output_dir,"analysis_report.txt"))
cat(paste(report,collapse="\n"),"\n")
if(!all(vapply(check_results,function(x)x$passed,logical(1))))stop("Outcome-free cluster validation failed")
source("R/county_adjusted.R")
output_dir <- parent_of_output
