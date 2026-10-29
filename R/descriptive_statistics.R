# Source-cohort descriptives from the designated August 31 CSV.
# Full 2012-2021 source cohort; regressions use their separately documented subsets.
cd_parent <- output_dir;output_dir <- result_path(output_dir,"descriptive")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
cd_raw <- read_stars_csv();cd_data <- harmonize_stars(cd_raw)
cd_data$EDUC_ATTAIN <- attainment_stars(cd_raw)
cd_data$AGE_GROUP <- cut(cd_data$AGE,breaks=c(-Inf,17,24,34,44,54,64,74,84,Inf),
  labels=c("<18","18-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"))
cd_data$COUNTY <- clean_text(cd_raw$DEM_RES03CountyCorrected)
cd_data$MONTH <- factor(cd_raw$DEATH_4MonthNum,levels=1:12,labels=month.name)
cd_data$WEEKDAY <- factor(clean_text(cd_raw$DEATH_5Day),levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))
cd_data$REPORTED_EDUCATION <- clean_text(cd_raw$DEM_EDU1EducationAchieved)
cd_data$REPORTED_MARITAL_STATUS <- clean_text(cd_raw$DEM_MAR1MaritalStatus)
cd_data$REPORTED_RACE <- clean_text(cd_raw$DEM_RACE1)
cd_county <- read_table(file.path(cd_parent,"clustering/geographic/county_k3_assignments.csv"))
cd_data$COUNTY_CLUSTER <- factor(cd_county$CLUSTER[match(cd_data$COUNTY,cd_county$COUNTY)],levels=1:3)
cd_mon <- read_table(file.path(cd_parent,"clustering/monthly/month_k3_assignments.csv"))
cd_day <- read_table(file.path(cd_parent,"clustering/weekday/weekday_k3_assignments.csv"))
cd_data$MONTH_CLUSTER <- factor(cd_mon$CLUSTER[match(cd_raw$DEATH_4MonthNum,cd_mon$DMONTH)],levels=1:3)
cd_data$WEEKDAY_CLUSTER <- factor(cd_day$CLUSTER[match(as.character(cd_data$WEEKDAY),cd_day$DDAY)],levels=1:3)
cd_vars <- c("SEX","AGE_GROUP","RACE","EDUC_ATTAIN","VETSTAT","FARM","DIVSEP",
 "REPORTED_EDUCATION","REPORTED_MARITAL_STATUS","REPORTED_RACE","COUNTY","MONTH","WEEKDAY",
 "COUNTY_CLUSTER","MONTH_CLUSTER","WEEKDAY_CLUSTER","COD","DYEAR")
cd_tables <- list();cd_missing <- list();cd_yearmiss <- list()
for(v in cd_vars) {
 x <- as.character(cd_data[[v]]);x[is.na(x) | x=="" | x=="NA"] <- NA_character_
 known <- sum(!is.na(x));missing <- sum(is.na(x));N <- length(x)
 cd_missing[[v]] <- data.frame(variable=v,n=N,known=known,missing=missing,percent_missing=100*missing/N)
 for(y in 2012:2021) {
  xx <- x[cd_data$DYEAR==y]
  cd_yearmiss[[paste(v,y)]] <- data.frame(variable=v,year=y,n=length(xx),known=sum(!is.na(xx)),
   missing=sum(is.na(xx)),percent_missing=100*mean(is.na(xx)))
 }
 values <- x;values[is.na(values)] <- "Missing"
 freq <- if(is.factor(cd_data[[v]]))table(factor(values,levels=c(levels(cd_data[[v]]),"Missing"))) else table(values)
 tab <- data.frame(variable=v,category=names(freq),n=as.integer(freq),
   percent_total=100*as.integer(freq)/N,percent_known=ifelse(names(freq)=="Missing" | known==0,NA_real_,100*as.integer(freq)/known),
   total_records=N,known_records=known,missing_records=missing)
 cd_tables[[v]] <- tab
 write_table(tab,paste0("frequency_",v,".csv"))
}
write_table(do.call(rbind,cd_tables),"demographic_frequencies.csv")
write_table(do.call(rbind,cd_missing),"variable_missingness.csv")
write_table(do.call(rbind,cd_yearmiss),"missingness_by_year.csv")
write_table(data.frame(n=nrow(cd_data),age_known=sum(!is.na(cd_data$AGE)),age_missing=sum(is.na(cd_data$AGE)),
 age_mean=mean(cd_data$AGE,na.rm=TRUE),age_sd=sd(cd_data$AGE,na.rm=TRUE),age_median=median(cd_data$AGE,na.rm=TRUE),
 age_min=min(cd_data$AGE,na.rm=TRUE),age_max=max(cd_data$AGE,na.rm=TRUE),
 under18=sum(cd_data$AGE<18,na.rm=TRUE),adults=sum(cd_data$AGE>=18,na.rm=TRUE)),"sample_description.csv")
cd_method <- as.data.frame(table(year=cd_data$DYEAR,method=cd_data$COD))
write_table(cd_method,"means_by_year_counts.csv")
cd_sex_method <- as.data.frame(table(sex=cd_data$SEX,method=cd_data$COD))
cd_sex_method$percent_within_known_sex_method <- 100*cd_sex_method$Freq/ave(cd_sex_method$Freq,cd_sex_method$sex,FUN=sum)
write_table(cd_sex_method,"means_by_sex.csv")
write_table(data.frame(input="20230831_Montana_Suicide_Database.csv",
 md5=unname(tools::md5sum(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv")))),"input_checksum.csv")
cd_checks <- data.frame(check=c("full cohort 2764","years 2012-2021","month missing 1530","weekday missing 489",
 "county cluster matched 2650","all frequency totals match"),passed=c(nrow(cd_data)==2764,
 setequal(cd_data$DYEAR,2012:2021),sum(is.na(cd_data$MONTH))==1530,sum(is.na(cd_data$WEEKDAY))==489,
 sum(!is.na(cd_data$COUNTY_CLUSTER))==2650,all(vapply(cd_tables,function(t)sum(t$n)==2764,logical(1)))))
write_table(cd_checks,"validation_checks.csv");stopifnot(all(cd_checks$passed))
writeLines(c("Source-cohort descriptives: August 31 CSV, 2012-2021, all 2764 source records.",
 "Percent total includes missing values in the denominator. Percent known excludes missing values.",
 "Missingness includes unavailable harmonized codes; raw education, marital and race labels are also retained as separate distributions.",
 "No demographic association tests are added to these descriptive tables.",
 "These summaries describe the source cohort; county mortality and adult regressions have separately reported eligibility restrictions.",
 paste("Age mean:",round(mean(cd_data$AGE,na.rm=TRUE),2)),
 "Month missing: 1530; weekday missing: 489; county not matched: 114."),result_path(output_dir,"descriptive_report.txt"))
output_dir <- cd_parent
