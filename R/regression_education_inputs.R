# Shared resident-only adult education-regression inputs.
# Caller defines output_dir and loads common.R in its own environment.
analysis_years <- c(2012:2019,2021L)
ipums_path <- primary_path("IPUMSData/ipums_2012-2019_2021_with_replicate_weights.csv")
codebook_path <- primary_path("IPUMSData/ipums_2012-2019_2021_with_replicate_weights.cbk.txt")
ipums_raw <- read_table(ipums_path)
replicate_names <- paste0("REPWTP",1:80)
require_columns(ipums_raw,c("YEAR","STATEFIP","PERWT","EDUCD","REPWTP",replicate_names),"IPUMS extract")
stopifnot(setequal(unique(ipums_raw$YEAR),analysis_years),all(ipums_raw$STATEFIP==30),
  all(ipums_raw$REPWTP==1),all(ipums_raw$PERWT>0),
  all(is.finite(as.matrix(ipums_raw[replicate_names]))),
  grepl("REPWTP80",paste(readLines(codebook_path,warn=FALSE),collapse=" "),fixed=TRUE))
ipums <- harmonize_ipums(ipums_raw)
ipums$YEAR <- factor(ipums_raw$YEAR,levels=analysis_years)
ipums$EDUC_ATTAIN <- attainment_ipums(ipums_raw) # validates the observed EDUCD codes
partial_codes <- c(30,40,50,61,65,71)
ipums$EDUC_PARTIAL <- factor(ifelse(ipums_raw$EDUCD %in% c(0,1),NA_character_,
  ifelse(ipums_raw$EDUCD %in% partial_codes,"PARTIAL","OTHER")),levels=c("OTHER","PARTIAL"))
ipums$SOURCE_ROW <- seq_len(nrow(ipums))
raw <- read_stars_csv()
all_stars <- harmonize_stars(raw)
all_stars$EDUC_ATTAIN <- attainment_stars(raw)
all_stars$EDUC_PARTIAL <- factor(ifelse(is.na(all_stars$EDUC),NA_character_,
  ifelse(all_stars$EDUC=="PARTIAL","PARTIAL","OTHER")),levels=c("OTHER","PARTIAL"))
stars_residence <- raw$DEM_RES09State
write_table(regression_residence_audit(raw,analysis_years),"regression_residence_eligibility.csv")
eligible <- which(all_stars$DYEAR %in% analysis_years & montana_resident(raw))
stars <- all_stars[eligible,]
stars$YEAR <- factor(stars$DYEAR,levels=analysis_years)
stars$weight <- 1
stars$SOURCE_ROW <- eligible
shared <- intersect(names(ipums),names(stars))
ipums_stars <- rbind(ipums[shared],stars[shared])
ipums_stars$AGE_ADULT <- cut(ipums_stars$AGE,breaks=c(17,24,34,44,54,64,74,84,Inf),
  labels=c("18-24","25-34","35-44","45-54","55-64","65-74","75-84","85+"))
ipums_stars <- ipums_stars[!is.na(ipums_stars$AGE) & ipums_stars$AGE>=18,]
stopifnot(identical(is.na(ipums_stars$EDUC_PARTIAL),is.na(ipums_stars$EDUC_ATTAIN)))
