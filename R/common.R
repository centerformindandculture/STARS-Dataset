source("R/output_paths.R", local=TRUE)
# Shared helpers for the STARS supplement.
# Scripts are run from the supplement root; no personal absolute paths are used.
# The cluster package supplies silhouette and gap-statistic diagnostics.

options(stringsAsFactors = FALSE, digits = 12)
Sys.setenv(TZ = "UTC")
invisible(Sys.setlocale("LC_COLLATE", "C"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

read_table <- function(path) {
  read.csv(result_path(dirname(path), basename(path)), check.names = FALSE, na.strings = c("", "NA"))
}

write_table <- function(data, name, directory = output_dir) {
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  write.csv(data, result_path(directory, name), row.names = FALSE, na = "")
}

primary_path <- function(relative) {
  # Run from the supplement/repository root. Primary files live in the parent
  # STARS directory: ../R - Data Exploration/, ../IPUMSData/, and
  # A fresh clone requires these data folders
  # supplied alongside it; filenames and required columns are in docs/INPUTS.md.
  # STARS_SOURCE_DIR optionally overrides that parent directory.
  base <- Sys.getenv("STARS_SOURCE_DIR", unset = "..")
  if (!nzchar(base)) base <- ".."
  path <- file.path(base, relative)
  if (!file.exists(path)) {
    stop("Missing primary input: ", path,
         "\nPlace the expected data folders beside the cloned repository, or set ",
         "STARS_SOURCE_DIR to their parent. See docs/INPUTS.md.", call. = FALSE)
  }
  path
}

read_stars_csv <- function() {
  result <- read_table(primary_path("R - Data Exploration/20230831_Montana_Suicide_Database.csv"))
  for (column in names(result)) if (is.character(result[[column]]))
    result[[column]] <- clean_text(result[[column]])
  result
}

montana_resident <- function(raw) {
  require_columns(raw,"DEM_RES09State","STARS residence")
  state <- clean_text(raw$DEM_RES09State)
  !is.na(state) & state == "MONTANA"
}

regression_residence_audit <- function(raw, years) {
  state <- clean_text(raw$DEM_RES09State)
  status <- factor(ifelse(is.na(state),"UNKNOWN",
    ifelse(state=="MONTANA","MONTANA","NON_MONTANA")),
    levels=c("MONTANA","NON_MONTANA","UNKNOWN"))
  result <- as.data.frame(table(year=raw$DEATH_2Year,residence=status))
  names(result)[3] <- "records"
  result$eligible_year <- as.integer(as.character(result$year)) %in% years
  result$residence_eligible <- result$residence=="MONTANA"
  result
}

attainment_stars <- function(raw) {
  allowed <- c("NO HIGH SCHOOL","SOME HIGH SCHOOL","HIGH SCHOOL OR GED","SOME COLLEGE",
    "ASSOCIATES","BACHELORS","MASTERS","DOCTORATE")
  if (any(!is.na(raw$DEM_EDU1EducationAchieved) & !raw$DEM_EDU1EducationAchieved %in% allowed))
    stop("Unexpected STARS education category; review before recoding.")
  factor(recode_values(raw$DEM_EDU1EducationAchieved, list(
    LESS_THAN_ASSOCIATE = c("NO HIGH SCHOOL","SOME HIGH SCHOOL","HIGH SCHOOL OR GED","SOME COLLEGE"),
    ASSOCIATE_OR_HIGHER = c("ASSOCIATES","BACHELORS","MASTERS","DOCTORATE"))),
    levels=c("LESS_THAN_ASSOCIATE","ASSOCIATE_OR_HIGHER"))
}

attainment_ipums <- function(raw) {
  below <- c(2,11,12,14,15,16,17,22,23,25,26,30,40,50,61,63,64,65,71)
  degree <- c(81,101,114,115,116)
  if (!all(raw$EDUCD %in% c(0,1,below,degree)))
    stop("Unexpected EDUCD category; review the supplied codebook.")
  factor(recode_values(raw$EDUCD,list(LESS_THAN_ASSOCIATE=below,ASSOCIATE_OR_HIGHER=degree)),
    levels=c("LESS_THAN_ASSOCIATE","ASSOCIATE_OR_HIGHER"))
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

county_name <- function(x) {
  x <- clean_text(x)
  x[x == "Lewis & Clark" & !is.na(x)] <- "Lewis and Clark"
  x
}

age_group <- function(age) {
  cut(age, breaks = c(-Inf, 17, 24, 34, 44, 54, 64, Inf),
      labels = c("<18", "18-24", "25-34", "35-44", "45-54", "55-64", "65+"))
}

recode_values <- function(x, groups) {
  result <- rep(NA_character_, length(x))
  for (label in names(groups)) result[x %in% groups[[label]]] <- label
  result
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
}

harmonize_stars <- function(raw) {
  # Derived from RegressionOverview.Rmd.
  # Recode numeric sex/veteran values ONCE, before combining datasets.
  require_columns(raw, c("DEM_MAR1MaritalStatus", "DEM_RACE1", "DEM_RACE6Hispanic",
    "DEM_EDU1EducationAchieved", "DEM_GEN2GenderNum", "DEM_AGE1Age", "DEM_OCC5ArmedForcesEver",
    "DEM_OCC3RancherFarmer", "DEATH_2Year", "COD_Means"), "STARS data")
  marital <- clean_text(raw$DEM_MAR1MaritalStatus)
  race <- recode_values(raw$DEM_RACE1, list(
    WHITE = "WHITE", `AMERICAN INDIAN` = "AMERICAN INDIAN",
    OTHER = c("OTHER", "HISPANIC", "ASIAN/PACIFIC ISLANDER",
              "ASIAN/ PACIFIC ISLANDER", "BLACK")))
  race[which(raw$DEM_RACE6Hispanic == 1)] <- "OTHER"
  education <- recode_values(raw$DEM_EDU1EducationAchieved, list(
    COMPLETE = c("NO HIGH SCHOOL", "HIGH SCHOOL OR GED", "ASSOCIATES",
                 "BACHELORS", "MASTERS", "DOCTORATE"),
    PARTIAL = c("SOME HIGH SCHOOL", "SOME COLLEGE")))
  divsep <- recode_values(marital, list(
    OTHER = c("NEVER MARRIED", "SINGLE", "WIDOWED",
              "MARRIED (PREVIOUSLY DIVORCED, NOT SEPARATED)",
              "MARRIED (NOT PREVIOUSLY DIVORCED, NOT SEPARATED)"),
    `DIVORCED/SEPARATED` = c("DIVORCED", "MARRIED (SEPARATED)")))
  marital2 <- ifelse(divsep == "OTHER", "OTHER",
                     ifelse(marital == "DIVORCED", "DIVORCED", "SEPARATED"))
  data.frame(
    SEX = factor(raw$DEM_GEN2GenderNum, c(0, 1), c("MALE", "FEMALE")),
    AGE = raw$DEM_AGE1Age, AGE7 = age_group(raw$DEM_AGE1Age),
    RACE = factor(race, c("WHITE", "AMERICAN INDIAN", "OTHER")),
    EDUC = factor(education, c("COMPLETE", "PARTIAL")),
    VETSTAT = factor(raw$DEM_OCC5ArmedForcesEver, c(0, 1), c("NOT AF", "AF")),
    FARM = factor(raw$DEM_OCC3RancherFarmer, c(0, 1), c("OTHER", "FARMER")),
    DIVSEP = factor(divsep, c("OTHER", "DIVORCED/SEPARATED")),
    MARST2 = factor(marital2, c("OTHER", "DIVORCED", "SEPARATED")),
    DYEAR = raw$DEATH_2Year, COD = factor(clean_text(raw$COD_Means)),
    SUIC = 1L, weight = 1
  )
}

harmonize_ipums <- function(raw) {
  # The logistic module adds YEAR and sets annual case weights explicitly.
  require_columns(raw, c("SEX", "AGE", "MARST", "RACE", "HISPAN", "EDUCD", "VETSTAT", "PERWT", "OCC2010"), "IPUMS extract")
  race <- recode_values(raw$RACE, list(WHITE = 1, `AMERICAN INDIAN` = 3,
                                      OTHER = c(2, 4:9)))
  race[raw$HISPAN %in% 1:4] <- "OTHER"
  marital <- recode_values(raw$MARST, list(OTHER = c(1, 2, 5, 6),
                                          `DIVORCED/SEPARATED` = c(3, 4)))
  data.frame(
    SEX = factor(raw$SEX, c(1, 2), c("MALE", "FEMALE")),
    AGE = raw$AGE, AGE7 = age_group(raw$AGE),
    RACE = factor(race, c("WHITE", "AMERICAN INDIAN", "OTHER")),
    VETSTAT = factor(raw$VETSTAT, c(1, 2), c("NOT AF", "AF")),
    FARM = factor(ifelse(raw$OCC2010 == 205, "FARMER", "OTHER"), c("OTHER", "FARMER")),
    DIVSEP = factor(marital, c("OTHER", "DIVORCED/SEPARATED")),
    MARST2 = factor(recode_values(raw$MARST,
      list(OTHER = c(1, 2, 5, 6), DIVORCED = 4, SEPARATED = 3)),
      c("OTHER", "DIVORCED", "SEPARATED")),
    SUIC = 0L, weight = raw$PERWT
  )
}

coefficient_table <- function(model, name) {
  tab <- summary(model)$coefficients
  result <- data.frame(term = rownames(tab), coefficient = tab[, 1],
    std_error = tab[, 2], z = tab[, 3], p = tab[, 4], odds_ratio = exp(tab[, 1]))
  write_table(result, paste0(name, ".csv"))
  result
}
