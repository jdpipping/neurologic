################
### PACKAGES ###
################

# install.packages(c("readr", "tableone", "tidyverse"))
library(readr)
library(tableone)
library(tidyverse)


#################
### LOAD DATA ###
#################

# load master dataset
master_data = read_csv("data/clean/master_data.csv")

# check data structure
glimpse(master_data)
cat("Total participants:", nrow(master_data), "\n")

########################
### DEFINE EXPOSURES ###
########################

# create binary exposure indicators
master_data = master_data |>
  mutate(
    # stroke exposure: - MCQ160F - ever told you had a stroke
    stroke_exposed = case_when(
      MCQ160F == 1 ~ 1,    # yes to stroke
      MCQ160F == 2 ~ 0,    # no to stroke
      TRUE ~ NA_real_
    ),
    # tbi exposure: - csq240 - consciousness loss due to head injury
    tbi_exposed = case_when(
      CSQ240 == 1 ~ 1,    # yes to TBI
      CSQ240 == 2 ~ 0,    # no to TBI
      TRUE ~ NA_real_
    ),
  )

################################
### MISSINGNESS & IMPUTATION ###
################################

# create missingness indicators and impute covariates used for matching

# report baseline missingness
cat("\n=== COVARIATE MISSINGNESS (pre-imputation) ===\n")
covars_to_check = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR")
print(sapply(master_data[covars_to_check], function(x) sum(is.na(x))))

# create missingness indicators BEFORE imputation
missingness_indicators = list()

# check each variable and create indicator only if it had missingness
for (var in covars_to_check) {
  if (sum(is.na(master_data[[var]])) > 0) {
    missingness_indicators[[paste0(var, "_missing")]] = as.integer(is.na(master_data[[var]]))
  }
}

# add missingness indicators to dataset
if (length(missingness_indicators) > 0) {
  master_data = master_data |>
    mutate(!!!missingness_indicators)
}

# now impute missing values
master_data = master_data |>
  mutate(
    # numeric variables - impute with mean if missing
    RIDAGEYR = ifelse(is.na(RIDAGEYR), mean(RIDAGEYR, na.rm = TRUE), RIDAGEYR),
    INDFMPIR = ifelse(is.na(INDFMPIR), mean(INDFMPIR, na.rm = TRUE), INDFMPIR),
    # Categorical variables - add explicit Missing level if missing
    RIAGENDR = forcats::fct_explicit_na(factor(RIAGENDR), na_level = "Missing"),
    RIDRETH3 = forcats::fct_explicit_na(factor(RIDRETH3), na_level = "Missing"),
    DMDEDUC2 = forcats::fct_explicit_na(factor(DMDEDUC2), na_level = "Missing")
  )

cat("\n=== COVARIATE MISSINGNESS (indicators created) ===\n")
if (length(missingness_indicators) > 0) {
  print(colSums(master_data[names(missingness_indicators)]))
} else {
  cat("No missingness indicators created - no missing values found\n")
}


##############################
### EXPOSURE DISTRIBUTIONS ###
##############################

cat("\n=== EXPOSURE DISTRIBUTIONS ===\n")

# stroke distribution
cat("Stroke (MCQ160f) distribution:\n")
print(table(master_data$MCQ160f, useNA = "ifany"))
cat("Stroke exposed:", sum(master_data$stroke_exposed == 1, na.rm = TRUE), "\n")
cat("Stroke not exposed:", sum(master_data$stroke_exposed == 0, na.rm = TRUE), "\n")
cat("Stroke missing:", sum(is.na(master_data$stroke_exposed)), "\n")

# tbi distribution
cat("\nTBI (CSQ240) distribution:\n")
print(table(master_data$CSQ240, useNA = "ifany"))
cat("TBI exposed:", sum(master_data$tbi_exposed == 1, na.rm = TRUE), "\n")
cat("TBI not exposed:", sum(master_data$tbi_exposed == 0, na.rm = TRUE), "\n")
cat("TBI missing:", sum(is.na(master_data$tbi_exposed)), "\n")

# Control groups (derived from exposure variables)
cat("\nStroke control group (no stroke):", sum(master_data$stroke_exposed == 0, na.rm = TRUE), "\n")
cat("TBI control group (no TBI):", sum(master_data$tbi_exposed == 0, na.rm = TRUE), "\n")

########################
### CALCULATE RATIOS ###
########################

# calculate treated/control ratios
stroke_ratio = sum(master_data$stroke_exposed == 0, na.rm = TRUE) / 
                sum(master_data$stroke_exposed == 1, na.rm = TRUE)
tbi_ratio = sum(master_data$tbi_exposed == 0, na.rm = TRUE) / 
             sum(master_data$tbi_exposed == 1, na.rm = TRUE)

cat("\n=== TREATED/CONTROL RATIOS ===\n")
cat("Stroke treated/control ratio:", round(stroke_ratio, 3), "\n")
cat("TBI treated/control ratio:", round(tbi_ratio, 3), "\n")

###################################
### DEMOGRAPHIC CHARACTERISTICS ###
###################################

# create demographic summary tables
cat("\n=== DEMOGRAPHIC CHARACTERISTICS ===\n")

# overall demographics
overall_demo = CreateTableOne(data = master_data, 
                               vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR"),
                               factorVars = c("RIAGENDR", "RIDRETH3", "DMDEDUC2"))
print(overall_demo)

# by stroke status
stroke_demo = CreateTableOne(data = master_data, 
                            vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR"),
                            factorVars = c("RIAGENDR", "RIDRETH3", "DMDEDUC2"),
                            strata = "stroke_exposed")
print(stroke_demo)

# by tbi status
tbi_demo = CreateTableOne(data = master_data, 
                          vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR"),
                          factorVars = c("RIAGENDR", "RIDRETH3", "DMDEDUC2"),
                          strata = "tbi_exposed")
print(tbi_demo)


################################
### SAVE DATA WITH EXPOSURES ###
################################

# save master dataset with exposure variables
write_csv(master_data, "data/clean/matching_data.csv")

cat("\n=== DATA EXPLORATION COMPLETE ===\n")
cat("Master dataset with exposures saved: matching_data.csv\n")
cat("Ready for matching analysis\n")