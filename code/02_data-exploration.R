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

##############################
### DEFINE HEALTH OUTCOMES ###
##############################

# create health outcome variables for matching
master_data = master_data |>
  mutate(
    # alcohol abuse: ALQ141Q - number of days with 5+ drinks in past year
    # create binary indicator: 1 if any days with 5+ drinks in past year, 0 otherwise
    # Note: ALQ141Q is only asked if ALQ151 (ever had 5+ drinks) = 1
    # ALQ151: 1=Yes (ever had 5+ drinks), 2=No (never had 5+ drinks), 7=Refused, 9=Don't know
    # ALQ141Q: number of days (0-777), 999=Refused/Don't know
    alcohol_abuse = case_when(
      # If never had 5+ drinks (ALQ151 = 2), then no alcohol abuse
      ALQ151 == 2 ~ 0,
      # If ever had 5+ drinks and had some days in past year (ALQ141Q > 0, not 999)
      ALQ151 == 1 & ALQ141Q > 0 & ALQ141Q != 999 ~ 1,
      # If ever had 5+ drinks but 0 days in past year (ALQ141Q = 0)
      ALQ151 == 1 & ALQ141Q == 0 ~ 0,
      # If ever had 5+ drinks but ALQ141Q is missing or 999 (refused/don't know)
      ALQ151 == 1 & (is.na(ALQ141Q) | ALQ141Q == 999) ~ NA_real_,
      # Refused or don't know for ever had 5+ drinks question
      ALQ151 == 7 | ALQ151 == 9 ~ NA_real_,
      TRUE ~ NA_real_  # missing
    ),
    
    # smoking: SMQ040 - current smoking status
    # Note: SMQ040 is only asked if SMQ020 (ever smoked) = 1
    # If SMQ020 = 2 (never smoked), then SMQ040 is missing but should be "Not at all"
    # SMQ020: 1=Yes (ever smoked), 2=No (never smoked), 7=Refused, 9=Don't know
    # SMQ040: 1=Every day, 2=Some days, 3=Not at all, 7=Refused
    smoking_status = case_when(
      # If never smoked (SMQ020 = 2), then not a current smoker
      SMQ020 == 2 ~ "Not at all",
      # If ever smoked, use SMQ040 values
      SMQ020 == 1 & SMQ040 == 1 ~ "Every day",
      SMQ020 == 1 & SMQ040 == 2 ~ "Some days",
      SMQ020 == 1 & SMQ040 == 3 ~ "Not at all",
      SMQ020 == 1 & SMQ040 == 7 ~ "Refused",
      # Refused or don't know for ever smoked question
      SMQ020 == 7 ~ "Refused",
      SMQ020 == 9 ~ NA_character_,  # don't know
      TRUE ~ NA_character_           # missing
    ),
    
    # hypertension: BPQ020 - ever told you had hypertension
    # 1=Yes, 2=No, 9=Don't know
    hypertension = case_when(
      BPQ020 == 1 ~ 1,    # yes
      BPQ020 == 2 ~ 0,    # no
      TRUE ~ NA_real_     # missing or don't know
    ),
    
    # diabetes: DIQ010 - ever told you have diabetes
    # 1=Yes, 2=No, 3=Borderline, 7=Refused, 9=Don't know
    diabetes = case_when(
      DIQ010 == 1 ~ 1,    # yes
      DIQ010 == 2 ~ 0,    # no
      DIQ010 == 3 ~ 1,    # borderline (treat as yes)
      TRUE ~ NA_real_     # missing, refused, or don't know
    ),
    
    # stroke history: MCQ160F - ever told you had a stroke
    # (for TBI study only - already used as exposure for stroke study)
    # 1=Yes, 2=No, 7=Refused, 9=Don't know
    stroke_history = case_when(
      MCQ160F == 1 ~ 1,    # yes
      MCQ160F == 2 ~ 0,    # no
      TRUE ~ NA_real_      # missing, refused, or don't know
    )
  )

# convert smoking_status to factor
master_data = master_data |>
  mutate(
    smoking_status = factor(smoking_status, 
                            levels = c("Not at all", "Some days", "Every day", "Refused"))
  )

################################
### MISSINGNESS & IMPUTATION ###
################################

# create missingness indicators and impute covariates used for matching

# report baseline missingness
cat("\n=== COVARIATE MISSINGNESS (pre-imputation) ===\n")
covars_to_check = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
                    "alcohol_abuse", "smoking_status", "hypertension", "diabetes", "stroke_history")
print(sapply(master_data[covars_to_check], function(x) sum(is.na(x))))

# create missingness indicators BEFORE imputation
missingness_indicators = list()

# check each variable and create indicator only if it had missingness
# (check before imputation)
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

# calculate modes for binary health outcomes (for imputation)
get_mode = function(x) {
  tab = table(x, useNA = "no")
  if (length(tab) > 0) {
    return(as.numeric(names(which.max(tab))))
  } else {
    return(0)  # default to 0 if all missing
  }
}

alcohol_mode = get_mode(master_data$alcohol_abuse)
hypertension_mode = get_mode(master_data$hypertension)
diabetes_mode = get_mode(master_data$diabetes)
stroke_history_mode = get_mode(master_data$stroke_history)

# now impute missing values
master_data = master_data |>
  mutate(
    # numeric variables - impute with mean if missing
    RIDAGEYR = ifelse(is.na(RIDAGEYR), mean(RIDAGEYR, na.rm = TRUE), RIDAGEYR),
    INDFMPIR = ifelse(is.na(INDFMPIR), mean(INDFMPIR, na.rm = TRUE), INDFMPIR),
    # Categorical variables - add explicit Missing level if missing
    RIAGENDR = forcats::fct_explicit_na(factor(RIAGENDR), na_level = "Missing"),
    RIDRETH3 = forcats::fct_explicit_na(factor(RIDRETH3), na_level = "Missing"),
    DMDEDUC2 = forcats::fct_explicit_na(factor(DMDEDUC2), na_level = "Missing"),
    # Health outcomes - add explicit Missing level if missing
    smoking_status = forcats::fct_explicit_na(smoking_status, na_level = "Missing"),
    # Binary health outcomes - impute with mode (most common value) if missing
    alcohol_abuse = ifelse(is.na(alcohol_abuse), alcohol_mode, alcohol_abuse),
    hypertension = ifelse(is.na(hypertension), hypertension_mode, hypertension),
    diabetes = ifelse(is.na(diabetes), diabetes_mode, diabetes),
    stroke_history = ifelse(is.na(stroke_history), stroke_history_mode, stroke_history)
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
cat("Stroke (MCQ160F) distribution:\n")
print(table(master_data$MCQ160F, useNA = "ifany"))
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