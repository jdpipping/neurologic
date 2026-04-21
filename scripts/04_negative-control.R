################
### PACKAGES ###
################

# install.packages(c("survey", "tidyverse", "readr", "cobalt", "MatchIt"))
library(survey)
library(tidyverse)
library(readr)
library(cobalt)
library(MatchIt)

#################
### LOAD DATA ###
#################

# Load final matched datasets (1:4 matching, same as primary analysis)
stroke_data = read_csv("data/matched/stroke_1_4.csv", show_col_types = FALSE)

# Load matching data for revised TBI matching (with marijuana in PS model)
matching_data = read_csv("data/clean/matching_data.csv", show_col_types = FALSE)

# Load drug use data (DUQ200) - may not be in matched files if pipeline ran before drug merge
drug_data = read_csv("data/clean/drug_data.csv", show_col_types = FALSE)

# Merge DUQ200 into matched data if not present (join on SEQN + year for correct cycle)
if (!"DUQ200" %in% names(stroke_data)) {
  stroke_data = stroke_data |>
    left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
}

############################
### NEGATIVE CONTROL (Z) ###
############################

# DUQ200: "Ever used marijuana or hashish?"
# 1 = Yes -> Z = 1
# 2 = No -> Z = 0
# 7 = Refused, 9 = Don't know, NA = not administered -> Z = NA

create_marijuana_ever = function(df) {
  df |>
    mutate(
      marijuana_ever = case_when(
        DUQ200 == 1 ~ 1,   # Yes
        DUQ200 == 2 ~ 0,   # No
        TRUE ~ NA_real_    # Refused, Don't know, not administered
      )
    )
}

stroke_data = stroke_data |> create_marijuana_ever()

apply_caliper_penalty = function(distance_matrix, caliper, penalty = 10000) {
  distance_matrix[distance_matrix > caliper] = distance_matrix[distance_matrix > caliper] + penalty
  distance_matrix
}

build_original_tbi_match = function(matching_data, drug_data, ratio = 4) {
  tbi_matching_data = matching_data
  if (!"DUQ200" %in% names(tbi_matching_data)) {
    tbi_matching_data = tbi_matching_data |>
      left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
  }

  missingness_vars = c(
    "RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
    "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing",
    "diabetes_missing", "stroke_history_missing"
  )
  existing_missingness_vars = missingness_vars[missingness_vars %in% names(tbi_matching_data)]
  base_matching_vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2")
  design_vars = c("SDMVPSU", "SDMVSTRA")
  sample_weight_var = "WTINT2YR"
  health_outcome_vars = c("alcohol_abuse", "smoking_status", "hypertension", "diabetes")
  tbi_matching_vars = c(base_matching_vars, design_vars, sample_weight_var, health_outcome_vars, "stroke_history")
  tbi_missingness_vars = existing_missingness_vars[existing_missingness_vars %in%
    c(paste0(base_matching_vars, "_missing"),
      paste0(health_outcome_vars, "_missing"),
      "stroke_history_missing")]
  tbi_matching_vars_all = c(tbi_matching_vars, tbi_missingness_vars)
  match_formula_tbi = as.formula(paste("tbi_exposed ~", paste(tbi_matching_vars_all, collapse = " + ")))

  tbi_matching_data = tbi_matching_data |>
    filter(!is.na(tbi_exposed))

  ps_model_tbi = glm(match_formula_tbi, data = tbi_matching_data, family = binomial)
  tbi_matching_data$ps = predict(ps_model_tbi, type = "response")
  tbi_matching_data$logit_ps = predict(ps_model_tbi, type = "link")

  treated_logit = tbi_matching_data$logit_ps[tbi_matching_data$tbi_exposed == 1]
  control_logit = tbi_matching_data$logit_ps[tbi_matching_data$tbi_exposed == 0]
  ps_caliper = 0.2 * sqrt((sd(treated_logit, na.rm = TRUE)^2 + sd(control_logit, na.rm = TRUE)^2) / 2)

  treated_idx = which(tbi_matching_data$tbi_exposed == 1)
  control_idx = which(tbi_matching_data$tbi_exposed == 0)
  ps_dist = abs(outer(tbi_matching_data$logit_ps[treated_idx], tbi_matching_data$logit_ps[control_idx], "-"))
  ps_dist = apply_caliper_penalty(ps_dist, ps_caliper)

  match_obj = tryCatch(
    matchit(
      match_formula_tbi,
      data = tbi_matching_data,
      method = "optimal",
      distance = ps_dist,
      ratio = ratio
    ),
    error = function(e) NULL
  )

  if (is.null(match_obj)) {
    return(NULL)
  }

  match.data(match_obj) |>
    create_marijuana_ever()
}

tbi_data = build_original_tbi_match(matching_data, drug_data)
if (is.null(tbi_data)) {
  cat("Original TBI 1:4 matching failed under the penalized caliper setup; skipping original TBI negative-control diagnostic.\n")
}

###############################
### CREATE ANALYSIS WEIGHTS ###
###############################

# Same as primary analysis: w_analysis = (WTINT2YR/2) * weights

stroke_data = stroke_data |>
  mutate(
    w_nhanes = WTINT2YR / 2,
    w_analysis = w_nhanes * weights
  )

if (!is.null(tbi_data)) {
  tbi_data = tbi_data |>
    mutate(
      w_nhanes = WTINT2YR / 2,
      w_analysis = w_nhanes * weights
    )
}

###################################
### NEGATIVE CONTROL DIAGNOSTIC ###
###################################

# For each matched dataset:
# 1. Compute weighted SMD for Z between exposed and unexposed
# 2. Fit survey-weighted logistic regression: Z ~ A, report exp(alpha1) with 95% CI
# Restrict to respondents with observed Z (nonmissing marijuana_ever)

# Weighted SMD for binary variable
# SMD = (p1 - p0) / sqrt(((p1*(1-p1) + p0*(1-p0)) / 2))
# where p1 = weighted proportion in exposed, p0 = weighted proportion in unexposed
compute_weighted_smd = function(df, z_var, a_var, w_var) {
  df_complete = df |> filter(!is.na(.data[[z_var]]), !is.na(.data[[a_var]]))

  exposed = df_complete |> filter(.data[[a_var]] == 1)
  unexposed = df_complete |> filter(.data[[a_var]] == 0)

  w1 = sum(exposed[[w_var]])
  w0 = sum(unexposed[[w_var]])
  p1 = sum(exposed[[w_var]] * exposed[[z_var]]) / w1
  p0 = sum(unexposed[[w_var]] * unexposed[[z_var]]) / w0

  pooled_var = ((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  pooled_sd = sqrt(max(pooled_var, 1e-10))  # avoid division by zero
  smd = (p1 - p0) / pooled_sd

  list(
    smd = smd,
    p_exposed = p1,
    p_unexposed = p0,
    n_exposed = nrow(exposed),
    n_unexposed = nrow(unexposed),
    n_total = nrow(df_complete)
  )
}

# Run negative control diagnostic for one matched dataset
run_negative_control = function(df, exposure_var, dataset_name) {
  cat("\n=== Negative control:", dataset_name, "===\n")

  df_complete = df |>
    filter(!is.na(marijuana_ever), !is.na(.data[[exposure_var]]))

  cat("Matched sample with observed marijuana_ever:", nrow(df_complete), "\n")
  cat("  Exposed:", sum(df_complete[[exposure_var]] == 1), "\n")
  cat("  Unexposed:", sum(df_complete[[exposure_var]] == 0), "\n")

  if (nrow(df_complete) == 0 || sum(df_complete[[exposure_var]] == 1) == 0 || sum(df_complete[[exposure_var]] == 0) == 0) {
    cat("Insufficient data for negative control diagnostic\n")
    return(NULL)
  }

  # 1. Weighted SMD
  smd_result = compute_weighted_smd(df_complete, "marijuana_ever", exposure_var, "w_analysis")
  cat("\n1. Standardized mean difference (SMD):\n")
  cat("   SMD =", round(smd_result$smd, 4), "\n")
  cat("   Weighted proportion exposed:", round(smd_result$p_exposed, 4), "\n")
  cat("   Weighted proportion unexposed:", round(smd_result$p_unexposed, 4), "\n")
  cat("   Balance criterion: SMD < 0.10", if (abs(smd_result$smd) < 0.10) "OK" else "VIOLATION", "\n")

  # 2. Survey-weighted logistic regression
  design_vars = c("SDMVPSU", "SDMVSTRA")
  svy_design = svydesign(
    id = ~ SDMVPSU,
    strata = ~ SDMVSTRA,
    weights = ~ w_analysis,
    data = df_complete,
    nest = TRUE
  )

  formula_str = paste("marijuana_ever ~", exposure_var)
  svy_model = svyglm(
    formula = as.formula(formula_str),
    design = svy_design,
    family = quasibinomial(link = "logit")
  )

  coef_summary = summary(svy_model)$coefficients
  alpha1 = coef_summary[exposure_var, "Estimate"]
  se_alpha1 = coef_summary[exposure_var, "Std. Error"]

  or = exp(alpha1)
  or_ci_lower = exp(alpha1 - 1.96 * se_alpha1)
  or_ci_upper = exp(alpha1 + 1.96 * se_alpha1)

  cat("\n2. Survey-weighted logistic regression (marijuana_ever ~ exposure):\n")
  cat("   OR (95% CI) =", round(or, 3), "(", round(or_ci_lower, 3), ",", round(or_ci_upper, 3), ")\n")
  cat("   Balance criterion: OR close to 1", if (or_ci_lower <= 1 && or_ci_upper >= 1) "OK" else "Review", "\n")

  data.frame(
    dataset = dataset_name,
    n_total = smd_result$n_total,
    n_exposed = smd_result$n_exposed,
    n_unexposed = smd_result$n_unexposed,
    smd = smd_result$smd,
    p_exposed = smd_result$p_exposed,
    p_unexposed = smd_result$p_unexposed,
    or = or,
    or_ci_lower = or_ci_lower,
    or_ci_upper = or_ci_upper,
    balance_ok = abs(smd_result$smd) < 0.10 & or_ci_lower <= 1 & or_ci_upper >= 1,
    stringsAsFactors = FALSE
  )
}

#######################
### RUN DIAGNOSTICS ###
#######################

stroke_neg = run_negative_control(stroke_data, "stroke_exposed", "stroke_1_4")
tbi_neg = if (!is.null(tbi_data)) run_negative_control(tbi_data, "tbi_exposed", "tbi_1_4") else NULL

###############################################
### REVISED TBI MATCHING (marijuana in PS) ###
###############################################

# TBI had imbalance on marijuana ever-use. Re-run 1:4 matching with marijuana_ever in the PS model.

cat("\n=== REVISED TBI MATCHING (1:4 with marijuana ever-use in PS model) ===\n")

# Merge drug data if DUQ200 not in matching_data, then create marijuana_ever
if (!"DUQ200" %in% names(matching_data)) {
  matching_data = matching_data |>
    left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
}
matching_data = matching_data |>
  mutate(
    marijuana_ever = case_when(
      DUQ200 == 1 ~ 1,
      DUQ200 == 2 ~ 0,
      TRUE ~ NA_real_
    )
  )

# Restrict to TBI sample
tbi_matching_data = matching_data |>
  filter(!is.na(tbi_exposed))

# Add missingness indicator and impute marijuana_ever (for PS model)
marijuana_mode = as.numeric(names(which.max(table(tbi_matching_data$marijuana_ever, useNA = "no"))))
if (length(marijuana_mode) == 0) marijuana_mode = 0
tbi_matching_data = tbi_matching_data |>
  mutate(
    marijuana_ever_missing = as.integer(is.na(marijuana_ever)),
    marijuana_ever = ifelse(is.na(marijuana_ever), marijuana_mode, marijuana_ever)
  )

# TBI matching vars: original + marijuana_ever + marijuana_ever_missing
base_matching_vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2")
design_vars = c("SDMVPSU", "SDMVSTRA")
sample_weight_var = "WTINT2YR"
health_outcome_vars = c("alcohol_abuse", "smoking_status", "hypertension", "diabetes")
missingness_vars = c(
  "RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
  "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing",
  "diabetes_missing", "stroke_history_missing"
)
existing_missingness_vars = missingness_vars[missingness_vars %in% names(tbi_matching_data)]
tbi_missingness_vars_revised = existing_missingness_vars[existing_missingness_vars %in%
  c(paste0(base_matching_vars, "_missing"),
    paste0(health_outcome_vars, "_missing"),
    "stroke_history_missing")]
tbi_matching_vars_revised = c(
  base_matching_vars, design_vars, sample_weight_var,
  health_outcome_vars, "stroke_history", tbi_missingness_vars_revised,
  "marijuana_ever", "marijuana_ever_missing"
)

match_formula_tbi_revised = as.formula(paste(
  "tbi_exposed ~",
  paste(tbi_matching_vars_revised, collapse = " + ")
))

# Fit PS model
ps_model_tbi_revised = glm(match_formula_tbi_revised, data = tbi_matching_data, family = binomial)
tbi_matching_data$ps = predict(ps_model_tbi_revised, type = "response")
tbi_matching_data$logit_ps = predict(ps_model_tbi_revised, type = "link")

# Caliper: 0.2 SD of logit(PS)
treated_logit = tbi_matching_data$logit_ps[tbi_matching_data$tbi_exposed == 1]
control_logit = tbi_matching_data$logit_ps[tbi_matching_data$tbi_exposed == 0]
spool = sqrt((sd(treated_logit, na.rm = TRUE)^2 + sd(control_logit, na.rm = TRUE)^2) / 2)
ps_caliper = 0.2 * spool

# 1:4 matching
treated_idx = which(tbi_matching_data$tbi_exposed == 1)
control_idx = which(tbi_matching_data$tbi_exposed == 0)
ps_dist = abs(outer(tbi_matching_data$logit_ps[treated_idx], tbi_matching_data$logit_ps[control_idx], "-"))
ps_dist = apply_caliper_penalty(ps_dist, ps_caliper)

match_tbi_revised = matchit(
  match_formula_tbi_revised,
  data = tbi_matching_data,
  method = "optimal",
  distance = ps_dist,
  ratio = 4
)

tbi_revised_data = match.data(match_tbi_revised)
cat("Matched sample:", sum(tbi_revised_data$tbi_exposed == 1), "treated,",
    sum(tbi_revised_data$tbi_exposed == 0), "controls\n")

# Love plot
love_tbi_revised = love.plot(
  match_tbi_revised,
  binary = "std",
  title = "Balance Plot: TBI 1:4 (incl. marijuana use)",
  thresholds = c(m = 0.1)
)
if (!dir.exists("matching")) {
  dir.create("matching", recursive = TRUE)
}
ggsave("matching/love_tbi_1_4_revised.png", love_tbi_revised, width = 6, height = 4, dpi = 600)
cat("Love plot saved: matching/love_tbi_1_4_revised.png\n")

# Save revised matched dataset for protocol-specified primary outcome analysis
if (!dir.exists("data/matched")) {
  dir.create("data/matched", recursive = TRUE)
}
write_csv(tbi_revised_data, "data/matched/tbi_1_4_revised.csv")
cat("Revised matched dataset saved: data/matched/tbi_1_4_revised.csv\n")

####################
### SAVE RESULTS ###
####################

if (!dir.exists("results")) {
  dir.create("results")
}

neg_control_results = bind_rows(
  if (!is.null(stroke_neg)) stroke_neg else tibble(),
  if (!is.null(tbi_neg)) tbi_neg else tibble()
)

if (nrow(neg_control_results) > 0) {
  write_csv(neg_control_results, "results/negative_control_results.csv")
  cat("\n=== Negative control results saved to: results/negative_control_results.csv ===\n")
}

#####################
### SUMMARY TABLE ###
#####################

cat("\n=== NEGATIVE CONTROL BALANCE CHECK SUMMARY ===\n")
if (nrow(neg_control_results) > 0) {
  print(neg_control_results |>
    select(dataset, n_total, smd, or, or_ci_lower, or_ci_upper, balance_ok) |>
    mutate(across(c(smd, or, or_ci_lower, or_ci_upper), ~ round(.x, 4))))
}
cat("\nSuccessful balance: SMD < 0.10 and OR 95% CI includes 1\n")
cat("Report alongside main balance diagnostics (Sec. Negative control balance check)\n")
