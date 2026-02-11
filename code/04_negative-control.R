################
### PACKAGES ###
################

# install.packages(c("survey", "tidyverse", "readr"))
library(survey)
library(tidyverse)
library(readr)

#################
### LOAD DATA ###
#################

# Load final matched datasets (1:4 matching, same as primary analysis)
stroke_data = read_csv("data/matched/stroke_1_4.csv", show_col_types = FALSE)
tbi_data = read_csv("data/matched/tbi_1_4.csv", show_col_types = FALSE)

# Load drug use data (DUQ200) - may not be in matched files if pipeline ran before drug merge
drug_data = read_csv("data/clean/drug_data.csv", show_col_types = FALSE)

# Merge DUQ200 into matched data if not present (join on SEQN + year for correct cycle)
if (!"DUQ200" %in% names(stroke_data)) {
  stroke_data = stroke_data |>
    left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
}
if (!"DUQ200" %in% names(tbi_data)) {
  tbi_data = tbi_data |>
    left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
}

#############################
### NEGATIVE CONTROL (Z) ###
#############################

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
tbi_data = tbi_data |> create_marijuana_ever()

###############################
### CREATE ANALYSIS WEIGHTS ###
###############################

# Same as primary analysis: w_analysis = (WTINT2YR/2) * weights

stroke_data = stroke_data |>
  mutate(
    w_nhanes = WTINT2YR / 2,
    w_analysis = w_nhanes * weights
  )

tbi_data = tbi_data |>
  mutate(
    w_nhanes = WTINT2YR / 2,
    w_analysis = w_nhanes * weights
  )

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
tbi_neg = run_negative_control(tbi_data, "tbi_exposed", "tbi_1_4")

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
#####################s

cat("\n=== NEGATIVE CONTROL BALANCE CHECK SUMMARY ===\n")
if (nrow(neg_control_results) > 0) {
  print(neg_control_results |>
    select(dataset, n_total, smd, or, or_ci_lower, or_ci_upper, balance_ok) |>
    mutate(across(c(smd, or, or_ci_lower, or_ci_upper), ~ round(.x, 4))))
}
cat("\nSuccessful balance: SMD < 0.10 and OR 95% CI includes 1\n")
cat("Report alongside main balance diagnostics (Sec. Negative control balance check)\n")
