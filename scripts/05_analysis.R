################
### PACKAGES ###
################

# install.packages(c("survey", "tidyverse", "kableExtra", "ggplot2"))
library(survey)
library(tidyverse)
library(kableExtra)
library(ggplot2)

#################
### LOAD DATA ###
#################

# Protocol-specified matched datasets:
# - Stroke: selected 1:4 design
# - TBI: revised 1:4 design after negative-control imbalance check
stroke_match_path = "data/matched/stroke_1_4.csv"
tbi_match_path = "data/matched/tbi_1_4_revised.csv"

if (!file.exists(stroke_match_path)) {
  stop("Missing stroke matched dataset: data/matched/stroke_1_4.csv")
}
if (!file.exists(tbi_match_path)) {
  stop("Protocol requires revised TBI 1:4 matched data. Run scripts/04_negative-control.R to create data/matched/tbi_1_4_revised.csv")
}

stroke_data = read_csv(stroke_match_path, show_col_types = FALSE)
tbi_data = read_csv(tbi_match_path, show_col_types = FALSE)

cat("Stroke matched dataset:", nrow(stroke_data), "participants\n")
cat("TBI revised matched dataset:", nrow(tbi_data), "participants\n")

# Primary outcomes and exposures
primary_outcomes = c("usual_place", "any_insurance")
stroke_exposure = "stroke_exposed"
tbi_exposure = "tbi_exposed"
design_vars = c("SDMVPSU", "SDMVSTRA")

# Prespecified covariates for post-match residual imbalance check (|SMD| > 0.10)
stroke_covariates = c(
  "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2",
  "alcohol_abuse", "smoking_status", "hypertension", "diabetes",
  "RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
  "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing", "diabetes_missing"
)
tbi_covariates = c(
  "RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2",
  "alcohol_abuse", "smoking_status", "hypertension", "diabetes", "stroke_history",
  "RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
  "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing", "diabetes_missing", "stroke_history_missing",
  "marijuana_ever", "marijuana_ever_missing"
)

##############################
### HELPER FUNCTIONS       ###
##############################

weighted_mean = function(x, w) {
  sum(w * x) / sum(w)
}

weighted_var = function(x, w) {
  mu = weighted_mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

binary_smd = function(z, a, w) {
  p1 = weighted_mean(z[a == 1], w[a == 1])
  p0 = weighted_mean(z[a == 0], w[a == 0])
  pooled_var = (p1 * (1 - p1) + p0 * (1 - p0)) / 2
  if (!is.finite(pooled_var) || pooled_var <= 0) return(NA_real_)
  abs((p1 - p0) / sqrt(pooled_var))
}

compute_abs_smd = function(df, var, exposure_var, weight_var = "weights", continuous_vars = c("RIDAGEYR", "INDFMPIR")) {
  d = df |>
    filter(
      !is.na(.data[[var]]),
      !is.na(.data[[exposure_var]]),
      .data[[exposure_var]] %in% c(0, 1),
      !is.na(.data[[weight_var]]),
      is.finite(.data[[weight_var]]),
      .data[[weight_var]] > 0
    )

  if (nrow(d) == 0) return(NA_real_)

  a = d[[exposure_var]]
  w = d[[weight_var]]
  x = d[[var]]
  if (sum(a == 1) == 0 || sum(a == 0) == 0) return(NA_real_)

  is_continuous = is.numeric(x) && (var %in% continuous_vars)
  if (is_continuous) {
    m1 = weighted_mean(x[a == 1], w[a == 1])
    m0 = weighted_mean(x[a == 0], w[a == 0])
    v1 = weighted_var(x[a == 1], w[a == 1])
    v0 = weighted_var(x[a == 0], w[a == 0])
    pooled_sd = sqrt((v1 + v0) / 2)
    if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
    return(abs((m1 - m0) / pooled_sd))
  }

  xf = as.factor(x)
  levs = levels(xf)
  smds = map_dbl(levs, function(lev) {
    z = as.numeric(xf == lev)
    binary_smd(z, a, w)
  })
  if (all(is.na(smds))) return(NA_real_)
  max(smds, na.rm = TRUE)
}

find_imbalanced_covariates = function(df, exposure_var, candidate_vars, threshold = 0.10) {
  covars = candidate_vars[candidate_vars %in% names(df)]
  if (length(covars) == 0) {
    return(list(imbalanced = character(0), smd_table = tibble(covariate = character(), abs_smd = double())))
  }

  if (requireNamespace("cobalt", quietly = TRUE)) {
    balance_formula = as.formula(paste(exposure_var, "~", paste(covars, collapse = " + ")))
    bal_obj = cobalt::bal.tab(
      balance_formula,
      data = df,
      weights = df$weights,
      method = "weighting",
      estimand = "ATT",
      s.d.denom = "pooled",
      continuous = "std",
      binary = "std",
      quick = FALSE
    )
    bal_df = as.data.frame(bal_obj$Balance)
    bal_df$term = rownames(bal_obj$Balance)
    diff_col = if ("Diff.Adj" %in% names(bal_df)) "Diff.Adj" else "Diff.Un"
    bal_df = bal_df |>
      mutate(abs_smd = abs(.data[[diff_col]]))

    map_term_to_covariate = function(term) {
      hits = covars[term == covars | startsWith(term, paste0(covars, "_"))]
      if (length(hits) == 0) return(term)
      hits[which.max(nchar(hits))]
    }

    smd_table = bal_df |>
      mutate(covariate = map_chr(term, map_term_to_covariate)) |>
      group_by(covariate) |>
      summarize(
        abs_smd = if (all(is.na(abs_smd))) NA_real_ else max(abs_smd, na.rm = TRUE),
        .groups = "drop"
      )

    smd_table = tibble(covariate = covars) |>
      left_join(smd_table, by = "covariate") |>
      mutate(imbalanced = is.finite(abs_smd) & abs_smd > threshold)
  } else {
    # Fallback when cobalt is unavailable
    smd_table = tibble(
      covariate = covars,
      abs_smd = map_dbl(covars, ~ compute_abs_smd(df, .x, exposure_var, weight_var = "weights"))
    ) |>
      mutate(imbalanced = is.finite(abs_smd) & abs_smd > threshold)
  }

  list(
    imbalanced = smd_table |> filter(imbalanced) |> arrange(desc(abs_smd)) |> pull(covariate),
    smd_table = smd_table
  )
}

summarize_weight_diagnostics = function(df, dataset_name, weight_var = "w_analysis") {
  w = df[[weight_var]]
  w = w[is.finite(w) & !is.na(w) & w > 0]
  if (length(w) == 0) {
    return(tibble(
      dataset = dataset_name,
      n_total = nrow(df),
      n_positive_weights = 0,
      min_weight = NA_real_,
      p01_weight = NA_real_,
      median_weight = NA_real_,
      p99_weight = NA_real_,
      max_weight = NA_real_,
      mean_weight = NA_real_,
      sd_weight = NA_real_,
      ess = NA_real_,
      ess_ratio = NA_real_
    ))
  }

  ess_val = (sum(w)^2) / sum(w^2)
  tibble(
    dataset = dataset_name,
    n_total = nrow(df),
    n_positive_weights = length(w),
    min_weight = min(w),
    p01_weight = as.numeric(quantile(w, 0.01)),
    median_weight = median(w),
    p99_weight = as.numeric(quantile(w, 0.99)),
    max_weight = max(w),
    mean_weight = mean(w),
    sd_weight = sd(w),
    ess = ess_val,
    ess_ratio = ess_val / length(w)
  )
}

filter_protocol_eligible = function(df, exposure_var, cohort_label) {
  required_cols = c(exposure_var, "WTINT2YR", "weights", design_vars)
  missing_cols = setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0(cohort_label, " is missing required protocol variables: ", paste(missing_cols, collapse = ", ")))
  }

  n_before = nrow(df)
  d = df |>
    filter(
      !is.na(.data[[exposure_var]]),
      .data[[exposure_var]] %in% c(0, 1),
      !is.na(WTINT2YR), is.finite(WTINT2YR), WTINT2YR > 0,
      !is.na(weights), is.finite(weights), weights > 0,
      !is.na(.data[[design_vars[1]]]), !is.na(.data[[design_vars[2]]])
    )
  cat(cohort_label, "participants eligible for survey analysis:", nrow(d), "of", n_before, "\n")
  d
}

##############################
### CREATE OUTCOME VARIABLES #
##############################

# Primary outcome 1: Usual place to go for healthcare (HUQ030)
# 1 = Has a usual place (response = 1 or 3), 0 = No usual place (response = 2)
# Primary outcome 2: Any health insurance coverage (HIQ011)
# 1 = Any coverage (public or private), 0 = Uninsured (response = 2)

stroke_data = stroke_data |>
  mutate(
    usual_place = case_when(
      HUQ030 == 1 | HUQ030 == 3 ~ 1,  # Has a usual place
      HUQ030 == 2 ~ 0,                 # No usual place
      TRUE ~ NA_real_                  # Missing/refused/don't know
    ),
    any_insurance = case_when(
      HIQ011 == 1 ~ 1,    # Yes, covered
      HIQ011 == 2 ~ 0,     # No, not covered
      TRUE ~ NA_real_      # Missing/refused/don't know
    )
  )

tbi_data = tbi_data |>
  mutate(
    usual_place = case_when(
      HUQ030 == 1 | HUQ030 == 3 ~ 1,  # Has a usual place
      HUQ030 == 2 ~ 0,                 # No usual place
      TRUE ~ NA_real_                  # Missing/refused/don't know
    ),
    any_insurance = case_when(
      HIQ011 == 1 ~ 1,    # Yes, covered
      HIQ011 == 2 ~ 0,     # No, not covered
      TRUE ~ NA_real_      # Missing/refused/don't know
    )
  )

###############################
### CREATE ANALYSIS WEIGHTS ###
###############################

# Protocol: w_i = (WTINT2YR/2) * matching_weight
# Rescale NHANES weight by number of cycles (2 cycles pooled)
# Final analysis weight = NHANES weight * matching weight

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

#######################
### STROKE ANALYSIS ###
#######################

# Enforce protocol inclusion criteria for valid survey design-based outcome analysis
stroke_data = filter_protocol_eligible(stroke_data, stroke_exposure, "Stroke")
tbi_data = filter_protocol_eligible(tbi_data, tbi_exposure, "TBI")

# Post-match imbalance check for prespecified covariates (protocol robustness rule)
stroke_balance = find_imbalanced_covariates(stroke_data, stroke_exposure, stroke_covariates, threshold = 0.10)
tbi_balance = find_imbalanced_covariates(tbi_data, tbi_exposure, tbi_covariates, threshold = 0.10)

stroke_adjustment_vars = stroke_balance$imbalanced
tbi_adjustment_vars = tbi_balance$imbalanced

cat("\n=== POST-MATCH COVARIATE BALANCE CHECK (|SMD| > 0.10) ===\n")
cat("Stroke imbalanced covariates:", ifelse(length(stroke_adjustment_vars) == 0, "None", paste(stroke_adjustment_vars, collapse = ", ")), "\n")
cat("TBI imbalanced covariates:", ifelse(length(tbi_adjustment_vars) == 0, "None", paste(tbi_adjustment_vars, collapse = ", ")), "\n")

balance_smd_table = bind_rows(
  stroke_balance$smd_table |> mutate(dataset = "stroke_1_4"),
  tbi_balance$smd_table |> mutate(dataset = "tbi_1_4_revised")
) |>
  relocate(dataset, covariate, abs_smd, imbalanced)

cat("\n=== STROKE ANALYSIS ===\n")

stroke_results_list = list()

for (outcome in primary_outcomes) {
  cat("\n--- Analyzing outcome:", outcome, "---\n")
  
  # Filter to complete cases for this outcome
  stroke_complete = stroke_data |>
    filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[stroke_exposure]]),
      !is.na(w_analysis), is.finite(w_analysis), w_analysis > 0,
      !is.na(.data[[design_vars[1]]]), !is.na(.data[[design_vars[2]]])
    )
  
  cat("Complete cases:", nrow(stroke_complete), "\n")
  
  # Check if we have both exposure groups
  exposure_counts = table(stroke_complete[[stroke_exposure]], useNA = "ifany")
  cat("Exposure distribution:\n")
  print(exposure_counts)
  
  if (nrow(stroke_complete) == 0 || length(exposure_counts) < 2 || any(exposure_counts < 1)) {
    cat("Insufficient data - skipping\n")
    next
  }
  
  # Create survey design object
  id_formula = as.formula(paste("~", design_vars[1]))
  strata_formula = as.formula(paste("~", design_vars[2]))
  weights_formula = as.formula("~ w_analysis")
  
  svy_design = svydesign(
    id = id_formula,        # PSU
    strata = strata_formula, # Strata
    weights = weights_formula, # Analysis weights
    data = stroke_complete,
    nest = TRUE
  )
  
  # Create formula: outcome ~ exposure
  formula_str = paste(outcome, "~", stroke_exposure)
  formula_obj = as.formula(formula_str)
  
  # Fit survey-weighted logistic regression
  # Use quasibinomial() as specified in protocol
  svy_model = svyglm(
    formula = formula_obj,
    design = svy_design,
    family = quasibinomial(link = "logit")
  )
  
  # Extract coefficients
  coef_summary = summary(svy_model)
  coef_table = coef_summary$coefficients
  
  # Get exposure coefficient
  exposure_coef = coef_table[stroke_exposure, ]
  beta1 = exposure_coef["Estimate"]
  se_beta1 = exposure_coef["Std. Error"]
  p_unadj = exposure_coef["Pr(>|t|)"]

  # Protocol robustness check: if post-match imbalance remains (|SMD| > 0.10),
  # fit an adjusted outcome model including imbalanced prespecified covariates.
  adj_or = NA_real_
  adj_or_ci_lower = NA_real_
  adj_or_ci_upper = NA_real_
  adj_p_unadj = NA_real_
  adj_vars_for_model = character(0)
  stroke_adjustment_vars_present = stroke_adjustment_vars[stroke_adjustment_vars %in% names(stroke_complete)]
  if (length(stroke_adjustment_vars_present) > 0) {
    stroke_adjustment_vars_present = stroke_adjustment_vars_present[
      map_lgl(stroke_adjustment_vars_present, ~ any(!is.na(stroke_complete[[.x]])))
    ]
    if (length(stroke_adjustment_vars_present) > 0) {
      adjusted_formula = as.formula(
        paste(outcome, "~", stroke_exposure, "+", paste(stroke_adjustment_vars_present, collapse = " + "))
      )
      adjusted_model = tryCatch(
        svyglm(formula = adjusted_formula, design = svy_design, family = quasibinomial(link = "logit")),
        error = function(e) NULL
      )
      if (!is.null(adjusted_model)) {
        adjusted_coef = summary(adjusted_model)$coefficients
        if (stroke_exposure %in% rownames(adjusted_coef)) {
          beta1_adj = adjusted_coef[stroke_exposure, "Estimate"]
          se_beta1_adj = adjusted_coef[stroke_exposure, "Std. Error"]
          adj_or = exp(beta1_adj)
          adj_or_ci_lower = exp(beta1_adj - 1.96 * se_beta1_adj)
          adj_or_ci_upper = exp(beta1_adj + 1.96 * se_beta1_adj)
          adj_p_unadj = adjusted_coef[stroke_exposure, "Pr(>|t|)"]
          adj_vars_for_model = stroke_adjustment_vars_present
        }
      }
    }
  }
  
  # Compute OR and 95% CI on log-odds scale, then transform
  log_or = beta1
  log_or_se = se_beta1
  log_or_ci_lower = log_or - 1.96 * log_or_se
  log_or_ci_upper = log_or + 1.96 * log_or_se
  
  # Transform to OR scale
  or = exp(log_or)
  or_ci_lower = exp(log_or_ci_lower)
  or_ci_upper = exp(log_or_ci_upper)
  
  # Compute marginal predicted probabilities
  # Protocol: "marginal predicted probabilities under A=1 and A=0"
  # For a model with only exposure (no other covariates), predictions are constant
  # We compute them directly from model coefficients
  coef_names = names(coef(svy_model))
  beta0 = coef(svy_model)[1]  # intercept
  beta1_idx = which(coef_names == stroke_exposure)
  beta1_coef = coef(svy_model)[beta1_idx]
  
  # Marginal probabilities from model
  # P(Y=1|A=0) = expit(beta0)
  # P(Y=1|A=1) = expit(beta0 + beta1)
  prob_0 = plogis(beta0)
  prob_1 = plogis(beta0 + beta1_coef)
  
  # For SEs, use delta method
  vcov_model = vcov(svy_model)
  se_beta0 = sqrt(vcov_model[1, 1])
  cov_beta0_beta1 = vcov_model[1, beta1_idx]
  
  # SE for prob_0 = expit(beta0)
  prob_se_0 = prob_0 * (1 - prob_0) * se_beta0
  
  # SE for prob_1 = expit(beta0 + beta1) using delta method
  # Var(beta0 + beta1) = Var(beta0) + Var(beta1) + 2*Cov(beta0, beta1)
  var_sum = vcov_model[1, 1] + vcov_model[beta1_idx, beta1_idx] + 2 * cov_beta0_beta1
  se_sum = sqrt(var_sum)
  prob_se_1 = prob_1 * (1 - prob_1) * se_sum
  
  # Risk difference and CI
  rd = prob_1 - prob_0
  rd_se = sqrt(prob_se_0^2 + prob_se_1^2)  # Conservative approximation
  rd_ci_lower = rd - 1.96 * rd_se
  rd_ci_upper = rd + 1.96 * rd_se
  
  # Store results
  results = data.frame(
    outcome = outcome,
    exposure = stroke_exposure,
    n = nrow(stroke_complete),
    beta1 = beta1,
    se_beta1 = se_beta1,
    or = or,
    or_ci_lower = or_ci_lower,
    or_ci_upper = or_ci_upper,
    p_unadj = p_unadj,
    prob_exposed = prob_1,
    prob_unexposed = prob_0,
    rd = rd,
    rd_ci_lower = rd_ci_lower,
    rd_ci_upper = rd_ci_upper,
    adjusted_model_run = length(adj_vars_for_model) > 0,
    n_adjustment_covariates = length(adj_vars_for_model),
    adjustment_covariates = ifelse(length(adj_vars_for_model) > 0, paste(adj_vars_for_model, collapse = ";"), ""),
    or_adj = adj_or,
    or_adj_ci_lower = adj_or_ci_lower,
    or_adj_ci_upper = adj_or_ci_upper,
    p_unadj_adj_model = adj_p_unadj,
    ratio = 4,
    match_type = "stroke_1_4",
    stringsAsFactors = FALSE
  )
  
  stroke_results_list[[outcome]] = results
  
  cat("OR:", round(or, 3), 
      "95% CI: [", round(or_ci_lower, 3), ",", 
      round(or_ci_upper, 3), "]\n")
  cat("p-value (unadj):", format(p_unadj, scientific = TRUE, digits = 3), "\n")
  cat("Marginal probabilities - Exposed:", round(prob_1, 3), 
      "Unexposed:", round(prob_0, 3), "\n")
  cat("Risk difference:", round(rd, 4), 
      "95% CI: [", round(rd_ci_lower, 4), ",", round(rd_ci_upper, 4), "]\n")
  if (length(adj_vars_for_model) > 0) {
    cat("Adjusted robustness model OR:", round(adj_or, 3),
        "95% CI: [", round(adj_or_ci_lower, 3), ",", round(adj_or_ci_upper, 3), "]\n")
  }
}

# Combine stroke results
if (length(stroke_results_list) > 0) {
  stroke_results = bind_rows(stroke_results_list)
} else {
  stroke_results = data.frame()
}

####################
### TBI ANALYSIS ###
####################

cat("\n=== TBI ANALYSIS ===\n")

tbi_results_list = list()

for (outcome in primary_outcomes) {
  cat("\n--- Analyzing outcome:", outcome, "---\n")
  
  # Filter to complete cases for this outcome
  tbi_complete = tbi_data |>
    filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[tbi_exposure]]),
      !is.na(w_analysis), is.finite(w_analysis), w_analysis > 0,
      !is.na(.data[[design_vars[1]]]), !is.na(.data[[design_vars[2]]])
    )
  
  cat("Complete cases:", nrow(tbi_complete), "\n")
  
  # Check if we have both exposure groups
  exposure_counts = table(tbi_complete[[tbi_exposure]], useNA = "ifany")
  cat("Exposure distribution:\n")
  print(exposure_counts)
  
  if (nrow(tbi_complete) == 0 || length(exposure_counts) < 2 || any(exposure_counts < 1)) {
    cat("Insufficient data - skipping\n")
    next
  }
  
  # Create survey design object
  id_formula = as.formula(paste("~", design_vars[1]))
  strata_formula = as.formula(paste("~", design_vars[2]))
  weights_formula = as.formula("~ w_analysis")
  
  svy_design = svydesign(
    id = id_formula,        # PSU
    strata = strata_formula, # Strata
    weights = weights_formula, # Analysis weights
    data = tbi_complete,
    nest = TRUE
  )
  
  # Create formula: outcome ~ exposure
  formula_str = paste(outcome, "~", tbi_exposure)
  formula_obj = as.formula(formula_str)
  
  # Fit survey-weighted logistic regression
  svy_model = svyglm(
    formula = formula_obj,
    design = svy_design,
    family = quasibinomial(link = "logit")
  )
  
  # Extract coefficients
  coef_summary = summary(svy_model)
  coef_table = coef_summary$coefficients
  
  # Get exposure coefficient
  exposure_coef = coef_table[tbi_exposure, ]
  beta1 = exposure_coef["Estimate"]
  se_beta1 = exposure_coef["Std. Error"]
  p_unadj = exposure_coef["Pr(>|t|)"]

  # Protocol robustness check: if post-match imbalance remains (|SMD| > 0.10),
  # fit an adjusted outcome model including imbalanced prespecified covariates.
  adj_or = NA_real_
  adj_or_ci_lower = NA_real_
  adj_or_ci_upper = NA_real_
  adj_p_unadj = NA_real_
  adj_vars_for_model = character(0)
  tbi_adjustment_vars_present = tbi_adjustment_vars[tbi_adjustment_vars %in% names(tbi_complete)]
  if (length(tbi_adjustment_vars_present) > 0) {
    tbi_adjustment_vars_present = tbi_adjustment_vars_present[
      map_lgl(tbi_adjustment_vars_present, ~ any(!is.na(tbi_complete[[.x]])))
    ]
    if (length(tbi_adjustment_vars_present) > 0) {
      adjusted_formula = as.formula(
        paste(outcome, "~", tbi_exposure, "+", paste(tbi_adjustment_vars_present, collapse = " + "))
      )
      adjusted_model = tryCatch(
        svyglm(formula = adjusted_formula, design = svy_design, family = quasibinomial(link = "logit")),
        error = function(e) NULL
      )
      if (!is.null(adjusted_model)) {
        adjusted_coef = summary(adjusted_model)$coefficients
        if (tbi_exposure %in% rownames(adjusted_coef)) {
          beta1_adj = adjusted_coef[tbi_exposure, "Estimate"]
          se_beta1_adj = adjusted_coef[tbi_exposure, "Std. Error"]
          adj_or = exp(beta1_adj)
          adj_or_ci_lower = exp(beta1_adj - 1.96 * se_beta1_adj)
          adj_or_ci_upper = exp(beta1_adj + 1.96 * se_beta1_adj)
          adj_p_unadj = adjusted_coef[tbi_exposure, "Pr(>|t|)"]
          adj_vars_for_model = tbi_adjustment_vars_present
        }
      }
    }
  }
  
  # Compute OR and 95% CI
  log_or = beta1
  log_or_se = se_beta1
  log_or_ci_lower = log_or - 1.96 * log_or_se
  log_or_ci_upper = log_or + 1.96 * log_or_se
  
  or = exp(log_or)
  or_ci_lower = exp(log_or_ci_lower)
  or_ci_upper = exp(log_or_ci_upper)
  
  # Compute marginal predicted probabilities
  coef_names = names(coef(svy_model))
  beta0 = coef(svy_model)[1]  # intercept
  beta1_idx = which(coef_names == tbi_exposure)
  beta1_coef = coef(svy_model)[beta1_idx]
  
  prob_0 = plogis(beta0)
  prob_1 = plogis(beta0 + beta1_coef)
  
  # SEs using delta method
  vcov_model = vcov(svy_model)
  se_beta0 = sqrt(vcov_model[1, 1])
  cov_beta0_beta1 = vcov_model[1, beta1_idx]
  
  prob_se_0 = prob_0 * (1 - prob_0) * se_beta0
  
  var_sum = vcov_model[1, 1] + vcov_model[beta1_idx, beta1_idx] + 2 * cov_beta0_beta1
  se_sum = sqrt(var_sum)
  prob_se_1 = prob_1 * (1 - prob_1) * se_sum
  
  # Risk difference and CI
  rd = prob_1 - prob_0
  rd_se = sqrt(prob_se_0^2 + prob_se_1^2)
  rd_ci_lower = rd - 1.96 * rd_se
  rd_ci_upper = rd + 1.96 * rd_se
  
  # Store results
  results = data.frame(
    outcome = outcome,
    exposure = tbi_exposure,
    n = nrow(tbi_complete),
    beta1 = beta1,
    se_beta1 = se_beta1,
    or = or,
    or_ci_lower = or_ci_lower,
    or_ci_upper = or_ci_upper,
    p_unadj = p_unadj,
    prob_exposed = prob_1,
    prob_unexposed = prob_0,
    rd = rd,
    rd_ci_lower = rd_ci_lower,
    rd_ci_upper = rd_ci_upper,
    adjusted_model_run = length(adj_vars_for_model) > 0,
    n_adjustment_covariates = length(adj_vars_for_model),
    adjustment_covariates = ifelse(length(adj_vars_for_model) > 0, paste(adj_vars_for_model, collapse = ";"), ""),
    or_adj = adj_or,
    or_adj_ci_lower = adj_or_ci_lower,
    or_adj_ci_upper = adj_or_ci_upper,
    p_unadj_adj_model = adj_p_unadj,
    ratio = 4,
    match_type = "tbi_1_4_revised",
    stringsAsFactors = FALSE
  )
  
  tbi_results_list[[outcome]] = results
  
  cat("OR:", round(or, 3), 
      "95% CI: [", round(or_ci_lower, 3), ",", 
      round(or_ci_upper, 3), "]\n")
  cat("p-value (unadj):", format(p_unadj, scientific = TRUE, digits = 3), "\n")
  cat("Marginal probabilities - Exposed:", round(prob_1, 3), 
      "Unexposed:", round(prob_0, 3), "\n")
  cat("Risk difference:", round(rd, 4), 
      "95% CI: [", round(rd_ci_lower, 4), ",", round(rd_ci_upper, 4), "]\n")
  if (length(adj_vars_for_model) > 0) {
    cat("Adjusted robustness model OR:", round(adj_or, 3),
        "95% CI: [", round(adj_or_ci_lower, 3), ",", round(adj_or_ci_upper, 3), "]\n")
  }
}

# Combine TBI results
if (length(tbi_results_list) > 0) {
  tbi_results = bind_rows(tbi_results_list)
} else {
  tbi_results = data.frame()
}

#############################
### BONFERRONI ADJUSTMENT ###
#############################

cat("\n=== APPLYING BONFERRONI ADJUSTMENT ===\n")

# Apply Bonferroni separately for each exposure and outcome family
# Protocol: "Adjustments will be applied separately for each exposure (stroke; TBI) 
# and separately for the set of primary outcomes versus the set of secondary outcomes"

# Stroke primary outcomes
if (nrow(stroke_results) > 0) {
  stroke_results = stroke_results |>
    mutate(
      # Count number of tests in this family (primary outcomes)
      m_tests = n(),
      p_bonf = pmin(m_tests * p_unadj, 1)
    )
}

# TBI primary outcomes
if (nrow(tbi_results) > 0) {
  tbi_results = tbi_results |>
    mutate(
      # Count number of tests in this family (primary outcomes)
      m_tests = n(),
      p_bonf = pmin(m_tests * p_unadj, 1)
    )
}

#####################
### SAVE RESULTS ###
#####################

# Create output directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Save protocol diagnostics: residual covariate imbalance and weight summaries
if (nrow(balance_smd_table) > 0) {
  write_csv(balance_smd_table, "results/postmatch_balance_smd.csv")
  cat("Post-match SMD diagnostics saved to: results/postmatch_balance_smd.csv\n")
}

weight_diagnostics = bind_rows(
  summarize_weight_diagnostics(stroke_data, "stroke_1_4"),
  summarize_weight_diagnostics(tbi_data, "tbi_1_4_revised")
)
write_csv(weight_diagnostics, "results/weight_diagnostics.csv")
cat("Weight diagnostics saved to: results/weight_diagnostics.csv\n")

# Save results
if (nrow(stroke_results) > 0) {
  write_csv(stroke_results, "results/stroke_primary_results.csv")
  cat("\nStroke results saved to: results/stroke_primary_results.csv\n")
}

if (nrow(tbi_results) > 0) {
  write_csv(tbi_results, "results/tbi_primary_results.csv")
  cat("TBI results saved to: results/tbi_primary_results.csv\n")
}

# Combine all results for summary
if (nrow(stroke_results) > 0 || nrow(tbi_results) > 0) {
  all_results = bind_rows(
    stroke_results |> mutate(exposure_type = "stroke"),
    tbi_results |> mutate(exposure_type = "tbi")
  )
  write_csv(all_results, "results/results.csv")
  cat("All results saved to: results/results.csv\n")

  primary_full_table = all_results |>
    mutate(
      outcome_label = case_when(
        outcome == "usual_place" ~ "Usual place for healthcare",
        outcome == "any_insurance" ~ "Any health insurance coverage",
        TRUE ~ outcome
      ),
      exposure_label = case_when(
        exposure_type == "stroke" ~ "Stroke vs No Stroke",
        exposure_type == "tbi" ~ "TBI vs No TBI",
        TRUE ~ exposure_type
      )
    ) |>
    select(
      exposure_type, exposure_label, outcome, outcome_label, n,
      or, or_ci_lower, or_ci_upper, p_unadj, p_bonf,
      prob_exposed, prob_unexposed, rd, rd_ci_lower, rd_ci_upper,
      adjusted_model_run, n_adjustment_covariates, adjustment_covariates,
      or_adj, or_adj_ci_lower, or_adj_ci_upper, p_unadj_adj_model,
      ratio, match_type
    )

  write_csv(primary_full_table, "results/primary_full_results_table.csv")
  cat("Primary full table saved to: results/primary_full_results_table.csv\n")
}

#####################
### CREATE TABLES ###
#####################

cat("\n=== CREATING FORMATTED TABLES ===\n")

# Create formatted table for stroke results
if (nrow(stroke_results) > 0) {
  stroke_table = stroke_results |>
    mutate(
      Outcome = case_when(
        outcome == "usual_place" ~ "Usual place for healthcare",
        outcome == "any_insurance" ~ "Any health insurance coverage",
        TRUE ~ outcome
      ),
      `OR (95% CI)` = paste0(
        sprintf("%.2f", or), 
        " (", sprintf("%.2f", or_ci_lower), 
        ", ", sprintf("%.2f", or_ci_upper), ")"
      ),
      `P-value (unadj)` = ifelse(p_unadj < 0.001, "<0.001", sprintf("%.3f", p_unadj)),
      `P-value (Bonf adj)` = ifelse(p_bonf < 0.001, "<0.001", sprintf("%.3f", p_bonf)),
      `Risk difference (95% CI)` = paste0(
        sprintf("%.3f", rd), 
        " (", sprintf("%.3f", rd_ci_lower), 
        ", ", sprintf("%.3f", rd_ci_upper), ")"
      )
    ) |>
    select(Outcome, `OR (95% CI)`, `P-value (unadj)`, `P-value (Bonf adj)`, `Risk difference (95% CI)`, n)
  
  # Save as HTML table
  stroke_table_html = stroke_table |>
    kbl(caption = "Primary Analysis Results: Stroke vs. No Stroke (1:4 Matching)", 
        align = c("l", "c", "c", "c", "c", "c"),
        row.names = FALSE) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE,
                  position = "left") |>
    footnote(general = "OR: Odds ratio; CI: Confidence interval; Bonf adj: Bonferroni-adjusted")
  
  save_kable(stroke_table_html, file = "results/stroke_results_table.html")
  cat("Stroke results table saved to: results/stroke_results_table.html\n")
}

# Create formatted table for TBI results
if (nrow(tbi_results) > 0) {
  tbi_table = tbi_results |>
    mutate(
      Outcome = case_when(
        outcome == "usual_place" ~ "Usual place for healthcare",
        outcome == "any_insurance" ~ "Any health insurance coverage",
        TRUE ~ outcome
      ),
      `OR (95% CI)` = paste0(
        sprintf("%.2f", or), 
        " (", sprintf("%.2f", or_ci_lower), 
        ", ", sprintf("%.2f", or_ci_upper), ")"
      ),
      `P-value (unadj)` = ifelse(p_unadj < 0.001, "<0.001", sprintf("%.3f", p_unadj)),
      `P-value (Bonf adj)` = ifelse(p_bonf < 0.001, "<0.001", sprintf("%.3f", p_bonf)),
      `Risk difference (95% CI)` = paste0(
        sprintf("%.3f", rd), 
        " (", sprintf("%.3f", rd_ci_lower), 
        ", ", sprintf("%.3f", rd_ci_upper), ")"
      )
    ) |>
    select(Outcome, `OR (95% CI)`, `P-value (unadj)`, `P-value (Bonf adj)`, `Risk difference (95% CI)`, n)
  
  # Save as HTML table
  tbi_table_html = tbi_table |>
    kbl(caption = "Primary Analysis Results: TBI vs. No TBI (Revised 1:4 Matching)", 
        align = c("l", "c", "c", "c", "c", "c"),
        row.names = FALSE) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE,
                  position = "left") |>
    footnote(general = "OR: Odds ratio; CI: Confidence interval; Bonf adj: Bonferroni-adjusted")
  
  save_kable(tbi_table_html, file = "results/tbi_results_table.html")
  cat("TBI results table saved to: results/tbi_results_table.html\n")
}

####################
### CREATE PLOTS ###
####################

cat("\n=== CREATING PLOTS ===\n")

# Forest plot for stroke results
if (nrow(stroke_results) > 0) {
  stroke_plot_data = stroke_results |>
    mutate(
      Outcome = case_when(
        outcome == "usual_place" ~ "Usual place\nfor healthcare",
        outcome == "any_insurance" ~ "Any health\ninsurance",
        TRUE ~ outcome
      ),
      Outcome = factor(Outcome, levels = rev(c("Usual place\nfor healthcare", "Any health\ninsurance")))
    )
  
  stroke_forest = ggplot(stroke_plot_data, aes(x = or, y = Outcome)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = or_ci_lower, xmax = or_ci_upper, color = p_bonf < 0.05), 
                   height = 0.2, linewidth = 1.5) +
    geom_point(aes(fill = p_bonf < 0.05), size = 4.5, shape = 21, color = "white", stroke = 1.2) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
                  labels = c("0.5", "0.75", "1", "1.5", "2", "3", "4"),
                  limits = c(0.4, 4.5)) +
    scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2E86AB"),
                       labels = c("FALSE" = "Not significant", "TRUE" = "Significant (p < 0.05)"),
                       name = "Bonferroni-adjusted") +
    scale_fill_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2E86AB"),
                      guide = "none") +
    labs(
      title = "Stroke vs. No Stroke: Healthcare Access Outcomes",
      subtitle = "Odds ratios with 95% confidence intervals (1:4 matching)",
      x = "Odds Ratio (log scale)",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  
  # Create results directory if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
  }
  
  ggsave("results/stroke_forest_plot.png", stroke_forest, 
         width = 8, height = 4, dpi = 300)
  cat("Stroke forest plot saved to: results/stroke_forest_plot.png\n")
}

# Forest plot for TBI results
if (nrow(tbi_results) > 0) {
  tbi_plot_data = tbi_results |>
    mutate(
      Outcome = case_when(
        outcome == "usual_place" ~ "Usual place\nfor healthcare",
        outcome == "any_insurance" ~ "Any health\ninsurance",
        TRUE ~ outcome
      ),
      Outcome = factor(Outcome, levels = rev(c("Usual place\nfor healthcare", "Any health\ninsurance")))
    )
  
  tbi_forest = ggplot(tbi_plot_data, aes(x = or, y = Outcome)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = or_ci_lower, xmax = or_ci_upper, color = p_bonf < 0.05), 
                   height = 0.2, linewidth = 1.5) +
    geom_point(aes(fill = p_bonf < 0.05), size = 4.5, shape = 21, color = "white", stroke = 1.2) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
                  labels = c("0.5", "0.75", "1", "1.25", "1.5", "2"),
                  limits = c(0.6, 2.0)) +
    scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2E86AB"),
                       labels = c("FALSE" = "Not significant", "TRUE" = "Significant (p < 0.05)"),
                       name = "Bonferroni-adjusted") +
    scale_fill_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2E86AB"),
                      guide = "none") +
    labs(
      title = "TBI vs. No TBI: Healthcare Access Outcomes",
      subtitle = "Odds ratios with 95% confidence intervals (revised 1:4 matching)",
      x = "Odds Ratio (log scale)",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave("results/tbi_forest_plot.png", tbi_forest, 
         width = 8, height = 4, dpi = 300)
  cat("TBI forest plot saved to: results/tbi_forest_plot.png\n")
}

# Diagnostic plot: Distribution of analysis weights
if (nrow(stroke_data) > 0) {
  weight_diag_stroke = ggplot(stroke_data, aes(x = w_analysis)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
    scale_x_log10() +
    labs(
      title = "Distribution of Analysis Weights: Stroke Analysis",
      subtitle = "Combined NHANES and matching weights (log scale)",
      x = "Analysis Weight (log scale)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  ggsave("results/stroke_weights.png", weight_diag_stroke, 
         width = 8, height = 5, dpi = 300)
  cat("Stroke weights plot saved to: results/stroke_weights.png\n")
}

if (nrow(tbi_data) > 0) {
  weight_diag_tbi = ggplot(tbi_data, aes(x = w_analysis)) +
    geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7, color = "white") +
    scale_x_log10() +
    labs(
      title = "Distribution of Analysis Weights: TBI Analysis",
      subtitle = "Combined NHANES and matching weights (log scale)",
      x = "Analysis Weight (log scale)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  ggsave("results/tbi_weights.png", weight_diag_tbi, 
         width = 8, height = 5, dpi = 300)
  cat("TBI weights plot saved to: results/tbi_weights.png\n")
}

#######################
### SUMMARY TABLES ###
#######################

cat("\n=== SUMMARY OF RESULTS ===\n")

# Print summary for stroke
if (nrow(stroke_results) > 0) {
  cat("\n--- Stroke Results Summary ---\n")
  print(stroke_results |> 
        select(outcome, or, or_ci_lower, or_ci_upper, p_unadj, p_bonf, or_adj, or_adj_ci_lower, or_adj_ci_upper, p_unadj_adj_model) |>
        mutate(across(where(is.numeric), ~ round(.x, 4))))
}

# Print summary for TBI
if (nrow(tbi_results) > 0) {
  cat("\n--- TBI Results Summary ---\n")
  print(tbi_results |> 
        select(outcome, or, or_ci_lower, or_ci_upper, p_unadj, p_bonf, or_adj, or_adj_ci_lower, or_adj_ci_upper, p_unadj_adj_model) |>
        mutate(across(where(is.numeric), ~ round(.x, 4))))
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Review results/ directory for detailed output files\n")
cat("  - CSV files: stroke_primary_results.csv, tbi_primary_results.csv, results.csv, primary_full_results_table.csv\n")
cat("  - Diagnostics: postmatch_balance_smd.csv, weight_diagnostics.csv\n")
cat("  - HTML tables: stroke_results_table.html, tbi_results_table.html\n")
cat("  - Plots: stroke_forest_plot.png, tbi_forest_plot.png, stroke_weights.png, tbi_weights.png\n")
cat("Review matching/ directory for matching diagnostics (love plots)\n")
