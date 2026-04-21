################
### PACKAGES ###
################

# install.packages(c("survey", "tidyverse"))
library(survey)
library(tidyverse)

#################
### LOAD DATA ###
#################

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

stroke_exposure = "stroke_exposed"
tbi_exposure = "tbi_exposed"
design_vars = c("SDMVPSU", "SDMVSTRA")

# Prespecified covariates used for post-match residual imbalance robustness check
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

is_special_missing = function(x) {
  x %in% c(7, 9, 77, 99, 777, 999, 7777, 9999, 77777, 99999)
}

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

derive_secondary_outcomes = function(df) {
  df |>
    mutate(
      insured = case_when(
        HIQ011 == 1 ~ 1,
        HIQ011 == 2 ~ 0,
        TRUE ~ NA_real_
      ),
      # HIQ031* are checkbox-style variables in these files:
      # selected category has its code value; otherwise typically NA.
      private_cov = case_when(
        HIQ031A == 14 ~ 1,
        insured == 1 & is.na(HIQ031A) ~ 0,
        TRUE ~ NA_real_
      ),
      public_cov = case_when(
        HIQ031B == 15 | HIQ031D == 17 | HIQ031E == 18 |
          HIQ031F == 19 | HIQ031G == 20 | HIQ031H == 21 | HIQ031I == 22 ~ 1,
        insured == 1 ~ 0,
        TRUE ~ NA_real_
      ),
      public_coverage_any = case_when(
        insured == 1 & public_cov == 1 ~ 1,
        insured == 1 & public_cov == 0 ~ 0,
        TRUE ~ NA_real_
      ),
      hiq270_clean = case_when(
        HIQ270 == 1 ~ 1,
        HIQ270 == 2 ~ 0,
        is_special_missing(HIQ270) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      rx_coverage = case_when(
        insured == 1 ~ hiq270_clean,
        TRUE ~ NA_real_
      ),
      insurance_gap_past_year = case_when(
        HIQ210 == 1 ~ 1,
        HIQ210 == 2 ~ 0,
        is_special_missing(HIQ210) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      mental_health_visit = case_when(
        HUQ090 == 1 ~ 1,
        HUQ090 == 2 ~ 0,
        is_special_missing(HUQ090) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      any_overnight_hospital = case_when(
        HUQ071 == 1 ~ 1,
        HUQ071 == 2 ~ 0,
        is_special_missing(HUQ071) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      hud080_clean = case_when(
        is.na(HUD080) ~ NA_real_,
        is_special_missing(HUD080) ~ NA_real_,
        TRUE ~ as.numeric(HUD080)
      ),
      overnight_hospital_stays = case_when(
        any_overnight_hospital == 0 ~ 0,
        any_overnight_hospital == 1 & !is.na(hud080_clean) ~ hud080_clean,
        TRUE ~ NA_real_
      ),
      # NHANES renamed this field from HUQ040 to HUQ041 in 2013-2014.
      usual_source_type = dplyr::coalesce(HUQ040, HUQ041),
      source_office_clinic = case_when(
        usual_source_type %in% c(1, 2) ~ 1,
        usual_source_type %in% c(3, 4, 5) ~ 0,
        usual_source_type == 6 ~ NA_real_,
        is_special_missing(usual_source_type) ~ NA_real_,
        TRUE ~ NA_real_
      )
    )
}

compute_binary_absolute_effects = function(model, exposure_name) {
  coef_names = names(coef(model))
  beta0 = coef(model)[1]
  beta1_idx = which(coef_names == exposure_name)
  beta1_coef = coef(model)[beta1_idx]

  prob_0 = plogis(beta0)
  prob_1 = plogis(beta0 + beta1_coef)

  vcov_model = vcov(model)
  se_beta0 = sqrt(vcov_model[1, 1])
  cov_beta0_beta1 = vcov_model[1, beta1_idx]

  prob_se_0 = prob_0 * (1 - prob_0) * se_beta0
  var_sum = vcov_model[1, 1] + vcov_model[beta1_idx, beta1_idx] + 2 * cov_beta0_beta1
  se_sum = sqrt(var_sum)
  prob_se_1 = prob_1 * (1 - prob_1) * se_sum

  diff = prob_1 - prob_0
  diff_se = sqrt(prob_se_0^2 + prob_se_1^2)

  list(
    metric_exposed = prob_1,
    metric_unexposed = prob_0,
    metric_diff = diff,
    metric_diff_ci_lower = diff - 1.96 * diff_se,
    metric_diff_ci_upper = diff + 1.96 * diff_se
  )
}

compute_count_absolute_effects = function(model, exposure_name) {
  coef_names = names(coef(model))
  beta0 = coef(model)[1]
  beta1_idx = which(coef_names == exposure_name)
  beta1_coef = coef(model)[beta1_idx]

  mean_0 = exp(beta0)
  mean_1 = exp(beta0 + beta1_coef)

  vcov_model = vcov(model)
  se_beta0 = sqrt(vcov_model[1, 1])
  cov_beta0_beta1 = vcov_model[1, beta1_idx]

  mean_se_0 = mean_0 * se_beta0
  var_sum = vcov_model[1, 1] + vcov_model[beta1_idx, beta1_idx] + 2 * cov_beta0_beta1
  se_sum = sqrt(var_sum)
  mean_se_1 = mean_1 * se_sum

  diff = mean_1 - mean_0
  diff_se = sqrt(mean_se_0^2 + mean_se_1^2)

  list(
    metric_exposed = mean_1,
    metric_unexposed = mean_0,
    metric_diff = diff,
    metric_diff_ci_lower = diff - 1.96 * diff_se,
    metric_diff_ci_upper = diff + 1.96 * diff_se
  )
}

run_secondary_models = function(df, exposure_var, outcomes_tbl, adjustment_vars, match_type_label) {
  out_rows = list()

  for (i in seq_len(nrow(outcomes_tbl))) {
    outcome_name = outcomes_tbl$outcome[i]
    outcome_type = outcomes_tbl$type[i]
    outcome_label = outcomes_tbl$label[i]

    cat("\n--- Secondary outcome:", outcome_name, "---\n")

    d = df |>
      filter(
        !is.na(.data[[outcome_name]]),
        !is.na(.data[[exposure_var]]),
        !is.na(w_analysis), is.finite(w_analysis), w_analysis > 0,
        !is.na(.data[[design_vars[1]]]), !is.na(.data[[design_vars[2]]])
      )

    if (outcome_type == "count") {
      d = d |>
        filter(.data[[outcome_name]] >= 0, is.finite(.data[[outcome_name]]))
    } else {
      d = d |>
        filter(.data[[outcome_name]] %in% c(0, 1))
    }

    exposure_counts = table(d[[exposure_var]], useNA = "ifany")
    if (nrow(d) == 0 || length(exposure_counts) < 2 || any(exposure_counts < 1)) {
      cat("Insufficient data - skipping\n")
      next
    }

    svy_design = svydesign(
      id = as.formula(paste("~", design_vars[1])),
      strata = as.formula(paste("~", design_vars[2])),
      weights = ~ w_analysis,
      data = d,
      nest = TRUE
    )

    model_formula = as.formula(paste(outcome_name, "~", exposure_var))
    model_family = if (outcome_type == "count") quasipoisson(link = "log") else quasibinomial(link = "logit")

    model = tryCatch(
      svyglm(model_formula, design = svy_design, family = model_family),
      error = function(e) NULL
    )
    if (is.null(model)) {
      cat("Model failed - skipping\n")
      next
    }

    coef_table = summary(model)$coefficients
    if (!(exposure_var %in% rownames(coef_table))) {
      cat("Exposure coefficient missing - skipping\n")
      next
    }

    beta1 = coef_table[exposure_var, "Estimate"]
    se_beta1 = coef_table[exposure_var, "Std. Error"]
    p_unadj = coef_table[exposure_var, "Pr(>|t|)"]
    effect = exp(beta1)
    effect_ci_lower = exp(beta1 - 1.96 * se_beta1)
    effect_ci_upper = exp(beta1 + 1.96 * se_beta1)
    effect_label = if (outcome_type == "count") "RR" else "OR"

    abs_effects = if (outcome_type == "count") {
      compute_count_absolute_effects(model, exposure_var)
    } else {
      compute_binary_absolute_effects(model, exposure_var)
    }

    # Robustness check: adjusted model if residual imbalance remains
    adj_effect = NA_real_
    adj_effect_ci_lower = NA_real_
    adj_effect_ci_upper = NA_real_
    adj_p_unadj = NA_real_
    adj_vars_for_model = character(0)

    adjustment_vars_present = adjustment_vars[adjustment_vars %in% names(d)]
    if (length(adjustment_vars_present) > 0) {
      adjustment_vars_present = adjustment_vars_present[
        map_lgl(adjustment_vars_present, ~ any(!is.na(d[[.x]])))
      ]
      if (length(adjustment_vars_present) > 0) {
        adjusted_formula = as.formula(
          paste(outcome_name, "~", exposure_var, "+", paste(adjustment_vars_present, collapse = " + "))
        )
        adjusted_model = tryCatch(
          svyglm(adjusted_formula, design = svy_design, family = model_family),
          error = function(e) NULL
        )
        if (!is.null(adjusted_model)) {
          adjusted_coef = summary(adjusted_model)$coefficients
          if (exposure_var %in% rownames(adjusted_coef)) {
            beta1_adj = adjusted_coef[exposure_var, "Estimate"]
            se_beta1_adj = adjusted_coef[exposure_var, "Std. Error"]
            adj_effect = exp(beta1_adj)
            adj_effect_ci_lower = exp(beta1_adj - 1.96 * se_beta1_adj)
            adj_effect_ci_upper = exp(beta1_adj + 1.96 * se_beta1_adj)
            adj_p_unadj = adjusted_coef[exposure_var, "Pr(>|t|)"]
            adj_vars_for_model = adjustment_vars_present
          }
        }
      }
    }

    out_rows[[outcome_name]] = tibble(
      outcome = outcome_name,
      outcome_label = outcome_label,
      outcome_type = outcome_type,
      effect_label = effect_label,
      exposure = exposure_var,
      n = nrow(d),
      effect = effect,
      effect_ci_lower = effect_ci_lower,
      effect_ci_upper = effect_ci_upper,
      p_unadj = p_unadj,
      metric_exposed = abs_effects$metric_exposed,
      metric_unexposed = abs_effects$metric_unexposed,
      metric_diff = abs_effects$metric_diff,
      metric_diff_ci_lower = abs_effects$metric_diff_ci_lower,
      metric_diff_ci_upper = abs_effects$metric_diff_ci_upper,
      adjusted_model_run = length(adj_vars_for_model) > 0,
      n_adjustment_covariates = length(adj_vars_for_model),
      adjustment_covariates = ifelse(length(adj_vars_for_model) > 0, paste(adj_vars_for_model, collapse = ";"), ""),
      effect_adj = adj_effect,
      effect_adj_ci_lower = adj_effect_ci_lower,
      effect_adj_ci_upper = adj_effect_ci_upper,
      p_unadj_adj_model = adj_p_unadj,
      ratio = 4,
      match_type = match_type_label
    )

    cat(effect_label, ":", round(effect, 3), "95% CI [", round(effect_ci_lower, 3), ",", round(effect_ci_upper, 3), "]\n")
    cat("p-value (unadj):", format(p_unadj, scientific = TRUE, digits = 3), "\n")
    if (outcome_type == "count") {
      cat("Predicted mean - Exposed:", round(abs_effects$metric_exposed, 3),
          "Unexposed:", round(abs_effects$metric_unexposed, 3), "\n")
    } else {
      cat("Predicted probability - Exposed:", round(abs_effects$metric_exposed, 3),
          "Unexposed:", round(abs_effects$metric_unexposed, 3), "\n")
    }
  }

  if (length(out_rows) == 0) return(tibble())
  bind_rows(out_rows)
}

##############################
### DERIVE ANALYSIS FIELDS ###
##############################

stroke_data = derive_secondary_outcomes(stroke_data)
tbi_data = derive_secondary_outcomes(tbi_data)

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

stroke_data = filter_protocol_eligible(stroke_data, stroke_exposure, "Stroke")
tbi_data = filter_protocol_eligible(tbi_data, tbi_exposure, "TBI")

##############################
### BALANCE / ROBUSTNESS   ###
##############################

stroke_balance = find_imbalanced_covariates(stroke_data, stroke_exposure, stroke_covariates, threshold = 0.10)
tbi_balance = find_imbalanced_covariates(tbi_data, tbi_exposure, tbi_covariates, threshold = 0.10)

stroke_adjustment_vars = stroke_balance$imbalanced
tbi_adjustment_vars = tbi_balance$imbalanced

cat("\n=== POST-MATCH COVARIATE BALANCE CHECK (|SMD| > 0.10) ===\n")
cat("Stroke imbalanced covariates:", ifelse(length(stroke_adjustment_vars) == 0, "None", paste(stroke_adjustment_vars, collapse = ", ")), "\n")
cat("TBI imbalanced covariates:", ifelse(length(tbi_adjustment_vars) == 0, "None", paste(tbi_adjustment_vars, collapse = ", ")), "\n")

##############################
### SECONDARY OUTCOME SET  ###
##############################

secondary_outcomes = tribble(
  ~outcome, ~label, ~type,
  "public_coverage_any", "Any public/government insurance (among insured)", "binary",
  "rx_coverage", "Plan covers prescriptions", "binary",
  "insurance_gap_past_year", "Any time without insurance in past year", "binary",
  "mental_health_visit", "Seen a mental health professional in the past year", "binary",
  "any_overnight_hospital", "Any overnight hospital stay in the past year", "binary",
  "source_office_clinic", "Usual source of care is office/clinic (vs other settings)", "binary",
  "overnight_hospital_stays", "Number of overnight hospital stays in the past year", "count"
)

##############################
### RUN SECONDARY MODELS   ###
##############################

cat("\n=== STROKE SECONDARY ANALYSIS ===\n")
stroke_secondary = run_secondary_models(
  df = stroke_data,
  exposure_var = stroke_exposure,
  outcomes_tbl = secondary_outcomes,
  adjustment_vars = stroke_adjustment_vars,
  match_type_label = "stroke_1_4"
)

cat("\n=== TBI SECONDARY ANALYSIS ===\n")
tbi_secondary = run_secondary_models(
  df = tbi_data,
  exposure_var = tbi_exposure,
  outcomes_tbl = secondary_outcomes,
  adjustment_vars = tbi_adjustment_vars,
  match_type_label = "tbi_1_4_revised"
)

##############################
### BONFERRONI (SECONDARY) ###
##############################

if (nrow(stroke_secondary) > 0) {
  stroke_secondary = stroke_secondary |>
    mutate(
      m_tests = n(),
      p_bonf = pmin(m_tests * p_unadj, 1)
    )
}

if (nrow(tbi_secondary) > 0) {
  tbi_secondary = tbi_secondary |>
    mutate(
      m_tests = n(),
      p_bonf = pmin(m_tests * p_unadj, 1)
    )
}

##############################
### SAVE RESULTS           ###
##############################

if (!dir.exists("results")) {
  dir.create("results")
}

if (nrow(stroke_secondary) > 0) {
  write_csv(stroke_secondary, "results/stroke_secondary_results.csv")
  cat("\nStroke secondary results saved to: results/stroke_secondary_results.csv\n")
}

if (nrow(tbi_secondary) > 0) {
  write_csv(tbi_secondary, "results/tbi_secondary_results.csv")
  cat("TBI secondary results saved to: results/tbi_secondary_results.csv\n")
}

if (nrow(stroke_secondary) > 0 || nrow(tbi_secondary) > 0) {
  all_secondary = bind_rows(
    stroke_secondary |> mutate(exposure_type = "stroke"),
    tbi_secondary |> mutate(exposure_type = "tbi")
  )
  write_csv(all_secondary, "results/secondary_results.csv")
  cat("Combined secondary results saved to: results/secondary_results.csv\n")

  secondary_full_table = all_secondary |>
    mutate(
      exposure_label = case_when(
        exposure_type == "stroke" ~ "Stroke vs No Stroke",
        exposure_type == "tbi" ~ "TBI vs No TBI",
        TRUE ~ exposure_type
      )
    ) |>
    select(
      exposure_type, exposure_label, outcome, outcome_label, outcome_type, effect_label, n,
      effect, effect_ci_lower, effect_ci_upper, p_unadj, p_bonf,
      metric_exposed, metric_unexposed, metric_diff, metric_diff_ci_lower, metric_diff_ci_upper,
      adjusted_model_run, n_adjustment_covariates, adjustment_covariates,
      effect_adj, effect_adj_ci_lower, effect_adj_ci_upper, p_unadj_adj_model,
      ratio, match_type
    )

  write_csv(secondary_full_table, "results/secondary_full_results_table.csv")
  cat("Secondary full table saved to: results/secondary_full_results_table.csv\n")
}

secondary_balance_smd = bind_rows(
  stroke_balance$smd_table |> mutate(dataset = "stroke_1_4"),
  tbi_balance$smd_table |> mutate(dataset = "tbi_1_4_revised")
) |>
  relocate(dataset, covariate, abs_smd, imbalanced)
write_csv(secondary_balance_smd, "results/postmatch_balance_smd_secondary.csv")
cat("Secondary analysis balance diagnostics saved to: results/postmatch_balance_smd_secondary.csv\n")

##############################
### SUMMARY               ###
##############################

cat("\n=== SECONDARY ANALYSIS COMPLETE ===\n")
if (nrow(stroke_secondary) > 0) {
  cat("\n--- Stroke secondary summary ---\n")
  print(stroke_secondary |>
          select(
            outcome, outcome_type, effect_label, effect, effect_ci_lower, effect_ci_upper,
            p_unadj, p_bonf, effect_adj, effect_adj_ci_lower, effect_adj_ci_upper, p_unadj_adj_model
          ) |>
          mutate(across(where(is.numeric), ~ round(.x, 4))))
}
if (nrow(tbi_secondary) > 0) {
  cat("\n--- TBI secondary summary ---\n")
  print(tbi_secondary |>
          select(
            outcome, outcome_type, effect_label, effect, effect_ci_lower, effect_ci_upper,
            p_unadj, p_bonf, effect_adj, effect_adj_ci_lower, effect_adj_ci_upper, p_unadj_adj_model
          ) |>
          mutate(across(where(is.numeric), ~ round(.x, 4))))
}
