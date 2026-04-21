#!/usr/bin/env Rscript
# Summary: max |SMD| and # covariates > 0.05 for stroke 1:4, TBI 1:4, TBI 1:4 revised

library(tidyverse)
library(readr)
library(MatchIt)
library(cobalt)

matching_data = read_csv("data/clean/matching_data.csv", show_col_types = FALSE)
drug_data = read_csv("data/clean/drug_data.csv", show_col_types = FALSE)

# Matching variable definitions
base_matching_vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2")
design_vars = c("SDMVPSU", "SDMVSTRA")
sample_weight_var = "WTINT2YR"
health_outcome_vars = c("alcohol_abuse", "smoking_status", "hypertension", "diabetes")
missingness_vars = c("RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
                     "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing", 
                     "diabetes_missing", "stroke_history_missing")
existing_missingness = missingness_vars[missingness_vars %in% names(matching_data)]

stroke_matching = c(base_matching_vars, design_vars, sample_weight_var, health_outcome_vars,
  existing_missingness[existing_missingness %in% c(paste0(base_matching_vars, "_missing"), paste0(health_outcome_vars, "_missing"))])
tbi_matching = c(base_matching_vars, design_vars, sample_weight_var, health_outcome_vars, "stroke_history",
  existing_missingness[existing_missingness %in% c(paste0(base_matching_vars, "_missing"), paste0(health_outcome_vars, "_missing"), "stroke_history_missing")])

run_match_stroke = function() {
  d = matching_data |> filter(!is.na(stroke_exposed))
  f = as.formula(paste("stroke_exposed ~", paste(stroke_matching, collapse = " + ")))
  m = glm(f, data = d, family = binomial)
  d$logit_ps = predict(m, type = "link")
  tr = d$logit_ps[d$stroke_exposed == 1]
  ct = d$logit_ps[d$stroke_exposed == 0]
  cal = 0.2 * sqrt((sd(tr)^2 + sd(ct)^2) / 2)
  dist = abs(outer(d$logit_ps[d$stroke_exposed == 1], d$logit_ps[d$stroke_exposed == 0], "-"))
  dist[dist > cal] = dist[dist > cal] + 10000
  matchit(f, data = d, method = "optimal", distance = dist, ratio = 4)
}

run_match_tbi = function() {
  d = matching_data |> filter(!is.na(tbi_exposed))
  f = as.formula(paste("tbi_exposed ~", paste(tbi_matching, collapse = " + ")))
  m = glm(f, data = d, family = binomial)
  d$logit_ps = predict(m, type = "link")
  tr = d$logit_ps[d$tbi_exposed == 1]
  ct = d$logit_ps[d$tbi_exposed == 0]
  cal = 0.2 * sqrt((sd(tr)^2 + sd(ct)^2) / 2)
  dist = abs(outer(d$logit_ps[d$tbi_exposed == 1], d$logit_ps[d$tbi_exposed == 0], "-"))
  dist[dist > cal] = dist[dist > cal] + 10000
  matchit(f, data = d, method = "optimal", distance = dist, ratio = 4)
}

run_match_tbi_revised = function() {
  md = matching_data
  if (!"DUQ200" %in% names(md)) {
    md = md |> left_join(drug_data |> select(SEQN, year, DUQ200), by = c("SEQN", "year"))
  }
  d = md |>
    filter(!is.na(tbi_exposed)) |>
    mutate(marijuana_ever = case_when(DUQ200 == 1 ~ 1, DUQ200 == 2 ~ 0, TRUE ~ NA_real_))
  mm = as.numeric(names(which.max(table(d$marijuana_ever, useNA = "no"))))
  if (length(mm) == 0) mm = 0
  d = d |> mutate(
    marijuana_ever_missing = as.integer(is.na(marijuana_ever)),
    marijuana_ever = ifelse(is.na(marijuana_ever), mm, marijuana_ever)
  )
  tbi_revised_vars = c(tbi_matching, "marijuana_ever", "marijuana_ever_missing")
  f = as.formula(paste("tbi_exposed ~", paste(tbi_revised_vars, collapse = " + ")))
  m = glm(f, data = d, family = binomial)
  d$logit_ps = predict(m, type = "link")
  tr = d$logit_ps[d$tbi_exposed == 1]
  ct = d$logit_ps[d$tbi_exposed == 0]
  cal = 0.2 * sqrt((sd(tr)^2 + sd(ct)^2) / 2)
  dist = abs(outer(d$logit_ps[d$tbi_exposed == 1], d$logit_ps[d$tbi_exposed == 0], "-"))
  dist[dist > cal] = dist[dist > cal] + 10000
  matchit(f, data = d, method = "optimal", distance = dist, ratio = 4)
}

extract_bal = function(m, name) {
  b = bal.tab(m, binary = "std")
  df = as.data.frame(b$Balance)
  smd_col = intersect(names(df), c("Diff.Adj", "Std.Diff.Adj", "M.Adj"))[1]
  smd = abs(df[[smd_col]])
  data.frame(
    match = name,
    max_abs_smd = max(smd, na.rm = TRUE),
    n_above_005 = sum(smd > 0.05, na.rm = TRUE)
  )
}

cat("Running stroke 1:4...\n")
m_stroke = run_match_stroke()
cat("Running TBI 1:4...\n")
m_tbi = run_match_tbi()
cat("Running TBI 1:4 revised...\n")
m_tbi_rev = run_match_tbi_revised()

out = bind_rows(
  extract_bal(m_stroke, "Stroke 1:4"),
  extract_bal(m_tbi, "TBI 1:4"),
  extract_bal(m_tbi_rev, "TBI 1:4 (incl. marijuana use)")
)

cat("\n=== BALANCE SUMMARY (1:4 matching) ===\n")
print(out)
cat("\n")
write_csv(out, "results/balance_summary_1_4.csv")
cat("Saved: results/balance_summary_1_4.csv\n")
