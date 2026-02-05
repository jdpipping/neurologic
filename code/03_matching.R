################
### PACKAGES ###
################

# install.packages(c("cobalt", "iTOS", "MatchIt", "tableone"))
library(cobalt)
library(iTOS)
library(MatchIt)
library(readr)
library(tableone)
library(tidyverse)

# CHOSEN DESIGN: 1:4 matching for both stroke and TBI (see protocol Section 6 and README).
# Outcome analysis (04_analysis.R) uses: data/matched/stroke_1_4.csv, data/matched/tbi_1_4.csv.
# Love plots for the chosen design: matching/love_stroke_1_4.png, matching/love_tbi_1_4.png.

#################
### LOAD DATA ###
#################

# load master dataset with exposure variables (from 02_explore.R)
matching_data = read_csv("data/clean/matching_data.csv")

# check data structure
glimpse(matching_data)
cat("Total participants:", nrow(matching_data), "\n")

# verify exposure variables exist
cat("TBI exposure distribution:\n")
print(table(matching_data$tbi_exposed, useNA = "ifany"))
cat("\nStroke exposure distribution:\n")
print(table(matching_data$stroke_exposed, useNA = "ifany"))

##########################
### MATCHING VARIABLES ###
##########################

# define matching variables using original variable names and missingness indicators
missingness_vars = c("RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing",
                     "alcohol_abuse_missing", "smoking_status_missing", "hypertension_missing", 
                     "diabetes_missing", "stroke_history_missing")
existing_missingness_vars = missingness_vars[missingness_vars %in% names(matching_data)]

# base matching variables (demographics)
base_matching_vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2")

# NHANES design variables (strata and PSU, included as predictors in PS model)
design_vars = c("SDMVPSU", "SDMVSTRA")

# NHANES sample weight to include in PS model
sample_weight_var = "WTINT2YR"

# health outcome matching variables (for both studies)
health_outcome_vars = c("alcohol_abuse", "smoking_status", "hypertension", "diabetes")

# stroke study: demographics + design vars + sample weight + health outcomes (NO stroke_history since stroke is the exposure)
stroke_matching_vars = c(base_matching_vars, design_vars, sample_weight_var, health_outcome_vars)
stroke_missingness_vars = existing_missingness_vars[existing_missingness_vars %in% 
                                                     c(paste0(base_matching_vars, "_missing"),
                                                       paste0(health_outcome_vars, "_missing"))]
stroke_matching_vars_all = c(stroke_matching_vars, stroke_missingness_vars)

# TBI study: demographics + design vars + sample weight + health outcomes + stroke_history
tbi_matching_vars = c(base_matching_vars, design_vars, sample_weight_var, health_outcome_vars, "stroke_history")
tbi_missingness_vars = existing_missingness_vars[existing_missingness_vars %in% 
                                                 c(paste0(base_matching_vars, "_missing"),
                                                   paste0(health_outcome_vars, "_missing"),
                                                   "stroke_history_missing")]
tbi_matching_vars_all = c(tbi_matching_vars, tbi_missingness_vars)

cat("Stroke matching variables:", paste(stroke_matching_vars_all, collapse = ", "), "\n")
cat("TBI matching variables:", paste(tbi_matching_vars_all, collapse = ", "), "\n")

# create matching formulas
match_formula_stroke = as.formula(paste("stroke_exposed ~", 
                                        paste(stroke_matching_vars_all, collapse = " + ")))
match_formula_tbi = as.formula(paste("tbi_exposed ~", 
                                     paste(tbi_matching_vars_all, collapse = " + ")))

#######################
### STROKE MATCHING ###
#######################

cat("\n=== STROKE vs CONTROLS MATCHING ===\n")

# create stroke dataset
stroke_data = matching_data |>
  filter(!is.na(stroke_exposed))

cat("Stroke dataset:", nrow(stroke_data), "participants\n")
cat("Stroke patients:", sum(stroke_data$stroke_exposed == 1), "\n")
cat("Controls:", sum(stroke_data$stroke_exposed == 0), "\n")

# fit propensity score model for stroke
ps_model_stroke = glm(match_formula_stroke, data = stroke_data, family = binomial)
stroke_data$ps = predict(ps_model_stroke, type = "response")
stroke_data$logit_ps = predict(ps_model_stroke, type = "link")

# compute caliper: 0.2 SD of logit(propensity score)
treated_logit_stroke = stroke_data$logit_ps[stroke_data$stroke_exposed == 1]
control_logit_stroke = stroke_data$logit_ps[stroke_data$stroke_exposed == 0]
st_stroke = sd(treated_logit_stroke, na.rm = TRUE)
sc_stroke = sd(control_logit_stroke, na.rm = TRUE)
spool_stroke = sqrt((st_stroke^2 + sc_stroke^2) / 2)
ps_caliper_stroke = 0.2 * spool_stroke

cat("\nStroke PS Caliper (0.2 SD):", round(ps_caliper_stroke, 4), "\n")

# perform matching for different ratios
stroke_matches = list()

for (ratio in 1:6) {
  cat("\n--- Stroke 1 :", ratio, "matching ---\n")
  
  # create distance matrix: absolute difference in logit(PS)
  treated_idx = which(stroke_data$stroke_exposed == 1)
  control_idx = which(stroke_data$stroke_exposed == 0)
  n_treated = length(treated_idx)
  n_control = length(control_idx)
  
  # compute propensity score distance matrix (vectorized)
  # distance = absolute difference in logit(PS) between treated and control
  treated_logit = stroke_data$logit_ps[treated_idx]
  control_logit = stroke_data$logit_ps[control_idx]
  ps_dist = abs(outer(treated_logit, control_logit, "-"))
  
  # apply caliper penalty: add large value if distance > caliper (prevents matching)
  ps_dist[ps_dist > ps_caliper_stroke] = ps_dist[ps_dist > ps_caliper_stroke] + 1000
  
  # perform optimal matching with custom distance matrix
  match_obj = matchit(match_formula_stroke,
                      data = stroke_data,
                      method = "optimal",
                      distance = ps_dist,
                      ratio = ratio)
  
  # store results
  stroke_matches[[paste0("stroke_1_", ratio)]] = match_obj
  
  # extract matched data
  matched_data = match.data(match_obj)
  
  # check how many treated units were discarded due to caliper
  n_discarded = n_treated - sum(matched_data$stroke_exposed == 1)
  
  cat("Matched sample size:", nrow(matched_data), "\n")
  cat("Treated:", sum(matched_data$stroke_exposed == 1), "\n")
  cat("Controls:", sum(matched_data$stroke_exposed == 0), "\n")
  if (n_discarded > 0) {
    cat("Treated units discarded (outside caliper):", n_discarded, "\n")
  }
}

####################
### TBI MATCHING ###
####################

cat("\n=== TBI vs CONTROLS MATCHING ===\n")

# create TBI dataset
tbi_data = matching_data |>
  filter(!is.na(tbi_exposed))

cat("TBI dataset:", nrow(tbi_data), "participants\n")
cat("TBI patients:", sum(tbi_data$tbi_exposed == 1), "\n")
cat("Controls:", sum(tbi_data$tbi_exposed == 0), "\n")

# fit propensity score model for TBI
ps_model_tbi = glm(match_formula_tbi, data = tbi_data, family = binomial)
tbi_data$ps = predict(ps_model_tbi, type = "response")
tbi_data$logit_ps = predict(ps_model_tbi, type = "link")

# compute caliper: 0.2 SD of logit(propensity score)
treated_logit_tbi = tbi_data$logit_ps[tbi_data$tbi_exposed == 1]
control_logit_tbi = tbi_data$logit_ps[tbi_data$tbi_exposed == 0]
st_tbi = sd(treated_logit_tbi, na.rm = TRUE)
sc_tbi = sd(control_logit_tbi, na.rm = TRUE)
spool_tbi = sqrt((st_tbi^2 + sc_tbi^2) / 2)
ps_caliper_tbi = 0.2 * spool_tbi

cat("\nTBI PS Caliper (0.2 SD):", round(ps_caliper_tbi, 4), "\n")

# perform matching for different ratios
tbi_matches = list()

for (ratio in 1:6) {
  cat("\n--- TBI 1 :", ratio, "matching ---\n")
  
  # create distance matrix: absolute difference in logit(PS)
  treated_idx = which(tbi_data$tbi_exposed == 1)
  control_idx = which(tbi_data$tbi_exposed == 0)
  n_treated = length(treated_idx)
  n_control = length(control_idx)
  
  # compute propensity score distance matrix (vectorized)
  # distance = absolute difference in logit(PS) between treated and control
  treated_logit = tbi_data$logit_ps[treated_idx]
  control_logit = tbi_data$logit_ps[control_idx]
  ps_dist = abs(outer(treated_logit, control_logit, "-"))
  
  # apply caliper penalty: add large value if distance > caliper (prevents matching)
  ps_dist[ps_dist > ps_caliper_tbi] = ps_dist[ps_dist > ps_caliper_tbi] + 1000
  
  # perform optimal matching with custom distance matrix
  match_obj = matchit(match_formula_tbi,
                      data = tbi_data,
                      method = "optimal",
                      distance = ps_dist,
                      ratio = ratio)
  
  # store results
  tbi_matches[[paste0("tbi_1_", ratio)]] = match_obj
  
  # extract matched data
  matched_data = match.data(match_obj)
  
  # check how many treated units were discarded due to caliper
  n_discarded = n_treated - sum(matched_data$tbi_exposed == 1)
  
  cat("Matched sample size:", nrow(matched_data), "\n")
  cat("Treated:", sum(matched_data$tbi_exposed == 1), "\n")
  cat("Controls:", sum(matched_data$tbi_exposed == 0), "\n")
  if (n_discarded > 0) {
    cat("Treated units discarded (outside caliper):", n_discarded, "\n")
  }
}

###########################
### BALANCE DIAGNOSTICS ###
###########################

cat("\n=== BALANCE DIAGNOSTICS ===\n")

# function to run all balance diagnostics
run_balance_diagnostics = function(match_obj, match_name) {
  cat("\n---", match_name, "---\n")
  
  # 1. standardized mean differences (SMD)
  smd_results = bal.tab(match_obj, binary = "std")
  print(smd_results)
  
  # 2. love plot (title e.g. "Balance Plot: Stroke 1:4")
  exposure = sub("_.*", "", match_name)
  ratio = sub(".*_1_", "", match_name)
  display_exposure = if (exposure == "tbi") "TBI" else paste0(toupper(substring(exposure, 1, 1)), substring(exposure, 2))
  plot_title = paste("Balance Plot:", display_exposure, "1:", ratio)
  love_plot = love.plot(match_obj, binary = "std", title = plot_title)
  print(love_plot)
  # save love plot to file (6 in x 4 in)
  if (!dir.exists("matching")) {
    dir.create("matching", recursive = TRUE)
  }
  out_plot_path = file.path("matching", paste0("love_", match_name, ".png"))
  ggsave(out_plot_path, love_plot, width = 6, height = 4, dpi = 300)

  # 3. concise SMD summary across thresholds
  bal_df = as.data.frame(smd_results$Balance)
  # cobalt names SMD-after column variably; try common possibilities
  smd_after_cols = intersect(colnames(bal_df), c("Diff.Adj", "Std.Diff.Adj", "M.Adj"))
  if (length(smd_after_cols) == 0) {
    warning("Could not find SMD-after column in balance table; skipping summary")
    return(invisible(NULL))
  }
  smd_after = abs(bal_df[[smd_after_cols[1]]])
  thr = c(0.05, 0.10, 0.20)
  counts_over = sapply(thr, function(t) sum(smd_after > t, na.rm = TRUE))
  max_smd = suppressWarnings(max(smd_after, na.rm = TRUE))
  pct_over_0_1 = mean(smd_after > 0.10, na.rm = TRUE)

  cat("SMD violations -> >0.05:", counts_over[1], ", >0.10:", counts_over[2], ", >0.20:", counts_over[3], "\n")
  cat("Max SMD:", round(max_smd, 3), "; % covariates > 0.10:", round(100*pct_over_0_1, 1), "%\n")
}

# run diagnostics for all stroke matches
for (i in seq_along(stroke_matches)) {
  run_balance_diagnostics(stroke_matches[[i]], names(stroke_matches)[i])
}

# run diagnostics for all TBI matches  
for (i in seq_along(tbi_matches)) {
  run_balance_diagnostics(tbi_matches[[i]], names(tbi_matches)[i])
}

############################
### SAVE MATCHED SAMPLES ###
############################

# save all stroke matched datasets
for (nm in names(stroke_matches)) {
  md = match.data(stroke_matches[[nm]])
  out_path = file.path("data/matched", paste0(nm, ".csv"))
  write_csv(md, out_path)
}

# save all TBI matched datasets
for (nm in names(tbi_matches)) {
  md = match.data(tbi_matches[[nm]])
  out_path = file.path("data/matched", paste0(nm, ".csv"))
  write_csv(md, out_path)
}

cat("\n=== MATCHING COMPLETE ===\n")
cat("Saved all matched datasets to data/matched/ (stroke_1_R and tbi_1_R)\n")
cat("Review balance diagnostics to select optimal ratios\n")