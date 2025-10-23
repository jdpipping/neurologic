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
missingness_vars = c("RIDAGEYR_missing", "RIAGENDR_missing", "RIDRETH3_missing", "INDFMPIR_missing", "DMDEDUC2_missing")
existing_missingness_vars = missingness_vars[missingness_vars %in% names(matching_data)]

# base matching variables
base_matching_vars = c("RIDAGEYR", "RIAGENDR", "RIDRETH3", "INDFMPIR", "DMDEDUC2")

# combine base variables with existing missingness indicators
matching_vars = c(base_matching_vars, existing_missingness_vars)

cat("Matching variables:", paste(matching_vars, collapse = ", "), "\n")

# create matching formula
match_formula_stroke = as.formula(paste("stroke_exposed ~", 
                                        paste(matching_vars, collapse = " + ")))
match_formula_tbi = as.formula(paste("tbi_exposed ~", 
                                      paste(matching_vars, collapse = " + ")))

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

# perform matching for different ratios
stroke_matches = list()

for (ratio in 1:6) {
  cat("\n--- Stroke 1 :", ratio, "matching ---\n")
  
  # perform matching
  match_obj = matchit(match_formula_stroke,
                      data = stroke_data,
                      method = "optimal",
                      ratio = ratio)
  
  # store results
  stroke_matches[[paste0("stroke_1_", ratio)]] = match_obj
  
  # extract matched data
  matched_data = match.data(match_obj)
  
  cat("Matched sample size:", nrow(matched_data), "\n")
  cat("Treated:", sum(matched_data$stroke_exposed == 1), "\n")
  cat("Controls:", sum(matched_data$stroke_exposed == 0), "\n")
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

# perform matching for different ratios
tbi_matches = list()

for (ratio in 1:6) {
  cat("\n--- TBI 1 :", ratio, "matching ---\n")
  
  # perform matching
  match_obj = matchit(match_formula_tbi,
                      data = tbi_data,
                      method = "optimal",
                      ratio = ratio)
  
  # store results
  tbi_matches[[paste0("tbi_1_", ratio)]] = match_obj
  
  # extract matched data
  matched_data = match.data(match_obj)
  
  cat("Matched sample size:", nrow(matched_data), "\n")
  cat("Treated:", sum(matched_data$tbi_exposed == 1), "\n")
  cat("Controls:", sum(matched_data$tbi_exposed == 0), "\n")
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
  
  # 2. love plot
  love_plot = love.plot(match_obj, binary = "std", 
                        title = paste("Balance Plot:", match_name))
  print(love_plot)
  # save love plot to file
  out_plot_path = file.path("plots", paste0("love_", match_name, ".png"))
  ggsave(out_plot_path, love_plot, width = 8, height = 6, dpi = 300)

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