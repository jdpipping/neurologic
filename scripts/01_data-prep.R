################
### PACKAGES ###
################

# install.packages("haven")
library(haven)
# install.packages("tidyverse)
library(tidyverse)


#################
### LOAD DATA ###
#################

# DEMOGRAPHICS #

# load
demo_2011_2012 = read_xpt("data/nhanes/demographics/DEMO_G.xpt.txt")
demo_2013_2014 = read_xpt("data/nhanes/demographics/DEMO_H.xpt.txt")
# preview
head(demo_2011_2012)
head(demo_2013_2014)

# HIQ #

# load
hiq_2011_2012 = read_xpt("data/nhanes/hiq-huq/HIQ_G.xpt.txt")
hiq_2013_2014 = read_xpt("data/nhanes/hiq-huq/HIQ_H.xpt.txt")
# preview
head(hiq_2011_2012)
head(hiq_2013_2014)

# HUQ #

# load
huq_2011_2012 = read_xpt("data/nhanes/hiq-huq/HUQ_G.xpt.txt")
huq_2013_2014 = read_xpt("data/nhanes/hiq-huq/HUQ_H.xpt.txt")
# preview
head(huq_2011_2012)
head(huq_2013_2014)

# STROKE #

# load
stroke_2011_2012 = read_xpt("data/nhanes/stroke/MCQ_G.xpt.txt")
stroke_2013_2014 = read_xpt("data/nhanes/stroke/MCQ_H.xpt.txt")
# preview
head(stroke_2011_2012)
head(stroke_2013_2014)

# TBI #

# load
tbi_2011_2012 = read_xpt("data/nhanes/tbi/CSQ_G.xpt.txt")
tbi_2013_2014 = read_xpt("data/nhanes/tbi/CSQ_H.xpt.txt")
# preview
head(tbi_2011_2012)
head(tbi_2013_2014)

# ALCOHOL #

# load
alcohol_2011_2012 = read_xpt("data/nhanes/alcohol/ALQ_G.xpt.txt")
alcohol_2013_2014 = read_xpt("data/nhanes/alcohol/ALQ_H.xpt.txt")
# preview
head(alcohol_2011_2012)
head(alcohol_2013_2014)

# SMOKING #

# load
smoking_2011_2012 = read_xpt("data/nhanes/smoking/SMQ_G.xpt.txt")
smoking_2013_2014 = read_xpt("data/nhanes/smoking/SMQ_H.xpt.txt")
# preview
head(smoking_2011_2012)
head(smoking_2013_2014)

# HYPERTENSION #

# load
hypertension_2011_2012 = read_xpt("data/nhanes/hypertension/BPQ_G.xpt.txt")
hypertension_2013_2014 = read_xpt("data/nhanes/hypertension/BPQ_H.xpt.txt")
# preview
head(hypertension_2011_2012)
head(hypertension_2013_2014)

# DIABETES #

# load
diabetes_2011_2012 = read_xpt("data/nhanes/diabetes/DIQ_G.xpt.txt")
diabetes_2013_2014 = read_xpt("data/nhanes/diabetes/DIQ_H.xpt.txt")
# preview
head(diabetes_2011_2012)
head(diabetes_2013_2014)

# DRUG USE #

# load
drug_2011_2012 = read_xpt("data/nhanes/drug-use/DUQ_G.xpt.txt")
drug_2013_2014 = read_xpt("data/nhanes/drug-use/DUQ_H.xpt.txt")
# preview
head(drug_2011_2012)
head(drug_2013_2014)


###################
### MERGE YEARS ###
###################

# DEMOGRAPHICS #

# merge
demo_2011_2012 = demo_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
demo_2013_2014 = demo_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
demo_data = bind_rows(demo_2011_2012, demo_2013_2014)
# preview
head(demo_data)

# HIQ #

# merge
hiq_2011_2012 = hiq_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
hiq_2013_2014 = hiq_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
hiq_data = bind_rows(hiq_2011_2012, hiq_2013_2014)
# preview
head(hiq_data)

# HUQ #

# merge
huq_2011_2012 = huq_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
huq_2013_2014 = huq_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
huq_data = bind_rows(huq_2011_2012, huq_2013_2014)
# preview
head(huq_data)

# STROKE #

# merge
stroke_2011_2012 = stroke_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
stroke_2013_2014 = stroke_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
stroke_data = bind_rows(stroke_2011_2012, stroke_2013_2014)
# preview
head(stroke_data)

# TBI #

# merge
tbi_2011_2012 = tbi_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
tbi_2013_2014 = tbi_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
tbi_data = bind_rows(tbi_2011_2012, tbi_2013_2014)
# preview
head(tbi_data)

# ALCOHOL #

# merge
alcohol_2011_2012 = alcohol_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
alcohol_2013_2014 = alcohol_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
alcohol_data = bind_rows(alcohol_2011_2012, alcohol_2013_2014)
# preview
head(alcohol_data)

# SMOKING #

# merge
smoking_2011_2012 = smoking_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
smoking_2013_2014 = smoking_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
smoking_data = bind_rows(smoking_2011_2012, smoking_2013_2014)
# preview
head(smoking_data)

# HYPERTENSION #

# merge
hypertension_2011_2012 = hypertension_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
hypertension_2013_2014 = hypertension_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
hypertension_data = bind_rows(hypertension_2011_2012, hypertension_2013_2014)
# preview
head(hypertension_data)

# DIABETES #

# merge
diabetes_2011_2012 = diabetes_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
diabetes_2013_2014 = diabetes_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
diabetes_data = bind_rows(diabetes_2011_2012, diabetes_2013_2014)
# preview
head(diabetes_data)

# DRUG USE #

# merge
drug_2011_2012 = drug_2011_2012 |>
  mutate(year = as.factor("2011-2012"))
drug_2013_2014 = drug_2013_2014 |>
  mutate(year = as.factor("2013-2014"))
drug_data = bind_rows(drug_2011_2012, drug_2013_2014)
# preview
head(drug_data)


######################
### MASTER DATASET ###
######################

# merge all components
master_data <- demo_data |>
  # health insurance (100% overlap)
  left_join(hiq_data, by = c("SEQN", "year")) |>
  # healthcare utilization (100% overlap)
  left_join(huq_data, by = c("SEQN", "year")) |>
  # stroke data (20+, 96% overlap)
  left_join(stroke_data, by = c("SEQN", "year")) |>
  # tbi data (40+, 37% overlap)
  left_join(tbi_data, by = c("SEQN", "year")) |>
  # health outcomes
  left_join(alcohol_data, by = c("SEQN", "year")) |>
  left_join(smoking_data, by = c("SEQN", "year")) |>
  left_join(hypertension_data, by = c("SEQN", "year")) |>
  left_join(diabetes_data, by = c("SEQN", "year")) |>
  # drug use
  left_join(drug_data, by = c("SEQN", "year"))

#################
### SAVE DATA ###
#################

# save individual datasets to file
write_csv(demo_data, "data/clean/demo_data.csv")
write_csv(hiq_data, "data/clean/hiq_data.csv")
write_csv(huq_data, "data/clean/huq_data.csv")
write_csv(stroke_data, "data/clean/stroke_data.csv")
write_csv(tbi_data, "data/clean/tbi_data.csv")
write_csv(drug_data, "data/clean/drug_data.csv")

# save master dataset to file
write_csv(master_data, "data/clean/master_data.csv")