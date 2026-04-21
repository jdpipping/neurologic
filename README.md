# Neurologic Conditions and Healthcare Access

This repository contains joint work with UPenn Professors Dylan Small, Randel Swanson, and Andrea Schneider; Stanford Professor Mike Baiocchi; and Penn student Ashil Srivastava on the effect of neurologic conditions on healthcare access and utilization.


## Data Description

This study utilizes data from the **National Health and Nutrition Examination Survey (NHANES)** from 2011-2012 and 2013-2014 cycles. NHANES is a comprehensive survey that collects health and healthcare data for thousands of participants, providing nationally representative estimates of the non-institutionalized U.S. population.

### Dataset Overview

**Study Population**: Adults aged 20+ (stroke questions) and 40+ (TBI questions) from the general U.S. population
**Sample Design**: Complex, multistage probability sampling design
**Data Collection**: In-person interviews, physical examinations, and laboratory tests
**Years**: Combined 2011-2012 and 2013-2014 cycles (both TBI and stroke data collected)

### Data Categories

#### 1. **Demographics** (`demographics/`)
- Age, sex, race/ethnicity, education, income
- Socioeconomic status indicators
- Sample weights for nationally representative estimates

#### 2. **Health Insurance & Healthcare Utilization** (`hiq-huq/`)
- **HIQ (Health Insurance Questionnaire)**: Insurance coverage status and type
- **HUQ (Health Utilization Questionnaire)**: Healthcare access patterns and utilization

#### 3. **Stroke** (`stroke/`)
- **MCQ (Medical Conditions Questionnaire)**: Self-reported stroke history
- Primary exposure variable: "Have you ever been told you had a stroke?" (Yes/No)

#### 4. **Traumatic Brain Injury** (`tbi/`)
- **CSQ (Consumer Survey Questionnaire)**: TBI history and severity
- Primary exposure variable: "Have you ever had a loss of consciousness because of a head injury?" (Yes/No)

#### 5. **Health Outcomes for Matching**
- **Alcohol** (`alcohol/`): ALQ141Q - Number of days with 5+ drinks in past year
- **Smoking** (`smoking/`): SMQ040 - Current smoking status (every day, some days, not at all)
- **Hypertension** (`hypertension/`): BPQ020 - Ever told you had hypertension
- **Diabetes** (`diabetes/`): DIQ010 - Ever told you have diabetes
- **Stroke History** (`stroke/`): MCQ160F - Ever told you had a stroke (used as confounder for TBI study only)

### Key Variables

**Exposures:**
- TBI: Loss of consciousness due to head injury (Yes/No)
- Stroke: Self-reported stroke diagnosis (Yes/No)

**Primary Outcomes:**
- Routine place for healthcare (Yes/No)
- Healthcare coverage status (Yes/No)

**Secondary Outcomes:**
- Type of insurance (government/private)
- Type of healthcare facility
- Time since last healthcare visit

**Covariates (for matching):**
- Demographics: Age (RIDAGEYR), sex (RIAGENDR), race/ethnicity (RIDRETH3)
- Socioeconomic: Income (INDFMPIR), education level (DMDEDUC2)
- Health Outcomes: Alcohol abuse (alcohol_abuse), smoking status (smoking_status), hypertension (hypertension), diabetes (diabetes)
- Additional (TBI study only): Stroke history (stroke_history)

### Study Limitations

- **Age restrictions**: TBI questions limited to adults 40+, stroke questions to adults 20+
- **Exclusion criteria**: Institutionalized individuals, active duty military, homeless population
- **Self-reported data**: Relies on participant recall for medical history


## Analysis Pipeline

### Step 1: Data Preparation (`scripts/01_data-prep.R`)

**Purpose**: Load and merge NHANES datasets from 2011-2012 and 2013-2014 cycles

**Process**:
1. Load 9 datasets:
   - Demographics (DEMO_G/H.xpt): Age, sex, race/ethnicity, education, income
   - Health Insurance (HIQ_G/H.xpt): Insurance coverage status and type
   - Healthcare Utilization (HUQ_G/H.xpt): Healthcare access patterns
   - Stroke (MCQ_G/H.xpt): Self-reported stroke history
   - TBI (CSQ_G/H.xpt): Traumatic brain injury history
   - Alcohol (ALQ_G/H.xpt): Alcohol consumption patterns
   - Smoking (SMQ_G/H.xpt): Smoking status
   - Hypertension (BPQ_G/H.xpt): Hypertension diagnosis
   - Diabetes (DIQ_G/H.xpt): Diabetes diagnosis
   - Drug Use (DUQ_G/H.xpt): Marijuana ever-use (for negative control diagnostic)

2. Merge survey cycles (2011-2012 + 2013-2014) for each dataset

3. Create master dataset by merging all datasets using `SEQN` (participant ID) and `year`

4. Save outputs:
   - `data/clean/master_data.csv` (19,931 participants, 324 variables)
   - Individual datasets: `demo_data.csv`, `hiq_data.csv`, `huq_data.csv`, `stroke_data.csv`, `tbi_data.csv`

### Step 2: Data Exploration (`scripts/02_data-exploration.R`)

**Purpose**: Create exposure and health outcome variables, handle missingness, prepare data for matching

**Process**:

1. **Define Exposure Variables**:
   - `stroke_exposed`: Binary (1 = stroke, 0 = no stroke) from MCQ160F
     - Sample: 431 exposed, 10,888 controls, 8,612 missing (age < 20)
   - `tbi_exposed`: Binary (1 = TBI, 0 = no TBI) from CSQ240
     - Sample: 948 exposed, 6,451 controls, 12,532 missing (age < 40)

2. **Define Health Outcome Variables** (for matching):
   - `alcohol_abuse`: Binary (1 = any days with 5+ drinks in past year, 0 = none)
     - Handles skip pattern: ALQ151 (ever had 5+ drinks) → ALQ141Q (days in past year)
   - `smoking_status`: Factor ("Every day", "Some days", "Not at all", "Refused", "Missing")
     - Handles skip pattern: SMQ020 (ever smoked) → SMQ040 (current status)
   - `hypertension`: Binary (1 = diagnosed, 0 = not diagnosed) from BPQ020
   - `diabetes`: Binary (1 = diagnosed/borderline, 0 = not diagnosed) from DIQ010
   - `stroke_history`: Binary (1 = diagnosed, 0 = not diagnosed) from MCQ160F
     - Used ONLY for TBI study (stroke is the exposure in stroke study)

3. **Handle Missingness**:
   - Create missingness indicators BEFORE imputation
   - Impute missing values:
     - Numeric variables (RIDAGEYR, INDFMPIR): Impute with mean
     - Categorical variables (RIAGENDR, RIDRETH3, DMDEDUC2, smoking_status): Add "Missing" level
     - Binary variables (alcohol_abuse, hypertension, diabetes, stroke_history): Impute with mode

4. **Generate Descriptive Statistics**:
   - Exposure distributions
   - Demographic characteristics (overall and by exposure status)
   - Treated/control ratios

5. **Save Output**:
   - `data/clean/matching_data.csv` (19,931 participants, 338 variables)

### Step 3: Propensity Score Matching (`scripts/03_matching.R`)

**Purpose**: Match exposed and unexposed participants on demographics and health outcomes

**Process**:

1. **Define Matching Variables**:

   **Stroke Study (18 variables)**:
   - Demographics (5): RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, DMDEDUC2
   - NHANES Design Variables (2): SDMVPSU, SDMVSTRA (strata and PSU, included as predictors)
   - NHANES Sample Weight (1): WTINT2YR (interview weights, included as predictor following Dugoff et al.)
   - Health Outcomes (4): alcohol_abuse, smoking_status, hypertension, diabetes
   - Missingness Indicators (6): INDFMPIR_missing, DMDEDUC2_missing, alcohol_abuse_missing, smoking_status_missing, hypertension_missing, diabetes_missing
   - **Excludes**: stroke_history (stroke IS the exposure)

   **TBI Study (20 variables)**:
   - Demographics (5): RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, DMDEDUC2
   - NHANES Design Variables (2): SDMVPSU, SDMVSTRA (strata and PSU, included as predictors)
   - NHANES Sample Weight (1): WTINT2YR (interview weights, included as predictor following Dugoff et al.)
   - Health Outcomes (4): alcohol_abuse, smoking_status, hypertension, diabetes
   - Additional (1): stroke_history (included as confounder)
   - Missingness Indicators (7): All of the above + stroke_history_missing
   - **Includes**: stroke_history (as a confounder/health outcome)

   **Note**: Design variables (`SDMVPSU`, `SDMVSTRA`) and sample weights (`WTINT2YR`) are included as predictors in the unweighted logistic regression for propensity score estimation (Dugoff et al., 2014). This approach enhances external validity by accounting for NHANES complex survey design.

2. **Perform Matching**:
   - Method: Optimal matching (greedy algorithm)
   - Distance: Propensity score (logistic regression with SDMVPSU, SDMVSTRA, and WTINT2YR included as predictors)
   - Caliper: 0.2 SD of logit(propensity score) - manually implemented
     - Propensity score model fitted using logistic regression (includes WTINT2YR)
     - Caliper computed as 0.2 × pooled SD of logit(PS)
     - Custom distance matrix created: absolute difference in logit(PS)
     - Penalty (1000) added to distances > caliper (prevents matching)
     - Custom distance matrix used in optimal matching
   - Replacement: No replacement
   - Ratios tested: 1:1, 1:2, 1:3, 1:4, 1:5, 1:6
   - Stroke: 431 treated, 10,888 controls (25.3:1 ratio)
   - TBI: 948 treated, 6,451 controls (6.8:1 ratio)
   
   **Note**: A caliper of 0.2 SD is manually implemented (MatchIt's `caliper` argument doesn't work with `method = "optimal"`). The caliper restricts matches to within 0.2 standard deviations of the logit of the propensity score by adding a large penalty to distances exceeding the caliper. This improves balance but may discard some treated units that don't have controls within the caliper. This is the standard recommendation for propensity score matching (Rosenbaum & Rubin, 1985).

3. **Balance Diagnostics**:
   - Calculate Standardized Mean Differences (SMD) for all matching variables
   - Generate Love plots (visualize balance before and after matching)
   - Assess balance using thresholds:
     - SMD < 0.05: Excellent balance
     - SMD < 0.10: Good balance
     - SMD > 0.20: Poor balance

4. **Save Outputs**:
   - Matched datasets: `data/matched/stroke_1_1.csv` through `stroke_1_6.csv` (6 files)
   - Matched datasets: `data/matched/tbi_1_1.csv` through `tbi_1_6.csv` (6 files)
   - Love plots: `matching/love_stroke_1_1.png` through `love_stroke_1_6.png` (6 files)
   - Love plots: `matching/love_tbi_1_1.png` through `love_tbi_1_6.png` (6 files)
   - Love plot (revised TBI): `matching/love_tbi_1_4_revised.png` (TBI 1:4 with marijuana in PS, from Step 4)

### Step 4: Negative Control Balance Check (`scripts/04_negative-control.R`)

**Purpose**: Post-matching diagnostic for residual confounding using marijuana ever-use as a falsification variable (not in PS model)

**Process**:
1. **Negative control variable**: Marijuana ever-use (DUQ200: "Ever used marijuana or hashish?") — Z=1 if Yes, Z=0 if No; Refused/Don't know/not administered = missing
2. **Data**: DUQ from `data/nhanes/drug-use/`, merged via pipeline; assessed within matched samples (stroke_1_4, tbi_1_4) among respondents with observed Z
3. **Diagnostics** (per protocol):
   - Weighted SMD for Z between exposed and unexposed (using same analysis weights as primary: w_nhanes × matching weight)
   - Survey-weighted logistic regression: marijuana_ever ~ exposure; report OR (95% CI)
4. **Interpretation**: Success = SMD < 0.10 and OR close to 1 (95% CI includes 1)
5. **Output**: `results/negative_control_results.csv`
6. **Revised TBI matching** (when TBI shows imbalance): Re-run 1:4 TBI matching with marijuana_ever in the PS model; save `matching/love_tbi_1_4_revised.png`

### Step 5: Outcome Analysis (`scripts/05_analysis.R`)

**Purpose**: Estimate associations between exposure (stroke/TBI) and healthcare access/utilization outcomes using survey-weighted regression models

**Process**:

1. **Load Matched Data**:
   - Load 1:4 matched datasets (final matching ratio)
   - Stroke: `data/matched/stroke_1_4.csv` (2,155 participants)
   - TBI: `data/matched/tbi_1_4_revised.csv` (4,740 participants)

2. **Create Outcome Variables**:
   - `usual_place`: Binary (1 = has usual place for healthcare, 0 = no usual place)
     - Derived from HUQ030 (response = 1 or 3 → 1, response = 2 → 0)
   - `any_insurance`: Binary (1 = any health insurance coverage, 0 = uninsured)
     - Derived from HIQ011 (response = 1 → 1, response = 2 → 0)

3. **Create Analysis Weights**:
   - Rescale NHANES interview weight: `w_nhanes = WTINT2YR / 2` (accounts for pooling 2 cycles)
   - Final analysis weight: `w_analysis = w_nhanes * weights` (product of NHANES weight and matching weight)
   - Following Dugoff et al. (2014) for combining propensity score matching with survey weights

4. **Fit Survey-Weighted Regression Models**:
   - Model: Survey-weighted logistic regression using `svyglm()` with `quasibinomial()` family
   - Formula: `logit(Pr(Y=1)) = β₀ + β₁*A` (exposure-only model)
   - Survey design: Accounts for complex survey design (PSU, strata) using `svydesign()`
   - Variance estimation: Taylor series linearization (design-based)
   - Estimand: Population Average Treatment Effect on the Treated (PATT)

5. **Extract Results**:
   - Odds ratios (OR) and 95% confidence intervals
   - Marginal predicted probabilities under A=0 and A=1
   - Risk differences with 95% confidence intervals
   - Unadjusted p-values

6. **Multiple Testing Adjustment**:
   - Bonferroni adjustment applied separately for each exposure (stroke/TBI)
   - Adjustment applied separately for primary outcomes (2 tests per exposure)
   - Formula: `p_Bonf = min(m × p_unadj, 1)` where m = number of tests
   - Reports both unadjusted and Bonferroni-adjusted p-values

7. **Create Formatted Outputs**:
   - **CSV files**: 
     - `results/stroke_primary_results.csv` - Stroke analysis results
     - `results/tbi_primary_results.csv` - TBI analysis results
     - `results/results.csv` - Combined results (stroke + TBI)
     - `results/negative_control_results.csv` - Negative control balance check
   - **HTML tables** (kableExtra):
     - `results/stroke_results_table.html` - Formatted table with ORs, CIs, p-values
     - `results/tbi_results_table.html` - Formatted table with ORs, CIs, p-values
   - **Forest plots** (`results/`):
     - `stroke_forest_plot.png` - Forest plot with colored points/bars (red = not significant, blue = significant)
     - `tbi_forest_plot.png` - Forest plot with colored points/bars
   - **Diagnostic plots** (`results/`):
     - `stroke_weights.png` - Distribution of analysis weights (log scale)
     - `tbi_weights.png` - Distribution of analysis weights (log scale)

### Analysis Results (1:4 Matching)

**Stroke Analysis** (n = 2,151-2,155):
- **Usual place for healthcare**: OR = 1.70 (95% CI: 0.79, 3.66), p = 0.185 (Bonf adj: 0.369)
  - Marginal probabilities: Exposed = 96.1%, Unexposed = 93.7%
  - Risk difference = 2.6% (95% CI: -0.6%, 5.8%)
- **Any health insurance**: OR = 1.24 (95% CI: 0.84, 1.83), p = 0.291 (Bonf adj: 0.582)
  - Marginal probabilities: Exposed = 91.3%, Unexposed = 89.8%
  - Risk difference = 1.9% (95% CI: -1.6%, 5.3%)

**TBI Analysis** (n = 4,738-4,740):
- **Usual place for healthcare**: OR = 0.92 (95% CI: 0.66, 1.30), p = 0.654 (Bonf adj: 1.000)
  - Marginal probabilities: Exposed = 89.6%, Unexposed = 90.0%
  - Risk difference = -0.7% (95% CI: -3.8%, 2.4%)
- **Any health insurance**: OR = 0.93 (95% CI: 0.74, 1.17), p = 0.547 (Bonf adj: 1.000)
  - Marginal probabilities: Exposed = 85.6%, Unexposed = 86.5%
  - Risk difference = -0.9% (95% CI: -4.4%, 2.7%)

**Key Findings**:
- Stroke: Primary outcomes remain directionally positive but statistically inconclusive after Bonferroni adjustment
- TBI: Primary outcomes remain null in the revised 1:4 matched cohort
- All analyses account for complex survey design and use combined NHANES and matching weights (PATT estimand)

### Matching Results

**Stroke Study**:
- **Caliper**: 0.2652 SD of logit(propensity score)
- **1:1**: Max SMD = 0.221, 1 violation > 0.10 (8 covariates > 0.05)
- **1:2**: Max SMD = 0.115, 1 violation > 0.10 (5 covariates > 0.05)
- **1:3**: Max SMD = 0.086, 0 violations > 0.10 (4 covariates > 0.05)
- **1:4**: Max SMD = 0.057, 0 violations > 0.10 (2 covariates > 0.05)
- **1:5**: Max SMD = 0.048, 0 violations > 0.10 (0 covariates > 0.05)
- **1:6**: Max SMD = 0.062, 0 violations > 0.10 (3 covariates > 0.05)
- **All treated units matched** (no discarded units)

**TBI Study**:
- **Caliper**: 0.1408 SD of logit(propensity score)
- **1:1**: Max SMD = 0.295, 1 violation > 0.10 (4 covariates > 0.05)
- **1:2**: Max SMD = 0.207, 1 violation > 0.10 (1 covariate > 0.05)
- **1:3**: Max SMD = 0.101, 1 violation > 0.10 (2 covariates > 0.05)
- **1:4**: Max SMD = 0.058, 0 violations > 0.10 (5 covariates > 0.05)
- **1:5**: Max SMD = 0.115, 2 violations > 0.10 (7 covariates > 0.05)
- **1:6**: Max SMD = 0.181, 6 violations > 0.10 (15 covariates > 0.05)
- **All treated units matched** (no discarded units)

**Key Findings**:
- Penalized caliper implementation successful: All treated units matched without discarding
- Sample weight (WTINT2YR) included as predictor in PS model (Dugoff et al.)
- Stroke: 1:4 achieves good balance and remains the protocol-specified final ratio
- Original TBI 1:4 also achieves good covariate balance, but fails the negative-control check and is not used for final outcome models
- Revised TBI 1:4 matching is the final analytic cohort for TBI after adding marijuana ever-use to the PS model
- Smaller caliper for TBI (0.1408 vs 0.2652 for stroke) reflects less variability in propensity scores

**Balance summary (1:4 matching)** — run `scripts/balance_summary.R` for max |SMD| and covariates > 0.05:

| Match | Max \|SMD\| | Covariates > 0.05 |
|-------|------------|-------------------|
| Stroke 1:4 | 0.057 | 2 |
| TBI 1:4 | 0.058 | 5 |
| TBI 1:4 (incl. marijuana use) | 0.057 | 2 |

**Negative control**: Stroke passes (SMD = -0.021, OR 0.96, 95% CI 0.50 to 1.82); original TBI 1:4 fails (SMD = 0.141, OR 1.35, 95% CI 1.10 to 1.66). The final TBI outcome models therefore use the revised 1:4 cohort with marijuana ever-use included in the PS model.


## Directory Structure

```
.
├── scripts/
│   ├── 01_data-prep.R          # Data preparation: load and merge NHANES datasets
│   ├── 02_data-exploration.R   # Data exploration: create exposure/outcome variables, handle missingness
│   ├── 03_matching.R           # Propensity score matching: match exposed and unexposed participants
│   ├── 04_negative-control.R   # Negative control balance check + revised TBI matching
│   ├── 05_analysis.R           # Outcome analysis: survey-weighted regression models
│   ├── 06_secondary-analysis.R # Secondary healthcare access/utilization outcomes
│   └── balance_summary.R       # Balance summary for stroke/TBI 1:4 (all 3 versions)
├── data/
│   ├── clean/                  # Cleaned datasets
│   │   ├── master_data.csv     # Master dataset (19,931 participants, 324 variables)
│   │   ├── matching_data.csv   # Matching dataset (19,931 participants, 338 variables)
│   │   └── [individual datasets]
│   ├── matched/                # Matched datasets
│   │   ├── stroke_1_1.csv through stroke_1_6.csv (6 files)
│   │   └── tbi_1_1.csv through tbi_1_6.csv (6 files)
│   └── nhanes/                 # Raw NHANES data
│       ├── demographics/
│       ├── hiq-huq/
│       ├── stroke/
│       ├── tbi/
│       ├── alcohol/
│       ├── smoking/
│       ├── hypertension/
│       ├── diabetes/
│       └── drug-use/
├── matching/                   # Matching diagnostics
│   ├── love_stroke_1_1.png through love_stroke_1_6.png (6 files)
│   ├── love_tbi_1_1.png through love_tbi_1_6.png (6 files)
│   └── love_tbi_1_4_revised.png (TBI 1:4 incl. marijuana use)
├── results/                    # Analysis results
│   ├── results.csv             # Combined results (stroke + TBI)
│   ├── stroke_primary_results.csv
│   ├── tbi_primary_results.csv
│   ├── stroke_results_table.html
│   ├── tbi_results_table.html
│   ├── stroke_forest_plot.png
│   ├── tbi_forest_plot.png
│   ├── stroke_weights.png
│   ├── tbi_weights.png
│   ├── negative_control_results.csv
│   └── balance_summary_1_4.csv
├── protocol/
│   ├── protocol.pdf
│   ├── outline.pdf
│   └── references.bib
└── neurologic.Rproj
```
