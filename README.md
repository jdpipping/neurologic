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

#### 5. **Health Outcomes**
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

### Step 1: Data Preparation (`code/01_data-prep.R`)

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

2. Merge survey cycles (2011-2012 + 2013-2014) for each dataset

3. Create master dataset by merging all datasets using `SEQN` (participant ID) and `year`

4. Save outputs:
   - `data/clean/master_data.csv` (19,931 participants, 324 variables)
   - Individual datasets: `demo_data.csv`, `hiq_data.csv`, `huq_data.csv`, `stroke_data.csv`, `tbi_data.csv`

### Step 2: Data Exploration (`code/02_data-exploration.R`)

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

### Step 3: Propensity Score Matching (`code/03_matching.R`)

**Purpose**: Match exposed and unexposed participants on demographics and health outcomes

**Process**:

1. **Define Matching Variables**:

   **Stroke Study (15 variables)**:
   - Demographics (5): RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, DMDEDUC2
   - Health Outcomes (4): alcohol_abuse, smoking_status, hypertension, diabetes
   - Missingness Indicators (6): INDFMPIR_missing, DMDEDUC2_missing, alcohol_abuse_missing, smoking_status_missing, hypertension_missing, diabetes_missing
   - **Excludes**: stroke_history (stroke IS the exposure)

   **TBI Study (17 variables)**:
   - Demographics (5): RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, DMDEDUC2
   - Health Outcomes (4): alcohol_abuse, smoking_status, hypertension, diabetes
   - Additional (1): stroke_history (included as confounder)
   - Missingness Indicators (7): All of the above + stroke_history_missing
   - **Includes**: stroke_history (as a confounder/health outcome)

2. **Perform Matching**:
   - Method: Optimal matching (greedy algorithm)
   - Ratios tested: 1:1, 1:2, 1:3, 1:4, 1:5, 1:6
   - Stroke: 431 treated, 10,888 controls (25.3:1 ratio)
   - TBI: 948 treated, 6,451 controls (6.8:1 ratio)

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
   - Love plots: `plots/love_stroke_1_1.png` through `love_stroke_1_6.png` (6 files)
   - Love plots: `plots/love_tbi_1_1.png` through `love_tbi_1_6.png` (6 files)

### Matching Results

**Stroke Study**:
- **1:1**: Max SMD = 0.091, 0 violations > 0.10 ✓
- **1:2**: Max SMD = 0.049, 0 violations > 0.10 ✓
- **1:3**: Max SMD = 0.062, 0 violations > 0.10 ✓
- **1:4**: Max SMD = 0.047, 0 violations > 0.10 ✓
- **1:5**: Max SMD = 0.068, 0 violations > 0.10 ✓
- **1:6**: Max SMD = 0.138, 1 violation > 0.10 (3.8% covariates > 0.10)

**TBI Study**:
- **1:1**: Max SMD = 0.078, 0 violations > 0.10 ✓
- **1:2**: Max SMD = 0.054, 0 violations > 0.10 ✓
- **1:3**: Max SMD = 0.069, 0 violations > 0.10 ✓
- **1:4**: Max SMD = 0.119, 1 violation > 0.10 (3.6% covariates > 0.10) - Propensity score distance
  - All individual covariates < 0.10 (max = 0.0942 for alcohol_abuse)
- **1:5**: Max SMD = 0.204, 4 violations > 0.10 (14.3% covariates > 0.10)
- **1:6**: Max SMD = 0.309, 6 violations > 0.10 (21.4% covariates > 0.10)


## Directory Structure

```
.
├── code/
│   ├── 01_data-prep.R          # Data preparation: load and merge NHANES datasets
│   ├── 02_data-exploration.R   # Data exploration: create exposure/outcome variables, handle missingness
│   └── 03_matching.R           # Propensity score matching: match exposed and unexposed participants
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
│       └── diabetes/
├── plots/                      # Balance plots (Love plots)
│   ├── love_stroke_1_1.png through love_stroke_1_6.png (6 files)
│   └── love_tbi_1_1.png through love_tbi_1_6.png (6 files)
├── protocol/
│   └── outline.pdf
└── neurologic.Rproj
```