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

**Covariates:**
- Demographics: Age, sex, race/ethnicity
- Socioeconomic: Income, education level

### Study Limitations

- **Age restrictions**: TBI questions limited to adults 40+, stroke questions to adults 20+
- **Exclusion criteria**: Institutionalized individuals, active duty military, homeless population
- **Self-reported data**: Relies on participant recall for medical history


## Directory Structure

```
.
├── code/
│   └── read.R
├── data/
│   ├── clean/
│   └── nhanes/
│       ├── demographics/
│       ├── hiq-huq/
│       ├── stroke/
│       └── tbi/
├── protocol/
│   └── outline.pdf
└── neurologic.Rproj
```