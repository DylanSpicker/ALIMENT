<div align="center">
  
# 🥗 ALIMENT

**A machine-learning driven approach to estimating the distribution of usual dietary intake from repeated 24-hour recalls.**

</div>

---

> **Note:** This repository contains the codebase for both the ALIMENT estimation procedure and the data simulations used to validate the method against existing standard approaches.

## 📖 Background

Estimating the distribution of usual dietary intake for a population is a complex challenge, primarily because we typically only observe short-term, error-prone measurements (like 24-hour dietary recalls). **ALIMENT** provides a novel, machine-learning-based framework to more accurately estimate these underlying usual intake distributions, taking into account covariates and complex measurement errors.

The repository is organized into two main parts:
1. 🧮 **Estimators**: Code that implements the ALIMENT estimation procedure.
2. 🔬 **Simulations**: Code used to generate synthetic datasets for validating and testing our estimators against current state-of-the-art methods.

---

## 🧮 Estimators

> 🚧 *(Code for the estimators will be added here in the future.)*

---

## 🔬 Simulations

The `Simulation/` directory contains scripts that create synthetic dietary datasets used to test the estimators. We generate datasets under both **parametric** (assuming a specific statistical distribution) and **non-parametric** (using real-world data patterns) settings. 

### What the Scripts Do

| Simulation Type | Script Location | Description |
| :--- | :--- | :--- |
| **Parametric** | `Simulation/Parametric/parametric_dietary_data.R` | Simulates datasets using predefined mathematical distributions (like Gamma and Normal distributions) modeled after standard methods (e.g., the NCI method). It creates synthetic intakes for nutrients such as Vitamin A and Calcium. |
| **Non-Parametric (Demographic Free)** | `Simulation/Demographic Free/nhanes_nonparametric_demographic_free.R` | Creates a synthetic dataset by resampling actual National Health and Nutrition Examination Survey (NHANES) dietary data. It generates an unobserved "usual intake" and two days of observed intake for individuals, adding realistic random distortions. |
| **Non-Parametric (With Demographics)** | `Simulation/With Demographic/nhanes_nonparametric_sim.R` | Similar to the demographic-free approach, but it uses demographic characteristics (e.g., age, income, household size) to group similar individuals. Dietary patterns are then borrowed among individuals to create a realistic, covariate-dependent synthetic dataset. |


### 📦 Requirements & Inputs

To run the simulations, you will need **R** installed along with the following packages:
* `tidyverse`
* `haven`
* `labelled`

> [!WARNING]
> **Input Data Note:** The non-parametric simulations rely on raw NHANES and Food Patterns Equivalents Database (FPED) data files. Because these files are large and subject to their own data use agreements or distribution practices, **they are not included in this public repository**. 

If you have acquired these files, they should be placed in the respective directories as follows:

<details>
<summary><b>View Required File Structure</b></summary>
<br>

```text
ALIMENT/
├── NHANES/
│   ├── DEMO_J.xpt (Demographics)
│   ├── DR1TOT_J.xpt (Day 1 Dietary Data)
│   └── DR2TOT_J.xpt (Day 2 Dietary Data)
└── FPED/
    ├── fped_dr1tot_1718.sas7bdat (Day 1 FPED Data)
    └── fped_dr2tot_1718.sas7bdat (Day 2 FPED Data)
```

</details>

<br>

### 🚀 How to Run the Simulations

You can run any of the simulation scripts using an R console or an IDE like RStudio. For example, to run the demographic-free simulation from your terminal:

```bash
Rscript "Simulation/Demographic Free/nhanes_nonparametric_demographic_free.R"
```
> **Tip:** The scripts are designed to run with the root of the repository as your working directory, ensuring the script can find the required `NHANES/` and `FPED/` input folders.

### 📊 Expected Outputs

Running these scripts will generate new `.csv` files containing the simulated populations and their intake records. These synthetic datasets act as the benchmark ground truth that we then try to recover using the ALIMENT estimators.

* 📈 **Parametric Scripts:** Outputs many datasets into a `ParametricCSV/` folder. 
* 📉 **Non-Parametric Scripts:** Produces `simulated_nhanes_complete.csv` or `simulated_nhanes_complete_no_demographic.csv`, along with sets of true population data files in a `TrueDistCSV/` directory.

---
