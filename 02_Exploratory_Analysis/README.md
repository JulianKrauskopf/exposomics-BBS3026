# Session 2 — Exploratory Analysis of the Study Population

**Goal:** Describe the study population, explore measured environmental exposures and learn to interpret plots and correlations.

This practical follows the exploratory-analysis tasks from the Session 2 lecture. Work in [02_Exploratory_Analysis.R](02_Exploratory_Analysis.R): fill in the missing R commands (`...`) and run the script **one task at a time**.

## Before class

1. Download the latest course ZIP from the [BBS3026 repository](https://github.com/JulianKrauskopf/exposomics-BBS3026) using **Code → Download ZIP** and extract it.
2. Open the **course folder** in RStudio. Keep the folder structure; both CSVs are already provided in `data/`.
3. Run `source("R/setup_packages.R")` in the Console to install **dplyr, ggplot2 and corrplot** if needed.
4. Open `02_Exploratory_Analysis/02_Exploratory_Analysis.R`. Save a personal copy as `02_Exploratory_Analysis_YourName.R` in the same folder.

**Important:** The script reads `data/covar.csv` and `data/expo.csv` from the **course root**. In RStudio, use **Session → Set Working Directory → Choose Directory** and select the extracted `exposomics-BBS3026` folder if `getwd()` shows a different location. Do not open or resave the CSVs in Excel.

## Tasks

| Tasks | What you will learn |
|---|---|
| 1–3 | Packages, CSV import, `dim()`, `str()`, factors |
| 4–6 | `table()`, `filter()`, `group_by()`, `summarise()` |
| 7–10 | Exposure data, histograms, log10 scales, density plots |
| 11–12 | Chemical scatterplots and `corrplot()` |

### Study information

The two files contain **632 participants**. Cohort `E` represents **Italy** and `N` represents **Sweden**. Sex is coded `F` and `M`. Smoking status: **1 current, 2 former, 3 never**. Education and activity are coded from **1 (low) to 4 (high)**.

The measured exposure columns are **DDE, HCB, PCB, cadmium, lead**. Use one as your main compound, and compare a persistent organic pollutant (POP) with a heavy metal.

## Questions to discuss

- How do the participant characteristics differ between sexes and cohorts?
- Are exposure measurements normally distributed? What changes on a log10 scale?
- Which exposures are correlated? Why is correlation not causation?

**Deliverable:** Your completed R script and a short written interpretation of one distribution plot, one chemical comparison and one correlation result.

**Next:** Continue with [Session 3 — Exploratory Analysis II](../03_Exploratory_Analysis_II/README.md) for scatterplots against covariates, linear regression and research-question development.

*Note:* The log10 plots assume positive exposure values. Check the exposure data before applying a logarithm.
