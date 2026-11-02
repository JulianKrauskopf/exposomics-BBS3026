# Session 3 — Exploratory Analysis II: Regression and Research Question

Continue directly from [Session 2](../02_Exploratory_Analysis/README.md). You do **not** need to finish or rerun Session 2 first: the script reloads both CSV files, restores factors and joins participants by ID.

Open [03_Exploratory_Analysis_II.R](03_Exploratory_Analysis_II.R) and save a personal copy as `03_Exploratory_Analysis_II_YourName.R`. Use the **same course working directory** and **same selected compound** as in Session 2. Follow the `setwd()` examples at the start if needed.

## Learning goals

- Extend exploratory graphs with fitted trends and boxplots.
- Fit and interpret simple linear models with continuous and categorical explanatory variables.
- Compare simple models using model p-values and R-squared.
- Fit a multiple (multivariable) model, consider confounding and compare adjusted results.
- Make an illustrative prediction and discuss extrapolation.
- Fit and compare **multivariate multiple regression** models for two exposures.
- Present and discuss an individual research question.

## Tasks

| Tasks | Topic |
|---|---|
| 1 | Reload Session 2 data (no need to rerun earlier work) |
| 2–3 | Continuous/categorical covariate plots |
| 4–6 | Simple linear regression, reference groups, results table |
| 7–8 | Multiple regression, adjusted results, confounding |
| 9 | Predict exposure for a hypothetical participant |
| 10–11 | Multivariate regression and nested-model comparison |
| 12 | Prepare a 5-minute research-question presentation |

The dataset includes `data/covar.csv` and `data/expo.csv` in the repository. Session 3 uses **dplyr** and **ggplot2**, installed via [SETUP.md](../SETUP.md). No new data files are needed.

**Interpretation:** The script models log10-transformed exposures. Coefficients describe associations on that scale; back-transform predictions with `10^x`. More predictors do not automatically give a more scientifically valid model, and associations are not proof of causal effects.
