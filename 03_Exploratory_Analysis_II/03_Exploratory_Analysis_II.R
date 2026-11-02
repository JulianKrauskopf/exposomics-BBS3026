# BBS3026 — Session 3: Exploratory Analysis II
# Simple, multiple and multivariate linear regression
# Continue from Session 2. Fill in the gaps (...), one task at a time.
# Save a personal copy as 03_Exploratory_Analysis_II_YourName.R.

# GET STARTED — Use the SAME extracted course folder as in Session 2.
# If needed, uncomment ONE path and replace YOUR_USERNAME:
# Windows:
# setwd("C:/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
# macOS:
# setwd("/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
getwd()
file.exists("data/covar.csv")
file.exists("data/expo.csv")
# Both checks should print TRUE.

# TASK 1 — Continue with your Session 2 data
# This section recreates the data automatically, even after restarting RStudio.
library(dplyr)
library(ggplot2)
covar <- read.csv("data/covar.csv", row.names = 1)
expo <- read.csv("data/expo.csv", row.names = 1)
covar$smoking_status <- factor(covar$smoking_status)
covar$education <- factor(covar$education)
covar$activity <- factor(covar$activity)
covar$sex <- factor(covar$sex)
covar$cohort <- factor(covar$cohort)
common_ids <- intersect(rownames(covar), rownames(expo))
stopifnot(!anyDuplicated(rownames(covar)), !anyDuplicated(rownames(expo)))
study <- cbind(covar[common_ids, , drop = FALSE],
               expo[common_ids, , drop = FALSE])
stopifnot(nrow(study) == length(common_ids))
str(study)

# Choose the SAME chemical you explored during Session 2.
my_compound <- "DDE"  # Choose DDE, HCB, PCB, cadmium or lead.
stopifnot(my_compound %in% names(expo))
# Exposure measurements are positive; use log10 to reduce right-skewness.
stopifnot(all(study[[my_compound]] > 0, na.rm = TRUE))
study$log10_selected_exposure <- log10(study[[my_compound]])

# TASK 2 — Look for associations with continuous variables
# Continue the plots from the end of Session 2.
# Add a linear trend line and compare cohorts.
ggplot(study, aes(x = age, y = log10_selected_exposure)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "...", se = FALSE) +
  facet_wrap(~ cohort) +
  labs(x = "Age (years)", y = paste("log10", my_compound))

# Try BMI and alcohol by changing the x variable.
# QUESTION: Do the slopes seem similar in both cohorts?

# TASK 3 — Compare categorical variables
# Try smoking_status, sex, education, activity and cohort.
ggplot(study, aes(x = ..., y = log10_selected_exposure)) +
  geom_boxplot() +
  labs(y = paste("log10", my_compound))

# QUESTION: Which categories appear different?

# TASK 4 — Simple linear regression
# lm(response ~ explanatory_variable, data = ...)
# The slope describes a change in log10 exposure per unit increase in age.
model_age <- lm(log10_selected_exposure ~ ..., data = study)
summary(model_age)
# QUESTION: Report the age coefficient, p-value and R-squared.
# QUESTION: What is the intercept? What does a negative slope mean?

# TASK 5 — Simple models with categorical predictors
# R compares each category with a reference category.
model_smoking <- lm(log10_selected_exposure ~ ..., data = study)
summary(model_smoking)
# QUESTION: Which group is the reference for smoking status?
# QUESTION: What do the other coefficients represent?

# TASK 6 — Compare explanatory variables
# Create one simple model for each variable. Complete the list and function.
predictors <- c("sex", "age", "bmi", "smoking_status",
                "education", "activity", "alcohol", "cohort")
# Each model can have MULTIPLE coefficients (e.g. education categories).
# Record the overall model p-value and R-squared instead of only one category coefficient.
simple_results <- lapply(predictors, function(v) {
  fit <- lm(reformulate(v, response = "log10_selected_exposure"), data = study)
  s <- summary(fit)
  f <- s$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  data.frame(variable = v, model_p_value = unname(p), r_squared = s$r.squared)
})
simple_table <- do.call(rbind, simple_results)
simple_table
# Hint: order(table$column) sorts rows.
simple_table[order(simple_table$..., decreasing = TRUE), ]
# QUESTION: Which explanatory variables have the highest R-squared?
# NOTE: Overall model p-values are not the same as individual coefficients.

# TASK 7 — Multiple linear regression
# One exposure outcome and several explanatory variables.
# Adjust for age, sex, BMI, smoking, education, activity, alcohol and cohort.
model_full <- lm(
  log10_selected_exposure ~ age + sex + bmi + smoking_status +
    education + activity + alcohol + ...,
  data = study
)
summary(model_full)
# QUESTION: Compare adjusted coefficients and R-squared with simple models.
# QUESTION: Which variables remain associated after adjustment?
# QUESTION: Why can confounding change regression coefficients?
# IMPORTANT: Choose covariates based on scientific reasoning, not p-values alone.

# TASK 8 — Compare with a smaller model
# Create a justified reduced model, for example age + sex + cohort.
model_reduced <- lm(
  log10_selected_exposure ~ age + ... + cohort,
  data = study
)
summary(model_reduced)
# QUESTION: Compare R-squared and adjusted R-squared across the models.
# QUESTION: Could excluding a confounder change the estimate?

# TASK 9 — Predict exposure for a hypothetical participant
# Use the reduced model above, and supply its predictor variables.
# Ensure factor levels match the original study data.
my_person <- data.frame(
  age = 25,
  sex = factor("F", levels = levels(study$sex)),
  cohort = factor("N", levels = levels(study$cohort))
)
# Adapt my_person if you changed the reduced-model predictors.
predicted_log10 <- predict(model_reduced, newdata = my_person)
predicted_exposure <- 10^predicted_log10
predicted_exposure
range(study[[my_compound]], na.rm = TRUE)
# QUESTION: Is your prediction inside the observed exposure range?
# QUESTION: Are the person's covariate values inside the study range?
# NOTE: Predictions outside observed ranges are extrapolations.

# TASK 10 — Multivariate multiple regression
# DIFFERENCE: multiple regression = several predictors, ONE response.
# Multivariate regression = several response variables in one model.
# Choose TWO different chemicals; here cadmium and PCB.
# Use log10 because the exposure measurements are right-skewed.
study$log10_cadmium <- log10(study$cadmium)
study$log10_PCB <- log10(study$PCB)
multi_model <- lm(
  cbind(log10_cadmium, log10_PCB) ~ age + sex + bmi +
    smoking_status + education + activity + alcohol + ...,
  data = study
)
summary(multi_model)
# QUESTION: Are the patterns the same for both outcomes?

# TASK 11 — Compare multivariate models
# Use the SAME outcome columns and sample for valid model comparison.
multi_reduced <- lm(
  cbind(log10_cadmium, log10_PCB) ~ age + sex + bmi +
    smoking_status + education + activity + cohort,
  data = study
)
anova(multi_reduced, multi_model)
# QUESTION: What does this comparison test?
# NOTE: This compares nested multivariate models jointly.

# TASK 12 — Develop your research question
# Choose a chemical exposure and a health outcome from your literature work.
# Prepare a 5-minute presentation:
# (1) brief background, (2) dataset and methods,
# (3) findings so far, (4) research question.
# QUESTION: Which covariates could confound your research question?
# QUESTION: Which limitations prevent causal conclusions?
