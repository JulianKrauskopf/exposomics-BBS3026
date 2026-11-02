# BBS3026 — Session 2: Exploratory Analysis
# Work through the tasks in order. Replace each ... or YOUR_* placeholder.
# Run one task at a time with Ctrl+Enter (Windows) / Cmd+Enter (macOS).
# Save your completed script as 02_Exploratory_Analysis_YourName.R.

# BEFORE YOU START — Set the working directory
# Download the GitHub ZIP, extract it in Downloads, and find the course folder.
# Replace YOUR_USERNAME and the folder name with the actual path on your computer.
# Run ONLY ONE of the setwd() lines below (remove the # from that line).
#
# Windows example:
# setwd("C:/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
#
# macOS example:
# setwd("/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
#
# If you moved the extracted folder elsewhere, use its actual location.
# The working directory should contain the "data" folder, NOT the Session 2 folder.
getwd()
file.exists("data/covar.csv")
file.exists("data/expo.csv")
# Both file.exists() checks should return TRUE before you continue.


# TASK 1 — Load the packages
# dplyr: filter, group_by, summarise; ggplot2: graphs; corrplot: correlations.
library(...)
library(ggplot2)
library(corrplot)

# TASK 2 — Import the participant data
# Open this course folder in RStudio before running the script.
# The first CSV column contains participant IDs.
getwd()
covar <- read.csv("data/covar.csv", row.names = 1)
dim(covar)
str(covar)

# QUESTION: What do the 632 rows represent? Which variables are continuous?

# TASK 3 — Identify categorical variables
# Smoking status: 1 = smoker; 2 = former smoker; 3 = never smoker.
# Education and activity: 1 = low; 4 = high.
# Change all three columns to factors.
covar$smoking_status <- factor(covar$smoking_status)
covar$education <- factor(covar$...)
covar$activity <- factor(covar$...)
str(covar)

# TASK 4 — Count participants using table()
# Sex: F = female, M = male. Cohort: E = Italy, N = Sweden.
nrow(covar)
table(covar$...)
table(covar$cohort)
table(covar$smoking_status)

# QUESTION: How many males, females, participants from Italy and Sweden?

# TASK 5 — Describe age and BMI
# Complete the summary table with number of participants, mean age, mean BMI.
# Hint: summarise(n = n(), mean_age = mean(age, na.rm = TRUE), ...)
covar %>%
  summarise(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            mean_bmi = mean(..., na.rm = TRUE))

# Summarise by sex and by cohort.
covar %>%
  group_by(...) %>%
  summarise(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            mean_bmi = mean(bmi, na.rm = TRUE),
            .groups = "drop")

covar %>%
  group_by(cohort) %>%
  summarise(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            mean_bmi = mean(bmi, na.rm = TRUE),
            .groups = "drop")

# TASK 6 — Filter and count categories
# What is the average age of female participants?
covar %>%
  filter(sex == "...") %>%
  summarise(mean_age = mean(age, na.rm = TRUE))

# Count smokers, former smokers and never smokers by sex.
covar %>%
  group_by(sex, smoking_status) %>%
  summarise(n = n(), .groups = "drop")

# Explore education and activity by cohort. Fill the missing column.
covar %>% count(cohort, education)
covar %>% count(cohort, ...)

# TASK 7 — Import and inspect the exposure data
expo <- read.csv("data/expo.csv", row.names = 1)
dim(expo)
str(expo)
names(expo)

# Exposure columns are DDE, HCB, PCB, cadmium and lead.
# All five are measured chemicals. Choose one to investigate.
my_compound <- "DDE"     # Replace with another column name if you prefer.

# Use [[ ]] to select a column by its name stored in my_compound.
range(expo[[my_compound]], na.rm = TRUE)
mean(expo[[...]], na.rm = TRUE)

# TASK 8 — Visualise the exposure distribution
# Make a histogram. aes() maps a column to the x-axis.
ggplot(expo, aes(x = .data[[my_compound]])) +
  geom_histogram(bins = 30, fill = "steelblue") +
  labs(x = my_compound, y = "Number of participants")

# QUESTION: Is the distribution symmetric or right-skewed?

# Make a histogram with a log10 x-axis. Only positive values can be log10-transformed.
ggplot(expo, aes(x = .data[[my_compound]])) +
  geom_histogram(bins = 30, fill = "steelblue") +
  scale_x_...() +
  labs(x = paste(my_compound, "(log10 scale)"), y = "Number of participants")

# TASK 9 — Join participant characteristics and exposures
# Match the two files by participant ID, never by row position alone.
stopifnot(!anyDuplicated(rownames(covar)),
          !anyDuplicated(rownames(expo)))
common_ids <- intersect(rownames(covar), rownames(expo))
study <- cbind(covar[common_ids, , drop = FALSE],
               expo[common_ids, , drop = FALSE])
dim(study)
stopifnot(nrow(study) == length(common_ids))

# TASK 10 — Compare exposure distributions across groups
# Draw a density plot, filled by cohort. Try sex too.
ggplot(study, aes(x = .data[[my_compound]], fill = ...)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  scale_x_log10() +
  labs(x = paste(my_compound, "(log10 scale)"))

# Repeat with fill = sex. What changes?
# Hint: change the line above, then run it again.

# TASK 11 — Scatterplots between chemicals
# Choose one POP (DDE, HCB, PCB) and one heavy metal (cadmium, lead).
pop <- "PCB"
metal <- "cadmium"

ggplot(study, aes(x = .data[[pop]], y = .data[[metal]])) +
  geom_...() +
  labs(x = pop, y = metal)

# Large values can hide patterns. Try log10 scaling of both axes.
ggplot(study, aes(x = .data[[pop]], y = .data[[metal]])) +
  geom_point() +
  scale_x_log10() +
  scale_y_...() +
  labs(x = pop, y = metal)

# TASK 12 — Explore correlation across the five measured exposures
# A correlation coefficient describes direction and strength of association.
# Correlation does not establish causation.
exposure_names <- c("DDE", "HCB", "PCB", "cadmium", "lead")
expo_log <- log10(expo[exposure_names])

# Fill in the matrix name and the input to cor().
correlation_matrix <- cor(..., use = "pairwise.complete.obs")
corrplot(..., method = "number")

# QUESTION: Which pair shows the strongest positive correlation?
# QUESTION: Are any pairs negatively correlated?

# TASK 13 — Continuous participant characteristics
# Plot age against your compound, separately for the two cohorts.
ggplot(study, aes(x = age, y = .data[[my_compound]])) +
  geom_point(alpha = 0.5) +
  scale_y_log10() +
  facet_wrap(~ ...) +
  labs(y = my_compound)

# Add a fitted straight line to the age plot.
ggplot(study, aes(x = age, y = .data[[my_compound]])) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "...", se = FALSE) +
  scale_y_log10() +
  facet_wrap(~ cohort)

# TRY IT: Repeat the scatterplot for BMI and alcohol consumption.
# What patterns do you see? Are they similar across cohorts?

# TASK 14 — Categorical participant characteristics
# Make a boxplot for smoking status. Try education and activity too.
ggplot(study, aes(x = smoking_status, y = .data[[my_compound]])) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Smoking status", y = my_compound)

# TASK 15 — Your first linear regression
# A simple model describes how exposure changes with age.
# The slope is the estimated change in exposure for one year of age.
# We model log10 exposure here because the distributions are skewed.
study$log10_selected_exposure <- log10(study[[my_compound]])

model_age <- lm(log10_selected_exposure ~ ..., data = study)
summary(model_age)

# QUESTION: What is the age coefficient? Is it positive or negative?
# QUESTION: What proportion of variation does R-squared describe?

# Fit a second model with a categorical explanatory variable.
model_smoking <- lm(log10_selected_exposure ~ ..., data = study)
summary(model_smoking)

# INTERPRETATION
# What is the reference smoking-status category?
# Are the coefficients differences between categories?
# Is there evidence of an association? Association is not causation.

# OPTIONAL CHALLENGE — Fit separate age models for Italy and Sweden.
# Compare the age coefficients and R-squared values.
