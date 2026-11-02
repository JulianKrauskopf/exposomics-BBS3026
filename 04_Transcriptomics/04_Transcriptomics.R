# BBS3026 — Session 4: Transcriptomics Analysis I
# Work directly in a personal copy: 04_Transcriptomics_YourName.R
# Complete the ... placeholders and run ONE TASK AT A TIME.
# This session continues your compound and disease from Sessions 1–3.

# BEFORE STARTING — Use the extracted course folder as working directory.
# Uncomment ONE of the examples and replace YOUR_USERNAME:
# Windows:
# setwd("C:/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
# macOS:
# setwd("/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
getwd()
file.exists("data/covar.csv")  # should be TRUE
file.exists("data/expo.csv")   # should be TRUE
file.exists("data/EGM_Transcriptomics_GeneSymbol_Full.Rdata")  # should be TRUE
# Download the large Rdata file separately; see SETUP.md.

# TASK 1 — Load the packages and shared Open Targets functions
library(dplyr)
library(ggplot2)
library(limma)
source("R/project_helpers.R")

# TASK 2 — Continue with the participant and exposure data
covar <- read.csv("data/covar.csv", row.names = 1)
expo  <- read.csv("data/expo.csv", row.names = 1)
covar$sex <- factor(covar$sex)
covar$cohort <- factor(covar$cohort)
covar$education <- factor(covar$education)
covar$activity <- factor(covar$activity)
covar$smoking_status <- factor(covar$smoking_status)

# Select the SAME compound you investigated in Sessions 2 and 3.
my_compound <- "PCB"  # Choose DDE, HCB, PCB, cadmium or lead.
stopifnot(my_compound %in% names(expo))
study <- prepare_exposure(covar, expo, my_compound)
# This helper aligns IDs and creates log2_exposure for regression.
str(study)
nrow(study)

# TASK 3 — Load and explore transcriptomics data
# The Rdata file should contain a numeric expression matrix named 'omics'.
omics <- load_expression()
dim(omics)
str(omics[, seq_len(min(3, ncol(omics))), drop = FALSE])
# QUESTION: Are genes stored in rows or columns?
# QUESTION: How many genes and how many participants?

# TASK 4 — Explore normalisation
# Visualise distributions before and after quantile normalisation.
# The course expression data are already on the intended expression scale:
# do NOT apply log2() again unless the source data documentation requires it.
plotDensities(omics, legend = FALSE)
omics_norm <- normalizeBetweenArrays(..., method = "quantile")
plotDensities(omics_norm, legend = FALSE)
# QUESTION: What changed in the density distributions?

# TASK 5 — Find your disease in Open Targets
# Open Targets replaces the manual OMIM/Excel gene retrieval step.
# Choose a disease linked to your research question.
# Search interactively and copy an appropriate disease_id from the results.
disease_search <- ot_search_diseases("asthma")  # CHANGE to your outcome.
head(disease_search, 10)
my_disease_id <- "YOUR_DISEASE_ID" # e.g. copy an ID from disease_search
# QUESTION: Does the selected disease label match your research question?

# TASK 6 — Define the gene set BEFORE looking at exposure p-values
# Open Targets scores describe disease–target association evidence,
# NOT the effect of your compound on gene expression.
# Choose the cutoff first (0.2 is the course starting point).
min_disease_score <- 0.2
disease_genes <- get_disease_genes(
  disease_id = my_disease_id,
  min_score = min_disease_score
)
head(disease_genes[, c("gene_symbol", "association_score")])
nrow(disease_genes)
# QUESTION: How many genes meet your Open Targets score cutoff?

# TASK 7 — Compare gene-set sizes (EXPLORATORY EXTENSION)
# What happens when you choose a more stringent cutoff (e.g. 0.4)?
# This comparison is exploratory; keep your original cutoff unchanged
# when interpreting your predefined disease-focused model.
genes_strict <- get_disease_genes(
  disease_id = my_disease_id,
  min_score = ...
)
nrow(genes_strict)
# OPTIONAL: Try a different disease and compare its number of genes.
# NOTE: Larger gene sets do not automatically provide stronger evidence.

# TASK 8 — Match the disease genes to measured transcripts
# Not every gene in Open Targets is measured in this dataset.
predefined_genes <- define_disease_gene_set(
  expression = omics_norm,
  genes = ...
)
nrow(predefined_genes)
head(predefined_genes)
# QUESTION: How many Open Targets genes are represented in our matrix?

# TASK 9 — Match transcriptomics samples with the participant data
# Use a common identifier and exactly the SAME sample order.
common_ids <- rownames(study)[rownames(study) %in% colnames(omics_norm)]
study_matched <- study[common_ids, , drop = FALSE]
omics_matched <- omics_norm[, common_ids, drop = FALSE]
identical(rownames(study_matched), colnames(omics_matched))
# The result should be TRUE. If not, stop and check the IDs.
# QUESTION: How many participants have expression AND exposure data?

# TASK 10 — Explore stratification
# Compare all participants with females, males or one cohort.
# Start with all participants; optionally choose "F" or "M".
my_sex <- "all"   # "all", "F" or "M"
if (my_sex == "all") {
  selected_ids <- rownames(study_matched)
} else {
  selected_ids <- rownames(study_matched)[study_matched$sex == my_sex]
}
analysis_data <- droplevels(study_matched[selected_ids, , drop = FALSE])
analysis_expression <- omics_matched[, selected_ids, drop = FALSE]
identical(rownames(analysis_data), colnames(analysis_expression))
nrow(analysis_data)
# QUESTION: Why must we subset both tables using the same IDs?

# TASK 11 — Build a regression design matrix
# Response: expression level of EACH gene.
# Main explanatory variable: log2 exposure of your selected compound.
# Adjust for age, sex and cohort where these variables vary in your group.
# For a sex-stratified analysis, sex is constant and must be omitted.
adjust_for <- c("age", "cohort")
if (my_sex == "all") adjust_for <- c(adjust_for, "sex")
# Use model.matrix to inspect the covariates before fitting limma.
design <- model.matrix(
  reformulate(c("log2_exposure", adjust_for)),
  data = analysis_data
)
head(design)
# QUESTION: Which column tests the exposure–gene-expression association?
# NOTE: Missing data and factor levels are checked in the next task.

# TASK 12 — Fit the FULL transcriptome using limma
# Use ALL measured genes to estimate variance consistently.
# The helper matches complete cases, fits lmFit and runs eBayes.
full_analysis <- run_limma_exposure(
  expression = analysis_expression,
  exposure_data = analysis_data,
  adjust_for = adjust_for,
  normalize = FALSE  # We already normalised above.
)
full_results <- full_analysis$results
head(full_results[, c("feature_id", "logFC", "P.Value", "adj.P.Val")])
nrow(full_results)
# QUESTION: What does logFC represent for a continuous log2 exposure?
# QUESTION: Why should we distinguish P.Value from adj.P.Val?

# TASK 13 — Explore the full-transcriptome results
# Count genes with nominal p < 0.05 and BH-adjusted p < 0.05.
sum(full_results$P.Value < 0.05, na.rm = TRUE)
sum(full_results$... < 0.05, na.rm = TRUE)
hist(full_results$P.Value, breaks = 30,
     main = "Full transcriptome: nominal p-values",
     xlab = "P-value")
plot_volcano(full_results,
             title = paste(my_compound, "— all measured genes"),
             highlight_column = "P.Value")
# QUESTION: Does nominal p < 0.05 guarantee a genuine biological effect?

# TASK 14 — Analyse your PREDEFINED disease gene set
# This uses the same full-transcriptome limma fit, but focuses on the genes
# selected from Open Targets BEFORE examining exposure p-values.
disease_results <- disease_focused_results(
  all_results = full_results,
  predefined_genes = ...
)
head(disease_results[, c("feature_id", "disease_gene_symbol",
                         "logFC", "P.Value", "FDR_disease_set",
                         "adj.P.Val", "association_score")], 10)
sum(disease_results$P.Value < 0.05, na.rm = TRUE)
sum(disease_results$FDR_disease_set < 0.05, na.rm = TRUE)
plot_volcano(disease_results,
             title = paste(my_compound, "— disease-focused gene set"),
             highlight_column = "P.Value")
# QUESTION: How do nominal p, disease-set FDR and whole-transcriptome FDR differ?

# TASK 15 — Export results
# Save both complete and disease-focused results for the next sessions.
dir.create("outputs", showWarnings = FALSE)
write.csv(full_results,
          file = "outputs/session4_full_transcriptome.csv",
          row.names = FALSE)
write.csv(disease_results,
          file = "outputs/session4_disease_focused.csv",
          row.names = FALSE)
# QUESTION: Which result columns would you include in a poster?

# TASK 16 — Discuss your research question
# Report: exposure, disease, Open Targets ID, selected cutoff,
# total disease-associated genes, measured disease-associated genes,
# number of nominal p < 0.05 hits and number with BH FDR < 0.05.
# OPTIONAL: Repeat the workflow with a different disease, cutoff or sex.
# Label all changes AFTER seeing p-values as exploratory comparisons.
# QUESTION: Why can post hoc gene selection bias your interpretation?
# QUESTION: What can we conclude about association versus causation?
