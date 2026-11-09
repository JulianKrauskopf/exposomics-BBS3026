# BBS3026 — Session 5: Transcriptomics Analysis II
# From gene-level results to Gene Ontology, pathways and GSEA
# Complete the ... placeholders; run ONE TASK AT A TIME.
# Save your personal copy as 05_Transcriptomics_II_YourName.R.

# BEFORE YOU START — Work in the SAME course folder as Session 4.
# Uncomment one example, and replace YOUR_USERNAME if required:
# Windows:
# setwd("C:/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
# macOS:
# setwd("/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
getwd()
file.exists("outputs/session4_full_transcriptome.csv")
file.exists("outputs/session4_disease_focused.csv")
# Both checks must be TRUE. Otherwise finish Session 4 Task 15 first.
# Run source("R/setup_packages.R") once if packages are missing.

# TASK 1 — Load the packages and Session 4 results
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
source("R/enrichment_helpers.R")

full_results <- read.csv("outputs/session4_full_transcriptome.csv")
disease_results <- read.csv("outputs/session4_disease_focused.csv")

dim(full_results)
dim(disease_results)
names(full_results)
names(disease_results)
head(full_results, 5)
head(disease_results, 5)
# QUESTION: Why does the disease-focused table contain fewer genes?
# QUESTION: Which columns occur only in the disease-focused table?

# TASK 2 — Understand the result columns
# feature_id = transcript/gene ID; logFC = fitted exposure coefficient;
# t = moderated t statistic; P.Value = nominal (unadjusted) p-value;
# adj.P.Val = BH-adjusted p-value over ALL measured genes.
# FDR_disease_set = BH-adjusted p-value within the predefined disease set.
# association_score = Open Targets DISEASE–gene evidence score.
# In Session 4, logFC is the change in expression per unit log2 exposure.
head(full_results[, c("feature_id", "logFC", "t", "P.Value",
                      "adj.P.Val")])
head(disease_results[, c("feature_id", "P.Value",
                         "FDR_disease_set", "adj.P.Val",
                         "association_score")])

# TASK 3 — Count nominal and adjusted discoveries
# Fill in the missing column name.
sum(full_results$P.Value < 0.05, na.rm = TRUE)
sum(full_results$... < 0.05, na.rm = TRUE)
sum(disease_results$P.Value < 0.05, na.rm = TRUE)
sum(disease_results$FDR_disease_set < 0.05, na.rm = TRUE)
sum(disease_results$adj.P.Val < 0.05, na.rm = TRUE)
# QUESTION: Why are there often more nominal hits than adjusted hits?
# QUESTION: Why do the two adjusted p-value columns differ?

# TASK 4 — Inspect and plot the p-value distributions
hist(full_results$P.Value, breaks = 30,
     main = "All measured genes", xlab = "Nominal p-value")
hist(disease_results$..., breaks = 20,
     main = "Disease-focused genes", xlab = "Nominal p-value")
# Scatterplot nominal p versus adjusted p (both are in [0,1]).
ggplot(full_results, aes(x = P.Value, y = adj.P.Val)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Nominal p-value", y = "Whole-transcriptome BH FDR")
# QUESTION: How does correction change the apparent significance?

# TASK 5 — Explore volcano plots
# A volcano plot displays gene coefficients against -log10(nominal p).
ggplot(full_results,
       aes(x = logFC, y = -log10(pmax(P.Value, .Machine$double.xmin)))) +
  geom_point(alpha = 0.4) +
  labs(x = "Exposure coefficient (expression per log2 exposure)",
       y = "-log10(nominal p)", title = "Full transcriptome")

# Add colored points highlighting nominal p < 0.05.
ggplot(disease_results,
       aes(x = logFC, y = -log10(pmax(P.Value, .Machine$double.xmin)),
           color = ...)) +
  geom_point(alpha = 0.7) +
  labs(title = "Predefined disease-related genes")
# Hint: color = P.Value < 0.05, or color = FDR_disease_set < 0.05.
# QUESTION: Does a large coefficient guarantee a small p-value?

# TASK 6 — Investigate individual genes
# Select a measured gene by feature_id from the full results.
# Explore the lowest nominal p-value genes, then choose one.
top_genes <- full_results[order(full_results$P.Value), ]
head(top_genes[, c("feature_id", "logFC", "P.Value", "adj.P.Val")])
my_gene <- "YOUR_GENE"  # Replace with one feature_id from the table.
gene_row <- subset(full_results, feature_id == my_gene)
gene_row
# QUESTION: Is this gene also included in the disease-focused list?
my_gene %in% disease_results$feature_id

# TASK 7 — Plot actual expression versus exposure
# This plot requires the large Session 4 transcriptomics Rdata file.
# Choose the SAME measured chemical and participant group as Session 4.
# Do not substitute a different exposure for the model in your result tables.
file.exists("data/EGM_Transcriptomics_GeneSymbol_Full.Rdata")
source("R/project_helpers.R")
covar <- read.csv("data/covar.csv", row.names = 1)
expo <- read.csv("data/expo.csv", row.names = 1)
my_compound <- "PCB"  # CHANGE to the compound used in Session 4.
my_sex <- "all"       # CHANGE if Session 4 used "F" or "M".
covar$sex <- factor(covar$sex)
covar$cohort <- factor(covar$cohort)
covar$education <- factor(covar$education)
covar$activity <- factor(covar$activity)
covar$smoking_status <- factor(covar$smoking_status)
study <- prepare_exposure(covar, expo, my_compound)
omics <- load_expression()
omics_norm <- limma::normalizeBetweenArrays(omics, method = "quantile")
stopifnot(my_gene %in% rownames(omics_norm))
matched_ids <- rownames(study)[rownames(study) %in% colnames(omics_norm)]
if (my_sex != "all") {
  matched_ids <- matched_ids[study[matched_ids, "sex"] == my_sex]
}
plot_data <- data.frame(
  log2_exposure = study[matched_ids, "log2_exposure"],
  expression = as.numeric(omics_norm[my_gene, matched_ids]),
  cohort = study[matched_ids, "cohort"]
)
ggplot(plot_data, aes(x = log2_exposure, y = expression)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ ...) +
  labs(x = paste("log2", my_compound), y = paste(my_gene, "expression"))
# QUESTION: Does the unadjusted scatterplot resemble the adjusted limma coefficient?
# Why might they differ? (Think about age, sex and cohort.)

# TASK 8 — Prepare gene identifiers for functional analysis
# Gene Ontology (GO) and Reactome analyses require Entrez IDs.
# The mapping helper handles SYMBOL or Ensembl identifiers.
full_mapped <- map_results_to_entrez(full_results)
disease_mapped <- full_mapped[
  full_mapped$feature_id %in% disease_results$feature_id, , drop = FALSE]
nrow(full_mapped)
nrow(disease_mapped)
# QUESTION: Why did some measured genes not map to an Entrez ID?

# TASK 9 — Choose the gene list and correct background (universe)
# For overrepresentation analysis (ORA), choose a list of "hits".
# Nominal p < 0.05 is an EXPLORATORY teaching choice, not confirmed hits.
# You can later try "adj.P.Val" and compare what changes.
selection_column <- "P.Value"
full_hits <- full_mapped$ENTREZID[
  full_mapped[[selection_column]] < 0.05]
full_universe <- full_mapped$ENTREZID
length(full_hits)
length(full_universe)
# QUESTION: Why must the background be the genes actually tested and mapped?

# TASK 10 — Gene Ontology overrepresentation (ORA)
# GO Biological Process: is a GO term represented more than expected
# among your selected hits, relative to the measured background?
if (length(full_hits) >= 1) {
  go_ora <- enrichGO(
    gene = full_hits,
    universe = full_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "...",                # Hint: "BP"
    pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1,
    readable = TRUE
  )
  head(as.data.frame(go_ora)[, c("Description", "GeneRatio",
                                  "BgRatio", "pvalue", "p.adjust")])
  if (nrow(as.data.frame(go_ora)) > 0) {
    barplot(go_ora, showCategory = 10)
  }
} else {
  message("No selected hits: GO ORA cannot be run.")
}
# QUESTION: What do GeneRatio, BgRatio and adjusted p-value mean?
# QUESTION: Are any GO terms below adjusted p < 0.05?

# TASK 11 — Reactome pathway overrepresentation
# Reactome describes curated biological pathways.
if (length(full_hits) >= 1) {
  pathway_ora <- enrichPathway(
    gene = full_hits,
    universe = full_universe,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1,
    readable = TRUE
  )
  head(as.data.frame(pathway_ora))
  if (nrow(as.data.frame(pathway_ora)) > 0) {
    dotplot(pathway_ora, showCategory = 10)
  }
} else {
  message("No selected hits: pathway ORA cannot be run.")
}
# QUESTION: Which pathways are overrepresented, if any?

# TASK 12 — Compare disease-focused versus full ORA
# Disease-focused ORA must use only MAPPED, MEASURED disease genes
# as its background. A small set may yield no informative terms.
disease_hits <- disease_mapped$ENTREZID[
  disease_mapped$P.Value < 0.05]
disease_universe <- disease_mapped$ENTREZID
length(disease_hits)
length(disease_universe)
if (length(disease_hits) >= 1) {
  disease_go_ora <- enrichGO(
    gene = disease_hits,
    universe = disease_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1, qvalueCutoff = 1,
    readable = TRUE
  )
  head(as.data.frame(disease_go_ora))
} else {
  message("No nominal disease-set hits: skip disease-focused ORA.")
}
# QUESTION: How does changing the gene universe change your interpretation?
# NOTE: Do not present disease-focused and full ORA as directly comparable
# if their universes and numbers of tested genes differ.

# TASK 13 — Gene Set Enrichment Analysis (GSEA)
# Unlike ORA, GSEA uses a RANKED LIST of ALL tested mapped genes,
# not a cutoff-selected list. Use moderated t from limma for direction.
gene_ranks <- full_mapped$t
names(gene_ranks) <- full_mapped$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = ...)
length(gene_ranks)
# A positive t statistic means a positive exposure-expression association.
set.seed(1)
go_gsea <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)
head(as.data.frame(go_gsea))
if (nrow(as.data.frame(go_gsea)) > 0) {
  dotplot(go_gsea, showCategory = 10)
}
# QUESTION: What does a positive versus negative NES mean?

# TASK 14 — Reactome pathway GSEA
set.seed(1)
pathway_gsea <- gsePathway(
  geneList = gene_ranks,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)
head(as.data.frame(pathway_gsea))
if (nrow(as.data.frame(pathway_gsea)) > 0) {
  dotplot(pathway_gsea, showCategory = 10)
}
# QUESTION: Do ORA and GSEA identify the same biological themes?
# QUESTION: Why can they disagree?

# TASK 15 — Save tables for your poster
dir.create("outputs", showWarnings = FALSE)
if (exists("go_ora")) {
  write.csv(as.data.frame(go_ora),
            "outputs/session5_GO_ORA_full.csv", row.names = FALSE)
}
if (exists("pathway_ora")) {
  write.csv(as.data.frame(pathway_ora),
            "outputs/session5_Reactome_ORA_full.csv", row.names = FALSE)
}
write.csv(as.data.frame(go_gsea),
          "outputs/session5_GO_GSEA.csv", row.names = FALSE)
write.csv(as.data.frame(pathway_gsea),
          "outputs/session5_Reactome_GSEA.csv", row.names = FALSE)

# FINAL REFLECTION
# 1. Which result is strongest after correction for multiple testing?
# 2. Which analyses use the full transcriptome, and which use your disease set?
# 3. What are the limitations of pathway enrichment and post hoc exploration?
# 4. Which result would you show in your research poster, and why?
