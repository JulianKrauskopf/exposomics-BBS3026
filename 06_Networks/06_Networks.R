# BBS3026 — Session 6: Pathway and Protein Interaction Networks
# Final practical: interpret network connections for your exposomics project.
# Work through the tasks one at a time and complete the ... placeholders.
# Save your own copy as 06_Networks_YourName.R.

# GET READY — Use the SAME course folder as in Sessions 4 and 5.
# Uncomment ONE example and replace YOUR_USERNAME:
# Windows:
# setwd("C:/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
# macOS:
# setwd("/Users/YOUR_USERNAME/Downloads/exposomics-BBS3026-main")
getwd()
file.exists("outputs/session4_full_transcriptome.csv")
file.exists("outputs/session4_disease_focused.csv")
# Both should return TRUE. Do not overwrite your personal Session 4 results.

# TASK 1 — Load packages, helper functions and your earlier results
library(dplyr)
library(ggplot2)
library(RCy3)
library(igraph)
source("R/project_helpers.R")
source("R/network_helpers.R")
full_results <- read.csv("outputs/session4_full_transcriptome.csv")
disease_results <- read.csv("outputs/session4_disease_focused.csv")
dim(full_results)
dim(disease_results)
head(disease_results[, c("feature_id", "logFC", "P.Value",
                         "FDR_disease_set", "adj.P.Val")])
# QUESTION: Which p-values reflect a disease-set versus all-gene correction?

# TASK 2 — Decide which gene set to show in your network
# Start with your PREDEFINED disease-related genes from Open Targets.
# Alternatively, use full_results and label that comparison exploratory.
network_results <- disease_results
# network_results <- full_results
network_annotations <- prepare_network_annotations(network_results)
head(network_annotations[, c("GeneName", "logFC", "P.Value",
                             "StudentNominal", "FDR_analysis")])
# QUESTION: What do GeneName and logFC mean?
# NOTE: logFC is an EXPOSURE-EXPRESSION regression coefficient, not
# a difference between two experimental treatment groups.

# TASK 3 — Select genes for a protein interaction network
# For teaching, select nominal p < 0.05 genes (UNCORRECTED).
# A selection based on nominal p is exploratory, not FDR confirmation.
selection_mode <- "nominal"
# If fewer than 3 map, try "exploratory_top" and label it explicitly.
network_genes <- select_network_genes(
  annotations = network_annotations,
  mode = selection_mode,
  n_top = 30,
  max_genes = 80
)
nrow(network_genes)
head(network_genes[, c("GeneName", "logFC", "P.Value")])
# QUESTION: Why would selecting genes by adjusted p < 0.05 give a
# smaller network? What are the risks of using nominal hits?

# TASK 4 — Connect RStudio to Cytoscape Desktop
# FIRST open the Cytoscape application and leave it RUNNING.
# Install Cytoscape and its apps using SETUP.md before this task.
source("R/check_cytoscape.R")
cytoscapePing()
# QUESTION: Is Cytoscape a separate program from RStudio?

# TASK 5 — Understand the STRING network settings
# Nodes = proteins; edges = known or predicted functional associations,
# which are not necessarily direct physical protein binding.
# A higher STRING confidence cutoff keeps fewer, stronger-evidence edges.
string_cutoff <- 0.9    # High: 0.9; medium: 0.7; exploratory: 0.4
# Start with 0.9. Later test whether conclusions change at 0.7.
# QUESTION: What do you expect to happen if you lower the cutoff?

# TASK 6 — Build a human protein interaction network
# This uses STRING in Cytoscape with NO extra proteins added.
# The helper imports the network, colors nodes by logFC, sizes by degree,
# and saves an image and a hub table.
ppi <- build_ppi_network(
  input_genes = network_genes,
  cutoff = ...
)
# Fill the gap with string_cutoff.
# If fewer than 3 mapped proteins or no edges exist, ppi may be NULL.
# If so, try a justified lower cutoff or exploratory_top selection.
if (!is.null(ppi)) {
  ppi$requested_genes
  ppi$matched_nodes
  head(ppi$stats)
}
# QUESTION: How many requested genes match proteins in STRING?
# QUESTION: How many nodes and edges can you see in Cytoscape?

# TASK 7 — Read the Cytoscape network
# The network created above has:
# - blue: negative exposure coefficient
# - orange: positive exposure coefficient
# - larger nodes: higher number of connections
# - green borders: nominal gene-level p < 0.05
# The colors do not prove that proteins or pathways are activated.
if (!is.null(ppi)) {
  getNodeCount(network = ppi$network)
  getEdgeCount(network = ppi$network)
}
# QUESTION: Is a highly connected protein necessarily most important
# biologically, or does network coverage influence its degree?

# TASK 8 — Identify hubs (degree and betweenness)
# Degree counts a node's connections.
# Betweenness measures how frequently a node lies on shortest paths.
if (!is.null(ppi)) {
  hub_table <- ppi$stats
  hub_table <- hub_table[order(-hub_table$...), ]
  head(hub_table[, c("Protein", "StudentDegree",
                     "StudentBetweenness")], 10)
}
# Hint: order by StudentDegree.
# QUESTION: Which protein has the most connections?
# QUESTION: Does the same protein have the highest betweenness?

# TASK 9 — Investigate the confidence threshold
# Rebuild the same network with a less stringent cutoff.
# WARNING: This makes ANOTHER network in your Cytoscape session.
cutoff_sensitivity <- 0.7
ppi_looser <- build_ppi_network(
  input_genes = network_genes,
  cutoff = cutoff_sensitivity
)
if (!is.null(ppi_looser)) {
  getNodeCount(network = ppi_looser$network)
  getEdgeCount(network = ppi_looser$network)
}
# QUESTION: What changed? Did the main hub change?
# NOTE: The helper reports statistics for the LARGEST connected component.

# TASK 10 — Compare your pathway findings with the network
# Retrieve your GO and Reactome results from Session 5.
# GSEA ranks all tested genes; ORA is based on a selected hit list.
go_file <- "outputs/session5_GO_GSEA.csv"
reactome_file <- "outputs/session5_Reactome_GSEA.csv"
if (file.exists(go_file)) {
  go_results <- read.csv(go_file)
  head(go_results[, intersect(c("ID", "Description", "NES",
                               "p.adjust"), names(go_results))])
}
if (file.exists(reactome_file)) {
  reactome_results <- read.csv(reactome_file)
  head(reactome_results[, intersect(c("ID", "Description", "NES",
                                     "p.adjust"), names(reactome_results))])
}
# QUESTION: Does your interaction network reflect any biological themes
# from the enriched GO terms or Reactome pathways?
# QUESTION: Why does an enriched pathway not guarantee connected proteins?

# TASK 11 — (OPTIONAL) Visualize a curated WikiPathways pathway
# Find a HUMAN pathway at https://www.wikipathways.org/
# Copy its WP identifier, for example WP179 (illustrative only).
# Select a pathway relevant to your research question, not one picked
# solely to obtain a striking image.
pathway_id <- "WP179"   # Replace with your chosen human WP ID.
# The full results supply broad expression annotation for the pathway.
pathway_annotations <- prepare_network_annotations(full_results)
# Once Cytoscape is running, uncomment this line:
# pathway_plot <- overlay_pathway(pathway_annotations, pathway_id)
# QUESTION: How is a curated pathway diagram different from a STRING PPI?
# QUESTION: How many pathway nodes matched measured genes?
# The pathway image is saved to outputs/pathway_<ID>.png.

# TASK 12 — Make an additional network figure
# To preserve the current network image, choose a distinct filename.
# If both network runs succeeded, choose the stricter network to export.
if (!is.null(ppi)) {
  exportImage("outputs/session6_PPI_high_confidence",
              type = "PNG", network = ppi$network)
}
# QUESTION: Which plot and legend would best explain your conclusions?

# TASK 13 — Save network statistics and the Cytoscape project
dir.create("outputs", showWarnings = FALSE)
if (!is.null(ppi)) {
  write.csv(ppi$stats,
            "outputs/session6_hub_statistics_high_confidence.csv",
            row.names = FALSE)
}
# Save everything currently open in Cytoscape in a .cys project:
saveSession("outputs/Session6_Exposomics_Networks")
# QUESTION: What is the difference between a PNG figure and .cys session?

# TASK 14 — Final interpretation and research poster
# Write down:
# 1. The compound and disease chosen in earlier sessions.
# 2. The gene-list selection rule (nominal p or exploratory top-ranked).
# 3. The STRING cutoff, organism, and number of extra proteins (zero).
# 4. Network nodes, edges, largest component and top-degree hub.
# 5. One biological theme supported by GO/Reactome, if any.
# 6. One limitation: confounding, annotation bias, missing proteins,
#    overlapping pathways, post hoc choices, or multiple testing.
# QUESTION: Does protein connectivity prove an exposure causes disease?
# QUESTION: What findings and caveats will you include in your poster?
