# BBS3026 — Install packages for the course.
# Open the course folder in RStudio, then run:
# source("R/setup_packages.R")

packages <- c("dplyr", "ggplot2", "corrplot", "httr", "jsonlite", "igraph")
missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# Transcriptomics (Session 4) and functional enrichment (Session 5).
# Bioconductor packages must be installed through BiocManager.
bioc_packages <- c("limma", "clusterProfiler", "org.Hs.eg.db",
                   "AnnotationDbi", "ReactomePA", "enrichplot", "RCy3")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package setup complete. You can start the course practicals.")
