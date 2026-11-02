# BBS3026 — Install packages for the course.
# Open the course folder in RStudio, then run:
# source("R/setup_packages.R")

packages <- c("dplyr", "ggplot2", "corrplot", "httr", "jsonlite")
missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# limma is a Bioconductor package (not installed from CRAN).
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

message("Package setup complete. You can start the course practicals.")
