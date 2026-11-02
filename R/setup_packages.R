# BBS3026 — Install packages for the course.
# Open the course folder in RStudio, then run:
# source("R/setup_packages.R")

packages <- c("dplyr", "ggplot2", "corrplot")
missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("Package setup complete. You can start Session 1.")
