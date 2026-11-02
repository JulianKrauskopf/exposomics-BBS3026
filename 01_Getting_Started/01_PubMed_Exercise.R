# BBS3026 — Session 1: Getting Started with R
# Open this script in RStudio and save a personal copy as
# 01_PubMed_Exercise_YourName.R. Complete the missing parts.
# Run your code one task at a time (Ctrl+Enter / Cmd+Enter).

# TASK 1 — Character vectors
# Choose one compound, a comparison compound, and up to six diseases.
# Replace the examples with your own choices.
# Add or remove diseases by editing the c(...) vector.

my_compound <- "YOUR COMPOUND"
other_compound <- "OTHER COMPOUND"
diseases <- c("DISEASE 1", "DISEASE 2", "DISEASE 3")

# Combine your selected compound with each disease.
compound_disease <- paste(my_compound, diseases, sep = " + ")
compound_disease


# TASK 2 — Numeric vectors
# Search PubMed: "compound" AND "disease".
# Replace NA with the actual PubMed counts.
# The numbers must correspond to the diseases in the same order.
# If you change the number of diseases, change both numeric vectors.

pubmed_my_compound <- c(NA, NA, NA)
pubmed_other_compound <- c(NA, NA, NA)


# TASK 3 — Create a data frame
# Fill in the missing vector names (replace ...).

my_pubmed_data <- data.frame(
  compound_disease = compound_disease,
  my_compound_hits = ...,
  other_compound_hits = ...
)

my_pubmed_data


# TASK 4 — Inspect the structure
# Which columns are character and which are numeric?

str(my_pubmed_data)


# TASK 5 — Filter rows with subset()
# Print only entries with more than 500 PubMed results
# for your selected compound.
# Hint: subset(dataframe, column > number)

many_publications <- subset(
  my_pubmed_data,
  my_compound_hits > ...
)

many_publications


# TASK 6 — Sort rows with order()
# Show the diseases with the most PubMed results first.
# Hint: order(x, decreasing = TRUE)

ordered_data <- my_pubmed_data[
  order(my_pubmed_data$..., decreasing = TRUE),
]

ordered_data


# TASK 7 — Make two barplots
# Each bar represents one disease.
# Use different colors for the two compounds.
# Hint: barplot(height = values, names.arg = labels, col = "blue")

# Barplot for your selected compound
barplot(
  height = my_pubmed_data$...,
  names.arg = diseases,
  col = "...",
  main = my_compound,
  ylab = "Number of PubMed results",
  las = 2
)

# Barplot for your comparison compound
barplot(
  height = my_pubmed_data$...,
  names.arg = diseases,
  col = "...",
  main = other_compound,
  ylab = "Number of PubMed results",
  las = 2
)


# THINK ABOUT IT
# 1. Which disease has the most publications for each compound?
# 2. How do the two barplots differ?
# 3. Does a higher PubMed count mean a greater health risk? Why or why not?
