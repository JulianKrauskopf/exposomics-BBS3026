# Exposomics — BBS3026

Welcome to **BBS3026: Exposomics**.

In this course, you will learn to use R to explore environmental exposures and think critically about scientific evidence linking chemicals to human health. You will develop your own research example based on a compound and related health outcomes.

## Start with Session 1

1. Follow the [setup guide](SETUP.md) to install R, RStudio and the required package(s).
2. Read [Session 1 — Getting Started](01_Getting_Started/README.md).
3. Open [01_PubMed_Exercise.R](01_Getting_Started/01_PubMed_Exercise.R) in RStudio.
4. Save a personal copy as `01_PubMed_Exercise_YourName.R` and complete the seven short exercises.

You will use PubMed search counts to practice **vectors, data frames, filtering, sorting and basic plotting**.

## Session 2 — Exploratory Analysis

Continue with [Session 2 instructions](02_Exploratory_Analysis/README.md) and the interactive [R exercise](02_Exploratory_Analysis/02_Exploratory_Analysis.R). The datasets `data/covar.csv` and `data/expo.csv` are included in this repository; no separate data download is needed.

## Session 3 — Exploratory Analysis II

Continue with [Session 3 instructions](03_Exploratory_Analysis_II/README.md) and the interactive [regression R script](03_Exploratory_Analysis_II/03_Exploratory_Analysis_II.R). It loads the same data as Session 2, so you can start directly after opening the new script.

## Session 4 — Transcriptomics Analysis I

Continue with [Session 4 instructions](04_Transcriptomics/README.md) and the [interactive Open Targets / limma R script](04_Transcriptomics/04_Transcriptomics.R). Download the separate transcriptomics `.Rdata` file as described in [SETUP.md](SETUP.md). This session compares disease-focused gene sets and full-transcriptome results.

## Session 5 — Transcriptomics Analysis II

Use [Session 5 instructions](05_Transcriptomics_II/README.md) and the [interactive result-interpretation and pathway-analysis R script](05_Transcriptomics_II/05_Transcriptomics_II.R). Continue from the two result CSVs saved in Session 4. Explore individual genes, Gene Ontology, Reactome pathways, overrepresentation analysis (ORA) and GSEA.

## Session 6 — Cytoscape and Protein Networks

Finish with [Session 6 instructions](06_Networks/README.md) and the [interactive Cytoscape / STRING network R script](06_Networks/06_Networks.R). Analyze protein interactions, explore hub proteins, compare confidence cutoffs and optionally overlay expression results on WikiPathways. **Install Cytoscape Desktop** and run the [setup check](R/check_cytoscape.R) before class.

**Note:** The number of PubMed publications is not a measure of how harmful a compound is.
