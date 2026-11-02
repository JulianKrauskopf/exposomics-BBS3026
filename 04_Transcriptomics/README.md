# Session 4 — Transcriptomics Analysis I

This session follows your Session 4 slides and original R script, but replaces the manual OMIM/Excel disease-gene list with **Open Targets** (as in the previous exposomics course). Work in [04_Transcriptomics.R](04_Transcriptomics.R): complete the missing code and run each task in order.

## Start here

1. Download or update the [BBS3026 course repository](https://github.com/JulianKrauskopf/exposomics-BBS3026).
2. Follow the updated [SETUP.md](../SETUP.md). In addition to the CSVs already included in the repository, Session 4 requires the separately downloaded `EGM_Transcriptomics_GeneSymbol_Full.Rdata` file in `data/`.
3. Install packages by running `source("R/setup_packages.R")` from the course root. You need an internet connection to retrieve Open Targets data.
4. Open `04_Transcriptomics/04_Transcriptomics.R`; save your own copy as `04_Transcriptomics_YourName.R`.
5. Use the **same exposure compound** you chose in Sessions 2 and 3, and choose the disease from your research question.

## Learning goals

| Tasks | Focus |
|---|---|
| 1–4 | Reopen data, examine expression matrix, quantile normalisation |
| 5–7 | Search Open Targets, choose a disease, set evidence-score cutoff, compare gene counts |
| 8–10 | Match genes and participants, explore stratification |
| 11–12 | Design matrix; fit genome-wide limma regression (expression ~ exposure + covariates) |
| 13–14 | Inspect all-gene results and compare with a predefined disease-gene family |
| 15–16 | Save findings and critically interpret your research question |

## Open Targets and statistical interpretation

Open Targets provides evidence connecting genes/targets to **diseases**; it does not identify genes altered by your selected exposure. Select the disease and score cutoff **before** viewing gene-expression results. A score of **0.2** is the teaching starting point, **not a validated biological threshold**.

You may explore other outcomes or stricter cutoffs to see how the number of eligible genes changes, but describe such comparisons as **exploratory**. The predefined analysis tests measured disease-associated genes; the full analysis tests all measured genes. The full transcriptome is fitted first to use limma's shared variance estimation.

For classroom exploration, report nominal **p < 0.05**, while always displaying both **disease-set BH FDR** and **whole-transcriptome BH FDR**. Nominal hits alone are not confirmation of gene regulation, biological relevance or disease causation.

The helper file `R/project_helpers.R` contains the same Open Targets search, retrieval/cache, gene matching and limma utilities as the earlier course.

## Save your work

Keep `outputs/`, the downloaded transcriptomics file, and student result files **locally**, not on public GitHub. Submit your completed personal R script and a short interpretation as instructed in class.
