# Session 5 — Transcriptomics Analysis II: From Genes to Pathways

Continue from [Session 4](../04_Transcriptomics/README.md). This is a **DataCamp-style plain R script**: students fill in the missing expressions, inspect results and answer short questions.

Open [05_Transcriptomics_II.R](05_Transcriptomics_II.R), save a copy as `05_Transcriptomics_II_YourName.R`, and run each task in order.

## Before class

1. Download the latest [BBS3026 course ZIP](https://github.com/JulianKrauskopf/exposomics-BBS3026) and extract it, **or update the existing course folder**.
2. Keep your own **Session 4 outputs** inside that course folder: `outputs/session4_full_transcriptome.csv` and `outputs/session4_disease_focused.csv`. They contain the full and disease-focused limma results from **your chosen compound and disease**. These personal output files are not stored in GitHub.
3. For Task 7 (individual-gene scatterplot), retain `data/EGM_Transcriptomics_GeneSymbol_Full.Rdata` and use the **same exposure and participant subgroup** as Session 4.
4. Run `source("R/setup_packages.R")` in the Console once to install the additional enrichment packages. Installation may take time and requires internet.
5. Set your working directory to the **course root**, where both `data/` and `outputs/` are located.

## Tasks

| Tasks | Topic |
|---|---|
| 1–3 | Import tables, understand result columns, nominal p versus BH FDR |
| 4–5 | p-value distributions, nominal–adjusted scatterplot and volcano plots |
| 6–7 | Choose a gene; plot expression against exposure and interpret adjustments |
| 8–9 | Convert gene symbols/Ensembl IDs to Entrez IDs, define hits and background |
| 10–12 | Gene Ontology and Reactome overrepresentation analysis (ORA) |
| 13–14 | Gene Ontology and Reactome gene set enrichment analysis (GSEA) |
| 15 | Export tables and interpret findings for a poster |

## Core concepts

**Nominal p-value (`P.Value`)**: unadjusted per-gene test result. **Whole-transcriptome FDR (`adj.P.Val`)**: Benjamini–Hochberg adjusted across all measured genes. **Disease-set FDR (`FDR_disease_set`)**: BH-adjusted within the preselected and measured Open Targets disease genes. Each addresses a different family of tests.

**Overrepresentation analysis (ORA)** starts from a cutoff-based hit list (here, nominal p < 0.05 for classroom exploration) and asks whether a functional category contains more hits than expected, using the **measured, mapped genes as the universe**. The selected hits are not necessarily confirmed discoveries.

**GSEA** uses the **entire ranked list** of measured, mapped genes, using limma's moderated t-statistics; it does not require choosing a nominal significance cutoff.

**Gene Ontology Biological Process** and **Reactome** are two different functional annotation systems. Their term coverage, gene mappings and available gene-set sizes differ, so the results may not agree.

## Interpretation and limitations

If very few genes pass your threshold, ORA may produce no informative terms. This is a valid result; do not change thresholds just to obtain significant pathways. Some gene identifiers will not map to annotation databases. Small disease-focused gene sets are often poorly suited to ORA. Adjusted p-values in enrichment tables correct **pathway/term tests**, not the original gene-level tests. Exploratory comparisons after seeing expression p-values must not be described as preregistered hypotheses.

**Important:** Session 5 interprets the Session 4 model already fitted; it does not refit that model. It does not demonstrate causal effects of exposures on disease.

Keep personal result tables, expression matrices, and analysis outputs locally rather than uploading them to a public GitHub repository.
