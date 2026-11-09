# Session 6 — Networks and Cytoscape (Final Practical)

This **plain R-script, DataCamp-style practical** continues directly from Sessions 4 and 5. Students will construct **human STRING protein interaction networks**, explore network hubs and confidence thresholds, and optionally overlay results on a **WikiPathways** diagram using **Cytoscape Desktop controlled from RStudio**.

Open [06_Networks.R](06_Networks.R), save a personal copy as `06_Networks_YourName.R`, and complete the gaps one task at a time.

## Before class — essential setup

1. Download the latest [BBS3026 course repository](https://github.com/JulianKrauskopf/exposomics-BBS3026) or update your extracted local folder.
2. Keep your **personal results** from Session 4: `outputs/session4_full_transcriptome.csv` and `outputs/session4_disease_focused.csv`. Optionally bring your Session 5 GO and Reactome tables.
3. In the RStudio Console, run `source("R/setup_packages.R")`. The setup includes `RCy3` and `igraph` for network analysis.
4. Install **[Cytoscape Desktop](https://cytoscape.org/download.html)** for Windows or macOS. It is a separate application, not an R package.
5. **Start Cytoscape Desktop** and wait for its window to open. Keep it running during the practical.
6. From the course root in RStudio, run `source("R/check_cytoscape.R")`. This checks the connection and installs the **stringApp** and **WikiPathways** apps if they are missing. You need internet access. If installation fails, open Cytoscape **Apps → App Manager** and install the apps manually; rerun the check afterward.

Detailed platform-specific instructions and example `setwd()` commands are also provided in [SETUP.md](../SETUP.md).

## Activities

| Tasks | Topic |
|---|---|
| 1–3 | Import Session 4 results; choose the disease-focused genes; select exploratory network candidates |
| 4–6 | Connect Cytoscape from R, choose STRING settings, build a protein interaction network |
| 7–9 | Interpret nodes and edges, calculate hubs and betweenness, compare confidence cutoffs |
| 10–11 | Relate networks to Session 5 GO/Reactome findings and optionally import a WikiPathways pathway |
| 12–13 | Export PNG figures, hub table and Cytoscape session (.cys) |
| 14 | Summarize findings and limitations for the final poster |

## What does the network mean?

STRING edges represent **known or predicted protein associations**; they are not all direct physical interactions. A network depends on **which proteins were supplied** and the chosen STRING confidence cutoff. Degree and betweenness are calculated for the **largest connected component** by the course helper.

The starting gene list uses **nominal p < 0.05** within the predefined disease-related genes for demonstration, not FDR-confirmed discoveries. When too few genes are available, a top-ranked exploratory gene list or a lower STRING confidence threshold can demonstrate the software, but must be labeled exploratory.

**Read the node styling carefully:** orange/blue represent positive/negative **exposure-expression regression coefficients**, not guaranteed protein activation/inhibition. Green borders indicate **nominal** associations.

## Files and privacy

The R utilities `R/network_helpers.R` and `R/check_cytoscape.R` come from the established course workflow. All network images, statistics and the `.cys` project are saved in your local `outputs/` folder. Do **not** upload your expression data, personal results or `.cys` project to public GitHub.

A missing interaction network, a disconnected graph or no significant pathways are valid findings. Network evidence alone does not establish exposure-induced disease.
