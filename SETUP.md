# Setup — R and RStudio

Complete these steps before starting the practical sessions.

## 1. Install R

Download and install R for your operating system from [CRAN](https://cran.r-project.org/).

## 2. Install RStudio Desktop

Download and install [RStudio Desktop](https://posit.co/download/rstudio-desktop/). Start RStudio after installing R.

## 3. Download the course files

On the [course GitHub page](https://github.com/JulianKrauskopf/exposomics-BBS3026), choose **Code → Download ZIP**. Extract the ZIP file and open the course folder in RStudio.

## 4. Install the course packages

In the RStudio **Console**, run:

```r
source("R/setup_packages.R")
```

You need an internet connection for installation. The setup installs **dplyr, ggplot2 and corrplot** for Session 2 and **httr, jsonlite and limma (Bioconductor)** for Session 4. For Session 5 it also installs **clusterProfiler, org.Hs.eg.db, AnnotationDbi, ReactomePA and enrichplot** from Bioconductor. The first installation may take several minutes. For Session 6 it also installs **igraph** (CRAN) and **RCy3** (Bioconductor) for Cytoscape integration. **Session 1 uses base R**, so it does not require packages.

The ZIP also contains both Session 2 datasets: `data/covar.csv` (participant characteristics) and `data/expo.csv` (exposure measurements). Keep them in the `data/` folder. Do not open or resave these CSVs in Excel.

## 5. Download the Session 4 transcriptomics dataset

For Session 4, download **[EGM_Transcriptomics_GeneSymbol_Full.Rdata (Dropbox)](https://www.dropbox.com/scl/fi/e3b3euudt0propoh4gsgy/EGM_Transcriptomics_GeneSymbol_Full.Rdata?rlkey=rtcakxm6s9ij3d5g35yghqo3w&dl=0)** and place it in the **existing `data/` folder**, alongside `covar.csv` and `expo.csv`.

Do **not** upload this transcriptomics file to GitHub. The Open Targets disease-gene retrieval in Session 4 requires internet access.

From the course folder, check that all three files are present:

```r
file.exists("data/covar.csv")
file.exists("data/expo.csv")
file.exists("data/EGM_Transcriptomics_GeneSymbol_Full.Rdata")
```

All three checks should return `TRUE` before Session 4.

## 6. Continue with Session 5

Session 5 reads **your own saved Session 4 results** from `outputs/session4_full_transcriptome.csv` and `outputs/session4_disease_focused.csv`. Complete **Task 15 of Session 4** first. These files are generated on your computer and are not included in GitHub.

If you downloaded a fresh copy of the repository, copy your personal `outputs/` folder and the transcriptomics `.Rdata` file into it. Keep the same working directory for Session 5. The [Session 5 instructions](05_Transcriptomics_II/README.md) explain gene plots, GO analysis, Reactome pathways, ORA and GSEA.

## 7. Install Cytoscape Desktop for Session 6

Cytoscape is a **separate desktop application**, not part of RStudio.

1. Download [Cytoscape Desktop](https://cytoscape.org/download.html) for **Windows or macOS**, run the installer and start the program.
2. Wait until Cytoscape is fully open. **Keep Cytoscape running** while you use RStudio.
3. From your extracted course folder, run the following in the RStudio Console:

```r
source("R/check_cytoscape.R")
```

This checks RStudio's connection to Cytoscape and automatically installs two Cytoscape apps when missing: **stringApp** (protein interaction networks) and **WikiPathways** (curated pathway diagrams). You need internet access for app installation and STRING queries.

If automatic app installation fails, in the Cytoscape window open **Apps → App Manager**, search for `stringApp` and `WikiPathways`, install them, and rerun the check. If R cannot connect, make sure Cytoscape has finished starting and check local firewall prompts.

The network analysis files are in [06_Networks](06_Networks/README.md). Keep your personal Session 4 and 5 results in the local `outputs/` folder. The network session and images should also stay local.

## 8. Check your installation

Type `2 + 2` in the Console and press Enter. R should return `[1] 4`.

If you encounter a problem, bring the error message to class.
