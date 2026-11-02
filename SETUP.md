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

You need an internet connection for installation. The setup installs **dplyr, ggplot2 and corrplot** for Session 2 and **httr, jsonlite and limma (Bioconductor)** for Session 4. **Session 1 uses base R**, so it does not require packages.

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

## 6. Check your installation

Type `2 + 2` in the Console and press Enter. R should return `[1] 4`.

If you encounter a problem, bring the error message to class.
