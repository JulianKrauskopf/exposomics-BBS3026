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

You need an internet connection for installation. The setup installs **dplyr, ggplot2 and corrplot** for Session 2. **Session 1 uses base R**, so it does not require packages.

The ZIP also contains both Session 2 datasets: `data/covar.csv` (participant characteristics) and `data/expo.csv` (exposure measurements). Keep them in the `data/` folder. Do not open or resave these CSVs in Excel.

## 5. Check your installation

Type `2 + 2` in the Console and press Enter. R should return `[1] 4`.

If you encounter a problem, bring the error message to class.
