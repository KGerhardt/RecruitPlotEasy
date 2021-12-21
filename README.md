# RecruitPlotEasy

A tool for interactive Recruitment Plot generation and viewing


# Requirements

- RStudio
- libgit2

# Installation

RecruitPlotEasy requires an existing installation of R and RStudio to work properly. If you do not already have these installed, you will have to install them first before following the numbered instructions below. Installations for both R and RStudio are available for all operating systems and instructions for installing both are available on their respective websites.

R can be found at https://www.r-project.org/

RStudio can be found at https://www.rstudio.com/ RecruitPlotEasy requires only the free RStudio Desktop to work.

To install RecruitPlotEasy:

1. Open an RStudio session
2. Ensure you have devtools installed:
   ```R
   if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
   ```
3. Install RecruitPlotEasy:
   ```R
   devtools::install_github("KGerhardt/RecruitPlotEasy")
   ```


# Use

Load the RecruitPlotEasy library and launch the GUI:
```R
library("RecruitPlotEasy")
RecruitPlotEasy()
```
