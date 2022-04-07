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
The installation command is also the way to update your version of RecruitPlotEasy to a newer version, as they become available.

# Use

Load the RecruitPlotEasy library and launch the GUI:
```R
library("RecruitPlotEasy")
RecruitPlotEasy()
```

# Update Python

Updating the version of the python script RecruitPlotEasy uses is a little less straightforward than updating the R code. To do so, open your RStudio and run the following commands:

```R
library(reticulate)
use_miniconda(condaenv = "recruitment_plots", required = T)
recplot_py <- import("RecruitPlotEasy")
py_install(packages = "RecruitPlotEasy", envname = "recruitment_plots", pip = T)
```
You should see a message similar to the following (version numbers may differ)

Collecting RecruitPlotEasy
  Downloading RecruitPlotEasy-2.0.6-py3-none-any.whl (29 kB)
Installing collected packages: RecruitPlotEasy
  Attempting uninstall: RecruitPlotEasy
    Found existing installation: RecruitPlotEasy 2.0.4
    Uninstalling RecruitPlotEasy-2.0.4:
      Successfully uninstalled RecruitPlotEasy-2.0.4
Successfully installed RecruitPlotEasy-2.0.6
