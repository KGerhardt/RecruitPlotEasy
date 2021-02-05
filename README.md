# RecruitPlotEasy

A tool for interactive Recruitment Plot generation and viewing


# Requirements

- RStudio
- libgit2

# Installation

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
