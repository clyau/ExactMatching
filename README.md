# ExactMatching

This repository contains the R code used for the simulated examples presented in [Glimm and Yau (2025)](https://arxiv.org/abs/2503.02850). The analysis are done in R version 4.4.0 on platform x86_64-w64-mingw32/x64, running under Windows 11 x64 (build 26100).

To perform the exact matching, R package `maicChecks` (version 0.2.0) is required. The package can be installed from CRAN repository. Additional details about the package can be find [here](https://clyau.github.io/maicChecks/).

This repository contains the following:

* `IllustrativeExamples`: related files for data presented in Section 2.1
    -   `illustrativeExamples.Rmd`: code used to simulate, match, and summarize the data
    -   `illustrative.RData`: simulated data and summary statistics
    -   `fig1.jpeg`, `fig2top.jpeg`, `fig2bott.jpeg`: Figures 1 and 2 in the manuscript
* `TenKSimulations`: related files for data presented in Section 2.3
    -   `TenKSimulations.Rmd`: code used to simulate, match, and summarize the data
    -   `fn_simdt.R`: R script containing the main function to simulate one dataset. This script is required by `TenKSimulations.Rmd`.
    -   `funs.R`: miscellaneous functions required by `TenKSimulations.Rmd` or `fn_simdt.R`
