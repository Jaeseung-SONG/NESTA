#!/usr/bin/env Rscript
Sys.setenv(R_LIBS_USER = "/home/js/R/x86_64-pc-linux-gnu-library/4.1")
source(file.path("/home/js/NESTA/simulation", "R", "utils.R"))
source(project_file("R/study_0628.R"))
run_study_0628()
