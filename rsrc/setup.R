
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: R PACKAGES SETUP
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# R PACKAGES SETUP

# INSTALL/LOAD GROUNDHOG PACKAGE (reproducibility solution for handling different versions of R and its packages)

# check if groundhog is installed and load it
if ("groundhog" %in% installed.packages()[,"Package"]) {
  library("groundhog")
} else {
  install.packages("groundhog")
  library("groundhog")
}

# define date of reference to load all packages (default works well with R 4.3.0 and groundhog 3.1.0, adjust the date to your needs)
groundhog.date <- "2023-05-06"

# guarantee version 3.1.0 of groundhog is being used (if you want a different version adjust the date)
groundhog::meta.groundhog(date = "2023-05-06")


# INSTALL/LOAD TINYTEX PACKAGE (A lightweight and easy-to-maintain LaTeX distribution used to generate pdfs with RMarkdown)
groundhog::groundhog.library("tinytex", groundhog.date)

# install tinytex version 2023.05 (only if it is not already installed)
if (tinytex::tinytex_root() == "") {
  tinytex::install_tinytex(version = "2023.05")
}


# ADDITIONAL PACKAGES
# INSTALL/LOAD PACKAGES USED IN THE TEMPLATE SCRIPTS (ADJUST THE LIST BELOW TO YOUR NEEDS)

# here - reproducibility solution for handling file paths
groundhog::groundhog.library("here", groundhog.date)

# tictoc - to register scripts time processing
groundhog::groundhog.library("tictoc", groundhog.date)

# rmarkdown - to use Rmarkdown files
groundhog::groundhog.library("rmarkdown", groundhog.date)

# sf - to manipulate spatial data (vector format)
groundhog::groundhog.library("sf", groundhog.date)

# tidyverse - an opinionated collection of R packages designed for data science that share an underlying design philosophy, grammar, and data structures
groundhog::groundhog.library("tidyverse", groundhog.date)

# sjlabelled - to label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("sjlabelled", groundhog.date)

# Hmisc - to use `describe` function to generate codebook
groundhog::groundhog.library("Hmisc", groundhog.date)

# conflicted - to provide an alternative conflict resolution strategy when multiple packages have functions with the same name
#  making every conflict an error and forcing you to choose which function to use
groundhog::groundhog.library("conflicted", groundhog.date)

# janitor - to use simple functions for examining and cleaning dirty data
groundhog::groundhog.library("janitor", groundhog.date)

# kableExtra - to create beautiful and flexible tables in R
groundhog::groundhog.library("kableExtra", groundhog.date)

# modelsummary - to create and customize tables of model output
groundhog::groundhog.library("modelsummary", groundhog.date)

# terra - to manipulate spatial data in raster format
groundhog::groundhog.library("terra", groundhog.date)

# fixest - provides a family of functions to perform estimations with multiple fixed-effects.
# This package is currently (Feb. 2020) the fastest software available to perform fixed-effects estimations.
groundhog::groundhog.library("fixest", groundhog.date)

# readxl - to read excel files
groundhog::groundhog.library("readxl", groundhog.date)

# markovchain - to calculate probability transition matrix
groundhog::groundhog.library("markovchain", groundhog.date)

# RColorBrewer - to add color palettes
groundhog::groundhog.library("RColorBrewer", groundhog.date)

# marignaleffects - to plot conditional marginal effects
groundhog::groundhog.library("marginaleffects", groundhog.date)

# ggpubr - to combine multiple ggplots
groundhog::groundhog.library("ggpubr", groundhog.date)


# equatiomatic - to convert lm models to latex
#groundhog::groundhog.library("equatiomatic", groundhog.date)

# # data.table - An alternative way to organize data sets for very, very fast operations. Useful for big data.
#groundhog::groundhog.library("data.table", groundhog.date)
# # declare preference to avoid conflicted error with tictoc::shift
#conflicted::conflicts_prefer(data.table::shift)

# # geobr - to read brazilian shapefiles
#groundhog::groundhog.library("geobr", groundhog.date)


# SOURCE FUNCTIONS

# ExportTimeProcessing - function to export scripts time processing, registered with ticotoc,to csv files
#source(here::here("code/_functions/ExportTimeProcessing.R"))

# END OF SCRIPT
