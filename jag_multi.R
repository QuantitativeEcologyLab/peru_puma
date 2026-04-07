setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(camtrapR)
library(dplyr)

#opening the data files
covs <- read.csv("covariates.csv", stringsAsFactors = FALSE)
jag  <- read.csv("jagu.csv", stringsAsFactors = FALSE)
prey <- read.csv("all_prey.csv", stringsAsFactors = FALSE)

unique(prey$Code)

#detection history for prey

