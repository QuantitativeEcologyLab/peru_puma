library(HDInterval)
library(mcmcOutput)
library(wiqid)

#Setting the working directory
setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Peru/puma_project")

#Covariates
prey_covs <- read.csv("covs_prey.csv", header = TRUE)
head(prey_covs)

#puma cam data
puma_data<-read.csv("cam_puma.csv", header = TRUE)
nrow(data)
unique(puma_data$Station)

#prey cam data
prey_data<-read.csv("camdata_prey.csv", header = TRUE)
nrow(prey_data)
table(prey_data$Station)

data()
