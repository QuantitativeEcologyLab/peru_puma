Interspecific Interactions of Pumas and Jaguars in the Peruvian Amazon
Camera trap analysis of spatial and temporal interactions between pumas (Puma concolor), jaguars (Panthera onca), and their prey species in Las Piedras, Madre de Dios, Peru.
Authors

Dwija Desai — Department of Biology, UBC Okanagan; Okanagan Institute for Biodiversity, Resilience, and Ecosystem Services, UBC; Department of Computer Science, Math, Physics, and Statistics, UBC; Hoja Nueva, Madre de Dios, Peru
Samantha Zwicker — Hoja Nueva, Madre de Dios, Peru
Michael J. Noonan — Department of Biology, UBC Okanagan; Okanagan Institute for Biodiversity, Resilience, and Ecosystem Services, UBC; Department of Computer Science, Math, Physics, and Statistics, UBC

Introduction
Pumas and jaguars are the two largest felids in the Americas, co-occurring across the Peruvian Amazon. Schoener's niche partitioning theory predicts that sympatric competitors must segregate along space, time, or diet to coexist. While studies in Central America and the Pantanal have documented such mechanisms, Amazonian baselines from intact forests remain scarce. This project investigates puma–jaguar coexistence in the low-disturbance rainforest of Las Piedras, Peru.
Research Questions

Do pumas and jaguars co-occur more or less than expected by chance?
Does prey presence influence predator occupancy differently for each felid?
Do the two predators partition activity time to reduce competition?

Study Area & Design
Five camera trapping grids (215 stations, 17,661 trap nights) were deployed across ~115 km² of low-disturbance rainforest in Las Piedras, Madre de Dios, Peru (2019–2023).
Methods

Spatial overlap: Kernel density estimation of detection locations
Co-occurrence & occupancy: Multi-species occupancy models (occuMulti from the unmarked package) with 5-day sampling occasions
Temporal overlap: Diel activity overlap using the coefficient of overlapping (Δ) with 999 bootstrap replicates
Independence filter: Consecutive detections of the same species at the same station within 30 minutes are treated as a single event

Prey Species Analysed
Agouti, Acouchi, Brazilian Rabbit, Brocket Deer, Collared Peccary, Giant Anteater, Opossum, Paca, Small Armadillo, and Tapir.
Repository Structure
FileDescriptionRAI.RCalculates Relative Abundance Index for all species after applying a 30-min independence filterjag_prey.RJaguar–prey multi-species occupancy models (marginal and conditional occupancy)puma_prey.RPuma–prey multi-species occupancy modelspuma_jag.RPuma–jaguar two-species occupancy model with site covariatestemporal_overlap.RTemporal activity overlap analysis using the overlap package (Dhat4 with bootstrap CIs)kde.RSpatial kernel density estimation and overlap maps for predator and prey detections
Key Findings

Pumas and jaguars showed a significant positive spatial association (estimate = 1.126 ± 0.393, p = 0.004), with both species showing higher occupancy when the other was present (~0.60) compared to when absent (~0.34).
Pumas showed significant positive spatial associations with agoutis, opossums, and Brazilian rabbits; jaguars showed a significant positive association only with agoutis.
Puma and jaguar temporal activity overlapped substantially (Δ = 0.84), with both species showing crepuscular and nocturnal activity peaks.
Positive spatial co-occurrence with high temporal overlap suggests coexistence is not driven by strong avoidance on either axis. Shared prey availability in intact forest likely reduces competitive pressure, consistent with niche partitioning theory which predicts relaxed segregation where resources are abundant.

Requirements
R Packages
rinstall.packages(c(
  "camtrapR", "unmarked", "overlap", "astroFns",
  "dplyr", "lubridate", "ggplot2", "patchwork",
  "MASS", "sf", "eks", "ggnewscale", "colorspace",
  "RColorBrewer", "rphylopic",
  "AICcmodavg", "MuMIn"
))
Data Files (not included)
The analysis expects the following CSVs in the working directory:

puma.csv — Puma detection records
jagu.csv — Jaguar detection records
all_prey.csv — Prey species detection records
puma_jag.csv — Combined puma and jaguar records
covariates.csv — Station-level covariates (coordinates, setup/retrieval dates, days operable, grid, trail)

Usage

Place all CSV data files in a single directory.
Update the setwd() path at the top of each script to point to your data directory.
Run scripts in this suggested order:

RAI.R — species-level detection summaries
puma_jag.R — puma–jaguar occupancy
jag_prey.R and puma_prey.R — predator–prey occupancy (run both before the combined marginal occupancy plot in jag_prey.R)
temporal_overlap.R — activity overlap coefficients
kde.R — spatial density maps



Affiliations

Department of Biology, The University of British Columbia Okanagan, Kelowna, BC, Canada
Okanagan Institute for Biodiversity, Resilience, and Ecosystem Services, The University of British Columbia, Kelowna, BC, Canada
Department of Computer Science, Math, Physics, and Statistics, The University of British Columbia Okanagan, Kelowna, BC, Canada
Hoja Nueva, Madre de Dios, Peru
