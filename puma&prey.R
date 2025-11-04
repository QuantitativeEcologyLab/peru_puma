#Setting the working directory
setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project")

#Loading the packages from library
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(camtrapR)
library(parallel)
library(secr)
library(astroFns)
library(overlap)
library(dplyr)
library(tidyverse)
library(usethis)
library(devtools)
devtools::install_github("r-lib/conflicted")
library(conflicted)
library(dplyr)
library(camtrapR)


#Covariates correlation matrix
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

#Camera operation matrix
preydatacameraOperation <- cameraOperation(prey_covs,
                                       stationCol = "Station", 
                                       setupCol = "Setup_Date", 
                                       retrievalCol = "Retrieval_Date",
                                       hasProblems = FALSE,
                                       allCamsOn = FALSE,
                                       camerasIndependent = FALSE,
                                       dateFormat = "%Y-%m-%d", 
                                       writecsv = FALSE, 
)

puma5 <- detectionHistory(recordTable       = puma_data,
                           camOp                = preydatacameraOperation,
                           stationCol           = "Station",
                           speciesCol           = "Code",
                           recordDateTimeCol    = "DateTimeOriginal",
                           species              = "PUMA",
                           occasionLength       = 5, #number of days per sampling occasion
                           day1                 = "station",
                           datesAsOccasionNames = FALSE,
                           includeEffort        = FALSE,
                           scaleEffort          = FALSE,
                           timeZone             = "America/Lima"
)


summary(puma5)
table(puma5$detection_history)
#0     1
#3000  293

#preparing covariates
nights<-scale(prey_covs$Days_Operable)
trail<-as.factor(prey_covs$Trail)
agou <- scale(prey_covs$AGOU_captr)
acou <- scale(prey_covs$ACOU_captr)
paca <- scale(prey_covs$PACA_captr)
brab <- scale(prey_covs$BRAB_captr)
opos <- scale(prey_covs$OPOS_captr)
sarm <- scale(prey_covs$SARM_captr)

summary(opos)
# Subset continuous covariates
continuous_covariates <- prey_covs[, c("AGOU_captr","ACOU_captr","PACA_captr", "BRAB_captr","OPOS_captr", "SARM_captr"  )]

# Correlation matrix
cor_matrix <- cor(continuous_covariates, use = "complete.obs")
print("Correlation Matrix:")
print(cor_matrix)

#make it work faster
Nthreads <- detectCores() 
setNumThreads(Nthreads)

UFO_pumaPrey <- unmarkedFrameOccu(y = puma5$detection_history, siteCovs = data.frame(nights, trail, agou, acou, paca, brab, opos, sarm))
summary(UFO_pumaPrey)

#Dredge the global model for all possible combinations
pumaPrey_dredge <- occu(~nights+trail~ agou+acou+paca+brab+opos+sarm, UFO_pumaPrey) 

#Count the number of dredged models
pumaPrey_dredged <- dredge(pumaPrey_dredge, rank = "AIC") 

totalmods_pumaPrey <- nrow(pumaPrey_dredged)
print(totalmods_pumaPrey)

#Keep top models (within 2 deltaAIC) & review the top model
print(pumaPrey_dredged[1:5,]) 

pumaPrey_modsaic <- get.models(pumaPrey_dredged, subset = delta < 2,)
pumaPrey_modsaic
pumaPrey_modsaic[[1]]

avgmodpumaPrey<-model.avg(pumaPrey_modsaic, fit=T)
summary(avgmodpumaPrey) 

write.csv(avgmodpumaPrey, file = "results/occu_results.csv") #needs fixing

pumaPrey_topmod <- occu(formula = ~nights + trail + 1 ~ agou + brab + 1, data = UFO_pumaPrey)
pumaPreyavocc <- predict(pumaPrey_topmod, type = 'state', newdata = UFO_pumaPrey)
pumaPreyavdet <- predict(pumaPrey_topmod, type = 'det', newdata = UFO_pumaPrey)
mean(pumaPreyavocc$Predicted)
#0.6207147
mean(pumaPreyavdet$Predicted)
#0.145357

#ACTIVITY PATTERNS

prey_covs <- read.csv("covs_prey.csv", header = TRUE)
pumaPrey_activity <- read.csv("activity.csv", header = TRUE)
puma<-pumaPrey_activity[pumaPrey_activity$Code == "PUMA",]
agou<-pumaPrey_activity[pumaPrey_activity$Code == "AGOU",]
acou<-pumaPrey_activity[pumaPrey_activity$Code == "ACOU",]
brab<-pumaPrey_activity[pumaPrey_activity$Code == "BRAB",]
opos<-pumaPrey_activity[pumaPrey_activity$Code == "OPOS",]
paca<-pumaPrey_activity[pumaPrey_activity$Code == "PACA",]
sarm<-pumaPrey_activity[pumaPrey_activity$Code == "SARM",]
summary(agou)


pumaPrey_activity<-  pumaPrey_activity %>% 
  rename(
    Species = Code)

PUMA<-  puma %>% 
  rename(
    Species = Code)

AGOU<-  agou %>% 
  rename(
    Species = Code)

ACOU<-  acou %>% 
  rename(
    Species = Code)

BRAB<-  brab %>% 
  rename(
    Species = Code)

OPOS<-  opos %>% 
  rename(
    Species = Code)

SARM<-  sarm %>% 
  rename(
    Species = Code)

PACA<-  paca %>% 
  rename(
    Species = Code)

PUMAact<- "PUMA"
AGOUact<- "AGOU"
ACOUact<- "ACOU"
BRABact<- "BRAB"
OPOSact<- "OPOS"
SARMact<- "SARM"
PACAct <- "PACA"

##Puma 
#activity map

#dropping unmatched stations
pumaPrey_activity <- pumaPrey_activity[pumaPrey_activity$Station %in% prey_covs$Station, ]

puma_map <- detectionMaps(CTtable      = prey_covs,
                     recordTable   = pumaPrey_activity,
                     Xcol          = "utm_x",
                     Ycol          = "utm_y",
                     stationCol    = "Station",
                     speciesCol    = "Species",
                     speciesToShow = "PUMA",     
                     printLabels   = FALSE,
                     richnessPlot  = FALSE,     
                     speciesPlots  = TRUE,      
                     addLegend     = TRUE
)

activityDensity(recordTable = pumaPrey_activity,
                recordDateTimeCol    = "DateTime",
                species     = "PUMA",
                writePNG          = TRUE,
                plotDirectory    = "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Peru/puma_project",
)

activityHistogram(recordTable = pumaPrey_activity,
                  recordDateTimeCol    = "DateTime",
                  species     = "PUMA",
                  writePNG          = TRUE,
                  plotDirectory    = "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Peru/puma_project",
)

activityRadial(recordTable       = pumaPrey_activity,
               species           = "PUMA",
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTime",
               plotR             = TRUE,
               writePNG          = TRUE,
               plotDirectory     = "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Peru/puma_project",
               lwd               = 3,
               rp.type           = "p",     # plot type = polygon
               poly.col          = gray(0.5, alpha = 0.5),  # optional. remove for no fill 
)

#Plots for Agouti
activityDensity(recordTable = pumaPrey_activity,
                recordDateTimeCol = "DateTime",
                species     = AGOUact,
                writePNG          = FALSE,
)
activityHistogram (recordTable = pumaPrey_activity,
                   recordDateTimeCol = "DateTime",
                   species     = AGOUact,
                   writePNG          = FALSE,
)
activityRadial(recordTable       = pumaPrey_activity,
               species           = AGOUact,
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTime",
               plotR             = TRUE,
               lwd               = 3,
               rp.type           = "p",     # plot type = polygon
               poly.col          = gray(0.5, alpha = 0.5),  # optional. remove for no fill 
)

#to make sure DateTime is in the correct format (i.e., as.POSIXct) and not characters
pumaPrey_activity$DateTime <- as.POSIXct(pumaPrey_activity$DateTime,
                                         format = "%Y-%m-%d %H:%M:%S")

#extract time (HH:MM:SS) as a character
pumaPrey_activity$TimeOnly <- format(pumaPrey_activity$DateTime, "%H:%M:%S")

#convert to radians
timeradians <- hms2rad(h = pumaPrey_activity$TimeOnly)

timeradians <- hms2rad(h = pumaPrey_activity$Time)
puma.tr <- timeradians[pumaPrey_activity$Species == 'PUMA']
agou.tr<- timeradians[pumaPrey_activity$Species == 'AGOU']
acou.tr<- timeradians[pumaPrey_activity$Species == 'ACOU']
brab.tr<- timeradians[pumaPrey_activity$Species == 'BRAB']
opos.tr<- timeradians[pumaPrey_activity$Species == 'OPOS']
sarm.tr<- timeradians[pumaPrey_activity$Species == 'SARM']
paca.tr<- timeradians[pumaPrey_activity$Species == 'PACA']

length(puma.tr)
length(agou.tr)
length(acou.tr)
length(brab.tr)
length(opos.tr)
length(sarm.tr)
length(paca.tr)

#puma vs agouti
overlapPlot(puma.tr, agou.tr, 
            xscale = 24, 
            xcenter = c("noon"),
            main = expression(""),
            linecol = c("black", "darkolivegreen"),
            olapcol = "antiquewhite3")
legend( 'bottomright', c("Pumas (n=404)", "Agoutis (n=2101)"), 
        lty=c(1, 2), col=c("black", "darkolivegreen"), bty='o', cex = 0.7)

overlap.pagou<-overlapEst(puma.tr, agou.tr,type="Dhat4")
overlap.pagou
# Do 999 smoothed bootstrap values:
bs.pagou <- bootstrap(puma.tr, agou.tr, 999, type="Dhat4")
mean(bs.pagou)
#0.4343352
bootCI(overlap.pagou, bs.pagou)['norm0', ]
#0.3676374 & 0.4539507 


#puma vs acouchi
overlapPlot(puma.tr, acou.tr, 
            xscale = 24, 
            xcenter = "noon",
            main = "Puma vs Acouchi",
            linecol = c("black", "orange3"),
            olapcol = "antiquewhite3")
legend('bottomright', c( "Pumas (n=404)", "Acouchi (n=384)"),
       lty = c(1, 2), col = c("black", "orange3"), bty = 'o', cex = 0.7)

overlap.pacou <- overlapEst(puma.tr, acou.tr, type="Dhat4")
bs.pacou <- bootstrap(puma.tr, acou.tr, 999, type="Dhat4")
mean(bs.pacou)
#0.5408821
bootCI(overlap.pacou, bs.pacou)['norm0', ]
#0.4444102 0.5427711 


#puma vs brazilian rabbit
overlapPlot(puma.tr, brab.tr, 
            xscale = 24, 
            xcenter = "noon",
            main = "Puma vs Brazilian Rabbit",
            linecol = c("black", "darkred"),
            olapcol = "antiquewhite3")
legend('bottomright', c( "Pumas (n=404)", "Brazilian Rabbit (n=201)"),
       lty = c(1, 2), col = c("black", "darkred"), bty = 'o', cex = 0.7)

overlap.pbrab <- overlapEst(puma.tr, brab.tr, type="Dhat4")
bs.pbrab <- bootstrap(puma.tr, brab.tr, 999, type="Dhat4")
mean(bs.pbrab)
#0.7823108
bootCI(overlap.pbrab, bs.pbrab)['norm0', ]
#0.7192357 0.8320992 


#puma vs opossum
overlapPlot(puma.tr, opos.tr, 
            xscale = 24, 
            xcenter = "noon",
            main = "Puma vs Opossum",
            linecol = c("black", "darkgreen"),
            olapcol = "antiquewhite3")
legend('bottomright', c( "Pumas (n=404)", "Opossum (n=885)"),
       lty = c(1, 2), col = c("black", "darkgreen"), bty = 'o', cex = 0.7)

overlap.popos <- overlapEst(puma.tr, opos.tr, type="Dhat4")
bs.popos <- bootstrap(puma.tr, opos.tr, 999, type="Dhat4")
mean(bs.popos)
#0.7299953
bootCI(overlap.popos, bs.popos)['norm0', ]
#0.6808597 0.7680517 

#puma vs small armadillo
overlapPlot(puma.tr, sarm.tr, 
            xscale = 24, 
            xcenter = "noon",
            main = "Puma vs Armadillo",
            linecol = c("black", "brown"),
            olapcol = "antiquewhite3")
legend('bottomright', c( "Pumas (n=404)", "Small Armadillo (n=428)"),
       lty = c(1, 2), col = c("black", "brown"), bty = 'o', cex = 0.7)

overlap.psarm <- overlapEst(puma.tr, sarm.tr, type="Dhat4")
bs.psarm <- bootstrap(puma.tr, sarm.tr, 999, type="Dhat4")
mean(bs.psarm)
#0.7352699
bootCI(overlap.psarm, bs.psarm)['norm0', ]
#0.6689982 0.7660586 


#puma vs paca
overlapPlot(puma.tr, paca.tr, 
            xscale = 24, 
            xcenter = "noon",
            main = "Puma vs Paca",
            linecol = c("black", "darkblue"),
            olapcol = "antiquewhite3")
legend('bottomright', c( "Pumas (n=404)", "Paca (n=1047)"),
       lty = c(1, 2), col = c("black", "darkblue"), bty = 'o', cex = 0.7)

overlap.ppaca <- overlapEst(puma.tr, paca.tr, type="Dhat4")
bs.ppaca <- bootstrap(puma.tr, paca.tr, 999, type="Dhat4")
mean(bs.ppaca)
#0.6976537
bootCI(overlap.ppaca, bs.ppaca)['norm0', ]
#0.6455619 0.7340520 

