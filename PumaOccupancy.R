## PUMA OCCUPANCY
setwd("/Users/clemenslukasser/Documents/Presentations & Sessions/Puma Occupancy student session")

install.packages("camtrapR")
install.packages("AICcmodavg")
install.packages("unmarked")
install.packages("MuMIn")


library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(camtrapR)

#Covariates correlation matrix
covs <- read.csv("covs_puma.csv", header = TRUE)
head(covs)


#cam data
data<-read.csv("cam_puma.csv", header = TRUE)
nrow(data)
datacameraOperation <- cameraOperation(covs, 
                                       stationCol = "Station", 
                                       setupCol = "Setup_Date", 
                                       retrievalCol = "Retrieval_Date",
                                       hasProblems = FALSE,
                                       allCamsOn = FALSE,
                                       camerasIndependent = FALSE,
                                       dateFormat = "%Y-%m-%d", 
                                       writecsv = FALSE, 
)

pumas5 <- detectionHistory(recordTable       = data,
                           camOp                = datacameraOperation,
                           stationCol           = "Station",
                           speciesCol           = "Code",
                           recordDateTimeCol    = "DateTimeOriginal",
                           species              = "PUMA",
                           occasionLength       = 5,
                           day1                 = "station",
                           datesAsOccasionNames = FALSE,
                           includeEffort        = FALSE,
                           scaleEffort          = FALSE,
                           timeZone             = "America/Lima"
)
summary(pumas5)
table(pumas5$detection_history)

#preparing covariates
grid <- as.factor(covs$Grid)
nights<-scale(covs$Days_Operable)
trail<-as.factor(covs$Trail)
hab<-as.factor(covs$Habitat)
season<-as.factor(covs$Season)
river<-scale(covs$RIVER)
NDVI<-scale(covs$NDVI_noriver_500)
dist.per <- scale(covs$DIST_500....)
dist.dist <- scale(covs$DIST_DIST)
prey <- scale(covs$PREY)
jag <- scale(covs$JAGU)

# Subset continuous covariates
continuous_covariates <- covs[, c("Days_Operable", "RIVER", "NDVI_noriver_500", "DIST_500....", "DIST_DIST", "PREY", "JAGU" )]

# Correlation matrix
cor_matrix <- cor(continuous_covariates, use = "complete.obs")
print("Correlation Matrix:")
print(cor_matrix)

#make it work faster
library(parallel)
Nthreads <- detectCores() 
setNumThreads(Nthreads)

UFO_puma <- unmarkedFrameOccu(y = pumas5$detection_history, siteCovs = data.frame(grid, nights, trail, hab, season, river, NDVI, dist.per, dist.dist, prey, jag))

summary(UFO_puma)

puma_dredge <- occu(~grid+nights+trail+season~ hab+river+season+NDVI+dist.per+dist.dist+prey+jag, UFO_puma) #Dredge the global model for all possible combinations

puma_dredged <- dredge(puma_dredge, rank = "AIC") #Count the number of dredged models

totalmods_puma <- nrow(puma_dredged)
print(totalmods_puma)
print(puma_dredged[1:5,]) #Keep top models (within 2 deltaAIC) & review the top model

puma_modsaic <- get.models(puma_dredged, subset = delta < 2,)
puma_modsaic
puma_modsaic[[1]]

avgmodpuma<-model.avg(puma_modsaic, fit=T)
summary(avgmodpuma) 

write.csv(avgmodpuma, file = "results/occu_results.csv")

puma_topmod <- occu(formula = ~nights + trail + 1 ~ hab + jag + prey + river + 1, data = UFO_puma)
pumaavocc <- predict(puma_topmod, type = 'state', newdata = UFO_puma)
pumaavdet <- predict(puma_topmod, type = 'det', newdata = UFO_puma)
mean(pumaavocc$Predicted)
#0.5962958
mean(pumaavdet$Predicted)
#0.1356275







##Activity pattern - if there's time :)

install.packages("tidyverse")
install.packages("dplyr")
install.packages("overlap")
install.packages("astroFns")
library(astroFns)
library(overlap)
library(dplyr)
library(tidyverse)
install.packages("devtools")
devtools::install_github("r-lib/conflicted")
library(conflicted)
library(dplyr)
library(camtrapR)

covs <- read.csv("covs_puma.csv", header = TRUE)
activity <- read.csv("activity.csv", header = TRUE)
Puma<-activity[activity$Code == "PUMA",]
Prey<-activity[activity$Code == "PREY",]
Jag<-activity[activity$Code == "JAGU",]
summary(Jag)
activity<-  activity %>% 
  rename(
    Species = Code)

Puma<-  Puma %>% 
  rename(
    Species = Code)

Prey<-  Prey %>% 
  rename(
    Species = Code)
Jag<-  Jag %>% 
  rename(
    Species = Code)

Pumact<- "PUMA"
Preyact<- "PREY"
Jagact <- "JAGU"

##Puma 
#activity map
map <- detectionMaps(CTtable      = covs,
                     recordTable   = activity,
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
activityDensity(recordTable = activity,
                species     = Pumact,
                writePNG          = TRUE,
                plotDirectory    = "/Users/clemenslukasser/Documents/research/Pumas/code/data/activity",
)
activityHistogram(recordTable = activity,
                  species     = Pumact,
                  writePNG          = TRUE,
                  plotDirectory    = "/Users/clemenslukasser/Documents/research/Pumas/code/data/activity",
)
activityRadial(recordTable       = activity,
               species           = Pumact,
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTimeOriginal",
               plotR             = TRUE,
               writePNG          = TRUE,
               plotDirectory    = "/Users/clemenslukasser/Documents/research/Pumas/code/data/activity",
               lwd               = 3,
               rp.type           = "p",     # plot type = polygon
               poly.col          = gray(0.5, alpha = 0.5),  # optional. remove for no fill 
)

#Prey
activityDensity(recordTable = activity,
                species     = Preyact,
                writePNG          = FALSE,
)
activityHistogram (recordTable = activity,
                   species     = Preyact,
                   writePNG          = FALSE,
)
activityRadial(recordTable       = activity,
               species           = Preyact,
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTimeOriginal",
               plotR             = TRUE,
               lwd               = 3,
               rp.type           = "p",     # plot type = polygon
               poly.col          = gray(0.5, alpha = 0.5),  # optional. remove for no fill 
)

timeradians <- hms2rad(h = activity$Time)
puma.tr <- timeradians[activity$Species == 'PUMA']
prey.tr<- timeradians[activity$Species == 'PREY']
jag.tr <- timeradians[activity$Species == 'JAGU']

#puma and prey overlap
overlapPlot(puma.tr, prey.tr, 
            xscale = 24, 
            xcenter = c("noon"),
            main = expression(""),
            linecol = c("black", "darkolivegreen"),
            olapcol = "antiquewhite3")
legend( 'bottomright', c("Pumas (n=399)", "Prey (n=5048)"), 
        lty=c(1, 2), col=c("black", "darkolivegreen"), bty='o', cex = 0.7)

length(puma.tr)
length(prey.tr)

overlap.pp<-overlapEst(puma.tr, prey.tr,type="Dhat4")
# Do 999 smoothed bootstrap values:
bs.pp <- bootstrap(puma.tr, prey.tr, 999, type="Dhat4")
mean(bs.pp)
#0.8800148
bootCI(overlap.pp, bs.pp)['norm0', ]
#0.8348561 0.9131086 

#puma and jag overlap
overlapPlot(puma.tr, jag.tr, 
            xscale = 24, 
            xcenter = c("noon"),
            main = expression(""),
            linecol = c("black", "brown2"),
            olapcol = "antiquewhite3")
legend( 'bottomright', c("Pumas (n=399)", "Jaguars (n=183)"), 
        lty=c(1, 2), col=c("black", "brown2"), bty='o', cex = 0.7)

length(puma.tr)
length(jag.tr)

overlap.pj<-overlapEst(puma.tr, jag.tr,type="Dhat4")
# Do 999 smoothed bootstrap values:
bs.pj <- bootstrap(puma.tr, jag.tr, 999, type="Dhat4")
mean(bs.pj)
#0.8302345
bootCI(overlap.pj, bs.pj)['norm0', ]
#0.7754905 0.9081660  