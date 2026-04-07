setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(camtrapR)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(dplyr)
library(lubridate)
library(parallel)
library(ggplot2)

#load covariates
covs <- read.csv("covariates.csv", header = TRUE)
jags <- read.csv("jagu.csv", header = TRUE)
prey <- read.csv("all_prey.csv", header = TRUE)

#count detections per prey species
prey_counts <- prey %>%
  count(Code) %>%        
  filter(n > 100)

prey_counts

prey_keep <- prey_counts$Code
length(prey_keep)  

#adjusting time and date for prey
prey$Hour <- ifelse(prey$AM_PM == "PM" & prey$Hour != 12, prey$Hour + 12, prey$Hour)
prey$Hour <- ifelse(prey$AM_PM == "AM" & prey$Hour == 12, 0, prey$Hour)

prey$DateTimeOriginal <- as.POSIXct(
  paste(prey$Year, prey$Month, prey$Day, prey$Hour, prey$Min),
  format = "%Y %m %d %H %M",
  tz = "America/Lima"
)

#adjusting time and date for jag
jags$Hour <- ifelse(jags$AM_PM == "PM" & jags$Hour != 12, jags$Hour + 12, jags$Hour)
jags$Hour <- ifelse(jags$AM_PM == "AM" & jags$Hour == 12, 0, jags$Hour)

jags$DateTimeOriginal <- as.POSIXct(
  paste(jags$Year, jags$Month, jags$Day, jags$Hour, jags$Min),
  format = "%Y %m %d %H %M",
  tz = "America/Lima"
)


#camera data
datacameraOperation <- cameraOperation(
  CTtable            = covs,
  stationCol         = "Station",
  setupCol           = "Setup_Date",
  retrievalCol       = "Retrieval_Date",
  dateFormat         = "%Y-%m-%d",
  hasProblems        = FALSE,
  allCamsOn          = FALSE,
  camerasIndependent = FALSE,
  writecsv           = FALSE
)

#detection history for jaguar
bad_rowsj <- c(33)
jags_data <- jags[-bad_rowsj, ]
nrow(jags_data)

jag_hist <- detectionHistory(
  recordTable        = jags_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "JAGU",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#detection history for all prey
#AGOU
bad_rows_agou <- c(3336, 3340, 3343, 3344, 3347, 3348, 3351, 3356, 3358, 3361, 3365, 3370, 3371, 3375, 3379, 3386, 3387, 3388, 3390, 3514, 3515, 3517, 3519, 3520, 3523, 3529, 3533, 3535, 3546, 3550, 3554, 3570, 3575, 3576, 3579, 3581, 3584, 3588, 3590, 3593, 3595, 3599, 3601, 3602, 3605, 3606, 3608)
agou_data <- prey[-bad_rows_agou, ]
nrow(agou_data)

agou_hist <- detectionHistory(
  recordTable        = agou_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "AGOU",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#ACOU
bad_rows_acou <- c(3341)
acou_data <- prey[-bad_rows_acou, ]
nrow(acou_data)

acou_hist <- detectionHistory(
  recordTable        = acou_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "ACOU",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#BRAB
bad_rows_brab <- c(3428,8161)
brab_data <- prey[-bad_rows_brab, ]
nrow(brab_data)

brab_hist <- detectionHistory(
  recordTable        = brab_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "BRAB",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#BROC
#subset species
broc_data <- prey[prey$Code == "BROC", ]

#clean station and camera operation names
broc_data$Station <- trimws(broc_data$Station)
rownames(datacameraOperation) <- trimws(rownames(datacameraOperation))

#remove bad stations
broc_data <- broc_data[
  !broc_data$Station %in% c("BG38","BG48","HOJA44"),
]

#remove missing timestamps
broc_data <- broc_data[
  !is.na(broc_data$DateTimeOriginal),
]

#run detection history
broc_hist <- detectionHistory(
  recordTable       = broc_data,
  camOp             = datacameraOperation,
  stationCol        = "Station",
  speciesCol        = "Code",
  recordDateTimeCol = "DateTimeOriginal",
  species           = "BROC",
  occasionLength    = 5,
  day1              = "station",
  includeEffort     = FALSE
)

#CPEC
cpec_data <- prey[prey$Code == "CPEC", ]

#remove bad stations
cpec_data <- cpec_data[
  !cpec_data$Station %in% c("HOJA44","HN2143","BG38","BG48"),
]

cpec_hist <- detectionHistory(
  recordTable        = cpec_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "CPEC",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#GANT
gant_hist <- detectionHistory(
  recordTable        = prey,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "GANT",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#OPOS
bad_rows_opos <- c(3402, 3470, 3525, 3527, 3531, 3534, 3553, 3571)
opos_data <- prey[-bad_rows_opos, ]
nrow(opos_data)

opos_hist <- detectionHistory(
  recordTable        = opos_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "OPOS",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#PACA
bad_rows_paca <- c(1461, 3337, 3342, 3355, 3363, 3376, 3378, 3405, 3521, 5679, 8040, 13539)
paca_data <- prey[-bad_rows_paca, ]
nrow(paca_data)

paca_hist <- detectionHistory(
  recordTable        = paca_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "PACA",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#SARM - review
sarm_data <- prey %>%
  dplyr::filter(Code == "SARM") %>%
  dplyr::mutate(
    DateTimeOriginal = ymd_hms(DateTimeOriginal, tz = "America/Lima")
  ) %>%
  dplyr::filter(!is.na(DateTimeOriginal))


sarm_hist <- detectionHistory(
  recordTable        = sarm_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "SARM",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#TAPI
bad_rows_tapi <- c(1347, 3338, 3350, 3362, 3364, 3367, 3374, 3377, 3385, 3391, 3403, 3407, 3528, 3536, 3537, 3540, 3545, 3551, 3555, 3558, 3560, 3563, 3573, 3577)
tapi_data <- prey[-bad_rows_tapi, ]
nrow(tapi_data)

tapi_hist <- detectionHistory(
  recordTable        = tapi_data,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "TAPI",
  occasionLength     = 5,
  day1               = "station",
  includeEffort      = FALSE
)

#combine for occuMulti
agou_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  agou   = as.matrix(agou_hist$detection_history)
)

#unmarkedFrameOccuMulti
agou_umf <- unmarkedFrameOccuMulti(
  y = agou_ylist
)

#formulas
agou_stateformulas <- c(
  "~1",  
  "~1",  
  "~1"   
)

agou_detformulas <- c("~1", "~1")

mod_jag_agou <- occuMulti(
  detformulas   = agou_detformulas,
  stateformulas = agou_stateformulas,
  data          = agou_umf
)

summary(mod_jag_agou)

#jaguar
jag_marginal_agou <- predict(mod_jag_agou, type = "state", species = "jaguar")
#agouti
agou_marginal <- predict(mod_jag_agou, type = "state", species = "agou")

jag_marginal_agou
agou_marginal

#plotting marginal occupancy 
agou_marginal <- rbind(
  jag_marginal_agou[1, ],
  agou_marginal[1, ]
)
agou_marginal$Species <- c("Jaguar", "Agouti")

plot(
  1:2,
  agou_marginal$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "",
  ylab = "Marginal occupancy (95% CI)"
)
axis(1, at = 1:2, labels = agou_marginal$Species)

top <- 0.05
for (i in 1:2) {
  segments(i, agou_marginal$lower[i], i, agou_marginal$upper[i])
  segments(i - top, agou_marginal$lower[i], i + top)
  segments(i - top, agou_marginal$upper[i], i + top)
}

#conditional occupancy
jag_given_agou <- predict(mod_jag_agou, type = "state", species = "jaguar", cond = "agou")
jag_no_agou <- predict(mod_jag_agou, type = "state", species = "jaguar", cond = "-agou")

jag_given_agou
jag_no_agou

#agou occupancy
agou_given_jag <- predict(mod_jag_agou, type = "state", species = "agou", cond = "jaguar")
agou_no_jag <- predict(mod_jag_agou, type = "state", species = "agou", cond = "-jaguar")

agou_given_jag
agou_no_jag

#plotting conditional occupancy
cond_agou <- rbind(
  jag_given_agou[1, ],
  jag_no_agou[1, ]
)

cond_agou$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_agou$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Agouti status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_agou$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_agou$lower[i], i, cond_agou$upper[i])
  segments(i - top, cond_agou$lower[i], i + top)
  segments(i - top, cond_agou$upper[i], i + top)
}

#plotting agou conditional occupancy
cond_agou_jag <- rbind(
  agou_given_jag[1, ],
  agou_no_jag[1, ]
)

cond_agou_jag$agouti_status <- c("Present", "Absent")

plot(
  1:2,
  cond_agou_jag$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Jaguar status",
  ylab = "Agouti occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_agou_jag$agouti_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_agou_jag$lower[i], i, cond_agou_jag$upper[i])
  segments(i - top, cond_agou_jag$lower[i], i + top)
  segments(i - top, cond_agou_jag$upper[i], i + top)
}

#jaguar and acouchi
acou_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  acou   = as.matrix(acou_hist$detection_history)
)

acou_umf <- unmarkedFrameOccuMulti(y = acou_ylist)
acou_stateformulas <- c("~1", "~1", "~1")  
acou_detformulas   <- c("~1", "~1")        
mod_jag_acou <- occuMulti(
  detformulas   = acou_detformulas,
  stateformulas = acou_stateformulas,
  data          = acou_umf
)

jag_marginal_acou <- predict(mod_jag_acou, type = "state", species = "jaguar")
acou_marginal      <- predict(mod_jag_acou, type = "state", species = "acou")

#plotting marginal occupancy 
acou_marginal <- rbind(
  jag_marginal_acou[1, ],
  acou_marginal[1, ]
)
acou_marginal$Species <- c("Jaguar", "Acouchi")

plot(
  1:2,
  acou_marginal$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "",
  ylab = "Marginal occupancy (95% CI)"
)
axis(1, at = 1:2, labels = acou_marginal$Species)

top <- 0.05
for (i in 1:2) {
  segments(i, acou_marginal$lower[i], i, acou_marginal$upper[i])
  segments(i - top, acou_marginal$lower[i], i + top)
  segments(i - top, acou_marginal$upper[i], i + top)
}

jag_given_acou <- predict(mod_jag_acou, type = "state", species = "jaguar", cond = "acou")
jag_no_acou    <- predict(mod_jag_acou, type = "state", species = "jaguar", cond = "-acou")

#plotting conditional occupancy
cond_acou <- rbind(
  jag_given_acou[1, ],
  jag_no_acou[1, ]
)

cond_acou$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_acou$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Acouchi status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_agou$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_acou$lower[i], i, cond_acou$upper[i])
  segments(i - top, cond_acou$lower[i], i + top)
  segments(i - top, cond_acou$upper[i], i + top)
}

#jaguar and brab
brab_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  brab   = as.matrix(brab_hist$detection_history)
)

brab_umf <- unmarkedFrameOccuMulti(y = brab_ylist)

brab_stateformulas <- c("~1", "~1","~1")
brab_detformulas   <- c("~1", "~1")

mod_jag_brab <- occuMulti(
  detformulas   = brab_detformulas,
  stateformulas = brab_stateformulas,
  data          = brab_umf
)
mod_jag_brab

#marginal
jag_marginal_brab <- predict(mod_jag_brab, type = "state", species = "jaguar")
brab_marginal     <- predict(mod_jag_brab, type = "state", species = "brab")
jag_marginal_brab
brab_marginal

#conditional
jag_given_brab <- predict(mod_jag_brab, type = "state", species = "jaguar", cond = "brab")
jag_no_brab    <- predict(mod_jag_brab, type = "state", species = "jaguar", cond = "-brab")

#plotting conditional occupancy
cond_brab <- rbind(
  jag_given_brab[1, ],
  jag_no_brab[1, ]
)

cond_brab$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_brab$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Brazilian Rabbit status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_brab$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_brab$lower[i], i, cond_brab$upper[i])
  segments(i - top, cond_brab$lower[i], i + top)
  segments(i - top, cond_brab$upper[i], i + top)
}

#jaguar and brocket deer
broc_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  broc   = as.matrix(broc_hist$detection_history)
)

broc_umf <- unmarkedFrameOccuMulti(y = broc_ylist)

broc_stateformulas <- c("~1", "~1","~1")
broc_detformulas   <- c("~1", "~1")

mod_jag_broc <- occuMulti(
  detformulas   = broc_detformulas,
  stateformulas = broc_stateformulas,
  data          = broc_umf
)

#marginal
jag_marginal_broc <- predict(mod_jag_broc, type = "state", species = "jaguar")
brab_marginal     <- predict(mod_jag_broc, type = "state", species = "broc")

#conditional
jag_given_broc <- predict(mod_jag_broc, type = "state", species = "jaguar", cond = "broc")
jag_no_broc<- predict(mod_jag_broc, type = "state", species = "jaguar", cond = "-broc")

#plotting conditional occupancy
cond_broc <- rbind(
  jag_given_broc[1, ],
  jag_no_broc[1, ]
)

cond_broc$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_broc$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Brocket deer status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_broc$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_broc$lower[i], i, cond_broc$upper[i])
  segments(i - top, cond_broc$lower[i], i + top)
  segments(i - top, cond_broc$upper[i], i + top)
}

#jaguar and cpec
cpec_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  cpec   = as.matrix(cpec_hist$detection_history)
)

cpec_umf <- unmarkedFrameOccuMulti(y = cpec_ylist)
cpec_stateformulas <- c("~1", "~1", "~1")
cpec_detformulas   <- c("~1", "~1")

mod_jag_cpec <- occuMulti(
  detformulas   = cpec_detformulas,
  stateformulas = cpec_stateformulas,
  data          = cpec_umf
)

jag_marginal_cpec <- predict(mod_jag_cpec, type = "state", species = "jaguar")
cpec_marginal     <- predict(mod_jag_cpec, type = "state", species = "cpec")

#conditional occupancy 
jag_given_cpec    <- predict(mod_jag_cpec, type = "state", species = "jaguar", cond = "cpec")
jag_no_cpec       <- predict(mod_jag_cpec, type = "state", species = "jaguar", cond = "-cpec")

#plotting conditional occupancy
cond_cpec <- rbind(
  jag_given_cpec[1, ],
  jag_no_cpec[1, ]
)

cond_cpec$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_cpec$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Collared Peccary status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_cpec$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_cpec$lower[i], i, cond_cpec$upper[i])
  segments(i - top, cond_cpec$lower[i], i + top)
  segments(i - top, cond_cpec$upper[i], i + top)
}

#jaguar and gant
gant_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  gant   = as.matrix(gant_hist$detection_history)
)

gant_umf <- unmarkedFrameOccuMulti(y = gant_ylist)
gant_stateformulas <- c("~1", "~1", "~1")
gant_detformulas   <- c("~1", "~1")

mod_jag_gant <- occuMulti(
  detformulas   = gant_detformulas,
  stateformulas = gant_stateformulas,
  data          = gant_umf
)

jag_marginal_gant <- predict(mod_jag_gant, type = "state", species = "jaguar")
gant_marginal     <- predict(mod_jag_gant, type = "state", species = "gant")

jag_given_gant    <- predict(mod_jag_gant, type = "state", species = "jaguar", cond = "gant")
jag_no_gant       <- predict(mod_jag_gant, type = "state", species = "jaguar", cond = "-gant")

#plotting conditional occupancy
cond_gant <- rbind(
  jag_given_gant[1, ],
  jag_no_gant[1, ]
)

cond_gant$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_gant$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Giant anteater status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_gant$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_gant$lower[i], i, cond_gant$upper[i])
  segments(i - top, cond_gant$lower[i], i + top)
  segments(i - top, cond_gant$upper[i], i + top)
}


#jaguar and opossum
opos_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  opos   = as.matrix(opos_hist$detection_history)
)

opos_umf <- unmarkedFrameOccuMulti(y = opos_ylist)
opos_stateformulas <- c("~1", "~1", "~1")
opos_detformulas   <- c("~1", "~1")

mod_jag_opos <- occuMulti(
  detformulas   = opos_detformulas,
  stateformulas = opos_stateformulas,
  data          = opos_umf
)

jag_marginal_opos <- predict(mod_jag_opos, type = "state", species = "jaguar")
opos_marginal     <- predict(mod_jag_opos, type = "state", species = "opos")

jag_given_opos    <- predict(mod_jag_opos, type = "state", species = "jaguar", cond = "opos")
jag_no_opos       <- predict(mod_jag_opos, type = "state", species = "jaguar", cond = "-opos")

#plotting conditional occupancy
cond_opos <- rbind(
  jag_given_opos[1, ],
  jag_no_opos[1, ]
)

cond_opos$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_opos$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Opossom status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_opos$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_opos$lower[i], i, cond_opos$upper[i])
  segments(i - top, cond_opos$lower[i], i + top)
  segments(i - top, cond_opos$upper[i], i + top)
}


#jaguar and paca
paca_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  paca   = as.matrix(paca_hist$detection_history)
)

paca_umf <- unmarkedFrameOccuMulti(y = paca_ylist)
paca_stateformulas <- c("~1", "~1", "~1")
paca_detformulas   <- c("~1", "~1")

mod_jag_paca <- occuMulti(
  detformulas   = paca_detformulas,
  stateformulas = paca_stateformulas,
  data          = paca_umf
)

jag_marginal_paca <- predict(mod_jag_paca, type = "state", species = "jaguar")
paca_marginal     <- predict(mod_jag_paca, type = "state", species = "paca")

jag_given_paca    <- predict(mod_jag_paca, type = "state", species = "jaguar", cond = "paca")
jag_no_paca       <- predict(mod_jag_paca, type = "state", species = "jaguar", cond = "-paca")

#plotting conditional occupancy
cond_paca <- rbind(
  jag_given_paca[1, ],
  jag_no_paca[1, ]
)

cond_paca$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_paca$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Paca status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_paca$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_paca$lower[i], i, cond_paca$upper[i])
  segments(i - top, cond_paca$lower[i], i + top)
  segments(i - top, cond_paca$upper[i], i + top)
}


#jaguar and small armadillo
sarm_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  sarm   = as.matrix(sarm_hist$detection_history)
)

sarm_umf <- unmarkedFrameOccuMulti(y = sarm_ylist)
sarm_stateformulas <- c("~1", "~1", "~1")
sarm_detformulas   <- c("~1", "~1")

mod_jag_sarm <- occuMulti(
  detformulas   = sarm_detformulas,
  stateformulas = sarm_stateformulas,
  data          = sarm_umf
)

jag_marginal_sarm <- predict(mod_jag_sarm, type = "state", species = "jaguar")
sarm_marginal     <- predict(mod_jag_sarm, type = "state", species = "sarm")

jag_given_sarm    <- predict(mod_jag_sarm, type = "state", species = "jaguar", cond = "sarm")
jag_no_sarm       <- predict(mod_jag_sarm, type = "state", species = "jaguar", cond = "-sarm")

#plotting conditional occupancy
cond_sarm <- rbind(
  jag_given_sarm[1, ],
  jag_no_sarm[1, ]
)

cond_sarm$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_sarm$Predicted,
  ylim = c(0,1),
  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Small Armadillo status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_sarm$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_sarm$lower[i], i, cond_sarm$upper[i])
  segments(i - top, cond_sarm$lower[i], i + top)
  segments(i - top, cond_sarm$upper[i], i + top)
}

#jaguar and tapir
tapi_ylist <- list(
  jaguar = as.matrix(jag_hist$detection_history),
  tapi   = as.matrix(tapi_hist$detection_history)
)

tapi_umf <- unmarkedFrameOccuMulti(y = tapi_ylist)
tapi_stateformulas <- c("~1", "~1", "~1")
tapi_detformulas   <- c("~1", "~1")

mod_jag_tapi <- occuMulti(
  detformulas   = tapi_detformulas,
  stateformulas = tapi_stateformulas,
  data          = tapi_umf
)

jag_marginal_tapi<- predict(mod_jag_tapi, type = "state", species = "jaguar")
tapi_marginal<- predict(mod_jag_tapi, type = "state", species = "tapi")

jag_given_tapi<- predict(mod_jag_tapi, type = "state", species = "jaguar", cond = "tapi")
jag_no_tapi<- predict(mod_jag_tapi, type = "state", species = "jaguar", cond = "-tapi")

#plotting conditional occupancy
cond_tapi <- rbind(
  jag_given_tapi[1, ],
  jag_no_tapi[1, ]
)

cond_tapi$Jaguar_status <- c("Present", "Absent")

plot(
  1:2,
  cond_tapi$Predicted,
  ylim = c(0,1),  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "Tapir status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_tapi$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_tapi$lower[i], i, cond_tapi$upper[i])
  segments(i - top, cond_tapi$lower[i], i + top)
  segments(i - top, cond_tapi$upper[i], i + top)
}

#Combined marginal occupancy plotting
#Requires running puma_prey.R first to produce puma_marginal_ * objects

marginal_all <- bind_rows(
  # Puma rows
  data.frame(Species = "puma_AGOU", Predicted = puma_marginal_agou$Predicted[1], lower = puma_marginal_agou$lower[1], upper = puma_marginal_agou$upper[1], Type = "Puma"),
  data.frame(Species = "puma_BRAB", Predicted = puma_marginal_brab$Predicted[1], lower = puma_marginal_brab$lower[1], upper = puma_marginal_brab$upper[1], Type = "Puma"),
  data.frame(Species = "puma_CPEC", Predicted = puma_marginal_cpec$Predicted[1], lower = puma_marginal_cpec$lower[1], upper = puma_marginal_cpec$upper[1], Type = "Puma"),
  data.frame(Species = "puma_GANT", Predicted = puma_marginal_gant$Predicted[1], lower = puma_marginal_gant$lower[1], upper = puma_marginal_gant$upper[1], Type = "Puma"),
  data.frame(Species = "puma_OPOS", Predicted = puma_marginal_opos$Predicted[1], lower = puma_marginal_opos$lower[1], upper = puma_marginal_opos$upper[1], Type = "Puma"),
  data.frame(Species = "puma_PACA", Predicted = puma_marginal_paca$Predicted[1], lower = puma_marginal_paca$lower[1], upper = puma_marginal_paca$upper[1], Type = "Puma"),
  data.frame(Species = "puma_SARM", Predicted = puma_marginal_sarm$Predicted[1], lower = puma_marginal_sarm$lower[1], upper = puma_marginal_sarm$upper[1], Type = "Puma"),
  data.frame(Species = "puma_TAPI", Predicted = puma_marginal_tapi$Predicted[1], lower = puma_marginal_tapi$lower[1], upper = puma_marginal_tapi$upper[1], Type = "Puma"),
  data.frame(Species = "puma_ACOU", Predicted = puma_marginal_acou$Predicted[1], lower = puma_marginal_acou$lower[1], upper = puma_marginal_acou$upper[1], Type = "Puma"),
  
  # Jaguar rows
  data.frame(Species = "Jaguar_AGOU", Predicted = jag_marginal_agou$Predicted[1], lower = jag_marginal_agou$lower[1], upper = jag_marginal_agou$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_BRAB", Predicted = jag_marginal_brab$Predicted[1], lower = jag_marginal_brab$lower[1], upper = jag_marginal_brab$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_CPEC", Predicted = jag_marginal_cpec$Predicted[1], lower = jag_marginal_cpec$lower[1], upper = jag_marginal_cpec$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_GANT", Predicted = jag_marginal_gant$Predicted[1], lower = jag_marginal_gant$lower[1], upper = jag_marginal_gant$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_OPOS", Predicted = jag_marginal_opos$Predicted[1], lower = jag_marginal_opos$lower[1], upper = jag_marginal_opos$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_PACA", Predicted = jag_marginal_paca$Predicted[1], lower = jag_marginal_paca$lower[1], upper = jag_marginal_paca$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_SARM", Predicted = jag_marginal_sarm$Predicted[1], lower = jag_marginal_sarm$lower[1], upper = jag_marginal_sarm$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_TAPI", Predicted = jag_marginal_tapi$Predicted[1], lower = jag_marginal_tapi$lower[1], upper = jag_marginal_tapi$upper[1], Type = "Jaguar"),
  data.frame(Species = "Jaguar_ACOU", Predicted = jag_marginal_acou$Predicted[1], lower = jag_marginal_acou$lower[1], upper = jag_marginal_acou$upper[1], Type = "Jaguar"),
  
  # Prey rows (one per species)
  data.frame(Species = "Agouti",           Predicted = agou_marginal$Predicted[1], lower = agou_marginal$lower[1], upper = agou_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Brazilian rabbit", Predicted = brab_marginal$Predicted[1], lower = brab_marginal$lower[1], upper = brab_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Collared peccary", Predicted = cpec_marginal$Predicted[1], lower = cpec_marginal$lower[1], upper = cpec_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Giant anteater",   Predicted = gant_marginal$Predicted[1], lower = gant_marginal$lower[1], upper = gant_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Opossum",          Predicted = opos_marginal$Predicted[1], lower = opos_marginal$lower[1], upper = opos_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Paca",             Predicted = paca_marginal$Predicted[1], lower = paca_marginal$lower[1], upper = paca_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Small armadillo",  Predicted = sarm_marginal$Predicted[1], lower = sarm_marginal$lower[1], upper = sarm_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Tapir",            Predicted = tapi_marginal$Predicted[1], lower = tapi_marginal$lower[1], upper = tapi_marginal$upper[1], Type = "Prey"),
  data.frame(Species = "Acouchi",          Predicted = acou_marginal$Predicted[1], lower = acou_marginal$lower[1], upper = acou_marginal$upper[1], Type = "Prey")
)

puma_avg <- marginal_all %>%
  filter(Type == "Puma") %>%
  summarise(Species = "Puma", Predicted = mean(Predicted), lower = mean(lower), upper = mean(upper))

jaguar_avg <- marginal_all %>%
  filter(Type == "Jaguar") %>%
  summarise(Species = "Jaguar", Predicted = mean(Predicted), lower = mean(lower), upper = mean(upper))

prey_only <- marginal_all %>%
  filter(Type == "Prey")

prey_order <- prey_only %>%
  arrange(desc(Predicted)) %>%
  pull(Species)

plot_data <- bind_rows(
  puma_avg %>% mutate(Type = "Puma"),
  jaguar_avg %>% mutate(Type = "Jaguar"),
  prey_only
)

plot_data$Species <- factor(plot_data$Species, levels = c("Puma", "Jaguar", prey_order))

ggplot(plot_data, aes(x = Species, y = Predicted, color = Type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Species",
    y = "Marginal occupancy (95% CI)",
    color = "Type"
  ) +
  scale_color_manual(values = c("Puma" = "#A65628", "Prey" = "#2D9E5F", "Jaguar" = "#E6AB02"))

ggsave("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/marginal_occupancy_all.png",
       width = 8, height = 6, dpi = 300, units = "in")

#plotting all conditional occupancy in one figure
library(dplyr)
library(ggplot2)

# Function to extract conditional jaguar occupancy
get_conditional_jag <- function(model, prey_code, prey_label) {
  
  jag_present <- predict(model,
                         type = "state",
                         species = "jaguar",
                         cond = prey_code)
  
  jag_absent  <- predict(model,
                         type = "state",
                         species = "jaguar",
                         cond = paste0("-", prey_code))
  
  df <- rbind(
    jag_present[1, ],
    jag_absent[1, ]
  )
  
  df$PreyStatus <- c("Present", "Absent")
  df$PreySpecies <- prey_label
  
  return(df[, c("PreySpecies","PreyStatus","Predicted","lower","upper")])
}

#building one combined dataframe
all_conditional_jag <- bind_rows(
  
  get_conditional_jag(mod_jag_agou, "agou", "Agouti"),
  get_conditional_jag(mod_jag_acou, "acou", "Acouchi"),
  get_conditional_jag(mod_jag_brab, "brab", "Brazilian Rabbit"),
  get_conditional_jag(mod_jag_broc, "broc", "Brocket Deer"),
  get_conditional_jag(mod_jag_cpec, "cpec", "Collared Peccary"),
  get_conditional_jag(mod_jag_gant, "gant", "Giant Anteater"),
  get_conditional_jag(mod_jag_opos, "opos", "Opossum"),
  get_conditional_jag(mod_jag_paca, "paca", "Paca"),
  get_conditional_jag(mod_jag_sarm, "sarm", "Small Armadillo"),
  get_conditional_jag(mod_jag_tapi, "tapi", "Tapir")
  
)

order_levels <- all_conditional_jag %>%
  filter(PreyStatus == "Present") %>%
  arrange(desc(Predicted)) %>%
  pull(PreySpecies)

all_conditional_jag$PreySpecies <- factor(
  all_conditional_jag$PreySpecies,
  levels = order_levels
)

ggplot(all_conditional_jag,
       aes(x = PreySpecies,
           y = Predicted,
           color = PreyStatus)) +
  
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  
  geom_errorbar(aes(ymin = lower,
                    ymax = upper),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  
  scale_color_manual(values = c("Present" = "#E6AB02",
                                "Absent"  = "#999999"),
                     breaks = c("Present", "Absent")) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  ylab("Jaguar conditional occupancy (95% CI)") +
  xlab("Prey species") +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  )

ggsave("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/conditional_jag_prey.png",
       width = 12, height = 5, dpi = 300, units = "in")

#FITTING MODELS
#agouti
mod_jag_agou <- occuMulti(
  detformulas   = agou_detformulas,
  stateformulas = agou_stateformulas,
  data          = agou_umf
)

summary(mod_jag_agou)

# acou
mod_jag_acou <- occuMulti(
  detformulas   = acou_detformulas,
  stateformulas = acou_stateformulas,
  data          = acou_umf
)
summary(mod_jag_acou)

# brab
mod_jag_brab <- occuMulti(
  detformulas   = brab_detformulas,
  stateformulas = brab_stateformulas,
  data          = brab_umf
)
summary(mod_jag_brab)

# cpec
mod_jag_cpec <- occuMulti(
  detformulas   = cpec_detformulas,
  stateformulas = cpec_stateformulas,
  data          = cpec_umf
)
summary(mod_jag_cpec)

# broc
mod_jag_broc <- occuMulti(
  detformulas   = broc_detformulas,
  stateformulas = broc_stateformulas,
  data          = broc_umf
)
summary(mod_jag_broc)

# sarm
mod_jag_sarm <- occuMulti(
  detformulas   = sarm_detformulas,
  stateformulas = sarm_stateformulas,
  data          = sarm_umf
)
summary(mod_jag_sarm)

# opos
mod_jag_opos <- occuMulti(
  detformulas   = opos_detformulas,
  stateformulas = opos_stateformulas,
  data          = opos_umf
)
summary(mod_jag_opos)

# gant
mod_jag_gant <- occuMulti(
  detformulas   = gant_detformulas,
  stateformulas = gant_stateformulas,
  data          = gant_umf
)
summary(mod_jag_gant)

# tapi
mod_jag_tapi <- occuMulti(
  detformulas   = tapi_detformulas,
  stateformulas = tapi_stateformulas,
  data          = tapi_umf
)
summary(mod_jag_tapi)

# paca
mod_jag_paca <- occuMulti(
  detformulas   = paca_detformulas,
  stateformulas = paca_stateformulas,
  data          = paca_umf
)
summary(mod_jag_paca)

