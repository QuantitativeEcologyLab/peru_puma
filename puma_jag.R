setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(camtrapR)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(dplyr)
library(lubridate)
library(parallel)

#opening the data files
covs <- read.csv("covariates.csv", stringsAsFactors = FALSE)
puma_jag  <- read.csv("puma_jag.csv", stringsAsFactors = FALSE)

head(covs)
head(puma_jag)

#adjusting the correct date and time format
puma_jag$Hour <- ifelse(puma_jag$AM_PM == "PM" & puma_jag$Hour != 12, puma_jag$Hour + 12, puma_jag$Hour)
puma_jag$Hour <- ifelse(puma_jag$AM_PM == "AM" & puma_jag$Hour == 12, 0, puma_jag$Hour)

puma_jag$DateTimeOriginal <- as.POSIXct(
  paste(puma_jag$Year, puma_jag$Month, puma_jag$Day, puma_jag$Hour, puma_jag$Min),
  format = "%Y %m %d %H %M",
  tz = "America/Lima"
)

data <- puma_jag
nrow(data)

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

#puma detection history
bad_rows <- c(76,77,78,80,81,82,83,91)
data_clean <- data[-bad_rows, ]
nrow(data_clean)


pumas <- detectionHistory(
  recordTable        = data_clean,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "PUMA",
  occasionLength     = 5,
  day1               = "station",
  datesAsOccasionNames = FALSE,
  includeEffort      = FALSE,
  scaleEffort        = FALSE,
  timeZone           = "America/Lima"
)

summary(pumas)
table(pumas$detection_history)

#jaguar detection history
jbad_rows <- c(84)
jdata_clean <- data[-jbad_rows, ]

jags <- detectionHistory(
  recordTable        = jdata_clean,
  camOp              = datacameraOperation,
  stationCol         = "Station",
  speciesCol         = "Code",
  recordDateTimeCol  = "DateTimeOriginal",
  species            = "JAGU",
  occasionLength     = 5,
  day1               = "station",
  datesAsOccasionNames = FALSE,
  includeEffort      = FALSE,
  scaleEffort        = FALSE,
  timeZone           = "America/Lima"
)

summary(jags)
table(jags$detection_history)

#check dimensions of the detection data
ylist <- list(
  puma = as.matrix(pumas$detection_history),
  jaguar = as.matrix(jags$detection_history)
)

lapply(ylist, head)

#preparing covariates
grid <- as.factor(covs$Grid)
nights <- scale(covs$Days_Operable)
trail <- as.factor(covs$Trail)

sitecovs <- data.frame(
  grid,
  nights,
  trail
)

#build the unmarkedFrameOccuMulti
umf <- unmarkedFrameOccuMulti(
  y = ylist,
  siteCovs = sitecovs
)

#checking the fDesign
umf@fDesign
colnames(umf@fDesign)

stateformulas <- c(
  "~1",  # puma
  "~1",  # jaguar
  "~1"   # puma:jaguar interaction
)

detformulas <- c( "~1","~1")

#fitting the model
mod_null <- occuMulti(detformulas=detformulas, 
                      stateformulas=stateformulas, 
                      data=umf)
summary(mod_null)

#occupancy state probability
occ_prob <- predict(mod_null, type = "state")
occ_prob

head(occ_prob$Predicted)
rowSums(occ_prob$Predicted)

#marginal occupancy for each of the species
puma_marginal <- predict(
  mod_null,
  type = "state",
  species = "puma"
)

jag_marginal <- predict(
  mod_null,
  type = "state",
  species = "jaguar"
)

head(puma_marginal)
head(jag_marginal)

all_marginal <- rbind(
  puma_marginal[1, ],
  jag_marginal[1, ]
)

all_marginal$Species <- c("Puma", "Jaguar")

all_marginal

#plotting the marginal occupancy
plot(
  1:2,
  all_marginal$Predicted,
  ylim = c(0, 1),
  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  xaxt = "n",
  xlab = "",
  ylab = "Marginal occupancy (95% CI)"
)

axis(1, at = 1:2, labels = all_marginal$Species)

#confidence intervals
top <- 0.05
for (i in 1:2) {
  segments(i, all_marginal$lower[i], i, all_marginal$upper[i])
  segments(i - top, all_marginal$lower[i], i + top)
  segments(i - top, all_marginal$upper[i], i + top)
}

#conditional occupancy
#PUMA
puma_given_jag <- predict(
  mod_null,
  type = "state",
  species = "puma",
  cond = "jaguar"
)

head(puma_given_jag)

puma_no_jag <- predict(
  mod_null,
  type = "state",
  species = "puma",
  cond = "-jaguar"
)

head(puma_no_jag)

#JAGUAR
jag_given_puma <- predict(
  mod_null,
  type = "state",
  species = "jaguar",
  cond = "puma"
)
head(jag_given_puma)

jag_no_puma <- predict(
  mod_null,
  type = "state",
  species = "jaguar",
  cond = "-puma"
)
head(jag_no_puma)

cond_puma <- rbind(
  puma_given_jag[1, ],
  puma_no_jag[1, ]
)
cond_puma$Jaguar_status <- c("Present", "Absent")

cond_jag <- rbind(
  jag_given_puma[1, ],
  jag_no_puma[1, ]
)
cond_jag$Puma_status <- c("Present", "Absent")

#Plotting the conditional occupancy
library(RColorBrewer)

# Puma conditional on jaguar
png("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/conditional_puma_jag.png", 
    width = 12, height = 5, res = 300, units = "in")

plot(
  1:2,
  cond_puma$Predicted,
  ylim = c(0, 1),
  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  col = c("#A65628", "#999999"),
  xaxt = "n",
  xlab = "Jaguar status",
  ylab = "Puma occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_puma$Jaguar_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_puma$lower[i], i, cond_puma$upper[i], col = c("#A65628", "#999999")[i])
  segments(i - top, cond_puma$lower[i], i + top, col = c("#A65628", "#999999")[i])
  segments(i - top, cond_puma$upper[i], i + top, col = c("#A65628", "#999999")[i])
}

dev.off()

# Jaguar conditional on puma
png("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/conditional_jag_puma.png", 
    width = 12, height = 5, res = 300, units = "in")

plot(
  1:2,
  cond_jag$Predicted,
  ylim = c(0, 1),
  xlim = c(0.5, 2.5),
  pch = 19,
  cex = 1.5,
  col = c("#E6AB02", "#999999"),
  xaxt = "n",
  xlab = "Puma status",
  ylab = "Jaguar occupancy (95% CI)"
)

axis(1, at = 1:2, labels = cond_jag$Puma_status)

top <- 0.05
for (i in 1:2) {
  segments(i, cond_jag$lower[i], i, cond_jag$upper[i], col = c("#E6AB02", "#999999")[i])
  segments(i - top, cond_jag$lower[i], i + top, col = c("#E6AB02", "#999999")[i])
  segments(i - top, cond_jag$upper[i], i + top, col = c("#E6AB02", "#999999")[i])
}

dev.off()

#Adding puma and jaguar silhouettes

library(ggplot2)
library(dplyr)
library(rphylopic)

# get silhouettes by UUID (reliable, no search needed)
puma_img <- pick_phylopic(name = "Puma concolor", n = 1)
jag_img  <- pick_phylopic(name = "Panthera onca", n = 1)

puma_img <- get_phylopic(uuid = "3f8eff77-2868-4121-8d7d-a55ebdd49e04")
jag_img  <- get_phylopic(uuid = "08d6595a-aaff-4e13-ab6a-5458719194c6")

# --- Puma conditional on Jaguar ---

cond_puma_df <- data.frame(
  JaguarStatus = c("Present", "Absent"),
  Predicted    = c(puma_given_jag$Predicted[1], puma_no_jag$Predicted[1]),
  lower        = c(puma_given_jag$lower[1],     puma_no_jag$lower[1]),
  upper        = c(puma_given_jag$upper[1],     puma_no_jag$upper[1])
)

cond_puma_df$JaguarStatus <- factor(cond_puma_df$JaguarStatus, levels = c("Present", "Absent"))

ggplot(cond_puma_df, aes(x = JaguarStatus, y = Predicted, color = JaguarStatus)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_color_manual(values = c("Present" = "#A65628", "Absent" = "#999999")) +
  add_phylopic(img = puma_img, x = 2, y = 0.8, ysize = 0.20, fill = "#A65628", color = "black", alpha = 1)+
  coord_cartesian(ylim = c(0, 1)) + 
  xlab("Jaguar status") +
  ylab("Puma conditional occupancy (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(
  "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/conditional_puma_jag.png",
  width = 12, height = 5, dpi = 300, units = "in"
)
# --- Jaguar conditional on Puma ---

cond_jag_df <- data.frame(
  PumaStatus = c("Present", "Absent"),
  Predicted  = c(jag_given_puma$Predicted[1], jag_no_puma$Predicted[1]),
  lower      = c(jag_given_puma$lower[1],     jag_no_puma$lower[1]),
  upper      = c(jag_given_puma$upper[1],     jag_no_puma$upper[1])
)

cond_jag_df$PumaStatus <- factor(cond_jag_df$PumaStatus, levels = c("Present", "Absent"))

ggplot(cond_jag_df, aes(x = PumaStatus, y = Predicted, color = PumaStatus)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_color_manual(values = c("Present" = "#E6AB02", "Absent" = "#999999")) +
  add_phylopic(img = jag_img, x = 2, y = 0.8, ysize = 0.20, fill = "#E6AB02", color = "black", alpha = 1)+
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Puma status") +
  ylab("Jaguar conditional occupancy (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(
  "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/conditional_jag_puma.png",
  width = 12, height = 5, dpi = 300, units = "in"
)

library(ggplot2)
library(rphylopic)

# Puma silhouette
ggplot() +
  add_phylopic(img = puma_img, x = 0.5, y = 0.5, ysize = 0.5, fill = "#A65628", color = "black", alpha = 1) +
  xlim(-0.5, 1.5) + ylim(-0.5, 1.5) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(40, 40, 40, 40))

ggsave(
  "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/puma_silhouette.png",
  width = 6, height = 4, dpi = 300, bg = "transparent"
)

# Jaguar silhouette
ggplot() +
  add_phylopic(img = jag_img, x = 0.5, y = 0.5, ysize = 0.5, fill = "#E6AB02", color = "black", alpha = 1) +
  xlim(-0.5, 1.5) + ylim(-0.5, 1.5) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(40, 40, 40, 40))

ggsave(
  "/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures/jag_silhouette.png",
  width = 6, height = 4, dpi = 300, bg = "transparent"
)
