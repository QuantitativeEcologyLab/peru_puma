setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(overlap)
library(astroFns)
library(dplyr)
library(camtrapR)

#Color pallet
col_puma <- "#A65628"
col_jag  <- "#E6AB02"
col_prey <- "#2D9E5F"

puma_jag  <- read.csv("puma_jag.csv",  stringsAsFactors = FALSE)
prey_act  <- read.csv("all_prey.csv",  header = TRUE)

#Rename code to species
puma_jag <- puma_jag %>% rename(Species = Code)
prey_act  <- prey_act  %>% rename(Species = Code)

#Fixing AM/PM
puma_jag$Hour <- ifelse(puma_jag$AM_PM == "PM" & puma_jag$Hour != 12, puma_jag$Hour + 12, puma_jag$Hour)
puma_jag$Hour <- ifelse(puma_jag$AM_PM == "AM" & puma_jag$Hour == 12, 0, puma_jag$Hour)

prey_act$Hour <- ifelse(prey_act$AM_PM == "PM" & prey_act$Hour != 12, prey_act$Hour + 12, prey_act$Hour)
prey_act$Hour <- ifelse(prey_act$AM_PM == "AM" & prey_act$Hour == 12, 0, prey_act$Hour)

#Converting to Radians
pj_radians   <- hms2rad(h = puma_jag$Hour + puma_jag$Min / 60)
prey_radians <- hms2rad(h = prey_act$Hour  + prey_act$Min  / 60)

puma.tr <- pj_radians[puma_jag$Species == "PUMA"]
jag.tr  <- pj_radians[puma_jag$Species == "JAGU"]

agou.tr <- prey_radians[prey_act$Species == "AGOU"]
acou.tr <- prey_radians[prey_act$Species == "ACOU"]
brab.tr <- prey_radians[prey_act$Species == "BRAB"]
opos.tr <- prey_radians[prey_act$Species == "OPOS"]
sarm.tr <- prey_radians[prey_act$Species == "SARM"]
paca.tr <- prey_radians[prey_act$Species == "PACA"]
broc.tr <- prey_radians[prey_act$Species == "BROC"]
gant.tr <- prey_radians[prey_act$Species == "GANT"]
cpec.tr <- prey_radians[prey_act$Species == "CPEC"]
tapi.tr <- prey_radians[prey_act$Species == "TAPI"]

length(puma.tr)
length(jag.tr)
length(agou.tr)
length(acou.tr)
length(brab.tr)
length(opos.tr)
length(sarm.tr)
length(paca.tr)
length(broc.tr)
length(gant.tr)
length(cpec.tr)
length(tapi.tr)

#activity radial for puma and jag
puma_jag$DateTime <- format(as.POSIXct(
  paste(puma_jag$Year, puma_jag$Month, puma_jag$Day, puma_jag$Hour, puma_jag$Min),
  format = "%Y %m %d %H %M", tz = "America/Lima"), "%Y-%m-%d %H:%M:%S")

activityRadial(recordTable       = puma_jag,
               species           = "PUMA",
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTime",
               plotR             = TRUE,
               lwd               = 3,
               rp.type           = "p",
               poly.col          = adjustcolor(col_puma, alpha.f = 0.5))

activityRadial(recordTable       = puma_jag,
               species           = "JAGU",
               allSpecies        = FALSE,
               speciesCol        = "Species",
               recordDateTimeCol = "DateTime",
               plotR             = TRUE,
               lwd               = 3,
               rp.type           = "p",
               poly.col          = adjustcolor(col_jag, alpha.f = 0.5))

# PUMA vs JAGUAR
overlap.pj <- overlapEst(puma.tr, jag.tr, type = "Dhat4")
bs.pj      <- bootstrap(puma.tr, jag.tr, 999, type = "Dhat4")
mean(bs.pj) #0.8312821
bootCI(overlap.pj, bs.pj)["norm0", ] 
#    lower     upper 
#0.7799466 0.9089019 

#Temporal plot
overlapPlot(puma.tr, jag.tr,
            xscale = 24, xcenter = "noon",
            main = "Puma and Jaguar",
            linecol = c(col_puma, col_jag),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Δ = ", round(overlap.pj, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_jag, NA), bty = "o", cex = 0.7)

# PUMA vs PREY
# puma vs agouti
overlapPlot(puma.tr, agou.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Agouti",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Agouti (n=", length(agou.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)

overlap.pagou <- overlapEst(puma.tr, agou.tr, type = "Dhat4")
bs.pagou      <- bootstrap(puma.tr, agou.tr, 999, type = "Dhat4")
mean(bs.pagou)[1] #0.445505
#lower     upper 
#0.3783686 0.4640120 
bootCI(overlap.pagou, bs.pagou)["norm0", ]

# puma vs acouchi
overlapPlot(puma.tr, acou.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Acouchi",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Acouchi (n=", length(acou.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)

overlap.pacou <- overlapEst(puma.tr, acou.tr, type = "Dhat4")
bs.pacou      <- bootstrap(puma.tr, acou.tr, 999, type = "Dhat4")
mean(bs.pacou) #0.5404129
bootCI(overlap.pacou, bs.pacou)["norm0", ]
#lower     upper 
#0.4431405 0.5421804 

# puma vs brazilian rabbit
overlapPlot(puma.tr, brab.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Brazilian Rabbit",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Brazilian rabbit (n=", length(brab.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)

overlap.pbrab <- overlapEst(puma.tr, brab.tr, type = "Dhat4")
bs.pbrab      <- bootstrap(puma.tr, brab.tr, 999, type = "Dhat4")
mean(bs.pbrab) #0.7853953
bootCI(overlap.pbrab, bs.pbrab)["norm0", ]
#lower     upper 
#0.7193643 0.8329044 

# puma vs opossum
overlapPlot(puma.tr, opos.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Opossum",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Opossum (n=", length(opos.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)

overlap.popos <- overlapEst(puma.tr, opos.tr, type = "Dhat4")
bs.popos      <- bootstrap(puma.tr, opos.tr, 999, type = "Dhat4")
mean(bs.popos) #0.7292775
bootCI(overlap.popos, bs.popos)["norm0", ]
#lower     upper 
#0.6812759 0.7688991

# puma vs armadillo
overlapPlot(puma.tr, sarm.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Small Armadillo",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Armadillo (n=", length(sarm.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)

overlap.psarm <- overlapEst(puma.tr, sarm.tr, type = "Dhat4")
bs.psarm      <- bootstrap(puma.tr, sarm.tr, 999, type = "Dhat4")
mean(bs.psarm) #0.7340762
bootCI(overlap.psarm, bs.psarm)["norm0", ]
#lower     upper 
#0.6699102 0.7672661 

# puma vs paca
overlapPlot(puma.tr, paca.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Paca",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Paca (n=", length(paca.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)
overlap.ppaca <- overlapEst(puma.tr, paca.tr, type = "Dhat4")
bs.ppaca      <- bootstrap(puma.tr, paca.tr, 999, type = "Dhat4")
mean(bs.ppaca) #0.6980496
bootCI(overlap.ppaca, bs.ppaca)["norm0", ]
#   lower     upper 
#0.6468155 0.7345922 

# puma vs brocket deer
broc.tr <- broc.tr[!is.na(broc.tr)] #omitting the NAs

overlapPlot(puma.tr, broc.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Brocket Deer",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Brocket deer (n=", length(broc.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)
overlap.pbroc <- overlapEst(puma.tr, broc.tr, type = "Dhat4")
bs.pbroc      <- bootstrap(puma.tr, broc.tr, 999, type = "Dhat4")
mean(bs.pbroc) #0.7578507
bootCI(overlap.pbroc, bs.pbroc)["norm0", ]
#lower     upper 
#0.7005995 0.7908361 

# puma vs giant anteater
overlapPlot(puma.tr, gant.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Giant Anteater",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Giant anteater (n=", length(gant.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)
overlap.pgant <- overlapEst(puma.tr, gant.tr, type = "Dhat4")
bs.pgant      <- bootstrap(puma.tr, gant.tr, 999, type = "Dhat4")
mean(bs.pgant) #0.3312697
bootCI(overlap.pgant, bs.pgant)["norm0", ]
#lower     upper 
#0.2126667 0.3352999 

# puma vs collared peccary
overlapPlot(puma.tr, cpec.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Collared Peccary",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Collared peccary (n=", length(cpec.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)
overlap.pcpec <- overlapEst(puma.tr, cpec.tr, type = "Dhat4")
bs.pcpec      <- bootstrap(puma.tr, cpec.tr, 999, type = "Dhat4")
mean(bs.pcpec) 
#0.4194736
bootCI(overlap.pcpec, bs.pcpec)["norm0", ]
#lower     upper 
#0.3520783 0.4368882 

# puma vs tapir
overlapPlot(puma.tr, tapi.tr,
            xscale = 24, xcenter = "noon", main = "Puma and Tapir",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Tapir (n=", length(tapi.tr), ")")),
       lty = c(1, 2), col = c(col_puma, col_prey), bty = "o", cex = 0.7)
overlap.ptapi <- overlapEst(puma.tr, tapi.tr, type = "Dhat4")
bs.ptapi      <- bootstrap(puma.tr, tapi.tr, 999, type = "Dhat4")
mean(bs.ptapi)
#0.8365043
bootCI(overlap.ptapi, bs.ptapi)["norm0", ]
#lower     upper 
#0.7949754 0.8812114 


# JAGUAR vs PREY
# jaguar vs agouti
overlapPlot(jag.tr, agou.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Agouti",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Agouti (n=", length(agou.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jagou <- overlapEst(jag.tr, agou.tr, type = "Dhat4")
bs.jagou      <- bootstrap(jag.tr, agou.tr, 999, type = "Dhat4")
mean(bs.jagou) #0.5341967
bootCI(overlap.jagou, bs.jagou)["norm0", ]
#lower     upper 
#0.4715178 0.5821482 

# jaguar vs acouchi
overlapPlot(jag.tr, acou.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Acouchi",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Acouchi (n=", length(acou.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jacou <- overlapEst(jag.tr, acou.tr, type = "Dhat4")
bs.jacou      <- bootstrap(jag.tr, acou.tr, 999, type = "Dhat4")
mean(bs.jacou) #0.620499
bootCI(overlap.jacou, bs.jacou)["norm0", ]
#lower     upper 
#0.5304355 0.6623402 

# jaguar vs brazilian rabbit
overlapPlot(jag.tr, brab.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Brazilian Rabbit",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Brazilian rabbit (n=", length(brab.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jbrab <- overlapEst(jag.tr, brab.tr, type = "Dhat4")
bs.jbrab      <- bootstrap(jag.tr, brab.tr, 999, type = "Dhat4")
mean(bs.jbrab) #0.6754648
bootCI(overlap.jbrab, bs.jbrab)["norm0", ]
#lower     upper 
#0.5672042 0.7112970 

# jaguar vs opossum
overlapPlot(jag.tr, opos.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Opossum",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Opossum (n=", length(opos.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jopos <- overlapEst(jag.tr, opos.tr, type = "Dhat4")
bs.jopos      <- bootstrap(jag.tr, opos.tr, 999, type = "Dhat4")
mean(bs.jopos) #0.6126522
bootCI(overlap.jopos, bs.jopos)["norm0", ]
#lower     upper 
#0.5224213 0.6532575 

# jaguar vs armadillo
overlapPlot(jag.tr, sarm.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Small Armadillo",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Armadillo (n=", length(sarm.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jsarm <- overlapEst(jag.tr, sarm.tr, type = "Dhat4")
bs.jsarm      <- bootstrap(jag.tr, sarm.tr, 999, type = "Dhat4")
mean(bs.jsarm) #0.6124745
bootCI(overlap.jsarm, bs.jsarm)["norm0", ]
#lower     upper 
#0.5101238 0.6477693 

# jaguar vs paca
overlapPlot(jag.tr, paca.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Paca",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Paca (n=", length(paca.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jpaca <- overlapEst(jag.tr, paca.tr, type = "Dhat4")
bs.jpaca      <- bootstrap(jag.tr, paca.tr, 999, type = "Dhat4")
mean(bs.jpaca) #0.5837444
bootCI(overlap.jpaca, bs.jpaca)["norm0", ]
#lower     upper 
#0.4910511 0.6280919

# jaguar vs brocket deer
overlapPlot(jag.tr, broc.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Brocket deer",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Brocket deer (n=", length(broc.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jbroc <- overlapEst(jag.tr, broc.tr, type = "Dhat4")
bs.jbroc      <- bootstrap(jag.tr, broc.tr, 999, type = "Dhat4")
mean(bs.jbroc) #0.8141
bootCI(overlap.jbroc, bs.jbroc)["norm0", ]
#lower     upper 
#0.7598977 0.8786782

# jaguar vs giant anteater
overlapPlot(jag.tr, gant.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Giant Anteater",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Giant anteater (n=", length(gant.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jgant <- overlapEst(jag.tr, gant.tr, type = "Dhat4")
bs.jgant      <- bootstrap(jag.tr, gant.tr, 999, type = "Dhat4")
mean(bs.jgant) #0.4808301
bootCI(overlap.jgant, bs.jgant)["norm0", ]
#lower     upper 
#0.3588934 0.5122738 

# jaguar vs collared peccary
overlapPlot(jag.tr, cpec.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Collared Peccary",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Collared peccary (n=", length(cpec.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jcpec <- overlapEst(jag.tr, cpec.tr, type = "Dhat4")
bs.jcpec      <- bootstrap(jag.tr, cpec.tr, 999, type = "Dhat4")
mean(bs.jcpec) #0.5369604
bootCI(overlap.jcpec, bs.jcpec)["norm0", ]
#lower     upper 
#0.4638390 0.5936334

# jaguar vs tapir
overlapPlot(jag.tr, tapi.tr,
            xscale = 24, xcenter = "noon", main = "Jaguar and Tapir",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3")
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Tapir (n=", length(tapi.tr), ")")),
       lty = c(1, 2), col = c(col_jag, col_prey), bty = "o", cex = 0.7)
overlap.jtapi <- overlapEst(jag.tr, tapi.tr, type = "Dhat4")
bs.jtapi      <- bootstrap(jag.tr, tapi.tr, 999, type = "Dhat4")
mean(bs.jtapi) #0.7149723
bootCI(overlap.jtapi, bs.jtapi)["norm0", ]
#lower     upper 
#0.6396909 0.7714173

#MAKING A SUMMERY TABLE
results_summary <- data.frame(
  Predator = c(
    rep("Puma",   10),
    rep("Jaguar", 10),
    "Puma"
  ),
  Prey = c(
    "Agouti", "Acouchi", "Brocket rabbit", "Opossum", "Armadillo",
    "Paca", "Brocket deer", "Giant anteater", "Collared peccary", "Tapir",
    "Agouti", "Acouchi", "Brocket rabbit", "Opossum", "Armadillo",
    "Paca", "Brocket deer", "Giant anteater", "Collared peccary", "Tapir",
    "Jaguar"
  ),
  n_predator = c(
    rep(length(puma.tr), 10),
    rep(length(jag.tr),  10),
    length(puma.tr)
  ),
  n_prey = c(
    length(agou.tr), length(acou.tr), length(brab.tr), length(opos.tr), length(sarm.tr),
    length(paca.tr), length(broc.tr), length(gant.tr), length(cpec.tr), length(tapi.tr),
    length(agou.tr), length(acou.tr), length(brab.tr), length(opos.tr), length(sarm.tr),
    length(paca.tr), length(broc.tr), length(gant.tr), length(cpec.tr), length(tapi.tr),
    length(jag.tr)
  ),
  Dhat4 = c(
    overlap.pagou, overlap.pacou, overlap.pbrab, overlap.popos, overlap.psarm,
    overlap.ppaca, overlap.pbroc, overlap.pgant, overlap.pcpec, overlap.ptapi,
    overlap.jagou, overlap.jacou, overlap.jbrab, overlap.jopos, overlap.jsarm,
    overlap.jpaca, overlap.jbroc, overlap.jgant, overlap.jcpec, overlap.jtapi,
    overlap.pj
  ),
  Boot_mean = c(
    mean(bs.pagou), mean(bs.pacou), mean(bs.pbrab), mean(bs.popos), mean(bs.psarm),
    mean(bs.ppaca), mean(bs.pbroc), mean(bs.pgant), mean(bs.pcpec), mean(bs.ptapi),
    mean(bs.jagou), mean(bs.jacou), mean(bs.jbrab), mean(bs.jopos), mean(bs.jsarm),
    mean(bs.jpaca), mean(bs.jbroc), mean(bs.jgant), mean(bs.jcpec), mean(bs.jtapi),
    mean(bs.pj)
  ),
  CI_low = c(
    bootCI(overlap.pagou, bs.pagou)["norm0", "lower"],
    bootCI(overlap.pacou, bs.pacou)["norm0", "lower"],
    bootCI(overlap.pbrab, bs.pbrab)["norm0", "lower"],
    bootCI(overlap.popos, bs.popos)["norm0", "lower"],
    bootCI(overlap.psarm, bs.psarm)["norm0", "lower"],
    bootCI(overlap.ppaca, bs.ppaca)["norm0", "lower"],
    bootCI(overlap.pbroc, bs.pbroc)["norm0", "lower"],
    bootCI(overlap.pgant, bs.pgant)["norm0", "lower"],
    bootCI(overlap.pcpec, bs.pcpec)["norm0", "lower"],
    bootCI(overlap.ptapi, bs.ptapi)["norm0", "lower"],
    bootCI(overlap.jagou, bs.jagou)["norm0", "lower"],
    bootCI(overlap.jacou, bs.jacou)["norm0", "lower"],
    bootCI(overlap.jbrab, bs.jbrab)["norm0", "lower"],
    bootCI(overlap.jopos, bs.jopos)["norm0", "lower"],
    bootCI(overlap.jsarm, bs.jsarm)["norm0", "lower"],
    bootCI(overlap.jpaca, bs.jpaca)["norm0", "lower"],
    bootCI(overlap.jbroc, bs.jbroc)["norm0", "lower"],
    bootCI(overlap.jgant, bs.jgant)["norm0", "lower"],
    bootCI(overlap.jcpec, bs.jcpec)["norm0", "lower"],
    bootCI(overlap.jtapi, bs.jtapi)["norm0", "lower"],
    bootCI(overlap.pj,    bs.pj)   ["norm0", "lower"]
  ),
  CI_high = c(
    bootCI(overlap.pagou, bs.pagou)["norm0", "upper"],
    bootCI(overlap.pacou, bs.pacou)["norm0", "upper"],
    bootCI(overlap.pbrab, bs.pbrab)["norm0", "upper"],
    bootCI(overlap.popos, bs.popos)["norm0", "upper"],
    bootCI(overlap.psarm, bs.psarm)["norm0", "upper"],
    bootCI(overlap.ppaca, bs.ppaca)["norm0", "upper"],
    bootCI(overlap.pbroc, bs.pbroc)["norm0", "upper"],
    bootCI(overlap.pgant, bs.pgant)["norm0", "upper"],
    bootCI(overlap.pcpec, bs.pcpec)["norm0", "upper"],
    bootCI(overlap.ptapi, bs.ptapi)["norm0", "upper"],
    bootCI(overlap.jagou, bs.jagou)["norm0", "upper"],
    bootCI(overlap.jacou, bs.jacou)["norm0", "upper"],
    bootCI(overlap.jbrab, bs.jbrab)["norm0", "upper"],
    bootCI(overlap.jopos, bs.jopos)["norm0", "upper"],
    bootCI(overlap.jsarm, bs.jsarm)["norm0", "upper"],
    bootCI(overlap.jpaca, bs.jpaca)["norm0", "upper"],
    bootCI(overlap.jbroc, bs.jbroc)["norm0", "upper"],
    bootCI(overlap.jgant, bs.jgant)["norm0", "upper"],
    bootCI(overlap.jcpec, bs.jcpec)["norm0", "upper"],
    bootCI(overlap.jtapi, bs.jtapi)["norm0", "upper"],
    bootCI(overlap.pj,    bs.pj)   ["norm0", "upper"]
  )
)

#Saving for poster figures
setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures")

# puma vs jaguar
png("temporal_puma_jag.png", width = 6, height = 6, res = 300, units = "in")
overlapPlot(puma.tr, jag.tr,
            xscale = 24, xcenter = "noon",
            main = "Puma and Jaguar",
            linecol = c(col_puma, col_jag),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Δ = ", round(overlap.pj, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_jag, NA), bty = "o", cex = 0.7)
dev.off()

# puma vs agouti
png("temporal_puma_agou.png", width = 6, height = 6, res = 300, units = "in")
overlapPlot(puma.tr, agou.tr,
            xscale = 24, xcenter = "noon",
            main = "Puma and Agouti",
            linecol = c(col_puma, col_prey),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Agouti (n=", length(agou.tr), ")"),
         paste0("Δ = ", round(overlap.pagou, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)
dev.off()

# puma vs brazilian rabbit
png("temporal_puma_brab.png", width = 6, height = 6, res = 300, units = "in")
overlapPlot(puma.tr, brab.tr,
            xscale = 24, xcenter = "noon",
            main = "Puma and Brazilian Rabbit",
            linecol = c(col_puma, col_prey),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Brazilian rabbit (n=", length(brab.tr), ")"),
         paste0("Δ = ", round(overlap.pbrab, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)
dev.off()

# puma vs opossum
png("temporal_puma_opos.png", width = 6, height = 6, res = 300, units = "in")
overlapPlot(puma.tr, opos.tr,
            xscale = 24, xcenter = "noon",
            main = "Puma and Opossum",
            linecol = c(col_puma, col_prey),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Puma (n=", length(puma.tr), ")"),
         paste0("Opossum (n=", length(opos.tr), ")"),
         paste0("Δ = ", round(overlap.popos, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)
dev.off()

# jaguar vs agouti
png("temporal_jag_agou.png", width = 6, height = 6, res = 300, units = "in")
overlapPlot(jag.tr, agou.tr,
            xscale = 24, xcenter = "noon",
            main = "Jaguar and Agouti",
            linecol = c(col_jag, col_prey),
            olapcol = "antiquewhite3",
            xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright",
       c(paste0("Jaguar (n=", length(jag.tr), ")"),
         paste0("Agouti (n=", length(agou.tr), ")"),
         paste0("Δ = ", round(overlap.jagou, 2))),
       lty = c(1, 2, NA), col = c(col_jag, col_prey, NA), bty = "o", cex = 0.7)
dev.off()

png("figure5_temporal_overlap.png", width = 18, height = 12, res = 300, units = "in", bg = "transparent")

layout(matrix(c(1,2,3,4,5,0), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1))

# Figure 5: Combined temporal overlap panel
# Run this AFTER your main temporal overlap script so all objects are in memory

setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/final_figures")

png("figure5_temporal_overlap.png", width = 18, height = 12, res = 300, units = "in", bg = "transparent")

layout(matrix(c(1, 2, 3, 4, 5, 0), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1))

# (A) Puma vs Agouti
overlapPlot(puma.tr, agou.tr, xscale = 24, xcenter = "noon", main = "(A) Puma and Agouti",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3", xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright", c(paste0("Puma (n=", length(puma.tr), ")"), paste0("Agouti (n=", length(agou.tr), ")"), paste0("\u0394 = ", round(overlap.pagou, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)

# (B) Puma vs Opossum
overlapPlot(puma.tr, opos.tr, xscale = 24, xcenter = "noon", main = "(B) Puma and Opossum",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3", xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright", c(paste0("Puma (n=", length(puma.tr), ")"), paste0("Opossum (n=", length(opos.tr), ")"), paste0("\u0394 = ", round(overlap.popos, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)

# (C) Puma vs Brazilian Rabbit
overlapPlot(puma.tr, brab.tr, xscale = 24, xcenter = "noon", main = "(C) Puma and Brazilian Rabbit",
            linecol = c(col_puma, col_prey), olapcol = "antiquewhite3", xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright", c(paste0("Puma (n=", length(puma.tr), ")"), paste0("Brazilian rabbit (n=", length(brab.tr), ")"), paste0("\u0394 = ", round(overlap.pbrab, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_prey, NA), bty = "o", cex = 0.7)

# (D) Jaguar vs Agouti
overlapPlot(jag.tr, agou.tr, xscale = 24, xcenter = "noon", main = "(D) Jaguar and Agouti",
            linecol = c(col_jag, col_prey), olapcol = "antiquewhite3", xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright", c(paste0("Jaguar (n=", length(jag.tr), ")"), paste0("Agouti (n=", length(agou.tr), ")"), paste0("\u0394 = ", round(overlap.jagou, 2))),
       lty = c(1, 2, NA), col = c(col_jag, col_prey, NA), bty = "o", cex = 0.7)

# (E) Puma vs Jaguar
overlapPlot(puma.tr, jag.tr, xscale = 24, xcenter = "noon", main = "(E) Puma and Jaguar",
            linecol = c(col_puma, col_jag), olapcol = "antiquewhite3", xaxt = "n")
axis(1, at = seq(0, 24, by = 6), labels = paste0(seq(0, 24, by = 6), ":00"))
legend("bottomright", c(paste0("Puma (n=", length(puma.tr), ")"), paste0("Jaguar (n=", length(jag.tr), ")"), paste0("\u0394 = ", round(overlap.pj, 2))),
       lty = c(1, 2, NA), col = c(col_puma, col_jag, NA), bty = "o", cex = 0.7)

dev.off()
