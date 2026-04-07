library(terra)
library(sf)
library(dplyr)

setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/puma_r")

#loading the raster data
prey_19 <- rast("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data/rasters/prey19.tif")
prey_21 <- rast("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data/rasters/prey21.tif")
prey_23 <- rast("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data/rasters/prey23.tif")

#loading the covs csv files for camera locations
cameras <- read.csv("covs.csv")
#convert to proper CRS
cameras_sf <- st_as_sf(
  cameras, 
  coords = c("utm_y", "utm_x"), 
  crs = 5389
  )

head(cameras_sf$year)

cameras$year <- as.integer(cameras$x)
table(cameras$year)

cams_2023 <- cameras_sf[cameras_sf$year == 2023, ]
nrow(cams_2023)
#Extraction function
extract_abundance <- function(raster, cameras_sf, year_val) {
  
  cams_year <- cameras_sf[cameras_sf$year == year_val, ]
  
  buffers <- st_buffer(cams_year, dist = 500)
  
  vals <- terra::extract(
    raster,
    vect(buffers),
    fun = mean,
    na.rm = TRUE
  )
  
  data.frame(
    Station = cams_year$Station,
    year = year_val,
    prey_abundance = vals[,2]
  )
}

table(cameras_sf$year)

#Running by year
abund_2019 <- extract_abundance(prey_19, cameras_sf, 2019)
abund_2021 <- extract_abundance(prey_21, cameras_sf, 2021)
abund_2023 <- extract_abundance(prey_23, cameras_sf, 2023)




