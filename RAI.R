setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(camtrapR)
library(dplyr)
library(lubridate)

pumas<- read.csv("puma.csv", header = TRUE)
jags <- read.csv("jagu.csv", header = TRUE)
prey <- read.csv("all_prey.csv", header = TRUE)
covs <- read.csv("covariates.csv", header = TRUE)

head(covs)
colnames(covs)
head(pumas)
head(jags)
head(prey)
colnames(pumas)

# Force the column to numeric in all three before combining
# Non-numeric values will become NA — we'll check for those after
pumas <- pumas %>% mutate(Number_Individuals. = as.numeric(Number_Individuals.))
jags  <- jags  %>% mutate(Number_Individuals. = as.numeric(Number_Individuals.))
prey  <- prey  %>% mutate(Number_Individuals. = as.numeric(Number_Individuals.))

# Now combine
all_records <- bind_rows(pumas, jags, prey)

nrow(all_records)
table(all_records$Code)  # see all species and their raw detection counts

all_records <- all_records %>%
  mutate(Code = trimws(Code))  # strips leading/trailing whitespace from all codes

# Verify it's gone
table(all_records$Code)["SARM"]  # should now show one combined count of 430

all_records <- all_records %>%
  mutate(
    time_str = paste0(Hour, ":", sprintf("%02d", Min), " ", AM_PM),
    datetime = mdy_hm(paste0(Month, "/", Day, "/", Year, " ", time_str))
  )

# Check it worked
head(all_records %>% select(Station, Code, datetime))

all_records <- all_records %>%
  arrange(Code, Station, datetime) %>%
  group_by(Code, Station) %>%
  mutate(
    time_diff = as.numeric(difftime(datetime, lag(datetime), units = "mins")),
    independent = is.na(time_diff) | time_diff >= 30
  ) %>%
  filter(independent) %>%
  ungroup()

# Check how many detections remain after filtering
table(all_records$Code)

range(all_records$datetime)

#calculating total trap nights
total_trap_nights <- sum(covs$Days_Operable)
total_trap_nights

#RAI
rai_table <- all_records %>%
  group_by(Code) %>%
  summarise(independent_detections = n()) %>%
  mutate(
    total_trap_nights = total_trap_nights,
    RAI = round((independent_detections / total_trap_nights) * 100, 4)
  ) %>%
  arrange(desc(RAI))

print(rai_table, n = Inf)

all_records <- all_records %>%
  mutate(
    datetime = case_when(
      AM_PM == "0" ~ ymd_hm(paste0(Year, "-", Month, "-", Day, " ", Hour, ":", sprintf("%02d", Min))),
      TRUE ~ mdy_hm(paste0(Month, "/", Day, "/", Year, " ", Hour, ":", sprintf("%02d", Min), " ", AM_PM))
    )
  )

# Drop the one bad B44 row
all_records <- all_records %>% filter(!is.na(datetime))

# Verify no more NAs
sum(is.na(all_records$datetime))

# Step 3 - rerun independence filter
all_records <- all_records %>%
  arrange(Code, Station, datetime) %>%
  group_by(Code, Station) %>%
  mutate(
    time_diff = as.numeric(difftime(datetime, lag(datetime), units = "mins")),
    independent = is.na(time_diff) | time_diff >= 30
  ) %>%
  filter(independent) %>%
  ungroup()

# Step 4 - recalculate RAI
rai_table <- all_records %>%
  group_by(Code) %>%
  summarise(independent_detections = n()) %>%
  mutate(
    total_trap_nights = total_trap_nights,
    RAI = round((independent_detections / total_trap_nights) * 100, 4)
  ) %>%
  arrange(desc(RAI))

print(rai_table, n = Inf)

# Save as CSV
write.csv(rai_table, "RAI_results.csv", row.names = FALSE)
