setwd("/Users/dwijadesai/Library/CloudStorage/OneDrive-UBC/Honours_peru/puma_project/raw_data")

library(dplyr)
library(ggplot2)
library(MASS)
library(sf)
library(eks)
library(colorspace)
library(ggnewscale)

covs <- read.csv("covariates.csv", stringsAsFactors = FALSE)
puma  <- read.csv("puma.csv", stringsAsFactors = FALSE)
jag  <- read.csv("jagu.csv", stringsAsFactors = FALSE)
prey <- read.csv("all_prey.csv", header = TRUE)

#join with the covs
puma_full <- left_join(puma, covs, by = "Station")
jag_full  <- left_join(jag,  covs, by = "Station")
prey_full <- left_join(prey, covs, by = "Station")

#removing missing coords
puma_full <- puma_full %>% filter(!is.na(utm_x), !is.na(utm_y))
jag_full  <- jag_full  %>% filter(!is.na(utm_x), !is.na(utm_y))
prey_full <- prey_full %>% filter(!is.na(utm_x), !is.na(utm_y))

#setting up sites for each species
puma_sites <- puma_full %>% distinct(Station, utm_x, utm_y)
jag_sites  <- jag_full  %>% distinct(Station, utm_x, utm_y)
prey_sites <- prey_full %>% distinct(Station, utm_x, utm_y)

#Colour palette - ColorBrewer "Set1"
library(RColorBrewer)
display.brewer.all()
brewer.pal(8, "Dark2")

col_puma <- "#A65628"
col_jag  <- "#E6AB02"
col_prey <- "#2D9E5F"

#puma kde
plot_puma <- ggplot(puma_sites, aes(utm_x, utm_y)) +
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity") +
  scale_fill_manual(
    values = colorRampPalette(c("white", col_puma))(11),
    guide  = "none"
  ) +
  geom_density_2d(
    contour_var = "ndensity",
    breaks = 0.50,
    color = "grey30",
    linewidth = 0.6
  ) +
  geom_point(shape = 21, fill = col_puma, color = "white", size = 2.5) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Puma — Kernel Density")

#jaguar kde
plot_jag <- ggplot(jag_sites, aes(utm_x, utm_y)) +
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity") +
  scale_fill_manual(
    values = colorRampPalette(c("white", col_jag))(11),
    guide  = "none"
  ) +
  geom_density_2d(
    contour_var = "ndensity",
    breaks = 0.50,
    color = "grey30",
    linewidth = 0.6
  ) +
  geom_point(shape = 21, fill = col_jag, color = "white", size = 2.5) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Jaguar — Kernel Density")

#prey kde
plot_prey <- ggplot(prey_sites, aes(utm_x, utm_y)) +
  geom_density_2d_filled(alpha = 0.7, contour_var = "ndensity") +
  scale_fill_manual(
    values = colorRampPalette(c("white", col_prey))(11),
    guide  = "none"
  ) +
  geom_density_2d(
    contour_var = "ndensity",
    breaks = 0.50,
    color = "grey30",
    linewidth = 0.6
  ) +
  geom_point(shape = 21, fill = col_prey, color = "white", size = 2.5) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Prey — Kernel Density")

# ── Helper: build overlap site table ────────────────────────────────────────
make_overlap <- function(sites_a, label_a, col_a,
                         sites_b, label_b, col_b) {
  df <- full_join(
    sites_a %>% mutate(A = 1),
    sites_b %>% mutate(B = 1),
    by = c("Station", "utm_x", "utm_y")
  )
  df$A[is.na(df$A)] <- 0
  df$B[is.na(df$B)] <- 0
  df$Category <- case_when(
    df$A == 1 & df$B == 1 ~ "Both",
    df$A == 1             ~ label_a,
    df$B == 1             ~ label_b
  )
  df
}

# ── Helper: plot two-species overlap KDE ────────────────────────────────────
# Two separate density layers (one per species) + coloured points on top.
# Requires the ggnewscale package.
make_overlap_kde <- function(sites_a, label_a, col_a,
                             sites_b, label_b, col_b,
                             title_str) {
  
  overlap_df <- make_overlap(sites_a, label_a, col_a,
                             sites_b, label_b, col_b)
  
  pt_cols <- setNames(
    c(col_a, col_b, "black"),   # gold for "Both"
    c(label_a, label_b, "Both")
  )
  
  ggplot() +
    # Species A density
    geom_density_2d_filled(
      data = sites_a, aes(utm_x, utm_y),
      alpha = 0.5, contour_var = "ndensity"
    ) +
    scale_fill_manual(
      values = colorRampPalette(c("white", col_a))(11),
      guide  = "none"
    ) +
    geom_density_2d(
      data = sites_a, aes(utm_x, utm_y),
      contour_var = "ndensity", breaks = 0.05,
      color = col_a, linewidth = 0.5, linetype = "dashed"
    ) +
    geom_density_2d(
      data = sites_a, aes(utm_x, utm_y),
      contour_var = "ndensity", breaks = 0.50,
      color = col_a, linewidth = 0.8
    ) +
    ggnewscale::new_scale_fill() +
    geom_density_2d_filled(
      data = sites_b, aes(utm_x, utm_y),
      alpha = 0.5, contour_var = "ndensity"
    ) +
    scale_fill_manual(
      values = colorRampPalette(c("white", col_b))(11),
      guide  = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_density_2d(
      contour_var = "ndensity",
      breaks = 0.50,
      color = "grey30",
      linewidth = 0.6
    ) +
    geom_point(
      data = overlap_df,
      aes(utm_x, utm_y, fill = Category),
      shape = 21, color = "white", size = 2.5
    ) +
    scale_fill_manual(values = pt_cols, name = "Category") +
    coord_equal() +
    theme_minimal() +
    labs(title = title_str)
}

#jaguar and puma overlap
plot_jag_puma <- make_overlap_kde(
  jag_sites,  "Jaguar", col_jag,
  puma_sites, "Puma",   col_puma,
  "Jaguar vs Puma — Kernel Density"
)

#jaguar and prey overlap
plot_jag_prey <- make_overlap_kde(
  jag_sites,  "Jaguar", col_jag,
  prey_sites, "Prey",   col_prey,
  "Jaguar vs Prey — Kernel Density"
)

#puma and prey overlap
plot_puma_prey <- make_overlap_kde(
  puma_sites, "Puma", col_puma,
  prey_sites, "Prey", col_prey,
  "Puma vs Prey — Kernel Density"
)

#plotting all
plot_puma
plot_jag
plot_prey
plot_jag_puma
plot_jag_prey
plot_puma_prey

#presenting
library(patchwork)

combined <- (plot_jag_puma | plot_jag_prey | plot_puma_prey) +
  plot_annotation(tag_levels = "A")
combined

ggsave("kde_figures.png", combined, width = 18, height = 12, dpi = 300)



