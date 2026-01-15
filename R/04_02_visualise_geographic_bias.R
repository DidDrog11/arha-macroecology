# ---
# 04_02_visualise_geographic_bias.R
#
# Purpose: Visualise the outputs of the Subnational (ADM2) ZINB Model
# ---

# 1. Setup ----------------------------------------------------------------
pacman::p_load(tidyverse, here, terra, tidyterra, sf, rnaturalearth, 
               brms, marginaleffects, patchwork, biscale, cowplot, DHARMa)

if(!exists("model_adm2_zinb")) {
  model_adm2_zinb <- read_rds(here("output", "models", "brms_adm2_zinb.rds"))
}

# Load Spatial Data
gadm_adm2_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp")) 
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
gadm_adm2_proj <- project(gadm_adm2_proj, mollweide_crs)

# Load the data used in the model
model_data <- read_rds(here("output", "models", "model_adm2_data_full.rds"))

# 2. Surveillance Gap Mapping ---------------------------------------------
world_map <- ne_countries(scale = "medium", returnclass = "sv") |>
  project(mollweide_crs) |>
  makeValid()

# Calculate residuals Observed vs. Predicted
pred_counts <- posterior_epred(model_adm2_zinb, ndraws = 500) %>% apply(2, mean)

residual_data <- model_data |>
  mutate(predicted = pred_counts,
         # Log-Ratio Residual: Positive = Undersampled, Negative = Oversampled
         log_residual = log10(n_hosts + 1) - log10(predicted + 1)) |>
  select(GID_2, n_hosts, predicted, log_residual)

# Join to spatial polygon
map_data_resid <- gadm_adm2_proj |>
  right_join(residual_data, by = "GID_2")

# Global map
max_res <- max(abs(residual_data$log_residual), na.rm = TRUE)
bg_colour   <- "white"
land_colour <- "#E0E0E0" 
text_colour <- "black"

colour_cold <- "#0072B2"
colour_mid  <- "#E0E0E0"
colour_hot  <- "#D55E00"

p_map_residuals <- ggplot() +
  # The Ocean (Background)
  theme(panel.background = element_rect(fill = bg_colour, colour = NA),
        plot.background = element_rect(fill = bg_colour, colour = NA),
        legend.background = element_rect(fill = bg_colour, colour = NA),
        legend.text = element_text(colour = text_colour),
        legend.title = element_text(colour = text_colour),
        plot.title = element_text(colour = text_colour, face = "bold"),
        plot.subtitle = element_text(colour = "grey20"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  # The Land Base Layer (Grey for no-data areas)
  geom_spatvector(data = world_map, fill = land_colour, colour = NA) +
  # The Residuals
  geom_spatvector(data = map_data_resid, aes(fill = log_residual), colour = NA) +
  scale_fill_gradient2(name = "Surveillance Gap\n(log10[Obs/Pred])",
                       low = colour_cold, mid = land_colour, high = colour_hot,
                       midpoint = 0, 
                       limits = c(-2.5, 2.5),
                       oob = scales::squish) +
  labs(title = "Global Mismatch in Sampling Effort",
       subtitle = "Blue = Oversampled (Hotspot) | Orange = Undersampled (Coldspot)") +
  theme(legend.position = "bottom")

ggsave(here("output", "figures", "global_residual_map.png"), 
       plot = p_map_residuals, 
       width = 12, height = 7, dpi = 300, bg = "white")

# 3. Inset Maps -----------------------------------------------------------
sampled_centroids <- map_data_resid |>
  centroids() |>
  crds()

# K-means clustering to find dense sampling regions
set.seed(123)
kmeans_result <- kmeans(sampled_centroids, centers = 15, nstart = 25)
clusters <- tibble(x = sampled_centroids[,1], y = sampled_centroids[,2], cluster = kmeans_result$cluster)

# Define regions of interest
roi_list <- list()

# Create zoomed maps function
create_inset <- function(data, xlim, ylim, title) {
  ggplot() +
    theme(panel.background = element_rect(fill = bg_colour, colour = NA),
          plot.background = element_rect(fill = bg_colour, colour = NA),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.8),
          plot.title = element_text(colour = text_colour, size = 10, face = "bold", hjust = 0.5)) +
    geom_spatvector(data = world_map, fill = land_colour, colour = NA) +
    geom_spatvector(data = data, aes(fill = log_residual), colour = NA, show.legend = FALSE) +
    scale_fill_gradient2(low = colour_cold, mid = land_colour, high = colour_hot,
                         midpoint = 0,
                         limits = c(-2.5, 2.5),
                         oob = scales::squish) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title)
}

# Define coordinates based on the global map (currently placeholders)
inset_asia <- create_inset(data = map_data_resid, 
                           xlim = c(5e6, 1.0e7),  ylim = c(3e6, 6.5e6), 
                           title = "Central Asia\n(Han et al. Novel Hotspot)")
inset_us <- create_inset(data = map_data_resid, xlim = c(-1.1e7, -0.7e7), ylim = c(3.5e6, 6e6),
                         title = "US Midwest\n(Han et al. Hotspot)")
inset_wafr <- create_inset(data = map_data_resid, xlim = c(-2.5e6, 2e6), ylim = c(0, 3e6),
                           title = "West Africa\n(Lassa Fever Belt)")

# Assemble
layout <- (p_map_residuals / (inset_us | inset_wafr | inset_asia)) + 
  plot_layout(heights = c(2, 1)) & 
  theme(plot.background = element_rect(fill = bg_colour, colour = NA))

ggsave(here("output", "figures", "fig_2_a.png"), 
       layout, width = 12, height = 10, dpi = 300)

# 4. Hotspot-Coldspot Table -----------------------------------------------

# Identify the most extreme residuals
mean_res <- mean(residual_data$log_residual, na.rm = TRUE)
sd_res   <- sd(residual_data$log_residual, na.rm = TRUE)
thresh_hot  <- mean_res + (3 * sd_res)
thresh_cold <- mean_res - (3 * sd_res)

message("Thresholds | Hotspot: > ", round(thresh_hot, 2), " | Coldspot: < ", round(thresh_cold, 2))

# Extract and Classify Outliers
anomalies <- residual_data |>
  mutate(type = case_when(log_residual > thresh_hot ~ "Hotspot (Oversampled)",
                          log_residual < thresh_cold ~ "Coldspot (Undersampled)",
                          TRUE ~ "Normal")) |>
  filter(type != "Normal") |>
  left_join(as.data.frame(gadm_adm2_proj) |>
              select(GID_2, COUNTRY, NAME_1, NAME_2), by = "GID_2")

# Format the Table for Manuscript
extreme_table <- anomalies |>
  group_by(type) |>
  arrange(desc(abs(log_residual))) |> # Sort by magnitude of error
  slice_head(n = 10) |>               # Take top 10 most extreme of each class
  ungroup() |>
  mutate(Location = paste0(NAME_2, ", ", COUNTRY),
         Fold_Difference = paste0(round(10^abs(log_residual), 0), "x"),
         Observed = comma(n_hosts),
         Predicted = comma(round(predicted, 1)),
         Residual_Z = round((log_residual - mean_res) / sd_res, 1)) |>
  select(Type = type, Location, Observed, Predicted, Fold_Difference, Z_Score = Residual_Z)

print(extreme_table)

# Save Outputs
write_csv(extreme_table, here("output", "tables", "table_1_surveillance_extremes.csv"))
write_csv(anomalies, here("output", "tables", "supp_table_all_anomalies.csv"))

# 5. Marginal Effects -----------------------------------------------------
pred_labels <- c(
  "s_light" = "Night-time Light Intensity (std)",
  "s_access" = "Accessibility (Travel Time, std)",
  "s_richness" = "Host Species Richness (std)"
)

create_mfx_plot <- function(model, predictor, title, y_limit = NULL) {
  p <- plot_predictions(model, condition = predictor, type = "response") +
    theme_minimal(base_family = "sans") +
    labs(title = title, 
         y = "Predicted Count", 
         x = pred_labels[[predictor]]) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.title = element_text(size = 10))
  
  if (!is.null(y_limit)) {
    p <- p + coord_cartesian(ylim = c(0, y_limit))
  }
  
  return(p)
}

common_ylim <- max(predict(model_adm2_zinb, newdata = datagrid(model = model_adm2_zinb, s_light = seq(-2, 5, length.out=100)), type = "response")[,1])

# 1. Night Lights (Wealth/Infra)
p_light  <- create_mfx_plot(model_adm2_zinb, "s_light", "b) Effect of Night Lights", y_limit = common_ylim)
# 2. Accessibility (Travel Time)
p_access <- create_mfx_plot(model_adm2_zinb, "s_access", "c) Effect of Remoteness", y_limit = common_ylim)
# 3. Biodiversity (Richness)
p_rich   <- create_mfx_plot(model_adm2_zinb, "s_richness", "d) Effect of Host Richness", y_limit = common_ylim)

# Combine
mfx_grid <- (p_light + p_access + p_rich)
ggsave(here("output", "figures", "marginal_effects_bias.png"), mfx_grid, width = 12, height = 4)

# Combine for Figure 1
row_1 <- wrap_elements(p_map_residuals + labs(title = "a) Global Surveillance Mismatch"))

# Define the middle row (Insets)
row_2 <- (inset_us | inset_wafr | inset_asia)

# Define the bottom row (MFX)
row_3 <- (p_light | p_access | p_rich)

# Stitch them together
figure_1 <- row_1 / row_2 / row_3 +
  plot_layout(heights = c(2, 1, 1.2)) & 
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave(here("output", "figures", "figure_1_combined.png"), 
       figure_1, width = 12, height = 14, dpi = 300, bg = "white")


# 6. DHARMa Diagnostics ---------------------------------------------------
# Simulate residuals
preds_sim <- posterior_predict(model_adm2_zinb, ndraws = 250)
dharma_obj <- createDHARMa(simulatedResponse = t(preds_sim),
                           observedResponse = model_adm2_zinb$data$n_hosts,
                           fittedPredictedResponse = pred_counts,
                           integerResponse = TRUE)

# QQ Plot and Residual vs Predicted
png(here("output", "figures", "dharma_diagnostics.png"), width = 800, height = 400)
plot(dharma_obj)
dev.off()

# Check for Spatial Autocorrelation
coords <-  tibble(x = crds(gadm_adm2_proj)[, 1],
                  y = crds(gadm_adm2_proj)[, 2],
                  GID_2 = gadm_adm2_proj$GID_2)

# This test can be slow, so we run it on a subset if N is huge
if(nrow(coords) > 5000) {
  idx <- sample(1:nrow(coords), 5000)
  testSpatialAutocorrelation(sim_res, x = coords$x[idx], y = coords$y[idx])
} else {
  testSpatialAutocorrelation(sim_res, x = coords$x, y = coords$y)
}

