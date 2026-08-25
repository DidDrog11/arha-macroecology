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
pred_counts <- posterior_epred(model_adm2_zinb, ndraws = 500) |> apply(2, mean)

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
       subtitle = "Blue = Undersampled (Surveillance Coldspot) | Orange = Oversampled (Surveillance Hotspot)") +
  theme(legend.position = "bottom")

ggsave(here("output", "figures", "global_residual_map.png"), 
       plot = p_map_residuals, 
       width = 12, height = 7, dpi = 300, bg = "white")

# 3. Inset Maps -----------------------------------------------------------
# Convert data to unprojected WGS84 for undistorted insets
map_data_resid_4326 <- project(map_data_resid, "EPSG:4326")
world_map_4326 <- project(world_map, "EPSG:4326")

# Define Lat/Lon coordinates for the insets
insets_bbox <- list(
  us   = c(xmin = -120, xmax = -104, ymin = 31, ymax = 49),
  wafr = c(xmin = -18,  xmax = 15,  ymin = 0,  ymax = 20),
  asia = c(xmin = 45,   xmax = 90,  ymin = 27, ymax = 58)
)

col_us   <- "#CC79A7"
col_wafr <- "#B8860B"
col_asia <- "#009E73"

# Zoomed maps function
create_inset <- function(data_4326, bg_4326, bbox, title, box_colour) {
  ggplot() +
    theme(panel.background = element_rect(fill = bg_colour, colour = NA),
          plot.background = element_rect(fill = bg_colour, colour = NA),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = box_colour, fill = NA, linewidth = 2),
          plot.title = element_text(colour = text_colour, size = 10, face = "bold", hjust = 0.5)) +
    geom_spatvector(data = bg_4326, fill = land_colour, colour = NA) +
    geom_spatvector(data = data_4326, aes(fill = log_residual), colour = NA, show.legend = FALSE) +
    scale_fill_gradient2(low = colour_cold, mid = land_colour, high = colour_hot,
                         midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
    geom_spatvector(data = bg_4326, fill = NA, colour = "grey20", linewidth = 0.6) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
             ylim = c(bbox["ymin"], bbox["ymax"]), 
             expand = FALSE, crs = 4326) +
    labs(title = title)
}

# Generate the three insets
us_states <- gadm_adm2_proj |>
  filter(GID_0 == "USA") |>
  filter(!NAME_1 %in% c("Alaska", "Hawaii")) |>
  aggregate(by = "NAME_1") |>
  makeValid()

inset_us <- create_inset(map_data_resid_4326, world_map_4326, insets_bbox$us, 
                         "Intermountain West\n(Han et al. Hotspot)", col_us) +
  geom_spatvector(data = project(us_states, "EPSG:4326"), fill = NA, colour = "grey20", linewidth = 0.4) +
  coord_sf(xlim = c(insets_bbox$us["xmin"], insets_bbox$us["xmax"]),
           ylim = c(insets_bbox$us["ymin"], insets_bbox$us["ymax"]),
           expand = FALSE, crs = 4326)
inset_wafr <- create_inset(map_data_resid_4326, world_map_4326, insets_bbox$wafr, "West Africa\n(Lassa Fever Belt)", col_wafr)
inset_asia <- create_inset(map_data_resid_4326, world_map_4326, insets_bbox$asia, "Central Asia\n(Han et al. Novel Hotspot)", col_asia)

# Create bounding box polygons to overlay on the main Mollweide map
create_bbox_poly <- function(bbox) {
  pol_str <- sprintf("POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
                     bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymin"],
                     bbox["xmax"], bbox["ymax"], bbox["xmin"], bbox["ymax"],
                     bbox["xmin"], bbox["ymin"])
  vect(pol_str, crs = "EPSG:4326") |> project(mollweide_crs)
}

bbox_us_moll   <- create_bbox_poly(insets_bbox$us)
bbox_wafr_moll <- create_bbox_poly(insets_bbox$wafr)
bbox_asia_moll <- create_bbox_poly(insets_bbox$asia)

# Add the semi-transparent coloured anchors to the main map
p_map_residuals_anchored <- p_map_residuals +
  geom_spatvector(data = bbox_us_moll,   fill = scales::alpha(col_us, 0.2),   colour = col_us,   linewidth = 1) +
  geom_spatvector(data = bbox_wafr_moll, fill = scales::alpha(col_wafr, 0.2), colour = col_wafr, linewidth = 1) +
  geom_spatvector(data = bbox_asia_moll, fill = scales::alpha(col_asia, 0.2), colour = col_asia, linewidth = 1)

# Assemble using the new anchored main map
row_1 <- wrap_elements(p_map_residuals_anchored + labs(title = "a) Global Surveillance Mismatch"))
row_2 <- (inset_us | inset_wafr | inset_asia)

figure_1 <- row_1 / row_2 + 
  plot_layout(heights = c(2, 1)) & 
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(here("output", "figures", "figure_1_surveillance_map.png"), 
       figure_1, width = 12, height = 10, dpi = 300, bg = "white")

# 4. Hotspot-Coldspot Table -----------------------------------------------

# Identify the most extreme residuals
mean_res <- mean(residual_data$log_residual, na.rm = TRUE)
sd_res   <- sd(residual_data$log_residual, na.rm = TRUE)
thresh_hot  <- mean_res + (3 * sd_res)
thresh_cold <- mean_res - (3 * sd_res)

message("Thresholds | Hotspot: > ", round(thresh_hot, 2), " | Coldspot: < ", round(thresh_cold, 2))

realm_path <- here("data", "external", "One Earth Realms", "Realm2023.shp")
realms_vect <- vect(realm_path)

district_centroids_4326 <- centroids(gadm_adm2_proj) |> project("EPSG:4326")
district_centroids_4326$realm <- extract(realms_vect, district_centroids_4326)$BioGeoRelm

realm_lookup <- tibble(GID_2 = district_centroids_4326$GID_2, 
                       realm = district_centroids_4326$realm)

# 2. Rebuild anomalies with State/Province and realm attached
anomalies <- residual_data |>
  mutate(type = case_when(log_residual > thresh_hot ~ "Hotspot (Oversampled)",
                          log_residual < thresh_cold ~ "Coldspot (Undersampled)",
                          TRUE ~ "Normal")) |>
  filter(type != "Normal") |>
  left_join(as.data.frame(gadm_adm2_proj) |>
              select(GID_2, COUNTRY, NAME_1, NAME_2), by = "GID_2") |>
  left_join(realm_lookup, by = "GID_2") |>
  mutate(NAME_2 = if_else(str_detect(NAME_2, "^Unorganized"), 
                          paste0(NAME_1, " (unincorporated)"), NAME_2))

# 3. Fold-difference fix: plain comparison for true-zero coldspots, 
#    genuine fold-multiplier only where both sides are non-zero
extreme_table_realm <- anomalies |>
  mutate(Location = paste0(NAME_2, ", ", NAME_1, ", ", COUNTRY),
         Comparison = case_when(
           n_hosts == 0 ~ paste0("0 observed vs. ", comma(round(predicted, 0)), " predicted"),
           TRUE ~ paste0(round(10^abs(log_residual), 0), "x difference")
         ),
         Observed = comma(n_hosts),
         Predicted = comma(round(predicted, 1)),
         Residual_Z = round((log_residual - mean_res) / sd_res, 1)) |>
  select(Type = type, Realm = realm, Location, Observed, Predicted, Comparison, Z_Score = Residual_Z)

# 4. Per-realm stratified extremes, not pure global top-N
extreme_table_by_realm <- extreme_table_realm |>
  group_by(Type, Realm) |>
  arrange(desc(abs(Z_Score))) |>
  slice_head(n = 3) |>
  ungroup() |>
  arrange(Realm, Type, desc(abs(Z_Score)))

print(extreme_table_by_realm)
write_csv(extreme_table_by_realm, here("output", "tables", "supp_table_surveillance_extremes_by_realm.csv"))
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
p_light  <- create_mfx_plot(model_adm2_zinb, "s_light", "a) Effect of Night Lights", y_limit = common_ylim)
# 2. Accessibility (Travel Time)
p_access <- create_mfx_plot(model_adm2_zinb, "s_access", "b) Effect of Remoteness", y_limit = common_ylim)
# 3. Biodiversity (Richness)
p_rich   <- create_mfx_plot(model_adm2_zinb, "s_richness", "c) Effect of Host Richness", y_limit = common_ylim)

# Combine the Marginal Effects for Supplementary Figure
mfx_grid <- (p_light + p_access + p_rich) +
  plot_annotation(title = "Marginal Effects of Surveillance Bias")
ggsave(here("output", "figures", "supplementary_figure_s3.png"), mfx_grid, width = 12, height = 4, bg = "white")

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

centroid_lookup <- tibble(
  GID_2 = gadm_adm2_proj$GID_2,
  x = crds(centroids(gadm_adm2_proj))[, 1],
  y = crds(centroids(gadm_adm2_proj))[, 2]
)

# residual_data already has exactly one row per modelled district, 
# in the same order the model was fit on
coords <- residual_data |>
  select(GID_2) |>
  left_join(centroid_lookup, by = "GID_2")

set.seed(123)  # reproducibility
if (nrow(coords) > 5000) {
  idx <- sample(seq_len(nrow(coords)), 5000)
} else {
  idx <- seq_len(nrow(coords))
}

dharma_sub <- createDHARMa(
  simulatedResponse = t(preds_sim)[idx, ],
  observedResponse = model_adm2_zinb$data$n_hosts[idx],
  fittedPredictedResponse = pred_counts[idx],
  integerResponse = TRUE
)

testSpatialAutocorrelation(dharma_sub, x = coords$x[idx], y = coords$y[idx])

# 7. Trait Data Completeness by Biogeographic Realm ------------------------
# R1's "Eltonian shortfall" comment — check whether Afrotropical genera
# genuinely show higher trait-data missingness. Reuses realm_lookup
# (already built above) and adm2_genus_lookup (already computed).

trait_data <- read_rds(here("data", "processed", "trait_data_phylo_matched.rds"))

if (!exists("pol_vars")) {
  pol_vars <- c("adult_mass_g", "litter_size", "litters_per_year", "gestation_d",
                "weaning_d", "age_first_reproduction_d", "max_longevity_d")
}

# Genus -> realm via majority-occurrence, same rule already validated in
# 07_04's region-assignment fix (a genus spans many districts/realms;
# assign it to whichever realm accounts for the plurality of its overlaps).
genus_realm_majority <- adm2_genus_lookup |>
  left_join(realm_lookup, by = "GID_2") |>
  drop_na(realm) |>
  count(genus, realm) |>
  group_by(genus) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(genus, realm)

trait_data_realm <- trait_data |>
  left_join(genus_realm_majority, by = "genus")

# Proportion imputed per trait, Afrotropic vs rest of world
imputation_by_realm <- trait_data_realm |>
  mutate(realm_group = if_else(realm == "Afrotropic", "Afrotropic", "Rest of World")) |>
  select(realm_group, all_of(pol_vars)) |>
  pivot_longer(-realm_group, names_to = "trait", values_to = "value") |>
  group_by(realm_group, trait) |>
  summarise(pct_missing = round(100 * mean(is.na(value)), 1), n = n(), .groups = "drop") |>
  pivot_wider(names_from = realm_group, values_from = c(pct_missing, n))

imputation_by_realm_full <- trait_data_realm |>
  filter(!is.na(realm)) |>
  select(realm, all_of(pol_vars)) |>
  pivot_longer(-realm, names_to = "trait", values_to = "value") |>
  group_by(realm, trait) |>
  summarise(pct_missing = round(100 * mean(is.na(value)), 1), .groups = "drop")

# Full trait x realm grid
imputation_by_realm_full |>
  pivot_wider(names_from = trait, values_from = pct_missing) |>
  print(n = Inf)

# Single summary: mean % missing across all 7 traits, per realm
imputation_by_realm_full |>
  group_by(realm) |>
  summarise(mean_pct_missing = round(mean(pct_missing), 1),
            n_species = sum(trait_data_realm$realm == first(realm), na.rm = TRUE)) |>
  arrange(desc(mean_pct_missing))


# 8. Family-Specific Testing Skew by Region ---------------------------------
# R2's "number sampled" comment, extended: does testing effort concentrate
# on one viral family within a region, independent of the species-level
# log_effort covariate? E.g. is North America predominantly Hantaviridae-
# tested, and the Lassa belt predominantly Arenaviridae-tested?

arha_db$pathogen |> count(non_integer)

family_effort_by_country <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  left_join(arha_db$host |> select(host_record_id, iso3c), by = "host_record_id") |>
  drop_na(iso3c, number_tested) |>
  group_by(iso3c, pathogen_family) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = pathogen_family, values_from = n_tested, values_fill = 0) |>
  mutate(total_tested = Arenaviridae + Hantaviridae,
         pct_arenaviridae = round(100 * Arenaviridae / total_tested, 1),
         pct_hantaviridae = round(100 * Hantaviridae / total_tested, 1)) |>
  filter(total_tested >= 50) |>  
  arrange(desc(total_tested))

print(family_effort_by_country, n = 30)

## Direct check on the two named examples
family_effort_by_country |> filter(iso3c %in% c("USA", "NGA", "SLE", "LBR", "GIN"))

family_effort_by_country |>
  ggplot(aes(x = reorder(iso3c, pct_arenaviridae), y = pct_arenaviridae)) +
  geom_col() +
  coord_flip() +
  labs(title = "Family-specific testing skew by country",
       subtitle = "% of individuals tested for Arenaviridae vs. Hantaviridae (min. 50 tested)",
       x = NULL, y = "% Arenaviridae")
