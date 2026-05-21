# ==============================================================================
# 07_04_spatial_validation.R
# Purpose: Project model predictions onto geographic regions to validate against Han et al. (2015).
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, brms, countrycode, terra, tidyterra, sf, biscale, cowplot, rnaturalearth, exactextractr)

# 1. Load Inputs ----------------------------------------------------------
# Load the Model (N=49k)
fit_dyadic <- read_rds(here("output", "models", "brms_dyadic_N49k.rds")) 

# Load Host Traits and Tree
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))
host_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))
A <- ape::vcv.phylo(host_tree)[unique(fit_dyadic$data$tip_label), unique(fit_dyadic$data$tip_label)]

# Load Country/ISO data
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-01-09.rds"))
species_geo <- arha_db$host |>
  distinct(host_species, iso3c) |>
  drop_na(iso3c)

# 2. Produce Predictions --------------------------------------------------
# We want to predict risk based on Biology (Pace of Life), 
# ignoring Sampling Effort (set to max) and Viral Identity (set to zero/population level).
full_pred_frame <- host_traits |>
  select(tip_label, pace_of_life_pc1) |>
  drop_na(pace_of_life_pc1) |>
  mutate(log_effort = max(fit_dyadic$data$log_effort, na.rm = TRUE),
         pathogen_species_cleaned = "NewPathogen")

# Predict using Fixed Effects Only (re_formula = NA)
# This gives the "Biological Potential" of the host
preds <- posterior_epred(fit_dyadic, newdata = full_pred_frame, re_formula = NA)

# Summarise (Mean Probability per Species)
species_risk <- full_pred_frame |>
  mutate(pred_prob = colMeans(preds)) |>
  select(tip_label, pred_prob) |>
  mutate(host_species = str_replace_all(tip_label, "_", " "))

# 3. Define the Regions ---------------------------------------------------
get_isos <- function(region_name) {
  code_list <- countrycode::codelist
  
  isos <- code_list |>
    filter(un.regionsub.name == region_name | 
             region == region_name | 
             region23 == region_name) |>
    pull(iso3c) |>
    unique() |>
    na.omit()
  
  return(as.character(isos))
}

regions_list <- list("West Africa"     = get_isos("Western Africa"),
                     "North America"   = c("USA", "MEX"),
                     "Atlantic Forest" = c("BRA", "ARG", "PRY", "URY"),
                     "East Asia"       = get_isos("Eastern Asia"),
                     "Western Europe"  = get_isos("Western Europe"),
                     "Rest of World"   = "ROW")

risk_map_data <- species_geo |>
  inner_join(species_risk, by = "host_species") |> 
  drop_na(pred_prob) |>
  mutate(region = "Rest of World")

risk_map_clean <- risk_map_data |>
  arrange(region == "Rest of World") |> 
  # Keep only the first occurrence
  distinct(tip_label, .keep_all = TRUE)

# We assign species to a Hotspot if they occur in any country within that hotspot.
for (r in names(regions_list)) {
  if(r != "Rest of World") {
    risk_map_clean <- risk_map_clean |>
      mutate(region = if_else(iso3c %in% regions_list[[r]], r, region))
  }
}

# 4. Leave-One-Region-Out Spatial Validation ------------------------------
model_data_full <- fit_dyadic$data |>
  left_join(risk_map_clean |> select(tip_label, iso3c, region), by = "tip_label") |>
  mutate(cv_continent = countrycode(iso3c, "iso3c", "continent")) |>
  mutate(cv_continent = if_else(cv_continent %in% c("Americas"), "Americas", cv_continent)) |>
  drop_na(cv_continent)

loro_predictions <- list()
cv_folds <- unique(model_data_full$cv_continent)

for (target_fold in cv_folds) {
  message("Withholding Continent: ", target_fold)
  
  # Train on all other continents, test on the target continent
  train_data <- model_data_full |> filter(cv_continent != target_fold)
  test_data  <- model_data_full |> filter(cv_continent == target_fold)
  
  if(nrow(test_data) == 0) next
  
  # Prune tree matrix strictly to training subset
  A_train <- A[unique(train_data$tip_label), unique(train_data$tip_label)]
  
  formula_cv <- bf(is_reservoir ~ log_effort + pace_of_life_pc1 + 
                     (1 | gr(tip_label, cov = A_train)) + 
                     (1 | host_id) + 
                     (1 | pathogen_species_cleaned))
  
  # Reduced iterations for CV efficiency
  fit_cv <- brm(formula = formula_cv,
                data = train_data,
                data2 = list(A_train = A_train),
                family = bernoulli(link = "logit"),
                prior = c(prior(normal(0, 1), class = "b"),
                          prior(normal(-4, 2), class = "Intercept"),
                          prior(normal(0, 1), class = "sd")),
                chains = 4, cores = 4, iter = 2000, warmup = 1000, 
                control = list(adapt_delta = 0.98),
                file = here("output", "models", paste0("brms_loro_", str_replace_all(target_fold, " ", "_"), ".rds")))
  
  # Predict back onto the withheld region (out-of-sample biological potential)
  # Max effort, new pathogen, fixed effects only
  test_pred_frame <- test_data |>
    distinct(tip_label, pace_of_life_pc1) |>
    mutate(log_effort = max(train_data$log_effort),
           pathogen_species_cleaned = "NewPathogen")
  
  preds_cv <- posterior_epred(fit_cv, newdata = test_pred_frame, re_formula = NA)
  
  test_pred_frame$pred_prob_cv <- colMeans(preds_cv)
  test_pred_frame$withheld_region <- target_fold
  
  loro_predictions[[target_fold]] <- test_pred_frame
}

# Combine and correlate LORO predictions vs Full Model predictions
final_cv_df <- bind_rows(loro_predictions) |>
  left_join(species_risk |> rename(pred_prob_full = pred_prob), by = "tip_label")

cv_correlation <- cor.test(final_cv_df$pred_prob_cv, final_cv_df$pred_prob_full, method = "spearman", exact = FALSE)

# 5. Visualise Out-of-Sample Validation Boxplots --------------------------
# Format LORO predictions for plotting
species_risk_loro <- final_cv_df |> 
  select(tip_label, pred_prob = pred_prob_cv) |>
  mutate(host_species = str_replace_all(tip_label, "_", " "))

risk_map_loro_data <- species_geo |>
  inner_join(species_risk_loro, by = "host_species") |> 
  drop_na(pred_prob) |>
  mutate(region = "Rest of World")

# Ensure 1 species = 1 data point, prioritising hotspots
risk_map_loro_clean <- risk_map_loro_data |>
  arrange(region == "Rest of World") |> 
  distinct(tip_label, .keep_all = TRUE)

for (r in names(regions_list)) {
  if(r != "Rest of World") {
    risk_map_loro_clean <- risk_map_loro_clean |>
      mutate(region = if_else(iso3c %in% regions_list[[r]], r, region))
  }
}

risk_summary <- risk_map_loro_clean |>
  group_by(region) |>
  summarise(mean_risk = mean(pred_prob), 
            median_risk = median(pred_prob), 
            n_species = n_distinct(host_species)) |>
  arrange(desc(median_risk))

plot_data <- risk_map_loro_clean |>
  mutate(region = factor(region, levels = risk_summary$region))

p_geo_risk <- ggplot(plot_data, aes(x = region, y = pred_prob, fill = region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  labs(title = "Intrinsic Reservoir Potential by Region (Validation)",
       subtitle = "LORO cross-validated model predictions based on Host Biology only",
       y = "Predicted Reservoir Probability",
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

ggsave(here("output", "figures", "supplementary_figure_s7.png"), p_geo_risk, width = 8, height = 6, bg = "white")

# Validation stats
stats <- pairwise.t.test(risk_map_loro_clean$pred_prob, risk_map_loro_clean$region, p.adjust.method = "bonferroni")
print(stats)

# 6. Subnational Risk Maps ------------------------------------------------
iucn_path <- here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp")
iucn <- vect(iucn_path)
iucn_harmonised <- read_rds(here("data", "external", "iucn_harmonised.rds"))
iucn_subset <- iucn |>
  filter(sci_name %in% iucn_harmonised$scientific_name) |>
  left_join(iucn_harmonised |>
              select(scientific_name, gbif_id, gbif_name_simple),
            by = c("sci_name" = "scientific_name"))

pred_species <- species_risk |>
  left_join(host_traits |>
              select(gbif_id, species), by = c("host_species" = "species"))

iucn_in_pred <- iucn_subset[iucn_subset$gbif_id %in% pred_species$gbif_id, ]
iucn_in_pred <- makeValid(iucn_in_pred)
iucn_in_pred <- aggregate(iucn_in_pred, by = "gbif_id")

iucn_ready <- iucn_in_pred |>
  left_join(pred_species |> select(gbif_id, species = host_species, pred_prob), by = "gbif_id") |>
  project("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Create a raster template
r_template <- rast(extent = ext(-17996000, 17996000, -8998000, 8998000), # Mollweide bounds
                   res = 20000, # 20,000 meters = 20km
                   crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
iucn_proj <- project(iucn_ready, crs(r_template))

r_richness <- rasterize(iucn_proj, r_template, field = "gbif_id", fun = "count", background = 0)
names(r_richness) <- "richness"
r_hazard <- rasterize(iucn_proj, r_template, field = "pred_prob", fun = "mean", background = NA)
names(r_hazard) <- "hazard"

r_hazard[r_richness < 2] <- NA
r_richness[r_richness < 2] <- NA

# World borders for context
world_sf <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(crs(r_template))

r_stack <- c(r_richness, r_hazard)
df_map <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)

data_bivariate <- bi_class(df_map, 
                           x = richness, 
                           y = hazard, 
                           style = "quantile", 
                           dim = 3)

p_raster_map <- ggplot() +
  geom_sf(data = world_sf, fill = "grey20", colour = NA) +
  geom_tile(data = data_bivariate, aes(x = x, y = y, fill = bi_class)) +
  bi_scale_fill(pal = "DkBlue", dim = 3) +
  labs(title = "Global Landscape of Zoonotic Potential",
       subtitle = "Community Competence (Hazard) vs. Species Richness",
       caption = "Resolution: 20km. Areas with <2 modelled species masked.") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5),
        legend.position = "none")

legend_plot <- bi_legend(pal = "DkBlue",
                         dim = 3,
                         xlab = "Species Richness",
                         ylab = "Mean Hazard",
                         size = 8)

final_figure <- ggdraw() +
  draw_plot(p_raster_map, 0, 0, 1, 1) +
  draw_plot(legend_plot, 0.1, 0.15, 0.2, 0.2)

ggsave(here("output", "figures", "raster_bivariate.png"), 
       final_figure, width = 12, height = 8, bg = "white")

# 6. Subnational Hotspots -------------------------------------------------
# We aggregate the 20km pixels up to the administrative district level
gadm_adm2 <- vect(here("data", "gadm", "gadm_adm2_simplified.shp"))
gadm_adm2 <- project(gadm_adm2, crs(r_template))

# Zonal Statistics
if (!require("exactextractr")) install.packages("exactextractr")
library(exactextractr)

gadm_sf <- st_as_sf(gadm_adm2)

# Calculate mean hazard per district
district_hazard <- exact_extract(r_hazard, gadm_sf, 'mean')
gadm_sf$mean_intrinsic_hazard <- district_hazard
district_richness <- exact_extract(r_richness, gadm_sf, 'mean')
gadm_sf$mean_district_richness <- district_richness

# Save Table
table_out <- gadm_sf |>
  st_drop_geometry() |>
  select(GID_2, NAME_1, NAME_2, COUNTRY, mean_intrinsic_hazard, mean_district_richness) |>
  drop_na(mean_intrinsic_hazard, mean_district_richness)

write_csv(table_out, here("output", "tables", "District_Hazard_Scores.csv"))

# We apply the same 3x3 Quantile breaks used in the map
# "3-3" = High Richness, High Hazard (Purple)
# "1-3" = Low Richness, High Hazard (Pink)
# "3-1" = High Richness, Low Hazard (Blue)

district_classed <- bi_class(table_out, 
                             x = mean_district_richness, 
                             y = mean_intrinsic_hazard, 
                             style = "quantile", 
                             dim = 3) |>
  separate(bi_class, into = c("richness_cat", "hazard_cat"), sep = "-", remove = FALSE)

# High Hazard + High Richness
top_purple <- district_classed |>
  filter(bi_class == "3-3") |> 
  arrange(desc(mean_intrinsic_hazard)) |>
  distinct(COUNTRY, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  mutate(Descriptor = "High Hazard / High Richness")

# High Hazard + Low Richness
top_pink <- district_classed |>
  filter(bi_class == "1-3") |> 
  arrange(desc(mean_intrinsic_hazard)) |>
  distinct(COUNTRY, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  mutate(Descriptor = "High Hazard / Low Richness")

# Low Hazard + High Richness
top_blue <- district_classed |>
  filter(bi_class == "3-1") |> 
  arrange(mean_intrinsic_hazard) |>
  distinct(COUNTRY, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  mutate(Descriptor = "Low Hazard / High Richness")

# Combine into one Summary Table
summary_table <- bind_rows(top_pink, top_blue, top_purple) |>
  select(Category = bi_class, Descriptor, Country = COUNTRY, Province = NAME_1, District = NAME_2, Mean_Hazard = mean_intrinsic_hazard, Mean_Richness = mean_district_richness)

write_csv(summary_table, here("output", "tables", "district_subset_table.csv"))
