# ==============================================================================
# 07_04_spatial_validation.R
# Purpose: Project model predictions onto geographic regions to validate against Han et al. (2015).
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, ape, brms, countrycode, terra, tidyterra, sf, biscale, cowplot, rnaturalearth, exactextractr, pROC)

# 1. Load Inputs ----------------------------------------------------------
# Load the Model (N=49k)
fit_dyadic <- read_rds(here("output", "models", "brms_dyadic_N49k.rds")) 

# Load Host Traits and Tree
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))
host_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))
A <- ape::vcv.phylo(host_tree)[unique(fit_dyadic$data$tip_label), unique(fit_dyadic$data$tip_label)]

# Load Country/ISO data
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-08-17.rds"))
species_geo <- arha_db$host |>
  distinct(host_species, iso3c) |>
  drop_na(iso3c)

# 2. Produce Predictions --------------------------------------------------
# Anchor effort at the 95th percentile rather than the literal max — the
# max belongs to Peromyscus maniculatus, the most reactively-surveilled
# species in the dataset (Sin Nombre reservoir), so using it as the
# universal "biological potential" ceiling undermines the paper's own
# surveillance-bias argument. 95th percentile is far less sensitive to
# any single record while still representing near-maximal effort.
effort_anchor <- quantile(fit_dyadic$data$log_effort, 0.95, na.rm = TRUE)

full_pred_frame <- host_traits |>
  select(tip_label, pace_of_life_pc1) |>
  drop_na(pace_of_life_pc1) |>
  mutate(log_effort = effort_anchor,
         pathogen_species_cleaned = "NewPathogen")

preds <- posterior_epred(fit_dyadic, newdata = full_pred_frame, re_formula = NA)

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

# Majority-occurrence region assignment (replaces the last-wins loop,
# which let list-definition order arbitrarily decide ~19% of species).
# Used for both risk_map_clean below and risk_map_loro_clean in Section 5.
assign_region_majority <- function(species_geo, regions_list) {
  species_geo |>
    mutate(region_match = case_when(
      iso3c %in% regions_list[["West Africa"]]     ~ "West Africa",
      iso3c %in% regions_list[["North America"]]   ~ "North America",
      iso3c %in% regions_list[["Atlantic Forest"]] ~ "Atlantic Forest",
      iso3c %in% regions_list[["East Asia"]]       ~ "East Asia",
      iso3c %in% regions_list[["Western Europe"]]  ~ "Western Europe",
      TRUE ~ "Rest of World"
    )) |>
    count(host_species, region_match) |>
    group_by(host_species) |>
    slice_max(n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(host_species, region = region_match)
}

species_region_lookup <- assign_region_majority(species_geo, regions_list)

risk_map_clean <- species_geo |>
  inner_join(species_risk, by = "host_species") |>
  drop_na(pred_prob) |>
  left_join(species_region_lookup, by = "host_species") |>
  mutate(region = replace_na(region, "Rest of World")) |>
  distinct(tip_label, .keep_all = TRUE)

# 4. Leave-One-Region-Out Spatial Validation ------------------------------
cv_continent_majority <- species_geo |>
  mutate(continent = countrycode(iso3c, "iso3c", "continent")) |>
  mutate(continent = if_else(continent == "Americas", "Americas", continent)) |>
  drop_na(continent) |>
  count(host_species, continent) |>
  group_by(host_species) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(host_species, cv_continent = continent)

model_data_full <- fit_dyadic$data |>
  mutate(host_species = str_replace_all(tip_label, "_", " ")) |>
  left_join(cv_continent_majority, by = "host_species") |>
  drop_na(cv_continent)

loro_predictions <- list()
fit_cv_list <- list()
auc_results <- list()
cv_folds <- unique(model_data_full$cv_continent)

for (target_fold in cv_folds) {
  message("Withholding Continent: ", target_fold)
  
  train_data <- model_data_full |> filter(cv_continent != target_fold)
  test_data  <- model_data_full |> filter(cv_continent == target_fold)
  
  if (nrow(test_data) == 0) next
  
  A_train <- A[unique(train_data$tip_label), unique(train_data$tip_label)]
  
  formula_cv <- bf(is_reservoir ~ log_effort + pace_of_life_pc1 + 
                     (1 | gr(tip_label, cov = A_train)) + 
                     (1 | host_id) + 
                     (1 | pathogen_species_cleaned))
  
  fit_cv <- brm(formula = formula_cv,
                data = train_data,
                data2 = list(A_train = A_train),
                family = bernoulli(link = "logit"),
                prior = c(prior(normal(0, 1), class = "b"),
                          prior(normal(-4, 2), class = "Intercept"),
                          prior(normal(0, 1), class = "sd")),
                chains = 4, cores = 8, iter = 2000, warmup = 1000, 
                control = list(adapt_delta = 0.98),
                file = here("output", "models", paste0("brms_loro_", str_replace_all(target_fold, " ", "_"), ".rds")))
  
  # "Biological potential" prediction for the regional comparison —
  # anchored at 95th percentile of TRAINING effort, not the max.
  train_effort_anchor <- quantile(train_data$log_effort, 0.95, na.rm = TRUE)
  
  test_pred_frame <- test_data |>
    distinct(tip_label, pace_of_life_pc1) |>
    mutate(log_effort = train_effort_anchor,
           pathogen_species_cleaned = "NewPathogen")
  
  preds_cv <- posterior_epred(fit_cv, newdata = test_pred_frame, re_formula = NA)
  test_pred_frame$pred_prob_cv <- colMeans(preds_cv)
  test_pred_frame$withheld_region <- target_fold
  
  loro_predictions[[target_fold]] <- test_pred_frame
  fit_cv_list[[target_fold]] <- fit_cv
  
  # AUC against ground truth: REAL observed per-dyad effort, real is_reservoir.
  # Wrapped in tryCatch so a fold with only one outcome class (unlikely here,
  # but possible) can't halt an unattended run.
  preds_real_effort <- posterior_epred(fit_cv, newdata = test_data, re_formula = NA)
  test_data$pred_prob_realistic <- colMeans(preds_real_effort)
  
  roc_obj <- tryCatch(
    pROC::roc(test_data$is_reservoir, test_data$pred_prob_realistic, quiet = TRUE),
    error = function(e) NULL
  )
  
  auc_results[[target_fold]] <- tibble(
    fold = target_fold,
    auc = if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_,
    n_dyads = nrow(test_data),
    n_positive = sum(test_data$is_reservoir)
  )
}

cv_auc_summary <- bind_rows(auc_results)
print(cv_auc_summary)
write_csv(cv_auc_summary, here("output", "tables", "revision_table_loro_auc.csv"))

# Phylogenetic leakage check (descriptive only, per earlier agreement):
# genus overlap between each withheld fold and its training set.
phylo_overlap <- map_dfr(cv_folds, function(fold) {
  train_genera <- model_data_full |> filter(cv_continent != fold) |> 
    distinct(tip_label) |> mutate(genus = word(tip_label, 1, sep = "_")) |> pull(genus) |> unique()
  test_genera <- model_data_full |> filter(cv_continent == fold) |> 
    distinct(tip_label) |> mutate(genus = word(tip_label, 1, sep = "_")) |> pull(genus) |> unique()
  tibble(fold = fold,
         n_test_genera = length(test_genera),
         n_shared_with_training = length(intersect(test_genera, train_genera)),
         pct_genera_shared = round(100 * length(intersect(test_genera, train_genera)) / length(test_genera), 1))
})
print(phylo_overlap)
write_csv(phylo_overlap, here("output", "tables", "revision_table_phylo_overlap.csv"))

# Combine and correlate LORO predictions vs Full Model predictions
final_cv_df <- bind_rows(loro_predictions) |>
  left_join(species_risk |> rename(pred_prob_full = pred_prob), by = "tip_label")

cv_correlation <- cor.test(final_cv_df$pred_prob_cv, final_cv_df$pred_prob_full, method = "spearman", exact = FALSE)
print(cv_correlation)

# 5. Visualise Out-of-Sample Validation Boxplots --------------------------
species_risk_loro <- final_cv_df |> 
  select(tip_label, pred_prob = pred_prob_cv) |>
  mutate(host_species = str_replace_all(tip_label, "_", " "))

risk_map_loro_clean <- species_geo |>
  inner_join(species_risk_loro, by = "host_species") |>
  drop_na(pred_prob) |>
  left_join(species_region_lookup, by = "host_species") |>
  mutate(region = replace_na(region, "Rest of World")) |>
  distinct(tip_label, .keep_all = TRUE)  |>
  filter(region != "Western Europe")

risk_summary <- risk_map_loro_clean |>
  group_by(region) |>
  summarise(mean_risk = mean(pred_prob), 
            median_risk = median(pred_prob), 
            n_species = n_distinct(host_species)) |>
  arrange(desc(median_risk))
print(risk_summary)

plot_data <- risk_map_loro_clean |>
  mutate(region = factor(region, levels = risk_summary$region))

p_geo_risk <- ggplot(plot_data, aes(x = region, y = pred_prob, fill = region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  labs(title = "Out-of-sample host probability by region",
       subtitle = "Leave-one-region-out cross-validated predictions from host biology only",
       y = "Predicted host probability (per host-virus dyad)",
       x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

ggsave(here("output", "figures", "supplementary_figure_s7.png"), p_geo_risk, width = 8, height = 6, bg = "white")

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

r_hazard[r_richness < 5] <- NA
r_richness[r_richness < 5] <- NA

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
  labs(title = "Global landscape of community host probability",
       caption = "Resolution: 20km. Areas with <5 modelled species masked.") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5),
        legend.position = "none")

legend_plot <- bi_legend(pal = "DkBlue",
                         dim = 3,
                         xlab = "Species richness",
                         ylab = "Mean host probability\n(per host-virus dyad)",
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

# Zonal statistics with coverage tracking.
# exact_extract's 'mean' silently ignores NA cells and renormalises over
# the remainder, so a district can receive a value derived from a small
# fraction of its area. Small, island and coastal districts are most
# affected, since ocean and species-poor cells are NA. We therefore
# record the proportion of each district covered by valid cells and
# retain only those with >=50% coverage.

gadm_sf$mean_intrinsic_hazard <- exact_extract(r_hazard, gadm_sf, 'mean')
gadm_sf$mean_district_richness <- exact_extract(r_richness, gadm_sf, 'mean')

gadm_sf$hazard_coverage <- exact_extract(r_hazard, gadm_sf, function(values, coverage) {
  if (sum(coverage) == 0) return(NA_real_)
  sum(coverage[!is.na(values)]) / sum(coverage)
})

gadm_sf$n_valid_cells <- exact_extract(r_hazard, gadm_sf, function(values, coverage) {
  sum(!is.na(values))
})

# Report how many districts this excludes, and their characteristics
coverage_summary <- gadm_sf |>
  st_drop_geometry() |>
  drop_na(mean_intrinsic_hazard) |>
  summarise(n_total = n(),
            n_retained = sum(hazard_coverage >= 0.5, na.rm = TRUE),
            n_excluded = sum(hazard_coverage < 0.5, na.rm = TRUE),
            pct_excluded = round(100 * n_excluded / n_total, 1))
print(coverage_summary)

table_out <- gadm_sf |>
  st_drop_geometry() |>
  select(GID_2, NAME_1, NAME_2, COUNTRY, mean_intrinsic_hazard,
         mean_district_richness, hazard_coverage, n_valid_cells) |>
  drop_na(mean_intrinsic_hazard, mean_district_richness) |>
  filter(hazard_coverage >= 0.5)

write_csv(table_out, here("output", "tables", "District_Hazard_Scores.csv"))

# Classes are labelled x-y, i.e. richness-probability:
# "3-3" = High richness, high mean host probability
# "1-3" = Low richness, high mean host probability
# "3-1" = High richness, low mean host probability
# Note: the two axes are near-independent (R2 = 0.038); they represent
# separate components of community-level hazard, not an interacting gradient.

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
  mutate(Descriptor = "High richness / High host probability")

# High Hazard + Low Richness
top_pink <- district_classed |>
  filter(bi_class == "1-3") |> 
  arrange(desc(mean_intrinsic_hazard)) |>
  distinct(COUNTRY, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  mutate(Descriptor = "Low richness / High host probability")

# Low Hazard + High Richness
top_blue <- district_classed |>
  filter(bi_class == "3-1") |> 
  arrange(mean_intrinsic_hazard) |>
  distinct(COUNTRY, .keep_all = TRUE) |>
  slice_head(n = 5) |>
  mutate(Descriptor = "High richness / Low host probability")

# Combine into one Summary Table
summary_table <- bind_rows(top_pink, top_blue, top_purple) |>
  select(Category = bi_class, Descriptor, Country = COUNTRY, Province = NAME_1, District = NAME_2, Mean_Hazard = mean_intrinsic_hazard, Mean_Richness = mean_district_richness)

write_csv(summary_table, here("output", "tables", "district_subset_table.csv"))

# Checks for Review -------------------------------------------------------
# 1. Confirms/quantifies the Results text claim
print(risk_summary)
print(stats)   # specifically check West Africa's row/column in the pairwise comparison

# 2. Confirms/refutes the Fig. 3 caption claim
summary_table |> filter(Category == "3-3")

district_classed |> 
  mutate(iso3c = str_sub(GID_2, end = 3)) |>
  filter(bi_class == "3-3") |> 
  count(COUNTRY)

# Assign each 20km cell to a country (cells are already equal-area under Mollweide)
cell_pts <- vect(data_bivariate, geom = c("x", "y"), crs = crs(r_template))
cell_countries <- terra::extract(vect(world_sf |> select(admin)), cell_pts)

country_area_pct <- data_bivariate |>
  mutate(country = cell_countries$admin) |>
  drop_na(country) |>
  count(country, bi_class) |>
  group_by(country) |>
  mutate(total_cells = sum(n), pct_area = round(100 * n / total_cells, 1)) |>
  ungroup() |>
  filter(total_cells >= 20) |>   # guard against noisy % from tiny/sparse countries
  filter(bi_class == "3-3") |>
  arrange(desc(pct_area))

print(country_area_pct, n = 20)

district_classed |> group_by(bi_class) |>
  summarise(n = n(), mean_host_prob = round(mean(mean_intrinsic_hazard), 4),
            mean_richness = round(mean(mean_district_richness), 1), .groups = "drop")
summary(lm(mean_intrinsic_hazard ~ mean_district_richness, data = table_out))
summary_table

