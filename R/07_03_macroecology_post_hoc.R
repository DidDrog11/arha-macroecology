# ==============================================================================
# 07_03_post_hoc_investigations.R
# Purpose: Interrogate the model for Viral Family differences and Geographic Bias.
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, brms, tidybayes, countrycode, cowplot)

fit_dyadic <- read_rds(here("output", "models", "brms_dyadic_N19k.rds")) 
# fit_dyadic <- read_rds(here("output", "models", "brms_dyadic_N246k.rds"))

# Load Lookup Data
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-01-09.rds"))
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))

# 1. Viral Family Differences ---------------------------------------------
# Extract Random Effects for Pathogens
viral_rfx <- fit_dyadic %>%
  spread_draws(r_pathogen_species_cleaned[pathogen, term]) %>%
  median_qi() %>% # Get Median and CI for each virus
  ungroup()

# Add Family Metadata
viral_meta <- arha_db$pathogen %>%
  distinct(pathogen_species_cleaned, pathogen_family) %>%
  filter(!is.na(pathogen_species_cleaned))

plot_data_virus <- viral_rfx %>%
  mutate(pathogen = str_replace_all(pathogen, "\\.", " ")) |>
  left_join(viral_meta, by = c("pathogen" = "pathogen_species_cleaned")) %>%
  mutate(pathogen = fct_reorder(pathogen, r_pathogen_species_cleaned)) %>%
  drop_na(pathogen_family)

# Plot: Caterpillar Plot coloured by Family
p_virus <- ggplot(plot_data_virus, aes(y = pathogen, x = r_pathogen_species_cleaned, colour = pathogen_family)) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper), size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "A) Viral Detectability (Random Intercepts)",
       subtitle = "Deviation from global baseline probability",
       x = "Random Intercept (Log-Odds)",
       y = NULL, color = "Family") +
  theme_minimal_vgrid() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "none")

# Statistical Check: Is there a mean difference?
p_virus_box <- ggplot(plot_data_virus, aes(x = pathogen_family, y = r_pathogen_species_cleaned, fill = pathogen_family)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(title = "B) Family Comparison",
       subtitle = "Distribution of intercepts",
       x = NULL, y = "Random Intercept") +
  theme_minimal() +
  theme(legend.position = "none")

# Combine
fig_posthoc_virus <- plot_grid(p_virus, p_virus_box, rel_widths = c(1.5, 1))
ggsave(here("output", "figures", "posthoc_viral_families.png"), fig_posthoc_virus, width = 12, height = 8)

# 2. Geographic Residuals -------------------------------------------------
# Host-Level Predictions
model_frame <- fit_dyadic$data %>% 
  as_tibble() %>%
  mutate(row_id = row_number())

preds <- fitted(fit_dyadic, summary = TRUE) %>% 
  as_tibble() %>%
  mutate(row_id = row_number())

# Residual = Observed (0/1) - Predicted Probability
residuals_df <- model_frame %>%
  left_join(preds, by = "row_id") %>%
  mutate(residual = is_reservoir - Estimate) %>%
  group_by(tip_label) %>%
  summarise(mean_residual = mean(residual),
            n_pairs = n()) %>%
  ungroup()

# Attach Geography (Centroid Latitude)
# Pulled from host_traits
geo_residuals <- residuals_df %>%
  left_join(host_traits %>% select(tip_label, species), by = "tip_label") %>%
  left_join(arha_db$host %>% 
              group_by(host_species) %>% 
              summarise(lat = mean(latitude, na.rm=TRUE)) |>
              select(species = host_species, lat), 
            by = "species") %>%
  drop_na(lat)

# Plot: Residuals vs Latitude
p_resid_lat <- ggplot(geo_residuals, aes(x = lat, y = mean_residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", colour = "steelblue", fill = "steelblue", alpha = 0.1) +
  labs(title = "Model Residuals by Latitude",
       subtitle = "Positive = Model Under-predicts Risk (Missing Reservoirs?)",
       y = "Mean Residual (Obs - Pred)",
       x = "Latitude") +
  theme_minimal()

# Map of Residuals
# Link species to countries via ArHa DB
resid_map_data <- geo_residuals %>%
  left_join(arha_db$host %>% distinct(species = host_species, iso3c), by = "species") %>%
  drop_na(iso3c) %>%
  group_by(iso3c) %>%
  summarise(avg_country_resid = mean(mean_residual))

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  left_join(resid_map_data, by = c("iso_a3" = "iso3c"))

p_resid_map <- ggplot(world) +
  geom_sf(aes(fill = avg_country_resid), colour = "white", size = 0.1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Mean\nResidual") +
  labs(title = "Geographic Distribution of Unexplained Risk",
       subtitle = "Red = More Reservoirs observed than predicted by traits") +
  theme_void()
