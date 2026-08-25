# ==============================================================================
# 07_03_post_hoc_investigations.R
# Purpose: Interrogate the model for Viral Family differences and Geographic Bias.
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, brms, tidybayes, countrycode, cowplot, rnaturalearth)

fit_dyadic <- read_rds(here("output", "models", "brms_dyadic_N49k.rds"))

# Load Lookup Data
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-08-17.rds"))
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))

# 1. Viral Family Differences ---------------------------------------------
# Extract Random Effects for Pathogens
viral_rfx <- fit_dyadic |>
  spread_draws(r_pathogen_species_cleaned[pathogen, term]) |>
  median_qi() |>
  ungroup()

# Add Family Metadata
viral_meta <- arha_db$pathogen |>
  distinct(pathogen_species_cleaned, pathogen_family) |>
  filter(!is.na(pathogen_species_cleaned))

plot_data_virus <- viral_rfx |>
  mutate(pathogen = str_replace_all(pathogen, "\\.", " ")) |>
  left_join(viral_meta, by = c("pathogen" = "pathogen_species_cleaned")) |>
  mutate(pathogen = fct_reorder(pathogen, r_pathogen_species_cleaned)) |>
  drop_na(pathogen_family)

# Plot: Caterpillar Plot coloured by Family
p_virus <- ggplot(plot_data_virus, aes(y = pathogen, x = r_pathogen_species_cleaned, colour = pathogen_family)) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper), linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "A) Viral Detectability (Random Intercepts)",
       subtitle = "Deviation from global baseline probability",
       x = "Random Intercept (Log-Odds)",
       y = NULL, colour = "Family") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 11, colour = "grey30", margin = margin(b = 10)),
        axis.text.y = element_text(size = 6),
        panel.grid.major.y = element_blank(),
        legend.position = "none")

# Statistical Check: Is there a mean difference?
p_virus_box <- ggplot(plot_data_virus, aes(x = pathogen_family, y = r_pathogen_species_cleaned, fill = pathogen_family)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(title = "B) Family Comparison",
       subtitle = "Distribution of intercepts",
       x = NULL, y = "Random Intercept") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 11, colour = "grey30", margin = margin(b = 10)),
        legend.position = "none")

# Descriptive: per-species points, violin instead of box
p_virus_violin <- ggplot(plot_data_virus, aes(x = pathogen_family, y = r_pathogen_species_cleaned, fill = pathogen_family)) +
  geom_violin(alpha = 0.6, scale = "count") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(title = "B) Family Comparison", subtitle = "Distribution of per-species intercepts (violin width = N)",
       x = NULL, y = "Random Intercept") +
  theme_minimal(base_size = 12) + theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)))

# C) Statistical evidence: full per-draw propagated family difference
viral_draws_full <- fit_dyadic |>
  spread_draws(r_pathogen_species_cleaned[pathogen, term]) |>
  mutate(pathogen = str_replace_all(pathogen, "\\.", " ")) |>
  left_join(viral_meta, by = c("pathogen" = "pathogen_species_cleaned")) |>
  drop_na(pathogen_family)

family_diff_draws <- viral_draws_full |>
  group_by(.draw, pathogen_family) |>
  summarise(family_mean = mean(r_pathogen_species_cleaned), .groups = "drop") |>
  pivot_wider(names_from = pathogen_family, values_from = family_mean) |>
  mutate(diff = Arenaviridae - Hantaviridae)

p_family_diff <- ggplot(family_diff_draws, aes(x = diff, y = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  stat_halfeye(fill = "#CC79A7", alpha = 0.7, .width = c(0.5, 0.95)) +
  labs(title = "C) Family-Level Difference",
       subtitle = "Posterior: Arenaviridae − Hantaviridae",
       x = "Difference in Random Intercept (Log-Odds)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 11, colour = "grey30", margin = margin(b = 10)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Combine all three
right_col <- plot_grid(p_virus_violin,
                       p_family_diff,
                       ncol = 1)
fig_posthoc_virus <- plot_grid(p_virus,
                               right_col,
                               ncol = 2,
                               rel_widths = c(1.3, 1))
ggsave(here("output", "figures", "supplementary_figure_s5.png"), fig_posthoc_virus, width = 15, height = 8, bg = "white")

# Full (unsummarized) posterior draws, not median_qi()
viral_draws_full <- fit_dyadic |>
  spread_draws(r_pathogen_species_cleaned[pathogen, term]) |>
  mutate(pathogen = str_replace_all(pathogen, "\\.", " ")) |>
  left_join(viral_meta, by = c("pathogen" = "pathogen_species_cleaned")) |>
  drop_na(pathogen_family)

# Family-level mean per posterior draw, then the difference per draw
family_diff_draws <- viral_draws_full |>
  group_by(.draw, pathogen_family) |>
  summarise(family_mean = mean(r_pathogen_species_cleaned), .groups = "drop") |>
  pivot_wider(names_from = pathogen_family, values_from = family_mean) |>
  mutate(diff = Arenaviridae - Hantaviridae)

family_diff_draws |>
  summarise(median_diff = median(diff),
            lower = quantile(diff, 0.025),
            upper = quantile(diff, 0.975),
            pd = mean(diff > 0))

# The posterior of the family-level difference
p_family_diff <- ggplot(family_diff_draws, aes(x = diff)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  stat_halfeye(fill = "#CC79A7", alpha = 0.7, .width = c(0.5, 0.95)) +
  labs(title = "Arenaviridae − Hantaviridae baseline detectability",
       subtitle = paste0("Posterior difference in family-mean random intercept (pd = ", 
                         round(mean(family_diff_draws$diff < 0) * 100, 1), "%)"),
       x = "Difference (Log-Odds)", y = NULL) +
  theme_minimal()

# 2. Geographic Residuals -------------------------------------------------
# Host-Level Predictions
model_frame <- fit_dyadic$data |> as_tibble()
preds <- fitted(fit_dyadic, summary = TRUE) |> as_tibble()

# Residual = Observed (0/1) - Predicted Probability
residuals_df <- bind_cols(model_frame, preds) |>
  mutate(residual = is_reservoir - Estimate) |>
  group_by(tip_label) |>
  summarise(mean_residual = mean(residual),
            n_pairs = n()) |>
  ungroup()

# Attach Geography (Centroid Latitude)
geo_residuals <- residuals_df |>
  left_join(host_traits |> select(tip_label, species), by = "tip_label") |>
  left_join(arha_db$host |> 
              group_by(host_species) |> 
              summarise(lat = mean(latitude, na.rm = TRUE)) |>
              select(species = host_species, lat), 
            by = "species") |>
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
resid_map_data <- geo_residuals |>
  left_join(arha_db$host |> distinct(species = host_species, iso3c), by = "species") |>
  drop_na(iso3c) |>
  group_by(iso3c) |>
  summarise(avg_country_resid = mean(mean_residual))

world <- ne_countries(scale = "medium", returnclass = "sv") |>
  merge(resid_map_data, by.x = "iso_a3", by.y = "iso3c", all.x = TRUE)

p_resid_map <- ggplot() +
  geom_spatvector(data = world, aes(fill = avg_country_resid), colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       na.value = "grey90", name = "Mean\nResidual") +
  labs(title = "Geographic Distribution of Unexplained Risk",
       subtitle = "Red = More Reservoirs observed than predicted by traits") +
  theme_void()

ggsave(here("output", "figures", "posthoc_geographic_residuals.png"), p_resid_map, width = 10, height = 6, bg = "white")
