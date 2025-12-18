# ---
# 05_temporal_bias_analysis.R
#
# Purpose: 
# 1. Quantify temporal sampling effort (Individuals & Studies).
# 2. Track species accumulation curves over time.
# 3. Model surveillance trends by continent using GAMs.
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, lubridate, scales, ggrepel, patchwork, mgcv, countrycode)

# Create output directory
output_dir <- here("output", "analysis_1")
dir.create(output_dir, showWarnings = FALSE)

# Load Database
# Using the harmonised version from our workflow
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2025-12-02.rds"))

host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
citations_data <- arha_db$citations

# --- 2. Data Preparation ---

# A. Host Temporal Data
# We coalesce sample date with publication date to maximize data utility
host_temporal <- host_data %>%
  select(study_id, host_record_id, start_date, number_of_hosts, host_species, iso3c) %>%
  left_join(citations_data %>% select(study_id, publication_year), by = "study_id") %>%
  mutate(
    year = coalesce(year(start_date), publication_year),
    decade = floor(year / 10) * 10
  ) %>%
  # Filter for the modern era of virology
  filter(year >= 1960, year <= 2025)

# B. Pathogen Temporal Data
pathogen_temporal <- pathogen_data %>%
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) %>%
  select(host_record_id, pathogen_family, pathogen_species_cleaned) %>%
  inner_join(host_temporal %>% select(host_record_id, year), by = "host_record_id")


# --- 3. Analysis 1: Effort Over Time ---

# A. Aggregate
effort_year <- host_temporal %>%
  group_by(year) %>%
  summarise(
    n_hosts = sum(number_of_hosts, na.rm = TRUE),
    n_studies = n_distinct(study_id)
  )

effort_decade <- host_temporal %>%
  group_by(decade) %>%
  summarise(
    n_hosts = sum(number_of_hosts, na.rm = TRUE),
    n_studies = n_distinct(study_id)
  )

# B. Visualisation (Stacked Plot)
# Annotations for key historical context
events <- tibble(
  year = c(1969, 1993),
  label = c("Lassa described\n(1969)", "Sin Nombre\nOutbreak (1993)")
)

p_hosts <- ggplot(effort_year, aes(x = year, y = n_hosts)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_vline(data = events, aes(xintercept = year), linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = events, aes(y = max(effort_year$n_hosts)*0.8, label = label), 
                  nudge_x = 5, size = 3) +
  scale_y_continuous(labels = comma) +
  labs(title = "A) Total Rodent/Shrew Individuals Sampled", x = NULL, y = "Count") +
  theme_minimal()

p_studies <- ggplot(effort_year, aes(x = year, y = n_studies)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick", size = 1.5) +
  geom_vline(data = events, aes(xintercept = year), linetype = "dashed", alpha = 0.5) +
  labs(title = "B) Unique Studies Published", x = "Year", y = "Count") +
  theme_minimal()

p_combined <- p_hosts / p_studies
ggsave(here(output_dir, "temporal_effort_combined.png"), p_combined, width = 8, height = 8)


# --- 4. Analysis 2: Diversity Over Time ---

diversity_year <- host_temporal %>%
  group_by(year) %>%
  summarise(n_host_spp = n_distinct(host_species)) %>%
  left_join(
    pathogen_temporal %>%
      group_by(year, pathogen_family) %>%
      summarise(n = n_distinct(pathogen_species_cleaned), .groups="drop") %>%
      pivot_wider(names_from = pathogen_family, values_from = n, values_fill = 0),
    by = "year"
  ) %>%
  pivot_longer(cols = -year, names_to = "taxa", values_to = "count")

p_diversity <- ggplot(diversity_year, aes(x = year, y = count, color = taxa)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("n_host_spp" = "black", "Arenaviridae" = "#E69F00", "Hantaviridae" = "#56B4E9"),
    labels = c("Arenaviruses", "Hantaviruses", "Host Species")
  ) +
  labs(title = "Species Discovery & Surveillance Breadth", y = "Unique Species Sampled") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here(output_dir, "temporal_diversity.png"), p_diversity, width = 8, height = 6)


# --- 5. Analysis 3: Continental Trends (GAMs) ---

message("--- Modelling Continental Trends ---")

# Add Continents
effort_continent <- host_temporal %>%
  drop_na(iso3c) %>%
  mutate(continent = countrycode(iso3c, "iso3c", "continent")) %>%
  drop_na(continent) %>%
  filter(continent != "Oceania") %>% # Remove sparse data
  group_by(year, continent) %>%
  summarise(n_hosts = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

# Fit GAM (Negative Binomial to handle overdispersed count data)
# We model a smooth trend per continent
gam_trend <- gam(n_hosts ~ s(year, by = as.factor(continent), k = 10) + as.factor(continent),
                 data = effort_continent,
                 family = nb())

# Predict for Plotting (Manual Link Scale Calculation for correct CIs)
pred_grid <- expand.grid(
  year = 1960:2024,
  continent = unique(effort_continent$continent)
)

preds <- predict(gam_trend, newdata = pred_grid, type = "link", se.fit = TRUE)
inv_link <- gam_trend$family$linkinv

plot_data <- pred_grid %>%
  mutate(
    fit = inv_link(preds$fit),
    lwr = inv_link(preds$fit - 1.96 * preds$se.fit),
    upr = inv_link(preds$fit + 1.96 * preds$se.fit)
  )

p_trends <- ggplot(plot_data, aes(x = year, y = fit, color = continent, fill = continent)) +
  geom_point(data = effort_continent, aes(y = n_hosts), alpha = 0.3, size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~continent, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Surveillance Trends by Continent (GAM)",
       y = "Individuals Sampled (Predicted vs Observed)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(here(output_dir, "temporal_trends_gam.png"), p_trends, width = 10, height = 7)

message("✅ Temporal Bias Analysis Complete.")