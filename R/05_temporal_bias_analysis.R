# ---
# 05_temporal_bias_analysis.R
#
# Purpose: 
# 1. Quantify temporal sampling effort (Individuals & Studies).
# 2. Track species accumulation curves over time.
# 3. Model surveillance trends by continent using GAMs.
# ---

# 1. Setup ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, lubridate, scales, ggrepel, patchwork, mgcv, countrycode)
output_dir <- here("output", "figures")

# Load Database
# Using the harmonised version from our workflow
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-01-09.rds"))

host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
citations_data <- arha_db$citations


# 2. Data Prep ------------------------------------------------------------
# Host Temporal Data
# We coalesce sample date with publication date to maximise data utility
host_temporal <- host_data |>
  select(study_id, host_record_id, start_date, number_of_hosts, host_species, iso3c) |>
  left_join(citations_data |> select(study_id, publication_year), by = "study_id") |>
  mutate(year = coalesce(year(start_date), publication_year),
         decade = floor(year / 10) * 10) |>
  filter(year >= 1960, year <= 2025)

# Pathogen Temporal Data
pathogen_temporal <- pathogen_data |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  select(host_record_id, pathogen_family, pathogen_species_cleaned) |>
  inner_join(host_temporal %>% select(host_record_id, year), by = "host_record_id")

# 3. Global Effort and Accumulation ---------------------------------------
# Aggregate
effort_year <- host_temporal |>
  group_by(year) |>
  summarise(n_hosts = sum(number_of_hosts, na.rm = TRUE),
            n_studies = n_distinct(study_id))

accumulation_curve <- host_temporal |>
  arrange(year) |>
  distinct(host_species, .keep_all = TRUE) |>
  count(year, name = "new_species") |>
  mutate(cumulative_species = cumsum(new_species)) |>
  right_join(tibble(year = 1960:2025), by = "year") |>
  arrange(year) |>
  fill(cumulative_species, .direction = "down") |>
  mutate(cumulative_species = replace_na(cumulative_species, 0))

# Key Events for Annotation
events <- tibble(year = c(1969, 1993, 2014, 2019),
                 label = c("Lassa\nDescribed", "Sin Nombre\nOutbreak", "W. Africa\nEbola", "SARS-CoV-2\nDetected"))

# Plot A: Effort
p_s1_a <- ggplot(effort_year, aes(x = year, y = n_hosts)) +
  geom_col(fill = "grey40", alpha = 0.8) +
  geom_vline(data = events, aes(xintercept = year), linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = events, aes(y = max(effort_year$n_hosts)*0.9, label = label), 
                  size = 3, nudge_x = 2, direction = "y") +
  scale_y_continuous(labels = comma) +
  labs(title = "A) Global Sampling Effort", y = "Individuals Sampled", x = NULL) +
  theme_minimal()

# Plot B: Accumulation
p_s1_b <- ggplot(accumulation_curve, aes(x = year, y = cumulative_species)) +
  geom_line(color = "orange", linewidth = 1.2) +
  labs(title = "B) Host Species Accumulation", y = "Cumulative Host Species", x = "Year") +
  theme_minimal()

fig_s1 <- p_s1_a / p_s1_b
ggsave(here(output_dir, "temporal_effort.png"), fig_s1, width = 8, height = 8, bg = "white")

# 4. Continental Trends ---------------------------------------------------
effort_continent <- host_temporal |>
  drop_na(iso3c) |>
  mutate(continent = countrycode(iso3c, "iso3c", "continent")) |>
  drop_na(continent) |>
  filter(continent != "Oceania") |>
  group_by(year, continent) |>
  summarise(n_hosts = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

# Fit GAM (Negative Binomial to handle overdispersed count data)
# We model a smooth trend per continent
gam_trend <- gam(n_hosts ~ s(year, by = as.factor(continent), k = 10) + as.factor(continent),
                 data = effort_continent,
                 family = nb())

# Predict for Plotting (Manual Link Scale Calculation for correct CIs)
pred_grid <- expand.grid(year = 1960:2025,
                         continent = unique(effort_continent$continent))

preds <- predict(gam_trend, newdata = pred_grid, type = "link", se.fit = TRUE)
inv_link <- gam_trend$family$linkinv

max_observed <- effort_continent %>%
  group_by(continent) %>%
  summarise(max_obs = max(n_hosts))

plot_data <- pred_grid |>
  mutate(fit = inv_link(preds$fit),
         lwr = inv_link(preds$fit - 1.96 * preds$se.fit),
         upr = inv_link(preds$fit + 1.96 * preds$se.fit))

plot_data_clamped <- plot_data %>%
  left_join(max_observed, by = "continent") %>%
  mutate(upr = pmin(upr, max_obs * 2))

p_trends <- ggplot(plot_data_clamped, aes(x = year, y = fit)) +
  geom_point(data = effort_continent, aes(y = n_hosts), 
             alpha = 0.4, size = 1.2, colour = "grey50") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = continent), alpha = 0.2) +
  geom_line(aes(colour = continent), linewidth = 1) +
  facet_wrap(~continent, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Surveillance Trajectories by Continent",
       subtitle = "Modelled sampling effort (GAM, Negative Binomial)",
       y = "Individuals Sampled",
       x = "Year") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 11))

ggsave(here(output_dir, "temporal_trends_gam.png"), p_trends, width = 10, height = 7)
