# ---
# 06_genetic_bias_analysis.R
#
# Purpose: 
# 1. Assess the "Genetic Gap": Where do we have surveillance but no sequences?
# 2. Map sequencing completeness (Host vs Pathogen) globally.
# 3. Compare PCR-positive counts to available sequences (Taxonomic gaps).
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, terra, tidyterra, rnaturalearth, 
               cowplot, patchwork, scales, biscale)

# Create output directory
output_dir <- here("output", "analysis_1", "genetic_bias")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Database
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2025-12-02.rds"))

host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
sequence_data <- arha_db$sequence


# --- 2. Geographic Sequencing Gaps ---

message("--- Analyzing Geographic Sequencing Gaps ---")

# A. Summarise Total Surveillance per Country
# (Denominator: How many hosts did we sample?)
hosts_per_country <- host_data %>%
  drop_na(iso3c) %>%
  group_by(iso3c, country) %>%
  summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

# B. Summarise Sequences per Country
# We need to link sequences back to their geographic origin
sequences_geo <- sequence_data %>%
  left_join(host_data %>% select(host_record_id, iso3c_host = iso3c), by = "host_record_id") %>%
  mutate(iso3c = coalesce(iso3c, iso3c_host)) %>% # Use sequence ISO, fallback to host ISO
  drop_na(iso3c) %>%
  group_by(iso3c, sequence_type) %>%
  summarise(n_sequences = n_distinct(accession_primary), .groups = "drop") %>%
  pivot_wider(names_from = sequence_type, values_from = n_sequences, values_fill = 0) %>%
  rename(n_seq_host = Host, n_seq_pathogen = Pathogen)

# C. Merge and Calculate Completeness
geo_gap_data <- hosts_per_country %>%
  left_join(sequences_geo, by = "iso3c") %>%
  mutate(
    n_seq_host = replace_na(n_seq_host, 0),
    n_seq_pathogen = replace_na(n_seq_pathogen, 0),
    # Proportions (capped at 1 for visualization logic)
    prop_host = pmin(n_seq_host / total_hosts_sampled, 1),
    prop_pathogen = pmin(n_seq_pathogen / total_hosts_sampled, 1)
  )

# --- 3. Bivariate Map (Host vs Pathogen Completeness) ---

message("--- Generating Bivariate Map ---")

# Prepare Map Data
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  select(iso_a3) %>%
  left_join(geo_gap_data, by = c("iso_a3" = "iso3c"))

# Define Bivariate Classes
data_bivariate <- bi_class(world_sf, x = prop_host, y = prop_pathogen, 
                           style = "quantile", dim = 3)

# Create Map
map_bi <- ggplot() +
  geom_sf(data = data_bivariate, aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet", dim = 3) +
  coord_sf(crs = "+proj=moll") +
  labs(title = "Genetic Completeness: Host vs Pathogen",
       subtitle = "Darker Purple = High completeness for both\nBlue = Host only, Pink = Pathogen only") +
  theme_void()

# Create Legend
legend_bi <- bi_legend(pal = "DkViolet", dim = 3, xlab = "Host %", ylab = "Pathogen %")

# Combine
final_map <- ggdraw() +
  draw_plot(map_bi, 0, 0, 1, 1) +
  draw_plot(legend_bi, 0.1, 0.1, 0.2, 0.2)

ggsave(here(output_dir, "map_bivariate_sequencing.png"), final_map, width = 12, height = 7)


# --- 4. Taxonomic Sequencing Gaps (Dumbbell Plots) ---

message("--- Analyzing Taxonomic Sequencing Gaps ---")

# A. Host Taxonomic Gaps
# Compare Total Sampled vs Total Sequenced
host_gap_data <- host_data %>%
  group_by(host_species) %>%
  summarise(n_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
  left_join(
    sequence_data %>%
      filter(sequence_type == "Host") %>%
      left_join(host_data, by = "host_record_id") %>%
      group_by(host_species) %>%
      summarise(n_sequenced = n_distinct(accession_primary)),
    by = "host_species"
  ) %>%
  mutate(n_sequenced = replace_na(n_sequenced, 0),
         prop_seq = n_sequenced / n_sampled) %>%
  filter(n_sampled > 1000) %>% # Filter for major hosts
  mutate(host_label = paste0(host_species, " (N=", comma(n_sampled), ")"),
         host_label = fct_reorder(host_label, prop_seq))

p_host_gap <- ggplot(host_gap_data, aes(y = host_label)) +
  geom_segment(aes(x = 0, xend = prop_seq, y = host_label, yend = host_label), color = "grey") +
  geom_point(aes(x = prop_seq), color = "steelblue", size = 3) +
  scale_x_continuous(labels = percent) +
  labs(title = "Host Sequencing Completeness",
       subtitle = "Proportion of sampled individuals with available genetic data (>1000 sampled)",
       x = "Proportion Sequenced", y = NULL) +
  theme_minimal()

ggsave(here(output_dir, "host_sequencing_gap.png"), p_host_gap, width = 8, height = 10)


# B. Pathogen Taxonomic Gaps
# Compare PCR Positives vs Sequences
pathogen_gap_data <- pathogen_data %>%
  filter(number_positive > 0) %>%
  group_by(pathogen_species_cleaned, pathogen_family) %>%
  summarise(n_pcr_pos = sum(number_positive, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    sequence_data %>%
      filter(sequence_type == "Pathogen") %>%
      left_join(pathogen_data, by = "pathogen_record_id") %>%
      group_by(pathogen_species_cleaned) %>%
      summarise(n_seq = n_distinct(accession_primary)),
    by = "pathogen_species_cleaned"
  ) %>%
  mutate(n_seq = replace_na(n_seq, 0)) %>%
  filter(n_pcr_pos > 10) %>% # Filter for relevant pathogens
  mutate(pathogen_species_cleaned = fct_reorder(pathogen_species_cleaned, n_pcr_pos))

p_pathogen_dumbbell <- ggplot(pathogen_gap_data, aes(y = pathogen_species_cleaned)) +
  geom_segment(aes(x = n_seq, xend = n_pcr_pos, y = pathogen_species_cleaned, yend = pathogen_species_cleaned),
               color = "grey", linewidth = 0.8) +
  geom_point(aes(x = n_seq, color = "Sequences"), size = 3) +
  geom_point(aes(x = n_pcr_pos, color = "PCR+ Detections"), size = 3) +
  scale_x_log10(labels = comma) +
  scale_color_manual(values = c("PCR+ Detections" = "firebrick", "Sequences" = "steelblue")) +
  facet_wrap(~ pathogen_family, scales = "free_y") +
  labs(title = "Pathogen Genetic Gap",
       subtitle = "Discrepancy between detections (PCR) and available sequences (GenBank)",
       x = "Count (log scale)", y = NULL, color = "Metric") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(here(output_dir, "pathogen_sequencing_gap.png"), p_pathogen_dumbbell, width = 12, height = 8)

message("✅ Genetic Bias Analysis Complete.")