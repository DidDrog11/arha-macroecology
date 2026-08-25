# ==============================================================================
# 06_genetic_bias_analysis.R
#
# Purpose: 
# 1. Generate Supplementary Figure S2: Geographic Sequencing Gap (Bivariate Map).
# 2. Generate Supplementary Figure S3: Pathogen Genetic Disconnect (Dumbbell Plot).
# ==============================================================================

# 1. Setup ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, terra, tidyterra, here, sf, rnaturalearth, cowplot, biscale, scales, patchwork, countrycode)

# Output Directory
output_dir <- here("output", "figures")
dir.create(output_dir, showWarnings = FALSE)

# Load Database
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-08-17.rds"))

host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
sequence_data <- arha_db$sequence


# 2. Geographic Sampling Gap ----------------------------------------------
# Calculate Surveillance Effort (Denominator)
hosts_per_country <- host_data |>
  drop_na(iso3c) |>
  group_by(iso3c) |>
  summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

# Calculate Sequencing Effort (Numerator)
# Link sequences to geography via the host record
sequences_geo <- sequence_data |>
  left_join(host_data |> select(host_record_id, iso3c_host = iso3c), by = "host_record_id") |>
  mutate(iso3c = coalesce(iso3c, iso3c_host)) |>
  drop_na(iso3c) |>
  group_by(iso3c, sequence_type) |>
  summarise(n_sequences = n_distinct(accession_primary), .groups = "drop") |>
  pivot_wider(names_from = sequence_type, values_from = n_sequences, values_fill = 0) |>
  rename(n_seq_host = Host, n_seq_pathogen = Pathogen)

# Merge & Calculate Completeness Ratios
geo_gap_data <- hosts_per_country |>
  left_join(sequences_geo, by = "iso3c") |>
  mutate(n_seq_host = replace_na(n_seq_host, 0),
         n_seq_pathogen = replace_na(n_seq_pathogen, 0),
         prop_host = pmin(n_seq_host / total_hosts_sampled, 1),
         prop_pathogen = pmin(n_seq_pathogen / total_hosts_sampled, 1))

# Create Bivariate Map
world_vect <- ne_countries(scale = "medium", returnclass = "sv") |>
  filter(admin != "Antarctica") |>
  select(iso_a3) |>
  left_join(geo_gap_data, by = c("iso_a3" = "iso3c"))

world_sf <- sf::st_as_sf(world_vect) %>%
  drop_na(prop_host, prop_pathogen)

# Define Classes (Quantile breaks for balanced colours)
data_bivariate <- bi_class(world_sf,
                           x = prop_host, 
                           y = prop_pathogen, 
                           style = "fisher", 
                           dim = 4)

map_bi <- ggplot() +
  geom_sf(data = data_bivariate, aes(fill = bi_class), colour = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet2", dim = 4) +
  coord_sf(crs = "+proj=moll") +
  labs(title = "Genetic Completeness",
       subtitle = "Comparison of Host vs. Pathogen sequencing rates per country") +
  theme_void()

# Legend
legend_bi <- bi_legend(pal = "DkViolet2", dim = 4, xlab = "Host Seq %", ylab = "Virus Seq %")

# Assemble
fig_s2 <- ggdraw() +
  draw_plot(map_bi, 0, 0, 1, 1) +
  draw_plot(legend_bi, 0.1, 0.15, 0.2, 0.2)

ggsave(here(output_dir, "supplementary_figure_s9.png"), fig_s2, width = 10, height = 6, bg = "white")

# Which pathogen records have at least one pathogen sequence
sequenced_records <- arha_db$sequence |>
  filter(sequence_type == "Pathogen") |>
  distinct(pathogen_record_id) |>
  mutate(has_sequence = TRUE)

pathogen_gap_data <- pathogen_data |>
  filter(str_detect(assay, "PCR|Sequencing|Isolation|Culture")) |>
  filter(number_positive > 0) |>
  drop_na(pathogen_species_ncbi) |> # Removes the unclassified NA rows
  left_join(sequenced_records, by = "pathogen_record_id") |>
  mutate(has_sequence = replace_na(has_sequence, FALSE)) |>
  group_by(pathogen_species_ncbi, pathogen_family) |>
  summarise(
    n_pcr_pos = sum(number_positive, na.rm = TRUE),
    n_seq = sum(number_positive[has_sequence == TRUE], na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n_pcr_pos > 10) |> 
  filter(str_detect(pathogen_family, "Hanta|Arena")) |>
  mutate(pathogen_species_ncbi = fct_reorder(pathogen_species_ncbi, n_pcr_pos))

# Dumbbell Plot
fig_s10 <- ggplot(pathogen_gap_data, aes(y = pathogen_species_ncbi)) +
  geom_segment(aes(x = n_seq, xend = n_pcr_pos, y = pathogen_species_ncbi, yend = pathogen_species_ncbi),
               colour = "grey60", linewidth = 0.8) +
  geom_point(aes(x = n_seq, colour = "Sequences"), size = 3) +
  geom_point(aes(x = n_pcr_pos, colour = "PCR+ Detections"), size = 3) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = c("PCR+ Detections" = "#D55E00", "Sequences" = "#0072B2")) +
  facet_wrap(~ pathogen_family, scales = "free_y") +
  labs(title = "The Pathogen Genetic Gap",
       subtitle = "Discrepancy between PCR detections (Orange) and available GenBank sequences (Blue)",
       x = "Count of Individual Hosts (log scale)", y = NULL, colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        axis.text.y = element_text(size = 9))

ggsave(here(output_dir, "supplementary_figure_s10.png"), fig_s10, width = 10, height = 8, bg = "white")

# Join to pathogen and calculate totals
gap_total <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(str_detect(assay, "PCR|Sequencing|Isolation")) |> 
  filter(number_positive > 0) |>
  left_join(sequenced_records, by = "pathogen_record_id") |>
  mutate(has_sequence = replace_na(has_sequence, FALSE)) |>
  summarise(total_positives = sum(number_positive, na.rm = TRUE),
            # If a record has sequences, we assume all 'number_positive' individuals in that record were sequenced.
            # if it doesn't, 0 were sequenced.
            total_sequenced = sum(number_positive[has_sequence == TRUE], na.rm = TRUE),
            pct_sequenced = round((total_sequenced / total_positives) * 100, 1))

print(gap_total)

# Family Comparison
gap_family <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(str_detect(assay, "PCR|Sequencing|Isolation")) |> 
  filter(number_positive > 0) |>
  left_join(sequenced_records, by = "pathogen_record_id") |>
  mutate(has_sequence = replace_na(has_sequence, FALSE)) |>
  group_by(pathogen_family) |>
  summarise(total_positives = sum(number_positive, na.rm = TRUE),
            total_sequenced = sum(number_positive[has_sequence == TRUE], na.rm = TRUE),
            pct_sequenced = round((total_sequenced / total_positives) * 100, 1))

print(gap_family)

# Identify unique Host-Virus associations across all data (including serology)
all_dyads <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(number_positive > 0) |>
  left_join(arha_db$host |> select(host_record_id, host_species), by = "host_record_id") |>
  distinct(host_species, pathogen_species_cleaned)

# Identify unique Host-Virus associations that have at least one sequence
sequenced_dyads <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(number_positive > 0) |>
  left_join(sequenced_records, by = "pathogen_record_id") |>
  filter(has_sequence == TRUE) |>
  left_join(arha_db$host |> select(host_record_id, host_species), by = "host_record_id") |>
  distinct(host_species, pathogen_species_cleaned)

# Calculate the gap
total_associations <- nrow(all_dyads)
supported_by_seq <- nrow(sequenced_dyads)
pct_no_sequence <- round((1 - (supported_by_seq / total_associations)) * 100, 1)

regional_gap <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(str_detect(assay, "PCR|Sequencing|Isolation")) |> 
  filter(number_positive > 0) |>
  left_join(sequenced_records, by = "pathogen_record_id") |>
  mutate(has_sequence = replace_na(has_sequence, FALSE)) |>
  left_join(arha_db$host |> select(host_record_id, iso3c), by = "host_record_id") |>
  mutate(region = countrycode(iso3c, origin = "iso3c", destination = "un.regionsub.name")) |>
  group_by(region) |>
  summarise(total_positives = sum(number_positive, na.rm = TRUE),
            total_sequenced = sum(number_positive[has_sequence == TRUE], na.rm = TRUE)) |>
  filter(total_positives > 50) |> 
  mutate(pct_unsequenced = round((1 - (total_sequenced / total_positives)) * 100, 1)) |>
  arrange(desc(pct_unsequenced))

# 3. Revision Table 3: Regional genetic sequencing completeness -----------
revision_table_3 <- regional_gap |>
  rename(region = region,
         n_positive_hosts = total_positives,
         n_sequenced = total_sequenced,
         pct_unsequenced = pct_unsequenced) |>
  mutate(pct_sequenced = round(100 - pct_unsequenced, 1)) |>
  select(region, n_positive_hosts, n_sequenced, pct_sequenced, pct_unsequenced) |>
  arrange(desc(pct_unsequenced))

write_csv(revision_table_3, here("output", "tables", "revision_table_3_regional_sequencing_gap.csv"))