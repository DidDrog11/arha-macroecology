# ---
# 03_taxonomic_bias_analysis.R
#
# Purpose: 
# 1. Compare the ArHa sampled hosts against the global checklist (trait_data_all).
# 2. Visualise taxonomic and geographic biases.
# 3. Model predictors of sampling effort (GAMs).
# 4. Quantify Phylogenetic Bias (Pagel's Lambda).
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, mgcv, patchwork, phytools, ape)

# Create output directories
dir.create(here("output", "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# --- 2. Load Processed Data ---
# A. The Global Checklist (All relevant Rodents/Eulipotyphla)
global_traits <- read_rds(here("data", "processed", "trait_data.rds"))

# B. The Analytic Subset (Phylo-matched & Imputed)
host_traits_analytic <- read_rds(here("data", "analytic", "host_traits_final.rds"))

# C. The ArHa Database (For raw sampling counts)
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2025-12-02.rds"))
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
mammal_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))

# --- 3. Sampling Effort Summary ---
# 3.1. Calculate Effort per Species from ArHa
sampling_summary <- host_data |>
  drop_na(gbif_id, host_species) |>
  left_join(pathogen_data |> 
      group_by(host_record_id) |> 
      summarise(pathogen_samples = sum(number_tested, na.rm = TRUE),
                positive_samples = sum(number_positive, na.rm = TRUE)),
    by = "host_record_id") |>
  group_by(gbif_id) |>
  summarise(n_studies = n_distinct(study_id),
            n_pathogen_samples = sum(pathogen_samples, na.rm = TRUE),
            n_positive_samples = sum(positive_samples > 0, na.rm = TRUE))

# 3.2. Join Effort to Global Checklist
bias_data <- global_traits |>
  left_join(sampling_summary, by = "gbif_id") |>
  mutate(sampling_status = if_else(!is.na(n_studies), "Sampled", "Not Sampled"),
         log_sampling_effort = log1p(replace_na(n_pathogen_samples, 0)))

# --- 4. Visualising Taxonomic Bias ---
# 4.1. By Family
family_bias <- bias_data |>
  drop_na(family) |>
  count(order, family, sampling_status) |>
  group_by(family) |>
  mutate(total = sum(n, na.rm = TRUE)) |>
  rowwise() |>
  mutate(prop_sampled = n / total) |>
  ungroup() |>
  filter(total > 2) |> 
  arrange(order, desc(total)) |>
  mutate(family_label = fct_reorder(paste0(family, " (", total, ")"), total))

p_family <- ggplot(family_bias, aes(x = family_label, y = n, fill = sampling_status)) +
  geom_col(position = "fill") +
  coord_flip() +
  facet_grid(order ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Not Sampled" = "grey80", "Sampled" = "steelblue")) +
  labs(title = "Sampling Coverage by Family", x = NULL, y = "Proportion of Species") +
  theme_minimal()

ggsave(here("output", "figures", "taxonomic_bias_family.png"), p_family, width = 10, height = 12)

# 4.2. By IUCN Status
p_iucn <- bias_data |>
  # Use IUCN column from trait_data_all (ensure column name matches, likely 'geographic_range_km2' etc sourced from IUCN)
  # Note: You might need to check if 'iucn_status' exists in your trait_data_all. 
  # If not, it might have been dropped in data_integration. 
  # Assuming it's there based on your previous script logic:
  count(synanthropy_status, sampling_status) |>
  drop_na(synanthropy_status) |>
  ggplot(aes(x = synanthropy_status, y = n, fill = sampling_status)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Not Sampled" = "grey80", "Sampled" = "firebrick")) +
  labs(title = "Sampling Bias by Synanthropy", x = NULL, y = "Proportion") +
  theme_minimal()

ggsave(here("output", "analysis_1", "bias_synanthropy.png"), p_iucn, width = 8, height = 6)


# --- 5. Modelling Sampling Intensity (GAMs) ---

message("--- Fitting GAMs for Sampling Effort ---")

# Prepare Data for Modeling
# We use the analytic dataset (imputed) because we need complete predictors
model_data <- host_traits_analytic |>
  inner_join(sampling_summary, by = "gbif_id") |>
  mutate(log_mass = log10(adult_mass_g),
         log_range = log10(geographic_range_km2)) |>
  drop_na(log_mass, log_range)

# Model 1: Core Predictors (Mass + Range + Synanthropy)
gam_effort <- gam(n_pathogen_samples ~ 
                    s(log_mass, k = 5) + 
                    s(log_range, k = 5) + 
                    synanthropy_status,
                  data = model_data, 
                  family = nb()) # Negative Binomial for overdispersed counts

summary(gam_effort)
plot(gam_effort, pages = 1, scheme = 1)

write_rds(gam_effort, here("output", "analysis_1", "models", "gam_sampling_effort.rds"))


# --- 6. Quantifying Phylogenetic Bias (Pagel's Lambda) ---

message("--- Calculating Phylogenetic Signal in Sampling ---")

# We test if sampling effort clusters by lineage
# Ensure we match the tree tips
lambda_data <- model_data |>
  filter(tip_label %in% mammal_tree$tip.label) |>
  column_to_rownames("tip_label")

# Extract vector
effort_vector <- log1p(lambda_data$n_pathogen_samples)
names(effort_vector) <- rownames(lambda_data)

# Run Test
lambda_res <- phylosig(mammal_tree, effort_vector, method = "lambda", test = TRUE)

print(lambda_res)

# Save result text
capture.output(lambda_res, file = here("output", "analysis_1", "phylogenetic_bias_lambda.txt"))


# --- 7. Pathogen Side Bias (Counts) ---
# This remains largely similar to your original code

virus_counts <- pathogen_data |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  count(pathogen_family, pathogen_species_cleaned, name = "n_records") |>
  arrange(desc(n_records)) |>
  head(20)

p_virus <- ggplot(virus_counts, aes(x = reorder(pathogen_species_cleaned, n_records), y = n_records)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  facet_wrap(~pathogen_family, scales = "free_y") +
  labs(title = "Top 20 Sampled Pathogens", x = NULL, y = "Number of Records") +
  theme_bw()

ggsave(here("output", "analysis_1", "virus_sampling_counts.png"), p_virus, width = 10, height = 8)

message("✅ Taxonomic Bias Analysis Complete.")