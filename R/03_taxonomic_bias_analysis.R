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
target_orders <- c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA")

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
  filter(toupper(order) %in% target_orders) |>
  drop_na(family) |>
  count(order, family, sampling_status) |>
  pivot_wider(names_from = sampling_status, values_from = n, values_fill = 0) |>
  mutate(total = Sampled + `Not Sampled`,
         prop_sampled = Sampled / total,
         family_label = paste0(family, " (", total, ")")) |>
  filter(total > 5) |> 
  arrange(order, desc(prop_sampled)) |>
  mutate(family_label = fct_reorder(paste0(family, " (", total, ")"), total))

p_family <- family_bias |>
  pivot_longer(cols = c(Sampled, `Not Sampled`), names_to = "status", values_to = "count") |>
  mutate(family_label = fct_reorder(family_label, total)) |>
  ggplot(aes(x = family_label, y = count, fill = status)) +
  coord_flip() +
  geom_col() +
  scale_y_log10() +
  facet_wrap(~order, scales = "free") +
  scale_fill_manual(values = c("Not Sampled" = "grey90", "Sampled" = "steelblue")) +
  labs(title = "Taxonomic Sampling Bias by Family", x = NULL, y = "Log10 Species Count") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 8),
        legend.position = "bottom")

# 4.2. By synanthropy Status
p_syn <- bias_data |>
  count(synanthropy_status, sampling_status) |>
  mutate(synanthropy_status = replace_na(synanthropy_status, "Unknown")) |>
  drop_na(synanthropy_status) |>
  ggplot(aes(x = synanthropy_status, y = n, fill = sampling_status)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Not Sampled" = "grey80", "Sampled" = "firebrick")) +
  labs(title = "Sampling Bias by Synanthropy", x = NULL, y = "Proportion") +
  theme_minimal()


# --- 5. Modelling Sampling Intensity (GAMs) ---
# Analytic dataset (imputed) is used
model_data <- host_traits_analytic |>
  left_join(sampling_summary, by = "gbif_id") |>
  mutate(n_pathogen_samples = replace_na(n_pathogen_samples, 0), # Fill Zeros: If species isn't in sampling_summary, it was sampled 0 times
         log_mass = log10(adult_mass_g),
         log_range = log10(geographic_range_km2),
         synanthropy_status = as.character(synanthropy_status),
         synanthropy_status = replace_na(synanthropy_status, "Unknown"),
         synanthropy_status = as.factor(synanthropy_status)) |>
  drop_na(log_mass, log_range)

# Model 1: Core Predictors (Mass + Range + Synanthropy)
gam_effort <- gam(n_pathogen_samples ~ 
                    s(log_mass, k = 5) + 
                    s(log_range, k = 5) + 
                    synanthropy_status,
                  data = model_data, 
                  family = nb())  # Negative Binomial for overdispersed counts

summary(gam_effort)
plot(gam_effort, pages = 1, scheme = 1)

write_rds(gam_effort, here("output", "analysis_1", "models", "gam_sampling_effort.rds"))


# --- 6. Quantifying Phylogenetic Bias (Pagel's Lambda) ---
# We test if sampling effort clusters by lineage
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
capture.output(lambda_res, file = here("output", "phylogenetic_bias_lambda.txt"))

# --- 7. Pathogen Side Bias (Counts) ---
pathogens_per_host <- pathogen_data |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  left_join(host_data |> select(host_record_id, host_species), by = "host_record_id") |>
  drop_na(host_species, pathogen_species_cleaned) |>
  separate_rows(pathogen_species_cleaned, sep = ",\\s*") |>
  group_by(host_species) |>
  summarise(n_viruses_tested = n_distinct(pathogen_species_cleaned)) |>
  arrange(desc(n_viruses_tested))

p_breadth <- ggplot(pathogens_per_host, aes(x = n_viruses_tested)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Surveillance Breadth", 
       x = "Unique Viruses Tested per Host Species", y = "Count of Host Species") +
  theme_minimal()

# --- 8. Host Pathogen ---
top_hosts <- pathogens_per_host |> top_n(62, n_viruses_tested) |> pull(host_species)

heatmap_data <- pathogen_data |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  left_join(host_data |> select(host_record_id, host_species), by = "host_record_id") |>
  filter(host_species %in% top_hosts) |>
  drop_na(pathogen_species_cleaned) |>
  separate_rows(pathogen_species_cleaned, sep = ",\\s*") |>
  group_by(host_species, pathogen_species_cleaned, pathogen_family) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE)) |>
  filter(n_tested > 0)

p_heatmap <- ggplot(heatmap_data, aes(x = pathogen_species_cleaned, y = host_species, fill = log1p(n_tested))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", name = "log(N Tested)") +
  facet_grid(~ pathogen_family, scales = "free_x", space = "free_x") +
  labs(title = "Co-surveillance Landscape (Top Hosts)", x = "Pathogen", y = "Host") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
