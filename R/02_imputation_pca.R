# ---
# 02_imputation_pca.R
#
# Purpose:
# 1. Load aligned trait data and phylogeny.
# 2. Impute missing life-history traits using phylogenetic covariance (Rphylopars).
# 3. Perform PCA to generate the "Pace of Life" axis (PC1).
# 4. Quantify phylogenetic bias in sampling effort.
# ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, ape, Rphylopars, phytools)

# Create output directory
dir.create(here("data", "analytic"), showWarnings = FALSE)

# Load the aligned outputs from 01_data_integration.R
trait_data <- read_rds(here("data", "processed", "trait_data_phylo_matched.rds"))
mammal_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))

# Define Pace of Life variables
pol_vars <- c("adult_mass_g", "litter_size", "litters_per_year", 
              "gestation_d", "weaning_d", "age_first_reproduction_d", 
              "max_longevity_d")

# Prepare input: Log10 transform to satisfy Brownian Motion assumptions
impute_input <- trait_data |> 
  select(species = tip_label, all_of(pol_vars)) |> 
  mutate(across(all_of(pol_vars), log10))

# Run Rphylopars (Brownian Motion model including covariance)
# phylo_correlated = TRUE allows the model to use correlations between traits
p_impute <- phylopars(trait_data = impute_input, tree = mammal_tree, phylo_correlated = TRUE)

# Extract Imputed Data
imputed_log_data <- p_impute$anc_recon[1:nrow(impute_input), ] |> 
  as.data.frame() |> 
  rownames_to_column("tip_label") |> 
  as_tibble()

# PCA on the imputed data
pca_input <- imputed_log_data |> 
  select(-tip_label)

pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)

# Save PCA
write_rds(pca_res, here("output", "models", "pca_result.rds"))

# Extract first 3 PCs
pc_scores <- as.data.frame(pca_res$x[, 1:3]) |>
  rename(
    pace_of_life_pc1 = PC1,            # Fast (Neg) vs Slow (Pos)
    repro_strategy_pc2 = PC2,          # Large Litters (Neg) vs Frequent Litters (Pos)
    maternal_investment_pc3 = PC3      # Somatic Maintenance (Neg) vs Parental Care (Pos)
  )

loadings <- pca_res$rotation

if (loadings["adult_mass_g", "PC1"] < 0) {
  message("Flipping PC1 axis to align with 'Slow' life history...")
  pc_scores$pace_of_life_pc1 <- pc_scores$pace_of_life_pc1 * -1
  # Flip loadings just for consistency in printing
  loadings[, "PC1"] <- loadings[, "PC1"] * -1
}

print(loadings[, 1:3])

# Save Final Analytic Dataset
final_host_traits <- trait_data |> 
  left_join(imputed_log_data, by = "tip_label", suffix = c("", "_log_imputed")) |> 
  bind_cols(pc_scores) |> 
  select(tip_label, gbif_id, 
         pace_of_life_pc1, repro_strategy_pc2, maternal_investment_pc3, 
         everything())

write_rds(final_host_traits, here("data", "analytic", "host_traits_final.rds"))
