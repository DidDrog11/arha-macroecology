# ---
# 01_data_integration.R
#
# Purpose: 
# 1. Load the finalized ArHa database.
# 2. Load external macroecological trait datasets (COMBINE, EltonTraits, PanTHERIA, etc.).
# 3. Harmonize taxonomy between ArHa (MDD-based) and trait databases (often MSW3-based).
# 4. Merge traits and perform PCA to generate "Pace of Life" and "Metabolic" predictors.
# 5. Export a 'model_ready_hosts.rds' for the GLMMs.
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, janitor, taxize, missForest, broom)

# Create directories if they don't exist
dir.create(here("data", "external"), showWarnings = FALSE)
dir.create(here("data", "processed"), showWarnings = FALSE)

# --- 2. Load ArHa Database ---
# Assuming you have downloaded the RDS to data/database/ per your previous structure
arha_path <- here("data", "database", "Project_ArHa_database_2025-12-02.rds")

if (!file.exists(arha_path)) {
  stop("ArHa database not found. Please download it from the GitHub repo.")
}

arha_db <- read_rds(arha_path)
host_list <- arha_db$host |> 
  distinct(host_species, gbif_id) |> 
  filter(!is.na(host_species))

# --- 3. Load External Trait Datasets (Placeholders) ---

# A. COMBINE (Soria-Auza et al.) - Good for mass, litter size, lifespan
# URL: https://figshare.com/articles/dataset/COMBINE_archive/14589288
combine_data <- read_csv(here("data", "external", "combine_data.csv")) |> 
  clean_names() |> 
  select(scientific_name = i_ucn_2020_binomial, 
         adult_mass_g, 
         litter_size, 
         litters_per_year, 
         maximum_longevity_y,
         gestation_d,
         weaning_d)

# B. EltonTraits 1.0 - Good for Diet and Activity
# URL: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887
elton_data <- read_tsv(here("data", "external", "MamFuncDat.txt")) |> 
  clean_names() |> 
  select(scientific_name, 
         diet_inv, diet_vend, diet_vect, diet_vfish, diet_vunk, 
         diet_scav, diet_fruit, diet_nect, diet_seed, diet_planto, 
         activity_nocturnal, activity_crepuscular, activity_diurnal,
         forstrat_ground) |> 
  # Calculate Diet Breadth (Shannon Index or simple count of categories > 0)
  rowwise() |> 
  mutate(diet_breadth = sum(c_across(starts_with("diet_")) > 0)) |> 
  ungroup()

# C. AnAge / PanTHERIA (Optional backup for life history)
# pantheria_data <- read_tsv(...)

# --- 4. Taxonomic Harmonization (The Hard Part) ---
# ArHa uses MDD (modern). COMBINE/Elton use MSW3 (older). 
# Names will mismatch (e.g., synonyms, split species).

# Helper function to try and match names if direct join fails
# In a real run, you might use the 'taxize' package or a synonym dictionary
harmonize_names <- function(df_traits, df_hosts) {
  
  # 1. Try Direct Join
  matched <- df_hosts |> inner_join(df_traits, by = c("host_species" = "scientific_name"))
  
  # 2. Find Unmatched
  unmatched_hosts <- df_hosts |> anti_join(df_traits, by = c("host_species" = "scientific_name"))
  
  message(paste0("Direct matches: ", nrow(matched)))
  message(paste0("Unmatched hosts: ", nrow(unmatched_hosts)))
  
  # TODO: Insert manual dictionary or taxize::synonyms() loop here for the unmatched 
  # For now, we will proceed with direct matches for the template
  
  return(matched)
}

# Merge COMBINE
host_combine <- harmonize_names(combine_data, host_list)

# Merge Elton
host_traits_full <- host_combine |> 
  left_join(elton_data, by = c("host_species" = "scientific_name"))

# --- 5. Imputation (Optional but Recommended) ---
# Macroecological models hate missing data. missForest is robust for this.
# Only impute if missingness is < ~30-40%
cols_to_impute <- c("adult_mass_g", "litter_size", "litters_per_year", "maximum_longevity_y", "gestation_d")

# Check missingness
colMeans(is.na(host_traits_full[, cols_to_impute]))

# Simple RF Imputation (Uncomment to run)
# set.seed(123)
# imputed_data <- missForest(host_traits_full |> select(all_of(cols_to_impute)) |> as.matrix())$ximp
# host_traits_full <- bind_cols(host_traits_full |> select(-all_of(cols_to_impute)), as_tibble(imputed_data))


# --- 6. Feature Engineering: "Pace of Life" PCA ---
# We want a single axis representing "Fast" vs "Slow" life history to avoid collinearity
# Fast = Small mass, large litter, short gestation, short life
# Slow = Large mass, small litter, long gestation, long life

pca_data <- host_traits_full |> 
  select(adult_mass_g, litter_size, litters_per_year, maximum_longevity_y, gestation_d) |> 
  drop_na() |> # PCA cannot handle NAs
  mutate(across(everything(), log1p)) # Log transform traits (standard in metabolic ecology)

pca_fit <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# Inspect loadings to interpret PC1 (usually the Fast-Slow axis)
# print(pca_fit$rotation)

# Add PC1 back to the data
# Note: Check sign of PC1. If 'Fast' traits are negative, multiply by -1 to make interpretation intuitive (Higher PC1 = Faster)
host_traits_final <- host_traits_full |> 
  drop_na(adult_mass_g, litter_size, litters_per_year, maximum_longevity_y, gestation_d) |> 
  bind_cols(life_history_pc1 = pca_fit$x[, 1])

# --- 7. Save Model-Ready Data ---
write_rds(host_traits_final, here("data", "processed", "host_traits_macroecology.rds"))

message("Data integration complete. Traits ready for modeling.")