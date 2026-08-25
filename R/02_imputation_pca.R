# ---
# 02_imputation_pca.R
#
# Purpose:
# 1. Load aligned trait data and phylogeny.
# 2. Impute missing life-history traits using phylogenetic covariance (Rphylopars).
# 3. Perform PCA to generate the "Pace of Life" axis (PC1).
#
# ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, ape, Rphylopars, phytools)

dir.create(here("data", "analytic"), showWarnings = FALSE)
dir.create(here("output", "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

# 1. Load and filter ------------------------------------------------------
trait_data  <- read_rds(here("data", "processed", "trait_data_phylo_matched.rds"))
mammal_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))

target_orders <- c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA")
trait_data <- trait_data |> filter(toupper(order) %in% target_orders)

pol_vars <- c("adult_mass_g", "litter_size", "litters_per_year",
              "gestation_d", "weaning_d", "age_first_reproduction_d",
              "max_longevity_d")

# Align tree and data to the same species set before imputation
mammal_tree <- ape::keep.tip(mammal_tree, intersect(mammal_tree$tip.label, trait_data$tip_label))
trait_data  <- trait_data |> filter(tip_label %in% mammal_tree$tip.label)

stopifnot(setequal(trait_data$tip_label, mammal_tree$tip.label))
stopifnot(!any(duplicated(trait_data$tip_label)))
message("Species retained after filtering and tree alignment: ", nrow(trait_data))

# 2. Phylogenetic imputation ----------------------------------------------
impute_input <- trait_data |>
  select(species = tip_label, all_of(pol_vars)) |>
  mutate(across(all_of(pol_vars), log10))

p_impute <- phylopars(trait_data = impute_input, tree = mammal_tree, phylo_correlated = TRUE)

imputed_log_data <- p_impute$anc_recon |>
  as.data.frame() |>
  rownames_to_column("tip_label") |>
  as_tibble() |>
  filter(tip_label %in% trait_data$tip_label)

stopifnot(setequal(imputed_log_data$tip_label, trait_data$tip_label))
stopifnot(nrow(imputed_log_data) == nrow(trait_data))

# Check: imputed values should not exceed plausible biological ranges
imputed_ranges <- imputed_log_data |>
  summarise(across(all_of(pol_vars), list(min = ~10^min(.x), max = ~10^max(.x))))
print(glimpse(imputed_ranges))

# 3. PCA on imputed traits ------------------------------------------------
pca_input <- imputed_log_data |> select(all_of(pol_vars))
pca_res   <- prcomp(pca_input, center = TRUE, scale. = TRUE)

write_rds(pca_res, here("output", "models", "pca_result.rds"))

pc_scores <- as.data.frame(pca_res$x[, 1:3]) |>
  rename(pace_of_life_pc1 = PC1,
         repro_strategy_pc2 = PC2,
         maternal_investment_pc3 = PC3) |>
  mutate(tip_label = imputed_log_data$tip_label)

# Orient PC1 so that HIGH = slow (large mass, long gestation, late maturity)
loadings <- pca_res$rotation
if (loadings["adult_mass_g", "PC1"] < 0) {
  message("Flipping PC1 axis to align with 'Slow' life history...")
  pc_scores$pace_of_life_pc1 <- pc_scores$pace_of_life_pc1 * -1
  loadings[, "PC1"] <- loadings[, "PC1"] * -1
}

message("PC1 variance explained: ",
        round(summary(pca_res)$importance[2, 1] * 100, 1), "%")
print(loadings)

# Orientation check: slow-lived taxa should score high, fast-lived low
orientation_check <- pc_scores |>
  filter(str_detect(tip_label, "Castor_|Hystrix_|Chinchilla_|Crocidura_|Mus_musculus|Rattus_norvegicus")) |>
  arrange(desc(pace_of_life_pc1))
orientation_check |>
  head(30) |>
  print()

# 4. Assemble analytic dataset --------------------------------------------
final_host_traits <- trait_data |>
  left_join(imputed_log_data, by = "tip_label", suffix = c("", "_log_imputed")) |>
  left_join(pc_scores, by = "tip_label") |>
  select(tip_label, gbif_id, pace_of_life_pc1, repro_strategy_pc2,
         maternal_investment_pc3, everything())

stopifnot(nrow(final_host_traits) == nrow(trait_data))
stopifnot(!any(is.na(final_host_traits$pace_of_life_pc1)))

write_rds(final_host_traits, here("data", "analytic", "host_traits_final.rds"))

# 5. Complete-case PCA (empirical traits only, no imputation) -------------
complete_traits <- trait_data |>
  drop_na(all_of(pol_vars)) |>
  mutate(across(all_of(pol_vars), log10))

message("Complete-case species: ", nrow(complete_traits))

pca_raw <- prcomp(select(complete_traits, all_of(pol_vars)), center = TRUE, scale. = TRUE)

pc_raw_scores <- as.data.frame(pca_raw$x[, 1:3]) |>
  rename(pace_of_life_pc1_raw = PC1,
         repro_strategy_pc2_raw = PC2,
         maternal_investment_pc3_raw = PC3) |>
  mutate(tip_label = complete_traits$tip_label)

if (pca_raw$rotation["adult_mass_g", "PC1"] < 0) {
  message("Flipping empirical PC1 axis...")
  pc_raw_scores$pace_of_life_pc1_raw <- pc_raw_scores$pace_of_life_pc1_raw * -1
}

host_traits_sensitivity <- complete_traits |>
  select(tip_label, gbif_id) |>
  left_join(pc_raw_scores, by = "tip_label")

stopifnot(nrow(host_traits_sensitivity) == nrow(complete_traits))
stopifnot(!any(is.na(host_traits_sensitivity$pace_of_life_pc1_raw)))

write_rds(host_traits_sensitivity, here("data", "analytic", "host_traits_sensitivity_raw.rds"))

# 6. Revision Table 1: covariate definitions ------------------------------
trait_labels <- tibble(
  variable = pol_vars,
  order_idx = seq_along(pol_vars),
  trait_name = c("Adult body mass", "Litter size", "Litters per year",
                 "Gestation length", "Weaning age",
                 "Age at first reproduction", "Maximum longevity"),
  unit = c("g", "n", "n", "days", "days", "days", "days"),
  source = "COMBINE (Soria-Auza et al. 2021)")

revision_table_1 <- trait_data |>
  select(all_of(pol_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(n_total     = n(),
            n_missing   = sum(is.na(value)),
            pct_imputed = round(100 * n_missing / n_total, 1),
            min         = min(value, na.rm = TRUE),
            median      = median(value, na.rm = TRUE),
            mean        = mean(value, na.rm = TRUE),
            sd          = sd(value, na.rm = TRUE),
            max         = max(value, na.rm = TRUE),
            .groups = "drop") |>
  left_join(trait_labels, by = "variable") |>
  arrange(order_idx) |>
  mutate(across(c(min, median, mean, sd, max), \(x) signif(x, 4))) |>
  select(trait_name, unit, source, n_total, pct_imputed, min, median, mean, sd, max)

write_csv(revision_table_1, here("output", "tables", "revision_table_1_covariate_definitions.csv"))
print(revision_table_1)

trait_dists <- trait_data |>
  select(all_of(pol_vars)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  drop_na(value) |>
  mutate(value = log10(value)) |>
  group_by(variable) |>
  summarise(dist = list(value), .groups = "drop")

write_rds(trait_dists, here("output", "tables", "revision_table_1_distributions.rds"))
