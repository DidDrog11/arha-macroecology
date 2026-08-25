# ==============================================================================
# 07_macroecological_analysis.R 
# Purpose: 
# 1. Join 'host_traits_final' with ArHa pathogen data to define "Reservoir Status".
# 2. Calculate Sampling Effort (Offset).
# 3. Fit Phylogenetic GLMMs (Bernoulli) using brms.
# 4. Test Hypotheses: Pace of Life, Synanthropy, Phylogeny.
# ==============================================================================

# 1. Setup ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, brms, ape, performance, tidybayes)

# Load Inputs
# The Analytic Trait Data (Predictors)
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))

# The Matched Phylogeny
host_tree <- read.tree(here("data", "processed", "mammal_tree_matched.tre"))

# The Raw ArHa Database
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-08-17.rds"))

# 2. Data Preparation -----------------------------------------------------
valid_pathogens <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  filter(!is.na(pathogen_species_cleaned)) |>
  separate_rows(pathogen_species_cleaned, sep = "\\s*[|,;]\\s*") |> 
  distinct(pathogen_species_cleaned, pathogen_family)

valid_hosts <- host_traits$tip_label

dyad_frame <- expand.grid(tip_label = valid_hosts,
                          pathogen_species_cleaned = valid_pathogens$pathogen_species_cleaned,
                          stringsAsFactors = FALSE) |>
  as_tibble() |>
  left_join(valid_pathogens, by = "pathogen_species_cleaned")

# Define Observed Positive Pairs
# Active Infection (PCR/Isolation)
obs_active <- arha_db$pathogen |>
  filter(!is.na(pathogen_species_cleaned)) |>
  separate_rows(pathogen_species_cleaned, sep = "\\s*[|,;]\\s*") |>
  filter(str_detect(assay, "PCR|Sequencing|Isolation|Culture|Antigen")) |>
  filter(number_positive > 0) |>
  left_join(arha_db$host |>
              select(host_record_id, host_species), by = "host_record_id") |>
  distinct(host_species, pathogen_species_cleaned) |>
  mutate(status_active = 1)

# Prior Exposure (Serology)
obs_exposed <- arha_db$pathogen |>
  filter(!is.na(pathogen_species_cleaned)) |>
  separate_rows(pathogen_species_cleaned, sep = "\\s*[|,;]\\s*") |>
  filter(str_detect(assay, "Serology|Antibody|ELISA|IFA")) |>
  filter(number_positive > 0) |>
  left_join(arha_db$host |>
              select(host_record_id, host_species), by = "host_record_id") |>
  distinct(host_species, pathogen_species_cleaned) |>
  mutate(status_exposed = 1)

model_data <- dyad_frame |>
  mutate(host_join_name = str_replace_all(tip_label, "_", " ")) |>
  left_join(obs_active, by = c("host_join_name" = "host_species", "pathogen_species_cleaned")) |>
  left_join(obs_exposed, by = c("host_join_name" = "host_species", "pathogen_species_cleaned")) |>
  mutate(status_active = replace_na(status_active, 0),
         status_exposed = replace_na(status_exposed, 0),
         is_reservoir = pmax(status_active, status_exposed),
         is_reservoir_strict = status_active)

# 3. Add Predictors -------------------------------------------------------
# Host Sampling Effort
host_effort <- arha_db$host |>
  group_by(host_species) |>
  summarise(n_individuals = sum(number_of_hosts, na.rm = TRUE)) |>
  mutate(log_effort = log(n_individuals + 1))

model_data_full <- model_data |>
  left_join(host_traits, by = "tip_label") |>
  left_join(host_effort, by = c("host_join_name" = "host_species")) |>
  mutate(log_effort = replace_na(log_effort, 0)) |>
  mutate(host_id = tip_label) |>
  filter(log_effort > 0) |>
  select(is_reservoir, log_effort, pace_of_life_pc1, tip_label, host_id, pathogen_species_cleaned) |>
  drop_na()

model_data_sens_full <- model_data |>
  left_join(host_traits, by = "tip_label") |>
  left_join(host_effort, by = c("host_join_name" = "host_species")) |>
  mutate(log_effort = replace_na(log_effort, 0)) |>
  mutate(host_id = tip_label) |>
  filter(log_effort > 0) |>
  select(is_reservoir_strict, log_effort, pace_of_life_pc1, tip_label, host_id, pathogen_species_cleaned) |>
  drop_na()

model_data_synanthropy <- model_data |>
  left_join(host_traits, by = "tip_label") |>
  left_join(host_effort, by = c("host_join_name" = "host_species")) |>
  mutate(log_effort = replace_na(log_effort, 0),
         synanthropy_status = factor(synanthropy_status, 
                                     levels = c("Not Synanthropic", "Occasionally Synanthropic", "Totally Synanthropic"))) |>
  filter(log_effort > 0) |> 
  select(is_reservoir, log_effort, synanthropy_status, pace_of_life_pc1, tip_label, pathogen_species_cleaned) |>
  drop_na() |>
  mutate(host_id = tip_label)

# 4. Prepare Phylogeny ----------------------------------------------------
A <- ape::vcv.phylo(host_tree)
# Filter matrix to only hosts in the model
A <- A[unique(model_data$tip_label), unique(model_data$tip_label)]

# 5. Fit Dyadic Model -----------------------------------------------------
# Formula: 
# Status ~ Effort + Traits + (1|Host_Phylo) + (1|Host_ID) + (1|Virus_ID)
# We add (1|Host_ID) to account for non-phylogenetic host variance (overdispersion)
# We add (1|Virus_ID) because some viruses (like Lassa) are easier to find or more studied.
priors_sparse <- c(prior(normal(0, 1), class = "b"),        # Tighter prior on slopes
                   prior(normal(-4, 2), class = "Intercept"), # Low baseline probability (~logit(0.002) is approx -6)
                   prior(normal(0, 1), class = "sd")          # Random effects
                   )

formula_global <- bf(is_reservoir ~ log_effort + pace_of_life_pc1 + 
                       (1 | gr(tip_label, cov = A)) +  # Phylogenetic Signal
                       (1 | host_id) +                 # Non-Phylogenetic Noise
                       (1 | pathogen_species_cleaned))

fit_dyadic_a <- brm(formula = formula_global,
                  data = model_data_full,
                  data2 = list(A = A),
                  family = bernoulli(link = "logit"),
                  prior = priors_sparse,
                  chains = 8, 
                  cores = 8, 
                  iter = 2500,   
                  warmup = 1500,
                  refresh = 40,
                  normalize = FALSE,
                  control = list(adapt_delta = 0.98, max_treedepth = 12),
                  file = here("output", "models", "brms_dyadic_N49k.rds"))

formula_sens <- bf(is_reservoir_strict ~ log_effort + pace_of_life_pc1 + 
                     (1 | gr(tip_label, cov = A)) + 
                     (1 | host_id) +                  
                     (1 | pathogen_species_cleaned))

fit_dyadic_a_sens <- brm(formula = formula_sens,
                         data = model_data_sens_full,
                         data2 = list(A = A),
                         family = bernoulli(link = "logit"),
                         prior = priors_sparse,
                         chains = 8, 
                         cores = 8, 
                         iter = 2500,   
                         warmup = 1500,
                         refresh = 40,
                         normalize = FALSE,
                         control = list(adapt_delta = 0.98, max_treedepth = 12),
                         file = here("output", "models", "brms_dyadic_sens_strict.rds"))

formula_2 <- bf(is_reservoir ~ log_effort + synanthropy_status + pace_of_life_pc1 + 
                (1 | gr(tip_label, cov = A)) +  # Phylogenetic Signal
                (1 | host_id) +                 # Non-Phylogenetic Noise
                (1 | pathogen_species_cleaned))

fit_dyadic_b <- brm(formula = formula_2,
                  data = model_data_synanthropy,
                  data2 = list(A = A),
                  family = bernoulli(link = "logit"),
                  prior = priors_sparse,
                  chains = 8, 
                  cores = 8, 
                  iter = 2500,   
                  warmup = 1500,
                  refresh = 20,
                  normalize = FALSE,
                  control = list(adapt_delta = 0.98, max_treedepth = 12),
                  file = here("output", "models", "brms_dyadic_N19k.rds"))

# 6. Complete-Case Sensitivity Analysis (Non-Imputed) ---------------------
# Load the raw sensitivity traits
host_traits_raw <- read_rds(here("data", "analytic", "host_traits_sensitivity_raw.rds"))

# Prepare the data: 476 complete cases
model_data_raw_sens <- model_data |>
  inner_join(host_traits_raw, by = "tip_label") |>
  left_join(host_effort, by = c("host_join_name" = "host_species")) |>
  mutate(log_effort = replace_na(log_effort, 0),
         host_id = tip_label) |>
  filter(log_effort > 0) |>
  select(is_reservoir, log_effort, pace_of_life_pc1_raw, tip_label, host_id, pathogen_species_cleaned) |>
  drop_na()

# Subset the covariance matrix to complete-case species
A_raw <- ape::vcv.phylo(host_tree)
valid_raw_hosts <- unique(model_data_raw_sens$tip_label)
A_raw <- A_raw[valid_raw_hosts, valid_raw_hosts]

# Define the formula using the raw PCA axis
formula_3 <- bf(is_reservoir ~ log_effort + pace_of_life_pc1_raw + 
                    (1 | gr(tip_label, cov = A_raw)) + 
                    (1 | host_id) + 
                    (1 | pathogen_species_cleaned))

# Fit the sensitivity model
fit_dyadic_c <- brm(formula = formula_3,
                    data = model_data_raw_sens,
                    data2 = list(A_raw = A_raw),
                    family = bernoulli(link = "logit"),
                    prior = priors_sparse,
                    chains = 8, 
                    cores = 8, 
                    iter = 2500,   
                    warmup = 1500,
                    refresh = 20,
                    normalize = FALSE,
                    control = list(adapt_delta = 0.98, max_treedepth = 12),
                    file = here("output", "models", "brms_dyadic_sens_raw.rds"))

# Extract fixed effects to check the stability of the β coefficient and pd
fixef(fit_dyadic_c, probs = c(0.025, 0.975))
bayesplot::mcmc_areas(fit_dyadic_c, pars = "b_pace_of_life_pc1_raw", prob = 0.95)

# Extract variance components to ensure phylogeny isn't artificially inflated
VarCorr(fit_dyadic_c)

# 7. Export model summaries -----------------------------------------------
export_full_summary <- function(mod, file) {
  s <- summary(mod)
  fix <- as.data.frame(s$fixed) |> rownames_to_column("Parameter") |> mutate(Type = "Fixed Effects")
  rnd <- bind_rows(lapply(names(s$random), \(g) as.data.frame(s$random[[g]]) |> mutate(Parameter = paste0("sd_", g), Type = "Variance Components")))
  
  bind_rows(fix, rnd) |>
    select(Type, Parameter, Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS) |>
    write_csv(here("output", "tables", file))
}

# Export all four models
export_full_summary(fit_dyadic_a, "table_s2_full.csv")
export_full_summary(fit_dyadic_b, "table_s3_synanthropy.csv")
export_full_summary(fit_dyadic_a_sens, "table_s4_strict.csv")
export_full_summary(fit_dyadic_c, "table_s5_complete.csv")

# 8. Synanthropy Box-Plot -------------------------------------------------

p_syn_pace <- ggplot(model_data_synanthropy, aes(x = synanthropy_status, y = pace_of_life_pc1, fill = synanthropy_status)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  labs(title = "Pace of life by synanthropy status",
       subtitle = "ANOVA: F(2, 18903) = 77.4, p < 0.001",
       x = NULL, y = "Pace of Life (PC1)") +
  theme_minimal() + theme(legend.position = "none")

ggsave(here("output", "figures", "revision_fig_synanthropy_pace.png"), p_syn_pace, width = 6, height = 5, bg = "white")


# --- Revision Tables: Descriptive summary of modelled taxa ---------------
# R1 methodological comment 3: taxonomic and geographic context for the
# species and viruses included in the analytic datasets.

# Apply the same separate_rows split to the lookup before joining, so
# ambiguous multi-name labels (e.g. "Virus A | Virus B") match the
# split names used in model_data_full.
pathogen_lookup_split <- arha_db$pathogen |>
  filter(!is.na(pathogen_species_cleaned)) |>
  separate_rows(pathogen_species_cleaned, sep = "\\s*[|,;]\\s*") |>
  distinct(pathogen_species_cleaned, pathogen_family, pathogen_genus_ncbi)

detections_split <- arha_db$pathogen |>
  filter(!is.na(pathogen_species_cleaned), number_positive > 0) |>
  separate_rows(pathogen_species_cleaned, sep = "\\s*[|,;]\\s*") |>
  left_join(arha_db$host |> select(host_record_id, iso3c), by = "host_record_id") |>
  group_by(pathogen_species_cleaned) |>
  summarise(n_positive = sum(number_positive, na.rm = TRUE),
            n_countries = n_distinct(iso3c, na.rm = TRUE),
            .groups = "drop")

virus_summary <- model_data_full |>
  distinct(pathogen_species_cleaned) |>
  left_join(pathogen_lookup_split, by = "pathogen_species_cleaned") |>
  left_join(detections_split, by = "pathogen_species_cleaned") |>
  left_join(
    model_data_full |>
      filter(is_reservoir == 1) |>
      group_by(pathogen_species_cleaned) |>
      summarise(n_host_species = n_distinct(tip_label), .groups = "drop"),
    by = "pathogen_species_cleaned") |>
  mutate(across(c(n_positive, n_countries, n_host_species), ~replace_na(.x, 0))) |>
  arrange(pathogen_family, desc(n_host_species)) |>
  select(pathogen_species_cleaned, pathogen_family, pathogen_genus_ncbi,
         n_host_species, n_positive, n_countries)

write_csv(virus_summary, here("output", "tables", "revision_table_virus_summary.csv"))

# --- Host species summary ---
# 704 species is too many to print individually; summarise by family.
host_summary <- model_data_full |>
  distinct(tip_label, log_effort) |>
  mutate(host_species = str_replace_all(tip_label, "_", " ")) |>
  left_join(host_traits |> select(tip_label, family, pace_of_life_pc1), by = "tip_label") |>
  left_join(
    model_data_full |>
      filter(is_reservoir == 1) |>
      group_by(tip_label) |>
      summarise(n_virus_assoc = n_distinct(pathogen_species_cleaned), .groups = "drop"),
    by = "tip_label") |>
  mutate(n_virus_assoc = replace_na(n_virus_assoc, 0),
         n_individuals = round(exp(log_effort) - 1)) |>
  group_by(family) |>
  summarise(n_species = n(),
            n_species_with_assoc = sum(n_virus_assoc > 0),
            median_individuals_tested = median(n_individuals, na.rm = TRUE),
            median_pace_of_life = round(median(pace_of_life_pc1, na.rm = TRUE), 2),
            .groups = "drop") |>
  arrange(desc(n_species))

write_csv(host_summary, here("output", "tables", "revision_table_host_summary.csv"))

print(virus_summary)
print(host_summary, n = 40)
