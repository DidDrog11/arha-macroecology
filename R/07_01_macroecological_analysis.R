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
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-01-09.rds"))

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

formula <- bf(is_reservoir ~ log_effort + pace_of_life_pc1 + 
                (1 | gr(tip_label, cov = A)) +  # Phylogenetic Signal
                (1 | host_id) +                 # Non-Phylogenetic Noise
                (1 | pathogen_species_cleaned))

fit_dyadic_a <- brm(formula = formula,
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
fixef(fit_dyadic_3, probs = c(0.025, 0.975))
bayesplot::mcmc_areas(fit_dyadic_3, pars = "b_pace_of_life_pc1_raw", prob = 0.95)

# Extract variance components to ensure phylogeny isn't artificially inflated
VarCorr(fit_dyadic_3)

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
export_full_summary(fit_dyadic_3, "table_s5_complete.csv")