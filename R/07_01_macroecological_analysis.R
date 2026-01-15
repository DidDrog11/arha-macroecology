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
         is_reservoir = pmax(status_active, status_exposed))

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
