# ---
# 01_data_integration.R
#
# Purpose: 
# 1. Load the finalised ArHa database.
# 2. Load external macroecological trait datasets
# 3. Harmonise taxonomy between ArHa and trait databases.
# ---

# 1. Setup ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, janitor, terra, taxize, missForest, broom, ape, geiger, phangorn)

# Create directories if they don't exist
dir.create(here("data", "external"), showWarnings = FALSE)
dir.create(here("data", "processed"), showWarnings = FALSE)
dir.create(here("data", "database"), showWarnings = FALSE)

add_gbif_ids <- function(df, name_col) {
  
  # Extract unique names to minimise API calls
  unique_names <- df[[name_col]] |> unique() |> na.omit()
  
  message(paste0("Resolving ", length(unique_names), " unique names against GBIF backbone..."))
  
  # Use taxize to get GBIF IDs (batch processing is faster/safer)
  # db = "gbif" checks the GBIF backbone
  # rows = 1 forces the top hit (caution: assumes top hit is correct)
  ids <- get_gbifid_(unique_names, rows = 1)
  
  # Bind the results into a lookup table
  lookup_table <- bind_rows(ids, .id = "original_name")
  
  # Processing matches
  lookup_clean <- lookup_table |> 
    select(original_name, usagekey, scientificname, matchtype, status) |> 
    rename(gbif_id = usagekey, gbif_scientific_name = scientificname) |> 
    mutate(gbif_name_simple = word(gbif_scientific_name, 1, 2),
           original_name_simple = word(original_name, 1, 2),
           # 2. Flag if names don't match OR if matchtype is fuzzy
           name_mismatch = tolower(gbif_name_simple) != tolower(original_name_simple),
           review_needed = case_when(matchtype != "EXACT" ~ TRUE,      # GBIF says it's not an exact match
                                     name_mismatch ~ TRUE,             # The names look different (Synonym?)
                                     is.na(gbif_id) ~ TRUE,            # No ID found
                                     TRUE ~ FALSE))
  
  # Join back to original dataframe
  df_out <- df |> 
    left_join(lookup_clean, by = c(setNames("original_name", name_col)))
  
  # Report success rate
  success_count <- sum(!is.na(df_out$gbif_id))
  message(paste0("Resolved ", success_count, " out of ", nrow(df), " rows."))
  
  return(df_out)
}

target_orders <- c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA")


# 2. Load ArHa ------------------------------------------------------------
arha_path <- here("data", "database", "Project_ArHa_database_2026-01-09.rds")

if (!file.exists(arha_path)) {
  message("ArHa database not found. Downloading it from the GitHub repo.")
  arha_rds_url <- "https://raw.githubusercontent.com/DidDrog11/arenavirus_hantavirus/main/data/database/Project_ArHa_database_2026-01-09.rds"
  download.file(arha_rds_url, arha_path, mode = "wb") 
  arha_db <- read_rds(arha_path)
} else {
  arha_db <- read_rds(here(arha_path))
}

host_list <- arha_db$host |> 
  distinct(host_species, gbif_id) |> 
  filter(!is.na(host_species))


# 3. Load External Databases ----------------------------------------------
# A. COMBINE (Soria-Auza et al.) - Good for mass, litter size, lifespan
# URL: https://figshare.com/articles/dataset/COMBINE_a_Coalesced_Mammal_Database_of_Intrinsic_and_extrinsic_traits/13028255/4?file=27703263
# Non-imputed source used
combine_data <- read_csv(here("data", "external", "combine_data_reported.csv")) |> 
  clean_names() |>
  filter(toupper(order) %in% target_orders)

if (!file.exists(here("data", "processed", "combine_harmonised.rds"))) {
  combine_with_ids <- add_gbif_ids(combine_data, "iucn2020_binomial")
  
  combine_unmatched <- combine_with_ids |> 
    filter(is.na(gbif_id) | is.na(gbif_name_simple) | review_needed == TRUE)
  
  combine_harmonised <- combine_with_ids |>
    filter(!is.na(gbif_id) & !is.na(gbif_name_simple)) |> 
    select(
      # ID
      scientific_name = iucn2020_binomial, 
      
      # Allometry
      adult_mass_g, 
      
      # Pace of Life (Reproduction & Longevity)
      litter_size = litter_size_n, 
      litters_per_year = litters_per_year_n, 
      gestation_d = gestation_length_d,
      max_longevity_d = max_longevity_d,
      weaning_d = weaning_age_d,
      age_first_reproduction_d,
      
      # Ecological Opportunity
      diet_breadth = det_diet_breadth_n,
      trophic_level,
      habitat_breadth = habitat_breadth_n,
      
      # Transmission Potential
      pop_density = density_n_km2,
      social_group_size = social_group_n,
      
      # Behavioural Controls
      fossoriality,
      activity_cycle,
      
      # Taxonomy
      gbif_name_simple,
      gbif_id,
      name_mismatch,
      review_needed
    ) |>
    # Filter out duplicates if multiple Combine rows map to same GBIF ID (e.g. synonyms)
    # We take the mean of traits for combined synonyms
    group_by(gbif_id) |> 
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              across(where(is.character), \(x) first(na.omit(x))))
  
  write_rds(combine_harmonised, here("data", "processed", "combine_harmonised.rds"))
} else {
  combine_harmonised <- read_rds(here("data", "processed", "combine_harmonised.rds"))
}

# B. Ecke - Synanthropy
ecke_path <- here("data", "external", "rodent_synanthropy.xlsx")

if (!file.exists(here("data", "processed", "ecke_synanthropy_clean.rds"))) {
  
  ecke_raw <- read_xlsx(ecke_path) |> 
    clean_names() |> 
    select(
      scientific_name = rodents, 
      synanthropic,
      s_index,           # Continuous synanthropy score
      hunted,            # Human interface
      cyclicity,         # Population dynamics
      seasonality        # Population dynamics
    )
  
  # Define Hierarchy for Summary
  ecke_scored <- ecke_raw |> 
    mutate(
      # Synanthropy Rank (Already defined)
      syn_rank = case_when(str_detect(tolower(synanthropic), "t") ~ 3,
                           str_detect(tolower(synanthropic), "o") ~ 2,
                           str_detect(tolower(synanthropic), "n") ~ 1,
                           TRUE ~ 0),
      cyclicity_binary = if_else(str_detect(tolower(cyclicity), "y"), 1, 0),
      seasonality_binary = if_else(str_detect(tolower(seasonality), "y"), 1, 0),
      hunted_binary = if_else(str_detect(tolower(hunted), "y"), 1, 0),
      # Clean S-index (Ensure numeric)
      s_index = as.numeric(s_index)
    )
  
  # Aggregate by Species (Take the highest observed status)
  ecke_aggregated <- ecke_scored |> 
    group_by(scientific_name) |> 
    summarise(
      max_syn_rank = if(all(is.na(syn_rank))) NA_real_ else max(syn_rank, na.rm = TRUE),
      is_hunted = if(all(is.na(hunted_binary))) NA_real_ else max(hunted_binary, na.rm = TRUE),
      is_cyclical = if(all(is.na(cyclicity_binary))) NA_real_ else max(cyclicity_binary, na.rm = TRUE),
      is_seasonal = if(all(is.na(seasonality_binary))) NA_real_ else max(seasonality_binary, na.rm = TRUE),
      mean_s_index = mean(s_index, na.rm = TRUE)
    ) |> 
    mutate(
      synanthropy_status = case_when(
        max_syn_rank == 3 ~ "Totally Synanthropic",
        max_syn_rank == 2 ~ "Occasionally Synanthropic",
        max_syn_rank == 1 ~ "Not Synanthropic",
        TRUE ~ NA_character_ 
      )
    ) |> 
    filter(!is.na(synanthropy_status)) |>
    ungroup() |>
    mutate(scientific_name = str_replace_all(scientific_name, "_", " "))
  
  # Harmonise Taxonomy (GBIF)
  ecke_with_ids <- add_gbif_ids(ecke_aggregated, "scientific_name")
  
  ecke_unmatched <- ecke_with_ids |> 
    filter(is.na(gbif_id) | is.na(gbif_name_simple) | review_needed == TRUE)
  
  ecke_harmonised <- ecke_with_ids |>
    filter(!is.na(gbif_id) & !is.na(gbif_name_simple))
  
  # Save
  write_rds(ecke_harmonised, here("data", "processed", "ecke_synanthropy_clean.rds"))
  
} else {
  ecke_harmonised <- read_rds(here("data", "processed", "ecke_synanthropy_clean.rds"))
}

# IUCN ranges
# Downloaded from IUCN redlist for all terrestrial mammals
iucn_path <- here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp")

if(!file.exists(here("data", "processed", "iucn_harmonised.rds"))) {
  
  iucn_raw <- vect(iucn_path)
  
  iucn_aggregated <- iucn_raw[toupper(iucn_raw$order_) %in% target_orders, ] |>
    aggregate("sci_name")
  
  range_areas <- expanse(iucn_aggregated, unit = "km")
  
  iucn_data <- tibble(scientific_name = iucn_aggregated$sci_name,
                      geographic_range_km2 = range_areas)
  
  iucn_with_ids <- add_gbif_ids(iucn_data, "scientific_name")
  
  iucn_unmatched <- iucn_with_ids |> 
    filter(is.na(gbif_id) | is.na(gbif_name_simple) | review_needed == TRUE)
  
  iucn_harmonised <- iucn_with_ids |> filter(!is.na(gbif_id))
  
  write_rds(iucn_harmonised, here("data", "processed", "iucn_harmonised.rds"))
  
} else {
  
  iucn_harmonised <- read_rds(here("data", "processed", "iucn_harmonised.rds"))
  
}

# Matching unmatched or manual review taxa
perform_manual_match = FALSE
if(perform_manual_match) {
  
  prep_unmatched <- function(df, source_name) {
    df |> 
      select(scientific_name, gbif_id, matchtype, gbif_name_simple) |> 
      mutate(source = source_name)
  }
  
  all_unmatched <- bind_rows(
    prep_unmatched(iucn_unmatched, "IUCN"),
    prep_unmatched(ecke_unmatched, "Ecke"),
    prep_unmatched(combine_unmatched |> rename(scientific_name = iucn2020_binomial), "COMBINE")
  ) |> 
    distinct(scientific_name, .keep_all = TRUE)
  
  manual_ids <- get_gbifid(all_unmatched$scientific_name)
  
  # Bind results
  manual_lookup <- tibble(scientific_name = all_unmatched$scientific_name,
                          gbif_id = as.character(manual_ids)) |> 
    filter(!is.na(gbif_id)) |> 
    mutate(gbif_id = as.integer(gbif_id),
           gbif_scientific_name = scientific_name)
  
  combine_manual <- combine_unmatched |>
    select(-any_of(c("gbif_id", "gbif_scientific_name", "matchtype", "status", "gbif_name_simple", "original_name_simple", "name_mismatch", "review_needed"))) |>
    left_join(manual_lookup, by = c("iucn2020_binomial" = "scientific_name")) |> 
    select(
      scientific_name = iucn2020_binomial, 
      adult_mass_g, 
      litter_size = litter_size_n, 
      litters_per_year = litters_per_year_n, 
      gestation_d = gestation_length_d,
      max_longevity_d = max_longevity_d,
      weaning_d = weaning_age_d,
      age_first_reproduction_d,
      diet_breadth = det_diet_breadth_n,
      trophic_level,
      habitat_breadth = habitat_breadth_n,
      pop_density = density_n_km2,
      social_group_size = social_group_n,
      fossoriality,
      activity_cycle,
      gbif_name_simple = gbif_scientific_name,
      gbif_id) |>
    group_by(gbif_id) |> 
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              across(where(is.character), \(x) first(na.omit(x)))) |>
    drop_na(gbif_id)
  ecke_manual <- ecke_unmatched |>
    select(-any_of(c("gbif_id", "gbif_scientific_name", "matchtype", "status", "gbif_name_simple", "original_name_simple", "name_mismatch", "review_needed"))) |>
    left_join(manual_lookup, by = c("scientific_name")) |>
    rename(gbif_name_simple = gbif_scientific_name)
  iucn_manual <- iucn_unmatched |>
    select(-any_of(c("gbif_id", "gbif_scientific_name", "matchtype", "status", "gbif_name_simple", "original_name_simple", "name_mismatch", "review_needed"))) |>
    left_join(manual_lookup, by = c("scientific_name"))
  
  combine_final <- bind_rows(combine_harmonised, combine_manual) |> 
    distinct(gbif_id, .keep_all = TRUE) |> 
    select(
      gbif_id,
      adult_mass_g, litter_size, litters_per_year, 
      gestation_d, max_longevity_d, weaning_d, age_first_reproduction_d,
      diet_breadth, trophic_level, habitat_breadth, 
      pop_density, social_group_size, fossoriality, activity_cycle
    )
  ecke_final <- bind_rows(ecke_harmonised, ecke_manual) |> 
    distinct(gbif_id, .keep_all = TRUE) |> 
    select(
      gbif_id,
      synanthropy_status, 
      mean_s_index, 
      is_hunted, 
      is_cyclical, 
      is_seasonal
    )
  iucn_final <- bind_rows(iucn_harmonised, iucn_manual) |> 
    distinct(gbif_id, .keep_all = TRUE) |> 
    select(
      gbif_id, 
      geographic_range_km2
    )
  
  species_ids <- unique(c(arha_db$host |>
                            filter(!is.na(host_species)) |>
                            distinct(gbif_id) |>
                            pull(gbif_id),
                          unique(combine_final$gbif_id),
                          unique(ecke_final$gbif_id),
                          unique(iucn_final$gbif_id))) |> 
  na.omit() |> 
  as.numeric()
  
  tax_list <- classification(species_ids, db = "gbif")
  
  tax_hierarchy <- cbind(tax_list) |> 
    as_tibble() |> 
    select(gbif_id = query, kingdom, phylum, class, order, family, genus, species) |> 
    clean_names() |> 
    mutate(gbif_id = as.numeric(gbif_id))
  
  final_trait_data <- tax_hierarchy |>
    # Join datasets one by one
    left_join(combine_final, by = "gbif_id") |> 
    left_join(ecke_final, by = "gbif_id") |> 
    left_join(iucn_final, by = "gbif_id")
  
  write_rds(final_trait_data, here("data", "processed", "trait_data.rds"))
} else {
  
  trait_data <- read_rds(here("data", "processed", "trait_data.rds"))
  
}

# --- 4. Load Trees ---
tree_input <- read.nexus("data/external/output.nex")

mammal_tree_mcc <- maxCladeCred(tree_input, tree = TRUE)

tree_tips_df <- tibble(tip_label = mammal_tree_mcc$tip.label,
                       clean_name = str_replace_all(tip_label, "_", " "))

tree_tips_resolved <- add_gbif_ids(tree_tips_df, name_col = "clean_name") |>
  filter(!is.na(gbif_id)) |>
  group_by(gbif_id) |>
  slice(1) |>
  ungroup()

common_ids <- intersect(trait_data |>
                          drop_na(gbif_id) |>
                          pull(gbif_id), tree_tips_resolved$gbif_id)

message(paste("Found", length(common_ids), "overlapping species via GBIF IDs."))

trait_data_phylo <- trait_data |>
  filter(gbif_id %in% common_ids) |>
  left_join(tree_tips_resolved |> select(gbif_id, tip_label), by = "gbif_id") |>
  group_by(gbif_id) |>
  slice(1) |>
  ungroup() |>
  select(tip_label, gbif_id, everything())

tree_matched <- keep.tip(mammal_tree_mcc, trait_data_phylo$tip_label)

trait_data_phylo_sorted <- trait_data_phylo[match(tree_matched$tip.label, trait_data_phylo$tip_label), ]

# D. Final Verification & Save
if (all(trait_data_phylo_sorted$tip_label == tree_matched$tip.label)) {
  
  message(paste("Alignment Successful: N =", nrow(trait_data_phylo_sorted)))
  
  write_rds(trait_data_phylo_sorted, here("data", "processed", "trait_data_phylo_matched.rds"))
  
  write.tree(tree_matched, here("data", "processed", "mammal_tree_matched.tre"))
  
}

# Data for descriptive analysis of range biases ---------------------------
realm_path <- here("data", "external", "One Earth Realms", "Realm2023.shp")
# Accessed through https://www.oneearth.org/bioregions-2023/
realms <- vect(realm_path)
