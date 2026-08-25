# --- Snapshot headline numbers after the taxonomic-filter re-run --------
library(bayestestR)

snapshot_model <- function(fit, label, pc1_param = "pace_of_life_pc1") {
  fx <- fixef(fit, probs = c(0.025, 0.975))
  vc <- VarCorr(fit)
  tibble(
    model          = label,
    n_obs          = nobs(fit),
    n_species      = length(unique(fit$data$tip_label)),
    pc1_mean       = round(fx[pc1_param, "Estimate"], 4),
    pc1_lower      = round(fx[pc1_param, "Q2.5"], 4),
    pc1_upper      = round(fx[pc1_param, "Q97.5"], 4),
    pc1_pd         = round(as.numeric(p_direction(fit, parameters = pc1_param)$pd) * 100, 2),
    effort_mean    = round(fx["log_effort", "Estimate"], 4),
    sd_host        = round(vc$host_id$sd[, "Estimate"], 4),
    sd_phylo       = round(vc$tip_label$sd[, "Estimate"], 4),
    sd_virus       = round(vc$pathogen_species_cleaned$sd[, "Estimate"], 4)
  )
}

model_snapshot <- bind_rows(
  snapshot_model(fit_dyadic_a,      "Global (S2)"),
  snapshot_model(fit_dyadic_b,      "Synanthropy subset (S3)"),
  snapshot_model(fit_dyadic_a_sens, "Strict / active infection (S4)"),
  snapshot_model(fit_dyadic_c,      "Complete case (S5)", pc1_param = "pace_of_life_pc1_raw")
)

# Synanthropy coefficient (only in fit_dyadic_b)
syn_fx <- fixef(fit_dyadic_b, probs = c(0.025, 0.975))
syn_snapshot <- tibble(
  metric = c("Synanthropy (Total) mean", "Synanthropy (Total) lower", "Synanthropy (Total) upper",
             "Synanthropy (Occasional) mean"),
  value  = c(round(syn_fx["synanthropy_statusTotallySynanthropic", "Estimate"], 4),
             round(syn_fx["synanthropy_statusTotallySynanthropic", "Q2.5"], 4),
             round(syn_fx["synanthropy_statusTotallySynanthropic", "Q97.5"], 4),
             round(syn_fx["synanthropy_statusOccasionallySynanthropic", "Estimate"], 4))
)

# Data-level and validation numbers
other_snapshot <- tibble(
  metric = c("N species in trait data", "N species in global model",
             "Effort anchor (95th pctile log_effort)",
             "PC1 variance explained (%)",
             "LORO Spearman rho",
             "AUC Americas", "AUC Asia", "AUC Europe", "AUC Africa",
             "Genus overlap Americas (%)", "Genus overlap Asia (%)",
             "Genus overlap Europe (%)", "Genus overlap Africa (%)"),
  value  = c(nrow(host_traits),
             length(unique(fit_dyadic_a$data$tip_label)),
             round(quantile(fit_dyadic_a$data$log_effort, 0.95), 4),
             NA,   # fill from summary(pca_res) if not in scope
             round(as.numeric(cv_correlation$estimate), 4),
             cv_auc_summary$auc[cv_auc_summary$fold == "Americas"],
             cv_auc_summary$auc[cv_auc_summary$fold == "Asia"],
             cv_auc_summary$auc[cv_auc_summary$fold == "Europe"],
             cv_auc_summary$auc[cv_auc_summary$fold == "Africa"],
             phylo_overlap$pct_genera_shared[phylo_overlap$fold == "Americas"],
             phylo_overlap$pct_genera_shared[phylo_overlap$fold == "Asia"],
             phylo_overlap$pct_genera_shared[phylo_overlap$fold == "Europe"],
             phylo_overlap$pct_genera_shared[phylo_overlap$fold == "Africa"])
)

dir.create(here("output", "tables", "snapshot"), recursive = TRUE, showWarnings = FALSE)
write_csv(model_snapshot,  here("output", "tables", "snapshot", "post_filter_models.csv"))
write_csv(syn_snapshot,    here("output", "tables", "snapshot", "post_filter_synanthropy.csv"))
write_csv(other_snapshot,  here("output", "tables", "snapshot", "post_filter_other.csv"))

print(model_snapshot)
print(syn_snapshot)
print(other_snapshot)

# --- Additional pre-filter snapshots: regional comparison & spatial ------
# Run from the 07_04 session where risk_summary, stats, district_classed,
# summary_table, r_hazard and r_richness are in scope.

snap_dir <- here("output", "tables", "snapshot")
dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Regional LORO comparison (medians, means, N per region)
write_csv(risk_summary, file.path(snap_dir, "post_filter_risk_summary.csv"))

# 2. Pairwise comparison matrix — flatten to long form so it's diffable
pairwise_long <- as.data.frame(stats$p.value) |>
  rownames_to_column("region_a") |>
  pivot_longer(-region_a, names_to = "region_b", values_to = "p_value") |>
  drop_na(p_value) |>
  mutate(p_value = signif(p_value, 3))

write_csv(pairwise_long, file.path(snap_dir, "post_filter_pairwise.csv"))

# 3. District archetype table (the Fig. 3 / Table 1 examples)
write_csv(summary_table, file.path(snap_dir, "post_filter_district_archetypes.csv"))

# 4. Bivariate class summary — the Mean Host Probability / Richness values
#    quoted in the Results text for Classes 1-3, 3-1 and 3-3
class_summary <- district_classed |>
  group_by(bi_class) |>
  summarise(n_districts = n(),
            mean_host_prob = round(mean(mean_intrinsic_hazard, na.rm = TRUE), 4),
            mean_richness = round(mean(mean_district_richness, na.rm = TRUE), 2),
            .groups = "drop")

write_csv(class_summary, file.path(snap_dir, "post_filter_bivariate_classes.csv"))

# 5. Raster-level summary — global distribution of the two mapped quantities
raster_summary <- tibble(
  metric = c("Hazard: min", "Hazard: median", "Hazard: mean", "Hazard: max",
             "Richness: min", "Richness: median", "Richness: mean", "Richness: max",
             "N non-NA cells (hazard)"),
  value = c(round(as.numeric(global(r_hazard, "min", na.rm = TRUE)), 5),
            round(as.numeric(global(r_hazard, fun = median, na.rm = TRUE)), 5),
            round(as.numeric(global(r_hazard, "mean", na.rm = TRUE)), 5),
            round(as.numeric(global(r_hazard, "max", na.rm = TRUE)), 5),
            as.numeric(global(r_richness, "min", na.rm = TRUE)),
            as.numeric(global(r_richness, fun = median, na.rm = TRUE)),
            round(as.numeric(global(r_richness, "mean", na.rm = TRUE)), 2),
            as.numeric(global(r_richness, "max", na.rm = TRUE)),
            as.numeric(global(!is.na(r_hazard), "sum", na.rm = TRUE)))
)

write_csv(raster_summary, file.path(snap_dir, "post_filter_raster_summary.csv"))

# 6. Country-level Class 3-3 area percentages (if country_area_pct in scope)
if (exists("country_area_pct")) {
  write_csv(country_area_pct, file.path(snap_dir, "post_filter_country_area_pct.csv"))
}

print(class_summary)
print(raster_summary)

# Class summary — the values quoted in the Results bivariate paragraph
district_classed |> group_by(bi_class) |>
  summarise(n = n(),
            mean_host_prob = round(mean(mean_intrinsic_hazard), 4),
            mean_richness = round(mean(mean_district_richness), 1), .groups = "drop")

# Richness-probability relationship on the masked data
summary(lm(mean_intrinsic_hazard ~ mean_district_richness, data = table_out))

# Range of the mapped surface, for the caption
summary(table_out$mean_intrinsic_hazard)

# Confirm top_purple composition still holds after the mask change
summary_table |> filter(Category == "3-3")


check <- gadm_sf |> filter(NAME_2 %in% c("Amami", "Bengkulu")) |> vect()

plot_district <- function(poly, title) {
  # Buffer the extent slightly so surrounding cells are visible for context
  e <- terra::ext(poly)
  buffer <- 60000  # 60km, i.e. 3 cells
  e_buffered <- terra::ext(e[1] - buffer, e[2] + buffer, e[3] - buffer, e[4] + buffer)
  r_crop <- crop(r_hazard, e_buffered)
  rich_crop <- crop(r_richness, e_buffered)
  
  ggplot() +
    geom_spatraster(data = r_crop) +
    geom_spatvector(data = poly, fill = NA, colour = "red", linewidth = 1) +
    scale_fill_viridis_c(option = "magma", na.value = "grey90",
                         name = "Mean host\nprobability") +
    labs(title = title,
         subtitle = paste0("Cells shown: ", sum(!is.na(values(r_crop))),
                           " | District richness: ", 
                           round(mean(values(rich_crop), na.rm = TRUE), 1))) +
    theme_minimal()
}

p1 <- plot_district(check[1], "Bengkulu, Indonesia")
p2 <- plot_district(check[2], "Amami, Japan")

p1 + p2  # patchwork


# ---
# 08_reconcile_counts.R
# Establish species and dyad counts in sequence from a fresh session.
# Run after rm(list = ls()).
# ---

pacman::p_load(tidyverse, here, brms)

arha_db     <- read_rds(here("data", "database", "Project_ArHa_database_2026-08-17.rds"))
host_traits <- read_rds(here("data", "analytic", "host_traits_final.rds"))

fit_global <- read_rds(here("output", "models", "brms_dyadic_N49k.rds"))
fit_syn    <- read_rds(here("output", "models", "brms_dyadic_N19k.rds"))
fit_strict <- read_rds(here("output", "models", "brms_dyadic_sens_strict.rds"))
fit_case   <- read_rds(here("output", "models", "brms_dyadic_sens_raw.rds"))

host_species_of <- function(fit) unique(fit$data$tip_label) |> str_replace_all("_", " ")

# --- 1. Model dimensions -------------------------------------------------
model_dims <- tibble(
  model = c("Global", "Synanthropy subset", "Strict / active infection", "Complete case"),
  dyads = c(nobs(fit_global), nobs(fit_syn), nobs(fit_strict), nobs(fit_case)),
  species = c(length(host_species_of(fit_global)), length(host_species_of(fit_syn)),
              length(host_species_of(fit_strict)), length(host_species_of(fit_case))),
  pathogens = c(n_distinct(fit_global$data$pathogen_species_cleaned),
                n_distinct(fit_syn$data$pathogen_species_cleaned),
                n_distinct(fit_strict$data$pathogen_species_cleaned),
                n_distinct(fit_case$data$pathogen_species_cleaned))
)
print(model_dims)

# --- 2. Database counts --------------------------------------------------
species_assayed_any <- arha_db$pathogen |>
  left_join(arha_db$host |> select(host_record_id, host_species), by = "host_record_id") |>
  drop_na(host_species) |> distinct(host_species) |> pull()

species_assayed_target <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  left_join(arha_db$host |> select(host_record_id, host_species), by = "host_record_id") |>
  drop_na(host_species) |> distinct(host_species) |> pull()

n_assays_target <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |> nrow()

n_individuals_target <- arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) |>
  summarise(n = sum(number_tested, na.rm = TRUE)) |> pull()

counts <- tibble(
  quantity = c("Distinct species in host table",
               "Species with any assay",
               "Species assayed for the two target families",
               "Assay records, target families",
               "Individuals tested, target families",
               "Species in trait/phylogeny dataset (post-filter)"),
  n = c(arha_db$host |> distinct(host_species) |> drop_na() |> nrow(),
        length(species_assayed_any), length(species_assayed_target),
        n_assays_target, n_individuals_target, nrow(host_traits))
)
print(counts)

# --- 3. Reconciliation ---------------------------------------------------
species_modelled <- host_species_of(fit_global)

reconciliation <- tibble(
  category = c("Assayed for target families AND modelled",
               "Assayed for target families, NOT modelled",
               "Modelled, NOT assayed for target families (all-zero dyads)"),
  n = c(length(intersect(species_assayed_target, species_modelled)),
        length(setdiff(species_assayed_target, species_modelled)),
        length(setdiff(species_modelled, species_assayed_target)))
)
print(reconciliation)

# --- 4. Why species drop out ---------------------------------------------
trait_all <- read_rds(here("data", "processed", "trait_data_phylo_matched.rds")) |>
  mutate(host_species = str_replace_all(tip_label, "_", " "))

dropped_reasons <- tibble(host_species = setdiff(species_assayed_target, species_modelled)) |>
  left_join(trait_all |> select(host_species, order), by = "host_species") |>
  mutate(reason = case_when(
    is.na(order) ~ "No trait/phylogeny match",
    !toupper(order) %in% c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA") ~ "Outside target orders",
    TRUE ~ "Target order, incomplete data"))

dropped_reasons |> count(reason, sort = TRUE) |> print()

# --- 5. Positive detections lost -----------------------------------------
arha_db$pathogen |>
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae"), number_positive > 0) |>
  left_join(arha_db$host |> select(host_record_id, host_species), by = "host_record_id") |>
  filter(host_species %in% dropped_reasons$host_species) |>
  group_by(host_species) |>
  summarise(n_positive = sum(number_positive, na.rm = TRUE), .groups = "drop") |>
  left_join(dropped_reasons, by = "host_species") |>
  arrange(desc(n_positive)) |> print(n = Inf)
