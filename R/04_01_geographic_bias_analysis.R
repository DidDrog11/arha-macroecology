# ---
# 04_01_geographic_bias_analysis.R
#
# Purpose: 
# 1. Map surveillance effort against host species ranges (IUCN).
# 2. Quantify "Proportion of Range Sampled" using Administrative (ADM2) intersections.
# 3. Interpolate sampling gaps to identify neglected regions.
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, terra, tidyterra, sf, rnaturalearth, 
               metR, scales, ggrepel, countrycode, akima, WDI, brms, ggnewscale,
               marginaleffects)

# Load Processed Data
target_orders <- c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA")
host_traits <- read_rds(here("data", "processed", "trait_data.rds"))
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2026-01-09.rds"))
host_data <- arha_db$host

# Constants
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Ensure directories exist
dir.create(here("data", "gadm"), showWarnings = FALSE)
dir.create(here("output", "coverage_cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "species_coverage_maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "models"), recursive = TRUE, showWarnings = FALSE)

# --- 2. Prepare Spatial Data ---
lookup_path <- here("data", "external", "adm2_genus_lookup.rds")
cache_dir <- here("data", "temp", "spatial_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(lookup_path)) {
  
  # 1. Setup Parallel Workers
  if(!require("furrr")) install.packages("furrr")
  library(furrr)
  options(future.globals.maxSize = 8000 * 1024^2) 
  
  plan(multisession, workers = 8) 
  
  # A. Load & Clean IUCN Ranges
  message("Loading and cleaning IUCN data...")
  iucn_vect <- vect(here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
  iucn_subset <- iucn_vect[iucn_vect$sci_name %in% host_traits$species, ]
  
  iucn_proj <- project(iucn_subset, mollweide_crs)
  iucn_proj <- makeValid(iucn_proj) 
  iucn_proj <- aggregate(iucn_proj, by = "sci_name") 
  
  traits_clean <- host_traits |> 
    select(species, genus) |> 
    distinct(species, .keep_all = TRUE)
  
  iucn_df_safe <- values(iucn_proj) |>
    distinct(sci_name) |>
    left_join(traits_clean, by = c("sci_name" = "species"))
  
  # B. Load GADM
  gadm_path <- here("data", "gadm", "gadm_adm2_combined.shp")
  gadm_vect <- vect(gadm_path)
  gadm_proj <- project(gadm_vect, mollweide_crs)
  gadm_proj$area_km2 <- expanse(gadm_proj, unit = "km")
  
  # C. Wrap Objects 
  message("Wrapping objects (This takes a moment)...")
  gadm_packed <- wrap(gadm_proj)
  iucn_packed <- wrap(iucn_proj)
  
  countries <- unique(gadm_proj$GID_0)
  
  # D. Parallel Loop with Caching
  message("Starting intersection. Checkpoints saved to: data/temp/spatial_cache/")
  
  future_walk(countries, function(iso) {
    # 0. CHECKPOINT: Skip if already done
    outfile <- file.path(cache_dir, paste0("lookup_", iso, ".rds"))
    if (file.exists(outfile)) return(NULL)
    
    # 1. Unwrap
    library(terra) 
    library(dplyr)
    g_vect <- unwrap(gadm_packed)
    i_vect <- unwrap(iucn_packed)
    
    # 2. Filter Country
    adm_country <- g_vect[g_vect$GID_0 == iso, ]
    country_border <- aggregate(adm_country)
    country_border <- buffer(country_border, width = 0)
    
    # 3. Filter Species
    has_species <- relate(i_vect, country_border, "intersects")[,1]
    iucn_country <- i_vect[has_species, ]
    
    # Early exit if empty (save empty file to prevent re-running)
    if (nrow(iucn_country) == 0) {
      write_rds(tibble(), outfile)
      return(NULL)
    }
    
    # 4. Intersect
    inter_mat <- relate(adm_country, iucn_country, "intersects")
    hits <- which(inter_mat, arr.ind = TRUE)
    
    if (nrow(hits) == 0) {
      write_rds(tibble(), outfile)
      return(NULL)
    }
    
    # 5. Build & Save Result
    species_names <- iucn_country$sci_name[hits[,2]]
    
    genus_map <- tibble(sci_name = species_names) |>
      left_join(iucn_df_safe, by = "sci_name")
    
    res <- tibble(GID_2 = adm_country$GID_2[hits[,1]],
                  genus = genus_map$genus) |>
      distinct()
    
    write_rds(res, outfile)
    
    # Clean up RAM explicitly inside worker
    rm(g_vect, i_vect, inter_mat)
    gc()
    
  }, .options = furrr_options(seed = TRUE, scheduling = 1)) 
  
  # Merge all checkpoints
  lookup_files <- list.files(cache_dir, full.names = TRUE, pattern = "\\.rds$")
  lookup_long <- map_dfr(lookup_files, read_rds)
  
  write_rds(lookup_long, lookup_path)
  writeVector(gadm_proj, here("data", "gadm", "gadm_adm2_simplified.shp"), overwrite = TRUE) 
  } else { 
    lookup_long <- read_rds(lookup_path)
    gadm_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp")) 
    }

# --- 2b. Finish Remaining Countries (Serial Mode + Batching) ---

if(!file.exists(lookup_path)) {
  # 1. Some countries were very slow to process the following code conducts that for the remainder
  cache_files <- list.files(cache_dir, pattern = "\\.rds$")
  done_iso <- str_remove_all(cache_files, "lookup_|.rds")
  if(!exists("gadm_proj")) {
    gadm_vect <- vect(here("data", "gadm", "gadm_adm2_combined.shp"))
    gadm_proj <- project(gadm_vect, mollweide_crs)
  }
  all_iso  <- unique(gadm_proj$GID_0)
  remaining_iso <- setdiff(all_iso, done_iso)
  
  if (length(remaining_iso) > 0) {
    # Load & Simplify IUCN (Global)
    if(!exists("iucn_proj")) {
      message("Reloading IUCN data...")
      iucn_vect <- vect(here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
      iucn_subset <- iucn_vect[iucn_vect$sci_name %in% host_traits$species, ]
      iucn_proj <- project(iucn_subset, mollweide_crs)
      iucn_proj <- simplifyGeom(iucn_proj, tolerance = 100) # 100m tolerance
      iucn_proj <- makeValid(iucn_proj)
      iucn_proj <- aggregate(iucn_proj, by = "sci_name")
    }
    
    # Safe Lookup Table
    traits_clean <- host_traits |> 
      select(species, genus) |> 
      distinct(species, .keep_all = TRUE)
    
    iucn_df_safe <- values(iucn_proj) |>
      distinct(sci_name) |>
      left_join(traits_clean, by = c("sci_name" = "species"))
    
    walk(remaining_iso, function(iso) {
      message("Processing: ", iso, " ...")
      outfile <- file.path(cache_dir, paste0("lookup_", iso, ".rds"))
      # 1. Get Country Polygons
      adm_country <- gadm_proj[gadm_proj$GID_0 == iso, ]
      # Simplify
      adm_country <- simplifyGeom(adm_country, tolerance = 100)
      adm_country <- makeValid(adm_country)
      # 2. Filter Species (Fast Bounding Box)
      country_border <- aggregate(adm_country)
      # Use buffer(0) if aggregate fails or creates bad topology
      country_border <- buffer(country_border, width = 0)
      
      has_species <- relate(iucn_proj, country_border, "intersects")[,1]
      iucn_country <- iucn_proj[has_species, ]
      
      if (nrow(iucn_country) == 0) {
        write_rds(tibble(), outfile)
        return(NULL)
      }
      
      # 3. Batch Processing
      # If a country has > 200 districts (like Romania), process in chunks.
      n_districts <- nrow(adm_country)
      batch_size <- 50
      
      results_list <- list()
      
      # Create batches
      batches <- split(1:n_districts, ceiling(seq_along(1:n_districts)/batch_size))
      
      message(paste0("  -> Splitting into ", length(batches), " batches for stability..."))
      
      for(k in seq_along(batches)) {
        
        idx <- batches[[k]]
        adm_chunk <- adm_country[idx, ]
        
        inter_mat <- relate(adm_chunk, iucn_country, "intersects")
        hits <- which(inter_mat, arr.ind = TRUE)
        
        if(nrow(hits) > 0) {
              species_names <- iucn_country$sci_name[hits[,2]]
              genus_map <- tibble(sci_name = species_names) |>
                left_join(iucn_df_safe, by = "sci_name")
          
          chunk_res <- tibble(GID_2 = adm_chunk$GID_2[hits[,1]], # Extract ID directly from chunk
                              genus = genus_map$genus)
          
          results_list[[k]] <- chunk_res
        }
        if(k %% 5 == 0) message(paste0("     Batch ", k, "/", length(batches), " done."))
      }
      
      # 4. Combine Batches & Save
      if(length(results_list) > 0) {
        res <- bind_rows(results_list) |> distinct()
        write_rds(res, outfile)
        message("  -> Success: ", iso)
      } else {
        write_rds(tibble(), outfile)
        message("  -> Success (Empty): ", iso)
      }
      
      gc() # Clean RAM
    })
  }
  
  # Final Merge
  lookup_files <- list.files(cache_dir, full.names = TRUE, pattern = "\\.rds$")
  lookup_long <- map_dfr(lookup_files, read_rds)
  write_rds(lookup_long, lookup_path)
} else {
  lookup_long <- read_rds(lookup_path)
  gadm_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp"))
}

# Section A: Maps and Coverage --------------------------------------------

# --- Quantify Sampling Coverage (By Genus) ---
sampled_districts <- host_data |>
  filter(toupper(host_order) %in% target_orders) |>
  drop_na(gadm_adm2, host_genus) |>
  filter(gadm_adm2 != "") |>
  group_by(gadm_adm2, host_genus) |>
  summarise(n_individuals = sum(number_of_hosts, na.rm = TRUE), .groups = "drop") |>
  rename(GID_2 = gadm_adm2)

target_genera <- host_traits |>
  drop_na(genus) |>
  distinct(genus) |>
  pull(genus)

if(!"area_km2" %in% names(gadm_proj)) {
  gadm_proj$area_km2 <- expanse(gadm_proj, unit = "km")
  writeVector(gadm_proj, here("data", "gadm", "gadm_adm2_simplified.shp"))
}

range_areas <- lookup_long |>
  inner_join(as_tibble(gadm_proj), by = "GID_2") |>
  group_by(genus) |>
  summarise(range_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")

sampled_areas <- sampled_districts |>
  rename(genus = host_genus) |>
  inner_join(lookup_long, 
             by = c("GID_2", "genus")) |>
  inner_join(as_tibble(gadm_proj), by = "GID_2") |>
  group_by(genus) |>
  summarise(sampled_area_km2 = sum(area_km2, na.rm = TRUE),
            n_individuals_sampled = sum(n_individuals, na.rm = TRUE),
            .groups = "drop")

coverage_summary <- range_areas |>
  left_join(sampled_areas, by = "genus") |>
  mutate(sampled_area_km2 = replace_na(sampled_area_km2, 0),
         n_individuals_sampled = replace_na(n_individuals_sampled, 0),
         prop_sampled = sampled_area_km2 / range_area_km2)

plot_data <- coverage_summary |>
  filter(sampled_area_km2 > 0, range_area_km2 > 0) |>
  mutate(label_text = case_when(rank(desc(n_individuals_sampled)) <= 5 ~ genus,
                                rank(desc(prop_sampled)) <= 5 ~ genus,
                                rank(desc(range_area_km2)) <= 8 ~ genus, 
                                TRUE ~  ""),
         label_category = case_when(label_text == "" ~ "Other",
                                    rank(desc(n_individuals_sampled)) <= 5 ~ "Sample N",
                                    rank(desc(prop_sampled)) <= 5 ~ "Coverage",
                                    rank(desc(range_area_km2)) <= 8 ~ "Range Size",
                                    TRUE ~ "Other"),
         label_category = fct_relevel(label_category, "Sample N", "Coverage", "Range Size", "Other"))

# Interpolate (Akima)
interp_res <- akima::interp(x = log10(plot_data$sampled_area_km2),
                            y = log10(plot_data$range_area_km2),
                            z = plot_data$prop_sampled,
                            linear = TRUE, nx = 200, ny = 200)

interp_df <- as.data.frame(akima::interp2xyz(interp_res)) |>
  as_tibble() |>
  mutate(x = 10^x, y = 10^y)

# Plot
p_bias <- ggplot() +
  # The interpolated heat map
  geom_raster(data = interp_df, aes(x = x, y = y, fill = z)) +
  scale_fill_viridis_c(name = "Sampling\nCoverage", labels = percent, option = "magma") +
  # Contours
  geom_contour(data = interp_df, aes(x = x, y = y, z = z), colour = "white", alpha = 0.5) +
  # The points (Genera)
  geom_point(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2), 
             shape = 21, fill = "white", alpha = 0.7) +
  new_scale_fill() +
  geom_label_repel(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2, label = label_text, fill = label_category),
                   size = 3, colour = "white", fontface = "bold", max.overlaps = 30, box.padding = 0.5, min.segment.length = 0) +
  # Scales
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_fill_manual(name = "Label Criteria", values = c("Sample N" = "darkgreen", "Coverage" = "firebrick", "Range Size" = "steelblue","Other" = NA),
                    na.translate = FALSE) +

  labs(title = "Geographic Sampling Bias by Genus",
       subtitle = "Relationship between Range Size and Sampled Area",
       x = "Sampled Area (km²)", y = "Total Range Size (km²)") +
  theme_minimal()

ggsave(plot = p_bias, filename = here("output", "figures", "sampling_bias_by_genus.png"), width = 8, height = 8)

genus_district_counts <- sampled_districts |>
  group_by(host_genus) |>
  summarise(n_districts = n_distinct(GID_2)) |>
  rename(genus = host_genus)

genus_stats_full <- coverage_summary |>
  left_join(genus_district_counts, by = "genus") |>
  arrange(desc(range_area_km2))

example_genus <- genus_stats_full |> slice(1)

# Country level modelling
# Aggregation
country_summary <- host_data |>
  drop_na(iso3c) |>
  group_by(iso3c) |>
  summarise(total_hosts = sum(number_of_hosts, na.rm = TRUE),
            n_species = n_distinct(gbif_id),
            .groups = "drop") |>
  filter(total_hosts > 0)

# WDI Data Fetching
wdi_path <- here("data", "external", "wdi_gdp_pop.rds")

if (!file.exists(wdi_path)) {
  iso_map <- country_summary |>
    distinct(iso3c) |>
    mutate(iso2c = countrycode(iso3c, "iso3c", "iso2c")) |>
    drop_na(iso2c)
  
  wdi_all <- WDI(country = "all", 
                 indicator = c("gdp_per_capita" = "NY.GDP.PCAP.KD", 
                               "population" = "SP.POP.TOTL"),
                 start = 2010, end = 2020, extra = TRUE) |>
    as_tibble()
  
  wdi_clean <- wdi_all |>
    filter(iso2c %in% iso_map$iso2c) |>
    group_by(iso2c) |>
    summarise(gdp_per_capita = mean(gdp_per_capita, na.rm = TRUE),
              population = mean(population, na.rm = TRUE),
              .groups = "drop") |>
    right_join(iso_map, by = "iso2c")
  
  wdi_final <- wdi_clean |>
    mutate(gdp_per_capita = case_when(# French Overseas (Approx as ~France)
      iso3c %in% c("GUF", "MYT", "REU") & is.na(gdp_per_capita) ~ 37000,
      # Taiwan (IMF Estimate avg 2010-2020)
      iso3c == "TWN" & is.na(gdp_per_capita) ~ 25000,                    
      # North Korea (Bank of Korea estimate)
      iso3c == "PRK" & (is.na(gdp_per_capita) | is.nan(gdp_per_capita)) ~ 640,                      
      # Venezuela
      iso3c == "VEN" & is.na(gdp_per_capita) ~ 16000,                    
      TRUE ~ gdp_per_capita),
      population = case_when(iso3c == "GUF" & is.na(population) ~ 300000,
                             iso3c == "MYT" & is.na(population) ~ 270000,
                             iso3c == "REU" & is.na(population) ~ 850000,
                             iso3c == "TWN" & is.na(population) ~ 23500000,
                             TRUE ~ population)) |>
    filter(!is.na(gdp_per_capita))
  
  write_rds(wdi_final, wdi_path)
} else {
  wdi_final <- read_rds(wdi_path)
}

# Fit Country Model (brms)
country_model_data <- country_summary |>
  left_join(wdi_final, by = "iso3c") |>
  drop_na(gdp_per_capita, population) |>
  mutate(log_gdp = as.numeric(scale(log10(gdp_per_capita))),
         log_pop = log10(population))

# Model: Count ~ GDP + offset(Population)
# Does wealth predict sampling intensity after accounting for population size?
model_country_int <- brm(total_hosts ~ log_gdp + offset(log_pop),
                         data = country_model_data,
                         family = negbinomial(),
                         prior = c(prior(normal(0, 1), class = "b"), prior(normal(0, 5), class = "Intercept")),
                         chains = 4, iter = 2000, cores = 6, 
                         file = here("output", "models", "brms_country_intensity"))

# Section C: Subnational (ADM2) modelling ---------------------------------
# Prepare Raster Predictors (Zonal Stats)
# This step extracts mean Nightlights, Pop Density, etc. for every ADM2 polygon
adm2_pred_path <- here("data", "external", "adm2_predictors.rds")

if (!file.exists(adm2_pred_path)) {
  
  r_pop <- rast(here("data", "external", "rasters", "worldpop_2025_25km.tif"))
  r_light <- rast(here("data", "external", "rasters", "viirs_2024_25km.tif"))
  r_access <- rast(here("data", "external", "rasters", "travel_time_25km.tif"))
  
  if (!exists("iucn_proj")) {
    iucn_vect <- vect(here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
    iucn_subset <- iucn_vect[toupper(iucn_vect$order_) %in% target_orders, ] 
    iucn_proj_filtered <- project(iucn_subset, mollweide_crs)
  } else {
    # If already loaded, just filter
    iucn_proj_filtered <- iucn_proj[toupper(iucn_proj$order_) %in% target_orders, ]
  }
  
  r_template <- r_pop
  iucn_proj_filtered$richness_constant <- 1
  r_richness <- rasterize(iucn_proj_filtered, r_template, field = "richness_constant", fun = "sum", background = 0)
  writeRaster(r_richness, here("data", "external", "rasters", "host_richness_25km.tif"), overwrite = TRUE)
  
  # Calculate Zonal Stats
  area_km2 = expanse(gadm_proj, unit = "km")
  pop_dens = extract(r_pop, gadm_proj, fun = mean, na.rm=TRUE)[,2]
  night_light = extract(r_light, gadm_proj, fun = mean, na.rm=TRUE)[,2]
  travel_time = extract(r_access, gadm_proj, fun = mean, na.rm=TRUE)[,2]
  host_richness = extract(r_richness, gadm_proj, fun = mean, na.rm=TRUE)[,2]
  
  adm2_stats <- tibble(GID_0 = gadm_proj$GID_0,
                       GID_2 = gadm_proj$GID_2,
                       area_km2 = area_km2,
                       pop_dens = pop_dens,
                       night_light = night_light,
                       travel_time = travel_time,
                       host_richness = host_richness)
  
  write_rds(adm2_stats, adm2_pred_path)
} else {
  adm2_stats <- read_rds(adm2_pred_path)
}

# Prepare Response (Aggregated counts per ADM2)
adm2_response <- host_data |>
  drop_na(gadm_adm2) |>
  group_by(gadm_adm2) |>
  summarise(n_hosts = sum(number_of_hosts, na.rm=TRUE))

# Merge & Scale
model_adm2_data <- adm2_stats |>
  left_join(adm2_response, by = c("GID_2" = "gadm_adm2")) |>
  mutate(n_hosts = replace_na(n_hosts, 0)) |>
  drop_na(pop_dens, night_light, travel_time) |>
  mutate(s_pop = as.numeric(scale(log1p(pop_dens))),
         s_light = as.numeric(scale(log1p(night_light))),
         s_access = as.numeric(scale(log1p(travel_time))),
         s_richness = as.numeric(scale(host_richness)),
         log_area = log(area_km2))

write_rds(model_adm2_data, here("output", "models", "model_adm2_data_full.rds"))

# Fit ZINB Model (Zero-Inflated Negative Binomial)
# Formula: 
#   Count ~ Density + Light + Access + (1|Country)
#   Zero_Inflation ~ Density + Light + Access
model_adm2_zinb <- brm(bf(n_hosts ~ s_pop + s_light + s_access + s_richness + offset(log_area) + (1 | GID_0),
                          zi ~ s_pop + s_light + s_access + s_richness),
                       data = model_adm2_data,
                       family = zero_inflated_negbinomial(),
                       prior = c(prior(normal(0, 1), class = "b"), prior(student_t(3, 0, 2.5), class = "sd")),
                       chains = 4, 
                       cores = 4,
                       control = list(adapt_delta = 0.95),
                       file = here("output", "models", "brms_adm2_zinb"))

summary(model_adm2_zinb)

# Geographic Realms -------------------------------------------------------
realm_path <- here("data", "external", "One Earth Realms", "Realm2023.shp")
realms_vect <- vect(realm_path)

points <- arha_db$host |> 
  drop_na(longitude, latitude) |>
  group_by(longitude, latitude) |>
  summarise(number_of_hosts = sum(number_of_hosts, na.rm = TRUE)) |>
  vect(geom = c("longitude", "latitude"), crs = "EPSG:4326")

points$realm <- extract(realms_vect, points)$BioGeoRelm

as_tibble(points) |>
  group_by(realm) |>
  summarise(proportion_of_hosts = sum(number_of_hosts, na.rm = TRUE) / sum(points$number_of_hosts))
