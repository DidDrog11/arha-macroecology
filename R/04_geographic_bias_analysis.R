# ---
# 04_geographic_bias_analysis.R
#
# Purpose: 
# 1. Map surveillance effort against host species ranges (IUCN).
# 2. Quantify "Proportion of Range Sampled" using Administrative (ADM2) intersections.
# 3. Interpolate sampling gaps to identify neglected regions.
# ---

# --- 1. Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, terra, tidyterra, sf, rnaturalearth, 
               metR, scales, ggrepel, countrycode, akima, WDI, brms)

# Load Processed Data
target_orders <- c("RODENTIA", "EULIPOTYPHLA", "SORICOMORPHA", "ERINACEOMORPHA")
host_traits <- read_rds(here("data", "processed", "trait_data.rds"))
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2025-12-02.rds"))
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

if (!file.exists(lookup_path)) {
  
  message("--- Generating Spatial Lookup Table (This takes time) ---")
  
  # A. Load IUCN Ranges
  iucn_vect <- vect(here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
  
  # Filter to our target species
  iucn_subset <- iucn_vect[iucn_vect$sci_name %in% host_traits$species, ]
  iucn_proj <- project(iucn_subset, mollweide_crs)
  
  # Attach Genus from our clean trait data
  iucn_proj_df <- as.data.frame(iucn_proj) |>
    left_join(host_traits |> select(species, genus), by = c("sci_name" = "species"))
  
  # B. Load GADM (Administrative Boundaries)
  # This file is available at an associated zenodo repository
  gadm_path <- here("data", "gadm", "gadm_adm2_combined.shp")
  if(!file.exists(gadm_path)) stop("Please ensure 'gadm_adm2_combined.shp' is in data/gadm/")
  
  gadm_vect <- vect(gadm_path)
  gadm_proj <- project(gadm_vect, mollweide_crs)
  
  # Calculate Area of every District
  gadm_proj$area_km2 <- expanse(gadm_proj, unit = "km")
  
  # C. The Intersection Loop
  countries <- unique(gadm_proj$GID_0)
  
  lookup_long <- map_dfr(countries, function(iso) {
    
    # Subset Country
    adm_country <- gadm_proj[gadm_proj$GID_0 == iso, ]
    
    # Subset IUCN ranges that touch this country
    country_border <- aggregate(adm_country)
    iucn_country <- iucn_proj[relate(iucn_proj, country_border, "intersects"), ]
    
    if (nrow(iucn_country) == 0) return(NULL)
    
    # Intersect: Which District contains which Genus?
    inter_mat <- relate(adm_country, iucn_country, "intersects")
    
    # Convert to Dataframe
    hits <- which(inter_mat, arr.ind = TRUE)
    
    if(nrow(hits) == 0) return(NULL)
    
    tibble(GID_2 = adm_country$GID_2[hits[,1]],
           genus = iucn_proj_df$genus[hits[,2]]) |>
      distinct()
    
  }, .progress = TRUE)
  
  # Save
  write_rds(lookup_long, lookup_path)
  # Save simplified vectors for plotting later
  writeVector(gadm_proj, here("data", "gadm", "gadm_adm2_simplified.shp"), overwrite = TRUE)
  
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

# Review this section once new spatial matching is complete
####

if(!"area_km2" %in% names(gadm_proj)) {
  gadm_proj$area_km2 <- expanse(gadm_proj, unit = "km")
}
  
range_areas <- temp_lookup_long |>
  rename(genus = gbif_genus) |>
  inner_join(as_tibble(gadm_proj), by = "GID_2") |>
  group_by(genus) |>
  summarise(range_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")

sampled_areas <- sampled_districts |>
  rename(genus = host_genus) |>
  inner_join(temp_lookup_long |> rename(genus = gbif_genus), 
             by = c("GID_2", "genus")) %>%
  inner_join(as_tibble(gadm_proj), by = "GID_2") %>%
  group_by(genus) %>%
  summarise(sampled_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")

coverage_summary <- range_areas |>
  left_join(sampled_areas, by = "genus") |>
  mutate(sampled_area_km2 = replace_na(sampled_area_km2, 0),
         prop_sampled = sampled_area_km2 / range_area_km2)
####

plot_data <- coverage_summary |>
  filter(sampled_area_km2 > 0, range_area_km2 > 0) |>
  mutate(label = if_else(rank(desc(range_area_km2)) <= 10, genus, ""))

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
  geom_contour(data = interp_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
  # The points (Genera)
  geom_point(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2), 
             shape = 21, fill = "white", alpha = 0.7) +
  geom_label_repel(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2, label = label),
                   size = 3, max.overlaps = 20) +
  # Scales
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  labs(title = "Geographic Sampling Bias by Genus",
       subtitle = "Relationship between Range Size and Sampled Area",
       x = "Sampled Area (km²)", y = "Total Range Size (km²)") +
  theme_minimal()

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
                         chains = 4, iter = 2000, cores = 4, 
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
         log_area = log(area_km2))

# Fit ZINB Model (Zero-Inflated Negative Binomial)
# Formula: 
#   Count ~ Density + Light + Access + (1|Country)
#   Zero_Inflation ~ Density + Light + Access
model_adm2_zinb <- brm(bf(n_hosts ~ s_pop + s_light + s_access + offset(log_area) + (1 | GID_0),
                          zi ~ s_pop + s_light + s_access),
                       data = model_adm2_data,
                       family = zero_inflated_negbinomial(),
                       prior = c(prior(normal(0, 1), class = "b"), prior(student_t(3, 0, 2.5), class = "sd")),
                       chains = 4, iter = 2000, cores = 4,
                       control = list(adapt_delta = 0.95),
                       file = here("output", "models", "brms_adm2_zinb"))

summary(model_adm2_zinb)


# Section D: Visualisation of marginal effects ----------------------------
# Plot Conditional Effects (How does Light affect Sampling?)
p_effects <- plot(conditional_effects(model_adm2_zinb), ask = FALSE)

p_light <- plot_predictions(model_adm2_zinb, condition = "s_light", type = "response") +
  labs(title = "Effect of Night Lights on Sampling", y = "Predicted Hosts Sampled") +
  theme_minimal()

# Residual Map (Observed vs Predicted)
# Regions that are over/under sampled relative to their wealth
preds <- posterior_epred(model_adm2_zinb) |> apply(2, mean)
model_adm2_data$resid <- log10(model_adm2_data$n_hosts + 1) - log10(preds + 1)
