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
               metR, scales, ggrepel, countrycode)

# Load Processed Data
# We use the full descriptive dataset for the checklist
host_traits <- read_rds(here("data", "processed", "trait_data_all.rds"))
arha_db <- read_rds(here("data", "database", "Project_ArHa_database_2025-12-02.rds")) # Update date
host_data <- arha_db$host

# Constants
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Ensure directories exist
dir.create(here("data", "gadm"), showWarnings = FALSE)
dir.create(here("output", "analysis_1", "coverage_cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "analysis_1", "species_coverage_maps"), recursive = TRUE, showWarnings = FALSE)


# --- 2. Prepare Spatial Data (Lookup Tables) ---
# This section is cached because intersection takes a long time.

lookup_path <- here("data", "external", "adm2_genus_lookup.rds")

if (!file.exists(lookup_path)) {
  
  message("--- Generating Spatial Lookup Table (This takes time) ---")
  
  # A. Load IUCN Ranges
  # Note: Ensure you have the IUCN shapefile at this path
  iucn_vect <- vect(here("data", "external", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
  
  # Filter to our target species
  iucn_subset <- iucn_vect[iucn_vect$sci_name %in% host_traits$species, ]
  iucn_proj <- project(iucn_subset, mollweide_crs)
  
  # Attach Genus from our clean trait data (safer than parsing strings)
  # We join by scientific name
  iucn_proj_df <- as.data.frame(iucn_proj) %>%
    left_join(host_traits %>% select(species, genus), by = c("sci_name" = "species"))
  
  # B. Load GADM (Administrative Boundaries)
  # If you don't have this, use geodata::gadm() to download, but assuming you have the file:
  gadm_path <- here("data", "gadm", "gadm_adm2_combined.shp")
  if(!file.exists(gadm_path)) stop("Please ensure 'gadm_adm2_combined.shp' is in data/gadm/")
  
  gadm_vect <- vect(gadm_path)
  gadm_proj <- project(gadm_vect, mollweide_crs)
  
  # Calculate Area of every District
  gadm_proj$area_km2 <- expanse(gadm_proj, unit = "km")
  
  # C. The Intersection Loop (Country by Country to save RAM)
  countries <- unique(gadm_proj$GID_0)
  
  lookup_long <- map_dfr(countries, function(iso) {
    
    # 1. Subset Country
    adm_country <- gadm_proj[gadm_proj$GID_0 == iso, ]
    
    # 2. Subset IUCN ranges that touch this country
    # (Quick spatial filter using the country's bounding box)
    country_border <- aggregate(adm_country)
    iucn_country <- iucn_proj[relate(iucn_proj, country_border, "intersects"), ]
    
    if (nrow(iucn_country) == 0) return(NULL)
    
    # 3. Intersect: Which District contains which Genus?
    # returns matrix: Rows = Districts, Cols = IUCN Polygons
    inter_mat <- relate(adm_country, iucn_country, "intersects")
    
    # Convert to Dataframe
    # We find indices where intersection is TRUE
    hits <- which(inter_mat, arr.ind = TRUE)
    
    if(nrow(hits) == 0) return(NULL)
    
    tibble(
      GID_2 = adm_country$GID_2[hits[,1]],
      # Use the Genus we attached earlier
      genus = iucn_proj_df$genus[hits[,2]]
    ) %>%
      distinct() # Remove duplicates if a genus has multiple polygons in one district
    
  }, .progress = TRUE)
  
  # Save
  write_rds(lookup_long, lookup_path)
  # Save simplified vectors for plotting later
  writeVector(gadm_proj, here("data", "gadm", "gadm_adm2_simplified.shp"), overwrite=TRUE)
  
} else {
  lookup_long <- read_rds(lookup_path)
  gadm_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp"))
}


# --- 3. Quantify Sampling Coverage (By Genus) ---

message("--- Calculating Sampling Coverage ---")

# 3.1. Map Sampling Points to Districts
sampling_sf <- host_data %>%
  drop_na(latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(mollweide_crs)

# Extract GID_2 for every sampling point
# (Using terra::extract is fast for points-in-polygons)
sampling_vect <- vect(sampling_sf)
extracted_adm <- extract(gadm_proj, sampling_vect)

sampling_points_tagged <- sampling_sf %>%
  mutate(GID_2 = extracted_adm$GID_2) %>%
  as_tibble()

# 3.2. Summarize Sampling by District & Genus
sampled_districts <- sampling_points_tagged %>%
  drop_na(GID_2, host_genus) %>%
  group_by(GID_2, host_genus) %>%
  summarise(n_records = n(), .groups = "drop")

# 3.3. Compare Range vs Sampled
# We iterate through genera to calculate the "Proportion of Range Sampled"

target_genera <- unique(host_traits$genus)

coverage_summary <- map_dfr(target_genera, function(gen) {
  
  # A. Total Range (All Districts this genus inhabits)
  inhabited_gids <- lookup_long %>% filter(genus == gen) %>% pull(GID_2)
  
  if(length(inhabited_gids) == 0) return(NULL)
  
  # Get area of these districts
  inhabited_polys <- gadm_proj[gadm_proj$GID_2 %in% inhabited_gids, ]
  range_area <- sum(inhabited_polys$area_km2, na.rm=TRUE)
  
  # B. Sampled Range (Districts with records)
  sampled_gids <- sampled_districts %>% filter(host_genus == gen) %>% pull(GID_2)
  
  # Overlap: Districts that are both inhabited AND sampled
  confirmed_gids <- intersect(inhabited_gids, sampled_gids)
  
  if(length(confirmed_gids) == 0) {
    sampled_area <- 0
  } else {
    sampled_polys <- gadm_proj[gadm_proj$GID_2 %in% confirmed_gids, ]
    sampled_area <- sum(sampled_polys$area_km2, na.rm=TRUE)
  }
  
  tibble(
    genus = gen,
    range_area_km2 = range_area,
    sampled_area_km2 = sampled_area,
    prop_sampled = sampled_area / range_area
  )
  
}, .progress = TRUE)

write_csv(coverage_summary, here("output", "analysis_1", "geographic_coverage_summary.csv"))


# --- 4. Visualizing the Bias (Interpolation Plot) ---

message("--- Generating Interpolated Bias Plot ---")

plot_data <- coverage_summary %>%
  filter(sampled_area_km2 > 0, range_area_km2 > 0) %>%
  # Label the top 10 most widespread genera
  mutate(label = if_else(rank(desc(range_area_km2)) <= 10, genus, ""))

# Interpolate (Akima)
interp_res <- akima::interp(
  x = log10(plot_data$sampled_area_km2),
  y = log10(plot_data$range_area_km2),
  z = plot_data$prop_sampled,
  linear = TRUE, nx = 200, ny = 200
)

interp_df <- as.data.frame(akima::interp2xyz(interp_res)) %>%
  as_tibble() %>%
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
  geom_text_repel(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2, label = label),
                  size = 3, max.overlaps = 20) +
  # Scales
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  labs(title = "Geographic Sampling Bias by Genus",
       subtitle = "Relationship between Range Size and Sampled Area",
       x = "Sampled Area (km²)", y = "Total Range Size (km²)") +
  theme_minimal()

ggsave(here("output", "analysis_1", "geographic_bias_interpolated.png"), p_bias, width = 10, height = 8)

message("✅ Geographic Bias Analysis Complete.")