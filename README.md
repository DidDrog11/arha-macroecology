# Project ArHa: Global Macroecology of Arenavirus and Hantavirus Reservoirs

## Project Overview
This repository contains the analysis pipeline for **Project ArHa**, a macroecological study aimed at identifying the intrinsic biological traits that define reservoir competence for Arenaviruses and Hantaviruses. 

Unlike traditional risk maps that often reflect sampling effort, this project explicitly models and corrects for surveillance bias. By integrating a harmonised global database of host-pathogen associations with phylogenetic mixed models, we decouple sampling effort from biological competence to identify silent hotspots—regions with high predicted reservoir diversity but low historical surveillance.

## Repository Structure

The project is organised into a reproducible pipeline where outputs from one stage feed into the next.

### Data Directories
* **`data/database/`**: Contains the primary *Project ArHa* relational database (Host, Pathogen, Sequence tables).
* **`data/external/`**: Raw inputs from external sources (IUCN shapefiles, COMBINE trait database, Ecke synanthropy data, WorldPop/VIIRS rasters).
* **`data/processed/`**: Intermediate files (e.g., harmonised taxonomy tables, matched phylogenies).
* **`data/analytic/`**: The final, imputed, and aligned datasets used for statistical modelling.

### Output Directories
* **`output/models/`**: Saved `brms` and `mgcv` model objects (`.rds`).
* **`output/figures/`**: Plots.
* **`output/tables/`**: Summary tables (e.g., Surveillance Anomalies, District Risk Scores).

---

## Analysis Pipeline

The analysis is split into numbered scripts that should be run in sequential order.

### Data Cleaning & Integration
**`01_data_integration.R`**
* **Purpose:** Harmonises inconsistent taxonomy across disconnected datasets.
* **Method:** Uses the GBIF taxonomic backbone (`taxize`) to resolve synonyms between the ArHa database, IUCN Red List ranges, and trait databases (COMBINE, Ecke et al.). It aligns all species to a consensus phylogeny.

**`02_imputation_pca.R`**
* **Purpose:** Prepares the analytic dataset by addressing missing life-history data.
* **Method:** Performs phylogenetic imputation (`Rphylopars`) on missing traits (e.g., litter size, mass) using the evolutionary covariance between species. It then runs a PCA to derive a single "Pace of Life" axis (PC1) to avoid collinearity in downstream models.

### Quantification of Bias
**`03_taxonomic_bias_analysis.R`**
* **Purpose:** Quantifies gaps in taxonomic coverage.
* **Method:** Compares the ArHa sampled species against the global checklist of Rodentia/Eulipotyphla. Uses Generalized Additive Models (GAMs) to determine if larger-bodied or wider-ranging species are disproportionately over-sampled. Calculates phylogenetic signal (Pagel's Lambda) in sampling effort.

**`04_01_geographic_bias_analysis.R`** & **`04_02_visualise_geographic_bias.R`**
* **Purpose:** Models the drivers of global surveillance intensity and identifies neglected regions.
* **Method:** Performs spatial intersection of IUCN ranges with GADM (ADM2) districts (Genus level). Fits a Zero-Inflated Negative Binomial (ZINB) model to test if surveillance counts are driven by wealth (Nightlights), accessibility, or biodiversity. Maps the residuals to highlight "Coldspots" (under-sampled) and "Hotspots" (over-sampled).

**`05_temporal_bias_analysis.R`**
* **Purpose:** Tracks the history of reservoir discovery.
* **Method:** Generates species accumulation curves and models continental trends in sampling effort over time (1960–2025) using GAMs.

**`06_genetic_bias_analysis.R`**
* **Purpose:** Assesses the "Genetic Gap" between detection and sequencing.
* **Method:** Produces bivariate maps comparing host vs. pathogen sequencing rates per country and plots quantifying the disconnect between PCR-positive assays and available GenBank sequences.

### 3. Macroecological Modelling 
**`07_macroecological_analysis.R`**
* **Purpose:** Tests the "Pace of Life" and "Synanthropy" hypotheses.
* **Method:** Fits Phylogenetic Generalized Linear Mixed Models (Bernoulli GLMMs) using `brms`. 
    * **Model Structure:** `Status ~ Effort + Traits + (1|Phylogeny) + (1|Species)`.
    * **Note:** Explicitly accounts for phylogenetic non-independence and varying sampling effort.

**`07_02_visualise_macroecological_analysis.R`**
* **Purpose:** Visualises model inference.
* **Method:** Generates posterior coefficient plots and marginal effects curves (e.g., probability of reservoir status vs. Pace of Life PC1).

**`07_03_post_hoc_investigations.R`**
* **Purpose:** Model diagnostics and validation.
* **Method:** Examines random effects to see if specific viral families are easier to detect. Maps model residuals to check for systematic geographic bias (e.g., latitudinal gradients in model failure).

**`07_04_spatial_validation.R`** (Risk Projection)
* **Purpose:** Generates the "Intrinsic Risk" map.
* **Method:** Projects the fitted model back onto the globe. Crucially, it predicts risk assuming **Saturated Sampling Effort** (setting effort to maximum), effectively asking: *"Where would we find reservoirs if we looked everywhere?"* Produces 20km resolution raster maps and district-level risk tables.

---

## Key Dependencies
* **Data Manipulation:** `tidyverse`, `janitor`, `here`
* **Spatial Analysis:** `terra`, `sf`, `tidyterra`, `rnaturalearth`
* **Phylogenetics:** `ape`, `phytools`, `Rphylopars`, `ggtree`
* **Modeling:** `brms` (Bayesian inference), `mgcv` (GAMs)