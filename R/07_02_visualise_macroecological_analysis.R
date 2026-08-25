# ==============================================================================
# 07_02_visualise_macroecological_analysis.R
# Purpose: Define a flexible function to visualize brms reservoir models.
#          Generates Figure 3 (Coefficients & Marginal Effects).
# ==============================================================================

# 1. Setup ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, brms, tidybayes, patchwork, ggdist, cowplot)

# Output Directory
output_dir <- here("output", "figures")
dir.create(output_dir, showWarnings = FALSE)

# Visualisation Function --------------------------------------------------
visualise_reservoir_model <- function(fit_model, plot_subtitle, filename) {
  
  # Find Maximum Sampling Effort
  max_effort_val <- quantile(fit_model$data$log_effort, 0.95, na.rm = TRUE)
  
  # Panel A: Coefficients (Posterior Distributions)
  post_draws <- fit_model |>
    gather_draws(`b_.*`, regex = TRUE) |>
    filter(.variable != "b_Intercept") |>
    mutate(variable_clean = case_when(.variable == "b_log_effort" ~ "Sampling Effort",
                                      .variable == "b_pace_of_life_pc1" ~ "Pace of Life (Slow)",
                                      .variable == "b_synanthropy_statusOccasionallySynanthropic" ~ "Synanthropy (Occasional)",
                                      .variable == "b_synanthropy_statusTotallySynanthropic" ~ "Synanthropy (Total)",
                                      TRUE ~ .variable)) |>
    mutate(variable_clean = fct_reorder(variable_clean, .value))
  
  pd_stats <- post_draws |>
    group_by(variable_clean) |>
    summarise(prob_dir = max(mean(.value > 0), mean(.value < 0)),
              median_val = median(.value),
              label = paste0(round(prob_dir * 100, 1), "%"))
  
  p_coefs <- ggplot(post_draws, aes(y = variable_clean, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    stat_halfeye(fill = "#56B4E9", alpha = 0.8, .width = c(0.5, 0.95), normalize = "groups") +
    geom_text(data = pd_stats, aes(x = median_val, label = label), nudge_y = 0.5, size = 3, fontface = "italic", colour = "grey30") +
    scale_y_discrete(labels = scales::label_wrap(12)) +
    labs(title = "a) Predictors of Host Status",
         subtitle = "Posterior log-odds distributions",
         x = "Effect Size (Log-Odds)", y = NULL) +
    theme_minimal_vgrid() +
    theme(plot.title = element_text(face = "bold", margin = margin(b = 10)),
          axis.text.y = element_text(face = "bold"))
  
  # Panel B: Surveillance Bias
  cond_effort <- conditional_effects(fit_model, effects = "log_effort")
  eff_data <- cond_effort$log_effort
  
  p_effort <- ggplot(eff_data, aes(x = log_effort, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "grey20") +
    geom_line(colour = "black", linewidth = 1) +
    labs(title = "b) Surveillance Bias",
         y = "Probability of Host Status",
         x = "Individuals Sampled (Log Scale)") +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000) + 1),
                       labels = c("1", "10", "100", "1k", "10k")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", margin = margin(b = 10)))
  
  # Panel C: Intrinsic Biology
  grid_pace <- fit_model$data |>
    distinct(pace_of_life_pc1) |>
    mutate(log_effort = max_effort_val)
  
  if ("synanthropy_status" %in% names(fit_model$data)) {
    grid_pace <- grid_pace |> 
      mutate(synanthropy_status = factor("Not Synanthropic", levels = levels(fit_model$data$synanthropy_status)))
  }
  
  draws_pace <- fit_model |>
    add_epred_draws(newdata = grid_pace, re_formula = NA, ndraws = 1000)
  
  p_pace <- ggplot(data = draws_pace, aes(x = pace_of_life_pc1, y = .epred)) +
    stat_lineribbon(.width = c(.50, .80, .95), alpha = 0.25, fill = "#009E73") +
    stat_summary(fun = median, geom = "line", colour = "#007858", linewidth = 1) +
    labs(title = "c) Intrinsic Biology",
         y = "Probability of Host Status",
         x = "Pace of Life (PC1)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    annotate("text", x = min(draws_pace$pace_of_life_pc1), y = Inf, label = "Fast", 
             hjust = 0, vjust = 1.5, colour = "grey50", fontface = "italic") +
    annotate("text", x = max(draws_pace$pace_of_life_pc1), y = Inf, label = "Slow", 
             hjust = 1, vjust = 1.5, colour = "grey50", fontface = "italic") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", margin = margin(b = 10)),
          legend.position = "none")
  
  # Panel D & Assembly
  if ("synanthropy_status" %in% names(fit_model$data)) {
    grid_syn <- expand.grid(synanthropy_status = unique(fit_model$data$synanthropy_status),
                            log_effort = max_effort_val,
                            pace_of_life_pc1 = 0)
    draws_syn <- fit_model |>
      add_epred_draws(newdata = grid_syn, re_formula = NA, ndraws = 1000)
    
    p_syn <- ggplot(draws_syn, aes(x = synanthropy_status, y = .epred)) +
      stat_eye(fill = "#D55E00", alpha = 0.6, .width = c(0.5, 0.95), point_colour = "black") +
      labs(title = "d) Ecological Opportunity",
           y = "Probability of Host Status",
           x = NULL) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                         expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(ylim = c(0, 0.4)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", margin = margin(b = 10)),
            axis.text.x = element_text(angle = 15, hjust = 1))
    
    # Assembly
    final_plot <- (p_coefs + p_effort) / (p_pace + p_syn) + 
      plot_layout(heights = c(1.2, 1)) & 
      theme(plot.margin = margin(15, 15, 15, 15))
    
  } else {
    # Assembly
    final_plot <- p_coefs / (p_effort | p_pace) +
      plot_layout(heights = c(1.2, 1)) & 
      theme(plot.margin = margin(15, 15, 15, 15))
  }
  
  ggsave(filename, final_plot, width = 10, height = 12, bg = "white")
  return(final_plot)
}

# 3. Execute Visualisation ------------------------------------------------
if (file.exists(here("output", "models", "brms_dyadic_N19k.rds"))) {
  fit_mechanism <- read_rds(here("output", "models", "brms_dyadic_N19k.rds"))
  
  visualise_reservoir_model(
    fit_model = fit_mechanism,
    plot_subtitle = "Subset with Synanthropy Data (N = 19,180)",
    filename = here(output_dir, "dyadic_glmm_synanthropy.png")
  )
}

if (file.exists(here("output", "models", "brms_dyadic_N49k.rds"))) {
  fit_global <- read_rds(here("output", "models", "brms_dyadic_N49k.rds"))

  visualise_reservoir_model(
    fit_model = fit_global,
    plot_subtitle = "Global Test: Full Dataset (N = 49,280)",
    filename = here(output_dir, "dyadic_glmm_global.png")
  )
}

if (file.exists(here("output", "models", "brms_dyadic_sens_strict.rds"))) {
  fit_sens <- read_rds(here("output", "models", "brms_dyadic_sens_strict.rds"))
  
  visualise_reservoir_model(
    fit_model = fit_sens,
    plot_subtitle = "Sensitivity Test: Molecular & Isolation Only",
    filename = here(output_dir, "dyadic_glmm_sensitivity_strict.png")
  )
}

if (file.exists(here("output", "models", "brms_dyadic_sens_raw.rds"))) {
  fit_raw <- read_rds(here("output", "models", "brms_dyadic_sens_raw.rds"))
}

# --- Generate Posterior Comparison Plot (Full vs Strict vs Raw) ---
if (exists("fit_global") && exists("fit_sens") && exists("fit_raw")) {
  
  draws_full <- fit_global |>
    gather_draws(b_pace_of_life_pc1) |>
    mutate(model_type = "Full Model (Imputed traits, all evidence types)")
  
  draws_strict <- fit_sens |>
    gather_draws(b_pace_of_life_pc1) |>
    mutate(model_type = "Strict Model (Imputed traits, active infection only)")
  
  draws_raw <- fit_raw |>
    gather_draws(b_pace_of_life_pc1_raw) |>
    mutate(model_type = "Complete Cases (Empirical traits, all evidence types)")
  
  compare_draws <- bind_rows(draws_full, draws_strict, draws_raw) |>
    mutate(model_type = fct_relevel(model_type, 
                                    "Complete Cases (Empirical traits, all evidence types)", 
                                    "Strict Model (Imputed traits, active infection only)", 
                                    "Full Model (Imputed traits, all evidence types)"))
  
  pd_labels <- compare_draws |>
    group_by(model_type) |>
    summarise(pd = max(mean(.value > 0), mean(.value < 0)),
              median_val = median(.value),
              label = paste0("pd = ", round(pd * 100, 1), "%"))
  
  p_compare <- ggplot(compare_draws, aes(x = .value, y = model_type, fill = model_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1) +
    stat_halfeye(alpha = 0.7, .width = c(0.5, 0.95), point_colour = "black") +
    geom_text(data = pd_labels, aes(x = median_val, label = label), 
              nudge_y = 0.3, fontface = "italic", colour = "grey20") +
    scale_fill_manual(values = c("Full Model (Imputed traits, all evidence types)" = "#56B4E9", 
                                 "Strict Model (Imputed traits, active infection only)" = "#D55E00",
                                 "Complete Cases (Empirical traits, all evidence types)" = "#009E73")) +
    labs(title = "Sensitivity Analysis: Pace of Life Effect",
         x = "Effect Size (Log-Odds)",
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold", size = 10))
    
    ggsave(here(output_dir, "supplementary_figure_s6.png"), p_compare, width = 8, height = 5, bg = "white")
}


# Review Response ---------------------------------------------------------
# Synanthropy vs Pace of Life
if (exists("fit_mechanism")) {
  p_syn_pace <- ggplot(fit_mechanism$data, 
                       aes(x = synanthropy_status, y = pace_of_life_pc1, fill = synanthropy_status)) +
    geom_violin(alpha = 0.7, scale = "count", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.8, outlier.shape = NA) +
    labs(title = "Pace of life by synanthropy status",
         x = NULL, y = "Pace of Life (PC1)") +
    theme_minimal() + theme(legend.position = "none")
  
  ggsave(here(output_dir, "revision_fig_synanthropy_pace.png"), p_syn_pace, width = 6, height = 5, bg = "white")
}

# Visualise host-vs-phylogeny
if (exists("fit_global")) {
  var_draws <- fit_global |>
    gather_draws(sd_host_id__Intercept, sd_tip_label__Intercept, sd_pathogen_species_cleaned__Intercept) |>
    mutate(component = case_when(
      str_detect(.variable, "host_id") ~ "Host Identity",
      str_detect(.variable, "tip_label") ~ "Phylogeny",
      str_detect(.variable, "pathogen") ~ "Viral Species",
      TRUE ~ .variable
    ))
  
  p_variance <- ggplot(var_draws, aes(x = .value, y = fct_reorder(component, .value))) +
    stat_halfeye(fill = "#CC79A7", alpha = 0.7, .width = c(0.5, 0.95)) +
    labs(title = "Variance partitioning: Global model",
         subtitle = "Random effect standard deviations",
         x = "SD (log-odds scale)", y = NULL) +
    theme_minimal()
  
  ggsave(here(output_dir, "revision_fig_variance_partitioning.png"), p_variance, width = 7, height = 4, bg = "white")
}
