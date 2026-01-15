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
  
  # Setup Data & Conditions
  # Find Maximum Sampling Effort
  max_effort_val <- max(fit_model$data$log_effort, na.rm = TRUE)
  
  # Create a condition list where effort is fixed at max
  # We use this so the signal isn't hidden by low effort
  cond_max_effort <- data.frame(log_effort = max_effort_val)
  
  # Panel A: Coefficients (Posterior Distributions)
  post_draws <- fit_model |>
    gather_draws(`b_.*`, regex = TRUE) |>
    filter(.variable != "b_Intercept") |>
    mutate(variable_clean = case_when(
      .variable == "b_log_effort" ~ "Sampling Effort",
      .variable == "b_pace_of_life_pc1" ~ "Pace of Life (Slow)",
      .variable == "b_synanthropy_statusOccasionallySynanthropic" ~ "Synanthropy (Occasional)",
      .variable == "b_synanthropy_statusTotallySynanthropic" ~ "Synanthropy (Total)",
      TRUE ~ .variable
    )) |>
    mutate(variable_clean = fct_reorder(variable_clean, .value))
  
  # Calculate Probability of Direction (pd)
  pd_stats <- post_draws |>
    group_by(variable_clean) |>
    summarise(prob_dir = max(mean(.value > 0), mean(.value < 0)),
              median_val = median(.value),
              label = paste0(round(prob_dir * 100, 1), "%"))
  
  p_coefs <- ggplot(post_draws, aes(y = variable_clean, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    stat_halfeye(fill = "#56B4E9", alpha = 0.8, .width = c(0.5, 0.95), normalize = "groups") +
    geom_text(data = pd_stats, aes(x = median_val, label = label), 
              nudge_y = 0.5, size = 3, fontface = "italic", color = "grey30") +
    labs(title = "a) Predictors of Reservoir Status",
         subtitle = "Posterior log-odds distributions",
         x = "Effect Size (Log-Odds)", y = NULL) +
    theme_minimal_vgrid() +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))
  
  # Panel B: Surveillance Bias
  # Here we don't fix effort
  cond_effort <- conditional_effects(fit_model, effects = "log_effort")
  eff_data <- cond_effort$log_effort
  
  p_effort <- ggplot(eff_data, aes(x = log_effort, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "grey20") +
    geom_line(color = "black", linewidth = 1) +
    labs(title = "b) Surveillance Bias",
         y = "Probability of Detection",
         x = "Individuals Sampled (Log Scale)") +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000) + 1),
                       labels = c("1", "10", "100", "1k", "10k")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # Panel C: Intrinsic Biology (fixed at max effort)
  grid_pace <- fit_model$data |>
    distinct(pace_of_life_pc1) |>
    mutate(log_effort = max(fit_model$data$log_effort, na.rm = TRUE),      
           # Set Synanthropy to reference level (if it exists) or remove it from grid if not needed for this plot
           synanthropy_status = factor("NonSynanthropic", levels = levels(fit_model$data$synanthropy_status)))
  
  if(!"synanthropy_status" %in% names(fit_model$data)) {
    grid_pace <- select(grid_pace, -synanthropy_status)
  }
  
  # Generate Epred Draws
  # This gives us the full distribution for every point on the x-axis
  draws_pace <- fit_model |>
    add_epred_draws(newdata = grid_pace, re_formula = NA, ndraws = 1000)
  
  max_zoom <- quantile(draws_pace$.epred, 0.99)
  
  p_pace <- ggplot(data = draws_pace, aes(x = pace_of_life_pc1, y = .epred)) +
    stat_lineribbon(.width = c(.50, .80, .95), alpha = 0.25, fill = "#009E73") +
    stat_summary(fun = median, geom = "line", colour = "#007858", linewidth = 1) +
    labs(title = "c) Intrinsic Biology (at maximum sampling effort)",
         y = "Probability of Reservoir Status",
         x = "Pace of Life (PC1)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    annotate("text", x = min(draws_pace$pace_of_life_pc1), y = Inf, label = "Fast", 
             hjust = 0, vjust = 1.5, colour = "grey50", fontface = "italic") +
    annotate("text", x = max(draws_pace$pace_of_life_pc1), y = Inf, label = "Slow", 
             hjust = 1, vjust = 1.5, colour = "grey50", fontface = "italic") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
  
  # 5. Panel D: Synanthropy (fixed at max effort)
  if ("synanthropy_status" %in% names(fit_model$data)) {
    grid_syn <- expand.grid(synanthropy_status = unique(fit_model$data$synanthropy_status),
                            log_effort = max(fit_model$data$log_effort, na.rm = TRUE),
                            pace_of_life_pc1 = 0)
    draws_syn <- fit_model |>
      add_epred_draws(newdata = grid_syn, re_formula = NA, ndraws = 1000)
    
    p_syn <- ggplot(draws_syn, aes(x = synanthropy_status, y = .epred)) +
      stat_eye(fill = "#D55E00", alpha = 0.6, .width = c(0.5, 0.95), point_colour = "black") +
      labs(title = "d) Ecological Opportunity (at maximum sampling effort)",
           y = "Probability of Reservoir Status",
           x = NULL) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 15, hjust = 1))
    
    # Assemble 4-Panel
    final_plot <- (p_coefs + p_effort) /
      (p_pace + p_syn) + 
      plot_layout(heights = c(1.2, 1))
    
  } else {
    # Assemble 3-Panel
    final_plot <- p_coefs / (p_effort + p_pace)
  }
  
  # Save
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

