############################################################
# Publication figure for the simulation performance results
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

.analysis_script_dir <- local({
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      candidate <- frame$ofile
      if (is.character(candidate) && length(candidate) == 1L &&
          !is.na(candidate) && nzchar(candidate)) {
        return(candidate)
      }
      NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0L) {
    frame_file <- frame_files[[length(frame_files)]]
    return(dirname(normalizePath(
      frame_file,
      winslash = "/",
      mustWork = TRUE
    )))
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
})

# By default, read and write results beside this script. Set
# STTE_MDR_RESULT_DIR to use another directory without editing this file.
result_dir <- Sys.getenv(
  "STTE_MDR_RESULT_DIR",
  unset = .analysis_script_dir
)
result_dir <- normalizePath(
  result_dir,
  winslash = "/",
  mustWork = FALSE
)
result_file <- file.path(result_dir, "result500.rds")
if (!file.exists(result_file)) {
  stop(
    "Simulation result file not found: ",
    result_file,
    "\nPlace result500.rds beside this script or set STTE_MDR_RESULT_DIR ",
    "to its containing directory.",
    call. = FALSE
  )
}

simulation <- readRDS(result_file)
performance <- as.data.frame(simulation$performance_table)

method_labels <- c(
  ipcw_pool = "IPCW-Entry",
  dr_ipcw_pool = "DR-IPCW-Entry",
  mdr_pool_follow_balanced_cell_normalized = "MDR-Entry",
  dr_mdr_pool_follow_balanced_cell_normalized = "DR-MDR-Entry",
  ipcw_state_gcomp = "IPCW-State",
  dr_ipcw_state_gcomp = "DR-IPCW-State",
  mdr_state_gcomp_follow_balanced_cell_normalized = "MDR-State",
  dr_mdr_state_gcomp_follow_balanced_cell_normalized = "DR-MDR-State",
  ltmle_tmle_glm = "LTMLE-GLM",
  trial_emulation = "TrialEmulation"
)

method_order <- unname(method_labels)

scenario_labels <- c(
  base = "Base",
  small_size = "Small\nsample",
  strong_confounding = "Strong\nconfounding",
  poor_positivity = "Poor\npositivity",
  poor_adherence = "Poor\nadherence",
  rare_event = "Rare\nevent",
  nonlinear = "Nonlinear"
)

scenario_order <- unname(scenario_labels)

plot_data <- performance |>
  mutate(
    Method = factor(
      unname(method_labels[method]),
      levels = rev(method_order)
    ),
    Scenario = factor(
      unname(scenario_labels[scenario]),
      levels = scenario_order
    ),
    bias_pp = 100 * bias,
    rmse_pp = 100 * rmse
  )

write.csv(
  plot_data |>
    transmute(
      scenario,
      method,
      n_replications,
      n_success,
      failure_rate,
      psi_true,
      mean_psi_hat,
      bias,
      empirical_se,
      rmse
    ),
  file.path(result_dir, "table_simulation_performance.csv"),
  row.names = FALSE
)

best_bias <- plot_data |>
  group_by(Scenario) |>
  slice_min(abs(bias_pp), n = 1L, with_ties = FALSE) |>
  ungroup()

best_rmse <- plot_data |>
  group_by(Scenario) |>
  slice_min(rmse_pp, n = 1L, with_ties = FALSE) |>
  ungroup()

common_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, lineheight = 0.9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

bias_limit <- max(abs(plot_data$bias_pp), na.rm = TRUE)
bias_text_cutoff <- 0.56 * bias_limit

bias_plot <- ggplot(
  plot_data,
  aes(x = Scenario, y = Method, fill = bias_pp)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_tile(
    data = best_bias,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", bias_pp),
      color = abs(bias_pp) >= bias_text_cutoff
    ),
    size = 2.8
  ) +
  geom_hline(
    yintercept = c(2.5, 6.5),
    color = "grey30",
    linewidth = 0.45
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-bias_limit, bias_limit),
    name = "Percentage\npoints"
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "black"),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "A. Signed bias") +
  common_theme

rmse_limit <- max(plot_data$rmse_pp, na.rm = TRUE)
rmse_text_cutoff <- 0.58 * rmse_limit

rmse_plot <- ggplot(
  plot_data,
  aes(x = Scenario, y = Method, fill = rmse_pp)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_tile(
    data = best_rmse,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", rmse_pp),
      color = rmse_pp >= rmse_text_cutoff
    ),
    size = 2.8
  ) +
  geom_hline(
    yintercept = c(2.5, 6.5),
    color = "grey30",
    linewidth = 0.45
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#084594",
    limits = c(0, rmse_limit),
    name = "Percentage\npoints"
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "black"),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "B. Root mean squared error") +
  common_theme

performance_figure <- bias_plot / rmse_plot

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_performance_heatmap.pdf"
  ),
  plot = performance_figure,
  device = grDevices::cairo_pdf,
  width = 8.6,
  height = 8.2,
  units = "in"
)

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_performance_heatmap.png"
  ),
  plot = performance_figure,
  width = 8.6,
  height = 8.2,
  units = "in",
  dpi = 300
)

############################################################
# Appendix figure for weight stability
############################################################

# Entry, State, and doubly robust estimators within each family use the
# same underlying weights. Retain one representative per construction
# to avoid plotting duplicate values.
weight_method_labels <- c(
  ipcw_pool = "IPCW",
  mdr_pool_follow_balanced_cell_normalized = "MDR",
  trial_emulation = "TrialEmulation"
)

weight_order <- unname(weight_method_labels)

weight_data <- performance |>
  filter(method %in% names(weight_method_labels)) |>
  mutate(
    Weight = factor(
      unname(weight_method_labels[method]),
      levels = rev(weight_order)
    ),
    Scenario = factor(
      unname(scenario_labels[scenario]),
      levels = scenario_order
    )
  )

write.csv(
  weight_data |>
    transmute(
      scenario,
      weight_construction = as.character(Weight),
      mean_maximum_weight = mean_weight_max,
      mean_weight_cv = mean_weight_cv
    ),
  file.path(result_dir, "table_simulation_weight_diagnostics.csv"),
  row.names = FALSE
)

best_weight_max <- weight_data |>
  group_by(Scenario) |>
  slice_min(mean_weight_max, n = 1L, with_ties = FALSE) |>
  ungroup()

best_weight_cv <- weight_data |>
  group_by(Scenario) |>
  slice_min(mean_weight_cv, n = 1L, with_ties = FALSE) |>
  ungroup()

max_log_range <- range(log10(weight_data$mean_weight_max), na.rm = TRUE)
max_text_cutoff <- 10^(
  max_log_range[1] + 0.58 * diff(max_log_range)
)

weight_max_plot <- ggplot(
  weight_data,
  aes(x = Scenario, y = Weight, fill = mean_weight_max)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_tile(
    data = best_weight_max,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f", mean_weight_max),
      color = mean_weight_max >= max_text_cutoff
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#084594",
    trans = "log10",
    name = "Mean\nmaximum"
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "black"),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "A. Mean maximum weight") +
  common_theme +
  theme(axis.text.y = element_text(size = 9.5))

cv_limit <- max(weight_data$mean_weight_cv, na.rm = TRUE)
cv_text_cutoff <- 0.58 * cv_limit

weight_cv_plot <- ggplot(
  weight_data,
  aes(x = Scenario, y = Weight, fill = mean_weight_cv)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_tile(
    data = best_weight_cv,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", mean_weight_cv),
      color = mean_weight_cv >= cv_text_cutoff
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#FFF7EC",
    high = "#7F0000",
    limits = c(0, cv_limit),
    name = "Mean\nCV"
  ) +
  scale_color_manual(
    values = c(`TRUE` = "white", `FALSE` = "black"),
    guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "B. Mean coefficient of variation") +
  common_theme +
  theme(axis.text.y = element_text(size = 9.5))

weight_figure <- weight_max_plot / weight_cv_plot

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_diagnostics.pdf"
  ),
  plot = weight_figure,
  device = grDevices::cairo_pdf,
  width = 8.6,
  height = 4.4,
  units = "in"
)

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_diagnostics.png"
  ),
  plot = weight_figure,
  width = 8.6,
  height = 4.4,
  units = "in",
  dpi = 300
)

############################################################
# Appendix figure for the distribution of maximum weights
############################################################

# The simulation output retains replication-level weight summaries rather
# than the individual row-level weights. Show the Monte Carlo distribution
# of the replication-specific maximum for each unique weight construction.
weight_distribution_data <- as.data.frame(simulation$rep_results) |>
  filter(
    method %in% names(weight_method_labels),
    is.finite(weight_max),
    weight_max > 0
  ) |>
  mutate(
    Weight = factor(
      unname(weight_method_labels[method]),
      levels = weight_order
    ),
    Scenario = factor(
      unname(scenario_labels[scenario]),
      levels = scenario_order
    )
  )

write.csv(
  weight_distribution_data |>
    group_by(scenario, Weight) |>
    summarise(
      n_replications = n(),
      minimum = min(weight_max),
      first_quartile = quantile(weight_max, 0.25),
      median = median(weight_max),
      third_quartile = quantile(weight_max, 0.75),
      percentile_95 = quantile(weight_max, 0.95),
      maximum = max(weight_max),
      .groups = "drop"
    ) |>
    rename(weight_construction = Weight),
  file.path(result_dir, "table_simulation_weight_max_distribution.csv"),
  row.names = FALSE
)

weight_distribution_plot <- ggplot(
  weight_distribution_data,
  aes(x = Scenario, y = weight_max, fill = Weight)
) +
  geom_boxplot(
    width = 0.72,
    position = position_dodge(width = 0.78),
    linewidth = 0.38,
    outlier.alpha = 0.35,
    outlier.size = 0.65,
    outlier.stroke = 0.15
  ) +
  scale_fill_manual(
    values = c(
      IPCW = "#4C78A8",
      MDR = "#59A14F",
      TrialEmulation = "#E15759"
    ),
    drop = FALSE
  ) +
  scale_y_log10(
    breaks = c(3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000),
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    x = NULL,
    y = "Replication-specific maximum weight (log scale)",
    fill = "Weight construction"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8.8, lineheight = 0.9),
    axis.text.y = element_text(size = 8.8),
    axis.title.y = element_text(size = 9.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.32),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 8, 4, 4)
  )

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_max_distribution.pdf"
  ),
  plot = weight_distribution_plot,
  device = grDevices::cairo_pdf,
  width = 7.2,
  height = 4.4,
  units = "in"
)

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_max_distribution.png"
  ),
  plot = weight_distribution_plot,
  width = 7.2,
  height = 4.4,
  units = "in",
  dpi = 300
)

############################################################
# Appendix figure for the distribution of mean weights
############################################################

# Show the Monte Carlo distribution of each replication's mean analysis
# weight. MDR weights are normalized to have mean one within each
# strategy-follow-up cell and therefore collapse to one in this summary.
weight_mean_distribution_data <- as.data.frame(simulation$rep_results) |>
  filter(
    method %in% names(weight_method_labels),
    is.finite(weight_mean)
  ) |>
  mutate(
    Weight = factor(
      unname(weight_method_labels[method]),
      levels = weight_order
    ),
    Scenario = factor(
      unname(scenario_labels[scenario]),
      levels = scenario_order
    )
  )

write.csv(
  weight_mean_distribution_data |>
    group_by(scenario, Weight) |>
    summarise(
      n_replications = n(),
      minimum = min(weight_mean),
      first_quartile = quantile(weight_mean, 0.25),
      median = median(weight_mean),
      monte_carlo_mean = mean(weight_mean),
      third_quartile = quantile(weight_mean, 0.75),
      percentile_95 = quantile(weight_mean, 0.95),
      maximum = max(weight_mean),
      .groups = "drop"
    ) |>
    rename(weight_construction = Weight),
  file.path(result_dir, "table_simulation_weight_mean_distribution.csv"),
  row.names = FALSE
)

weight_mean_distribution_plot <- ggplot(
  weight_mean_distribution_data,
  aes(x = Scenario, y = weight_mean, fill = Weight)
) +
  geom_hline(
    yintercept = 1,
    color = "#666666",
    linewidth = 0.42,
    linetype = "dashed"
  ) +
  geom_boxplot(
    width = 0.72,
    position = position_dodge(width = 0.78),
    linewidth = 0.38,
    outlier.alpha = 0.40,
    outlier.size = 0.70,
    outlier.stroke = 0.15
  ) +
  scale_fill_manual(
    values = c(
      IPCW = "#4C78A8",
      MDR = "#59A14F",
      TrialEmulation = "#E15759"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = c(0.9, 1, 1.25, 1.5, 2, 2.5, 3),
    labels = scales::label_number(accuracy = 0.01),
    limits = c(0.84, 3.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Replication-specific mean weight",
    fill = "Weight construction"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8.8, lineheight = 0.9),
    axis.text.y = element_text(size = 8.8),
    axis.title.y = element_text(size = 9.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.32),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 8, 4, 4)
  )

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_mean_distribution.pdf"
  ),
  plot = weight_mean_distribution_plot,
  device = grDevices::cairo_pdf,
  width = 7.2,
  height = 4.4,
  units = "in"
)

ggsave(
  filename = file.path(
    result_dir,
    "figure_simulation_weight_mean_distribution.png"
  ),
  plot = weight_mean_distribution_plot,
  width = 7.2,
  height = 4.4,
  units = "in",
  dpi = 300
)
