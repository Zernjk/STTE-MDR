############################################################
# Scenario diagnostics for 500 observed-data replications
# per scenario, loaded from one self-contained RDS file.
############################################################

options(stringsAsFactors = FALSE)

.script_dir <- local({
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

required_packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Reuse the exact DGP definitions that generated the observed data.
source(file.path(.script_dir, "stte_per_protocol_simulation.R"))

scenario_order <- c(
  "base", "small_size", "strong_confounding", "poor_positivity",
  "poor_adherence", "rare_event", "nonlinear"
)
scenario_labels <- c(
  base = "Base",
  small_size = "Small sample",
  strong_confounding = "Strong confounding",
  poor_positivity = "Poor positivity",
  poor_adherence = "Poor adherence",
  rare_event = "Rare event",
  nonlinear = "Nonlinear"
)

observed_data_file <- file.path(.script_dir, "observed_data500.rds")
if (!file.exists(observed_data_file)) {
  stop(
    "Required self-contained observed-data file not found: ",
    observed_data_file,
    call. = FALSE
  )
}

observed_data500 <- readRDS(observed_data_file)
required_components <- c("metadata", "manifest", "observed_data")
missing_components <- setdiff(required_components, names(observed_data500))
if (length(missing_components) > 0L) {
  stop(
    "observed_data500.rds is missing components: ",
    paste(missing_components, collapse = ", "),
    call. = FALSE
  )
}

inventory <- observed_data500$manifest
inventory$scenario <- factor(inventory$scenario, levels = scenario_order)
inventory <- inventory %>%
  arrange(scenario, global_replication) %>%
  mutate(scenario = as.character(scenario))

if (nrow(inventory) != 3500L) {
  stop("Expected 3,500 datasets in the self-contained manifest.", call. = FALSE)
}
scenario_count_check <- inventory %>% count(scenario, name = "n_replications")
if (nrow(scenario_count_check) != 7L ||
    any(scenario_count_check$n_replications != 500L)) {
  stop("Expected exactly 500 replications in every scenario.", call. = FALSE)
}
if (!identical(names(observed_data500$observed_data), scenario_order)) {
  stop(
    "The observed_data component must contain the seven scenarios in paper order.",
    call. = FALSE
  )
}
for (scenario_name in scenario_order) {
  scenario_data <- observed_data500$observed_data[[scenario_name]]
  scenario_manifest <- inventory[inventory$scenario == scenario_name, , drop = FALSE]
  if (!is.list(scenario_data) || length(scenario_data) != 500L) {
    stop(
      "Expected 500 observed datasets for scenario '", scenario_name, "'.",
      call. = FALSE
    )
  }
  if (!all(scenario_manifest$replication_id %in% names(scenario_data))) {
    stop(
      "The manifest and observed datasets do not agree for scenario '",
      scenario_name,
      "'.",
      call. = FALSE
    )
  }
}

output_dir <- file.path(.script_dir, "scenario_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inventory_check <- inventory %>%
  count(
    batch_number, batch, generation_seed, scenario,
    name = "n_replications"
  ) %>%
  arrange(factor(scenario, levels = scenario_order), batch_number)
write.csv(
  inventory_check,
  file.path(output_dir, "data_inventory.csv"),
  row.names = FALSE
)

discrete_event_curve <- function(dat, max_month = 24L) {
  survival_probability <- 1
  result <- vector("list", max_month)
  for (month in seq_len(max_month)) {
    month_data <- dat[dat$k == month - 1L, , drop = FALSE]
    n_at_risk <- nrow(month_data)
    if (n_at_risk == 0L || is.na(survival_probability)) {
      event_hazard <- NA_real_
      survival_probability <- NA_real_
    } else {
      event_hazard <- sum(month_data$dY_next, na.rm = TRUE) / n_at_risk
      survival_probability <- survival_probability * (1 - event_hazard)
    }
    result[[month]] <- data.frame(
      month = month,
      n_at_risk = n_at_risk,
      event_hazard = event_hazard,
      cumulative_event_risk_pct = 100 * (1 - survival_probability)
    )
  }
  bind_rows(result)
}

lp_A <- function(A_prev, L_k, X1, X2, X3, k, params) {
  value <- params$alpha0 + params$alpha_Aprev * A_prev +
    params$alpha_L * L_k + params$alpha_x1 * X1 +
    params$alpha_x2 * X2 + params$alpha_x3 * X3 +
    params$alpha_time * k
  if (isTRUE(params$nonlinear)) {
    value <- value + params$alpha_L2 * L_k^2
  }
  value
}

lp_Y <- function(A_k, L_k, X1, X2, X3, k, params) {
  value <- params$beta0 + params$beta_A * A_k +
    params$beta_L * L_k + params$beta_x1 * X1 +
    params$beta_x2 * X2 + params$beta_x3 * X3 +
    params$beta_time * k
  if (isTRUE(params$nonlinear)) {
    value <- value + params$beta_L2 * L_k^2 +
      params$beta_A_L * A_k * L_k
  }
  value
}

analyze_replication <- function(inventory_row, dat) {
  dataset_label <- paste0(
    inventory_row$scenario,
    "/",
    inventory_row$replication_id
  )
  expected_columns <- c(
    "id", "k", "X1", "X2", "X3", "L_k", "A_prev", "A_k",
    "Y_start", "C_start", "dY_next", "dC_next", "Y_next", "C_next"
  )
  missing_columns <- setdiff(expected_columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns in ", dataset_label, ": ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  params <- define_parameters(inventory_row$scenario)
  pA <- compute_pA(
    dat$A_prev, dat$L_k, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  pY <- compute_pY(
    dat$A_k, dat$L_k, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  L_quartiles <- as.numeric(quantile(
    dat$L_k,
    c(0.25, 0.75),
    names = FALSE,
    na.rm = TRUE
  ))
  L_low <- L_quartiles[1]
  L_high <- L_quartiles[2]

  lpA_low <- lp_A(
    dat$A_prev, L_low, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  lpA_high <- lp_A(
    dat$A_prev, L_high, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  lpY_low <- lp_Y(
    dat$A_k, L_low, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  lpY_high <- lp_Y(
    dat$A_k, L_high, dat$X1, dat$X2, dat$X3, dat$k, params
  )

  pA_low <- compute_pA(
    dat$A_prev, L_low, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  pA_high <- compute_pA(
    dat$A_prev, L_high, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  pY_low <- compute_pY(
    dat$A_k, L_low, dat$X1, dat$X2, dat$X3, dat$k, params
  )
  pY_high <- compute_pY(
    dat$A_k, L_high, dat$X1, dat$X2, dat$X3, dat$k, params
  )

  curve <- discrete_event_curve(dat, params$K) %>%
    mutate(
      batch_number = inventory_row$batch_number,
      batch = inventory_row$batch,
      generation_seed = inventory_row$generation_seed,
      scenario = inventory_row$scenario,
      within_batch_replication = inventory_row$within_batch_replication,
      global_replication = inventory_row$global_replication,
      .before = 1
    )

  untreated <- dat$A_prev == 0
  treated <- dat$A_prev == 1
  risk_12 <- curve$cumulative_event_risk_pct[curve$month == 12L]
  risk_24 <- curve$cumulative_event_risk_pct[curve$month == 24L]

  metrics <- data.frame(
    batch_number = inventory_row$batch_number,
    batch = inventory_row$batch,
    generation_seed = inventory_row$generation_seed,
    scenario = inventory_row$scenario,
    within_batch_replication = inventory_row$within_batch_replication,
    global_replication = inventory_row$global_replication,
    observed_seed = inventory_row$observed_seed,
    n_participants = length(unique(dat$id)),
    person_months = nrow(dat),
    event_count_12 = sum(dat$dY_next[dat$k < 12L], na.rm = TRUE),
    event_count_24 = sum(dat$dY_next, na.rm = TRUE),
    censor_count_24 = sum(dat$dC_next, na.rm = TRUE),
    event_rate_per_1000_pm = 1000 * mean(dat$dY_next, na.rm = TRUE),
    event_risk_12_pct = risk_12,
    event_risk_24_pct = risk_24,
    initiation_pct = 100 * mean(dat$A_k[untreated] == 1, na.rm = TRUE),
    continuation_pct = 100 * mean(dat$A_k[treated] == 1, na.rm = TRUE),
    discontinuation_pct = 100 * mean(dat$A_k[treated] == 0, na.rm = TRUE),
    confounder_treatment_or_iqr = exp(mean(lpA_high - lpA_low)),
    confounder_event_or_iqr = exp(mean(lpY_high - lpY_low)),
    confounder_treatment_rd_pctpt = 100 * mean(pA_high - pA_low),
    confounder_event_rd_pctpt = 100 * mean(pY_high - pY_low),
    pA_mean = mean(pA),
    pA_q01 = as.numeric(quantile(pA, 0.01, names = FALSE)),
    pA_q05 = as.numeric(quantile(pA, 0.05, names = FALSE)),
    pA_q50 = as.numeric(quantile(pA, 0.50, names = FALSE)),
    pA_q95 = as.numeric(quantile(pA, 0.95, names = FALSE)),
    pA_q99 = as.numeric(quantile(pA, 0.99, names = FALSE)),
    pA_extreme_pct = 100 * mean(pA <= 0.05 | pA >= 0.95),
    pA_severe_pct = 100 * mean(pA <= 0.01 + 1e-12 | pA >= 0.99 - 1e-12),
    mean_overlap_probability = mean(pmin(pA, 1 - pA)),
    pY_mean = mean(pY),
    L_q25 = L_low,
    L_q75 = L_high
  )

  pA_hist <- hist(
    pA,
    breaks = seq(0, 1, by = 0.02),
    include.lowest = TRUE,
    right = TRUE,
    plot = FALSE
  )
  positivity_histogram <- data.frame(
    batch_number = inventory_row$batch_number,
    batch = inventory_row$batch,
    scenario = inventory_row$scenario,
    global_replication = inventory_row$global_replication,
    probability_midpoint = pA_hist$mids,
    percent_person_months = 100 * pA_hist$counts / sum(pA_hist$counts)
  )

  list(
    metrics = metrics,
    event_curve = curve,
    positivity_histogram = positivity_histogram
  )
}

all_results <- vector("list", nrow(inventory))
for (i in seq_len(nrow(inventory))) {
  if (i == 1L || i %% 100L == 0L || i == nrow(inventory)) {
    message("Analyzing dataset ", i, " of ", nrow(inventory), " ...")
  }
  inventory_row <- inventory[i, , drop = FALSE]
  dat <- observed_data500$observed_data[[inventory_row$scenario]][[
    inventory_row$replication_id
  ]]
  all_results[[i]] <- analyze_replication(inventory_row, dat)
}

replication_diagnostics <- bind_rows(lapply(all_results, `[[`, "metrics")) %>%
  mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  arrange(scenario, global_replication)
event_curves <- bind_rows(lapply(all_results, `[[`, "event_curve")) %>%
  mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  arrange(scenario, global_replication, month)
positivity_histograms <- bind_rows(
  lapply(all_results, `[[`, "positivity_histogram")
) %>%
  mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  arrange(scenario, global_replication, probability_midpoint)

write.csv(
  replication_diagnostics,
  file.path(output_dir, "replication_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  event_curves,
  file.path(output_dir, "replication_event_curves.csv"),
  row.names = FALSE
)
write.csv(
  positivity_histograms,
  file.path(output_dir, "replication_positivity_histograms.csv"),
  row.names = FALSE
)
saveRDS(
  replication_diagnostics,
  file.path(output_dir, "replication_diagnostics_500_replications.rds")
)

id_columns <- c(
  "batch_number", "batch", "generation_seed", "scenario",
  "within_batch_replication", "global_replication", "observed_seed"
)
metric_columns <- setdiff(names(replication_diagnostics), id_columns)
replication_long <- replication_diagnostics %>%
  pivot_longer(
    cols = all_of(metric_columns),
    names_to = "metric",
    values_to = "value"
  )

summary_diagnostics_long <- replication_long %>%
  group_by(scenario, metric) %>%
  summarise(
    n_replications = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q025 = as.numeric(quantile(value, 0.025, na.rm = TRUE)),
    q975 = as.numeric(quantile(value, 0.975, na.rm = TRUE)),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, metric)

batch_comparison <- replication_long %>%
  group_by(scenario, metric, batch_number, batch) %>%
  summarise(
    n_replications = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, metric, batch_number)

write.csv(
  summary_diagnostics_long,
  file.path(output_dir, "summary_diagnostics_long.csv"),
  row.names = FALSE
)
write.csv(
  batch_comparison,
  file.path(output_dir, "batch_comparison.csv"),
  row.names = FALSE
)

dgp_table <- bind_rows(lapply(scenario_order, function(scenario_name) {
  params <- define_parameters(scenario_name)
  data.frame(
    scenario = scenario_name,
    scenario_label = unname(scenario_labels[scenario_name]),
    requested_or_overridden_n = if (is.null(params$observed_sample_size)) {
      2000L
    } else {
      params$observed_sample_size
    },
    alpha_L_treatment = params$alpha_L,
    beta_L_event = params$beta_L,
    alpha_Aprev_persistence = params$alpha_Aprev,
    beta0_event_intercept = params$beta0,
    nonlinear = params$nonlinear,
    alpha_L2_treatment = params$alpha_L2,
    beta_L2_event = params$beta_L2,
    beta_A_by_L_event = params$beta_A_L,
    treatment_probability_clip = .dgp_prob_clip(params, "treatment"),
    event_probability_clip = .dgp_prob_clip(params, "outcome")
  )
}))
write.csv(
  dgp_table,
  file.path(output_dir, "table_scenario_dgp_parameters.csv"),
  row.names = FALSE
)

get_metric_summary <- function(metric_name) {
  summary_diagnostics_long %>%
    filter(metric == metric_name) %>%
    select(scenario, median, q025, q975, min, max)
}

formatted_table <- data.frame(
  scenario = factor(scenario_order, levels = scenario_order),
  scenario_label = unname(scenario_labels[scenario_order])
)
formatted_specs <- list(
  n_participants = list("Participants, median [min, max]", 0L, TRUE),
  confounder_treatment_or_iqr = list("L-treatment OR (Q3 vs Q1)", 2L, FALSE),
  confounder_event_or_iqr = list("L-event OR (Q3 vs Q1)", 2L, FALSE),
  pA_extreme_pct = list("Practical positivity violations, %", 1L, FALSE),
  initiation_pct = list("Initiation when previously untreated, %", 1L, FALSE),
  continuation_pct = list("Continuation when previously treated, %", 1L, FALSE),
  event_count_12 = list("Events by month 12, count", 0L, FALSE),
  event_risk_12_pct = list("12-month cumulative event risk, %", 1L, FALSE)
)

for (metric_name in names(formatted_specs)) {
  specification <- formatted_specs[[metric_name]]
  dat <- get_metric_summary(metric_name)
  number_format <- paste0("%.", specification[[2]], "f")
  lower <- if (isTRUE(specification[[3]])) dat$min else dat$q025
  upper <- if (isTRUE(specification[[3]])) dat$max else dat$q975
  dat$value <- paste0(
    sprintf(number_format, dat$median), " [",
    sprintf(number_format, lower), ", ",
    sprintf(number_format, upper), "]"
  )
  names(dat)[names(dat) == "value"] <- specification[[1]]
  formatted_table <- left_join(
    formatted_table,
    dat[c("scenario", specification[[1]])],
    by = "scenario"
  )
}
formatted_table$scenario <- as.character(formatted_table$scenario)
write.csv(
  formatted_table,
  file.path(output_dir, "table_scenario_diagnostics_formatted.csv"),
  row.names = FALSE
)

# Publication figures ---------------------------------------------------------
paper_theme <- theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    plot.tag = element_text(face = "bold"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black")
  )

plot_summary <- summary_diagnostics_long %>%
  mutate(scenario_label = factor(
    unname(scenario_labels[as.character(scenario)]),
    levels = rev(unname(scenario_labels[scenario_order]))
  ))

confounding_plot_data <- plot_summary %>%
  filter(metric %in% c(
    "confounder_treatment_or_iqr", "confounder_event_or_iqr"
  )) %>%
  mutate(pathway = ifelse(
    metric == "confounder_treatment_or_iqr",
    "L -> treatment",
    "L -> event"
  ))

p_confounding <- ggplot(
  confounding_plot_data,
  aes(
    x = median, y = scenario_label, xmin = q025, xmax = q975,
    color = pathway, shape = pathway
  )
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.4) +
  geom_pointrange(position = position_dodge(width = 0.55), linewidth = 0.45) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
    labels = label_number(accuracy = 1)
  ) +
  scale_color_manual(values = c(
    "L -> treatment" = "#0072B2",
    "L -> event" = "#D55E00"
  )) +
  labs(x = "DGP-implied odds ratio: L Q3 versus Q1 (log scale)", y = NULL) +
  paper_theme

single_metric_plot <- function(metric_name, x_label, color, x_limits = NULL) {
  dat <- plot_summary %>% filter(metric == metric_name)
  result <- ggplot(
    dat,
    aes(x = median, y = scenario_label, xmin = q025, xmax = q975)
  ) +
    geom_pointrange(color = color, linewidth = 0.45) +
    labs(x = x_label, y = NULL) +
    paper_theme +
    theme(legend.position = "none")
  if (!is.null(x_limits)) {
    result <- result + coord_cartesian(xlim = x_limits)
  }
  result
}

p_positivity <- single_metric_plot(
  "pA_extreme_pct",
  "Monthly observations with p(A = 1 | history) <= 0.05 or >= 0.95 (%)",
  "#CC79A7"
)
p_adherence <- single_metric_plot(
  "continuation_pct",
  "Treatment continuation when previously treated (%)",
  "#009E73",
  c(0, 100)
)
event_upper <- max(
  plot_summary$q975[plot_summary$metric == "event_risk_12_pct"],
  na.rm = TRUE
) * 1.05
p_event <- single_metric_plot(
  "event_risk_12_pct",
  "12-month cumulative event risk (%)",
  "#D55E00",
  c(0, event_upper)
)

main_figure <- (p_confounding | p_positivity) /
  (p_adherence | p_event) +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(output_dir, "figure_main_scenario_diagnostics.pdf"),
  main_figure,
  width = 11,
  height = 8.2,
  units = "in"
)
ggsave(
  file.path(output_dir, "figure_main_scenario_diagnostics.png"),
  main_figure,
  width = 11,
  height = 8.2,
  units = "in",
  dpi = 320
)

positivity_summary <- positivity_histograms %>%
  group_by(scenario, probability_midpoint) %>%
  summarise(
    median = median(percent_person_months),
    q025 = as.numeric(quantile(percent_person_months, 0.025)),
    q975 = as.numeric(quantile(percent_person_months, 0.975)),
    .groups = "drop"
  ) %>%
  mutate(scenario_label = factor(
    unname(scenario_labels[as.character(scenario)]),
    levels = unname(scenario_labels[scenario_order])
  ))

positivity_figure <- ggplot(
  positivity_summary,
  aes(x = probability_midpoint, y = median)
) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#56B4E9", alpha = 0.20) +
  geom_line(color = "#0072B2", linewidth = 0.55) +
  geom_vline(xintercept = c(0.05, 0.95), linetype = 2, color = "grey45", linewidth = 0.35) +
  facet_wrap(~scenario_label, ncol = 2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1), limits = c(0, 1)) +
  labs(
    x = "DGP-implied treatment-assignment probability",
    y = "Monthly observations per 0.02 probability bin (%)"
  ) +
  paper_theme +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "figure_positivity_distributions.pdf"),
  positivity_figure,
  width = 8.2,
  height = 9.2,
  units = "in"
)
ggsave(
  file.path(output_dir, "figure_positivity_distributions.png"),
  positivity_figure,
  width = 8.2,
  height = 9.2,
  units = "in",
  dpi = 320
)

event_curve_summary <- event_curves %>%
  group_by(scenario, month) %>%
  summarise(
    median = median(cumulative_event_risk_pct, na.rm = TRUE),
    q025 = as.numeric(quantile(cumulative_event_risk_pct, 0.025, na.rm = TRUE)),
    q975 = as.numeric(quantile(cumulative_event_risk_pct, 0.975, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(scenario_label = factor(
    unname(scenario_labels[as.character(scenario)]),
    levels = unname(scenario_labels[scenario_order])
  ))

event_curve_figure <- ggplot(
  event_curve_summary,
  aes(x = month, y = median)
) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#E69F00", alpha = 0.20) +
  geom_line(color = "#D55E00", linewidth = 0.65) +
  geom_vline(xintercept = 12, linetype = 2, color = "grey45", linewidth = 0.35) +
  facet_wrap(~scenario_label, ncol = 2) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(1, 24)) +
  labs(x = "Month", y = "Cumulative event risk (%)") +
  paper_theme +
  theme(legend.position = "none")
ggsave(
  file.path(output_dir, "figure_event_risk_curves.pdf"),
  event_curve_figure,
  width = 8.2,
  height = 9.2,
  units = "in"
)
ggsave(
  file.path(output_dir, "figure_event_risk_curves.png"),
  event_curve_figure,
  width = 8.2,
  height = 9.2,
  units = "in",
  dpi = 320
)

saveRDS(
  list(
    metadata = list(
      n_replications_per_scenario = 500L,
      n_seed_batches = 6L,
      interval = "empirical 2.5th to 97.5th percentiles across replications"
    ),
    replication_diagnostics = replication_diagnostics,
    summary_diagnostics_long = summary_diagnostics_long,
    batch_comparison = batch_comparison,
    event_curve_summary = event_curve_summary,
    positivity_summary = positivity_summary,
    formatted_table = formatted_table
  ),
  file.path(output_dir, "scenario_diagnostics_500_replications.rds")
)

message("Finished 500-replication diagnostics: ", output_dir)
