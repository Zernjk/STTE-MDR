############################################################
# Doubly robust helper functions for STTE estimators
#
# These helpers implement AIPW-style risk estimators:
#   plug-in g-computation risk + weighted residual correction.
#
# The pooled estimators use residuals computed from remaining-risk predictions,
#   Delta_{k+1} = dY_{k+1} + (1 - dY_{k+1}) V_{k+1} - V_k,
# where V_k is computed recursively from fitted pooled-logistic hazards.
# In that remaining-risk setup, the IPCW correction is a
# strategy baseline-population average:
#   sum_k sum_i w_{igk} Delta_{i,k+1} / n_baseline_g.
# The MDR correction estimates the target residual mean in each
# strategy-follow cell, then multiplies by the target cell mass
# n_target_gk / n_target_g0.
#
# The state estimators use the same target-risk definition as their non-DR
# state counterparts: target-cell mean hazards combined as
#   1 - prod_k(1 - h_gk).
# Their residual correction is therefore applied on the hazard scale within
# strategy-follow cells before recomputing the cumulative risk.
############################################################

.dr_required_columns <- function(data, columns, data_name = "data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.dr_failure <- function(
    failure_reason,
    risk_g1 = NA_real_,
    risk_g0 = NA_real_,
    n_target_missing_cells = NA_integer_,
    n_target_present_adherent_missing_cells = NA_integer_,
    n_target_without_observed_outcome_cells = NA_integer_) {
  list(
    psi_hat = NA_real_,
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    n_target_missing_cells = n_target_missing_cells,
    n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells,
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    failure_reason = failure_reason
  )
}

.dr_weight_summary <- function(data, weight_col) {
  if (!weight_col %in% names(data) || nrow(data) == 0L) {
    return(list(weight_mean = NA_real_, weight_max = NA_real_, weight_cv = NA_real_))
  }

  w <- as.numeric(data[[weight_col]])
  w <- w[is.finite(w)]
  if (length(w) == 0L) {
    return(list(weight_mean = NA_real_, weight_max = NA_real_, weight_cv = NA_real_))
  }

  weight_mean <- mean(w, na.rm = TRUE)
  list(
    weight_mean = weight_mean,
    weight_max = max(w, na.rm = TRUE),
    weight_cv = ifelse(
      abs(weight_mean) < 1e-12,
      NA_real_,
      stats::sd(w, na.rm = TRUE) / weight_mean
    )
  )
}

.dr_add_remaining_risk_predictions <- function(
    data,
    outcome_fit,
    params,
    fallback_hazard,
    path_cols,
    event_col = "dY_next",
    follow_col = "follow_k") {

  .dr_required_columns(data, c(path_cols, follow_col), "remaining-risk data")
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "DR remaining-risk prediction requires package 'data.table'.",
      call. = FALSE
    )
  }

  data <- as.data.frame(data)

  if (!"follow_f" %in% names(data)) {
    data$follow_f <- factor(data[[follow_col]], levels = 0:(params$tau - 1L))
  }

  data$dr_hazard <- .predict_prob_safe(
    outcome_fit,
    newdata = data,
    fallback = fallback_hazard,
    eps = 1e-8
  )
  data$dr_V_k <- NA_real_
  data$dr_V_next <- NA_real_
  data$dr_has_contiguous_next <- FALSE
  data$dr_missing_nonterminal_next <- FALSE

  if (nrow(data) == 0L) {
    return(data)
  }

  duplicate_keys <- .dr_duplicate_key_rows(data, c(path_cols, follow_col))
  if (nrow(duplicate_keys) > 0L) {
    stop(
      .dr_duplicate_key_failure_reason(duplicate_keys, "remaining-risk data"),
      call. = FALSE
    )
  }

  has_event_col <- event_col %in% names(data)

  dt <- data.table::as.data.table(data)
  dt$.dr_original_row <- seq_len(nrow(dt))
  dt$.dr_follow_numeric <- as.numeric(dt[[follow_col]])

  bad_follow <- !is.finite(dt$.dr_follow_numeric)
  if (any(bad_follow)) {
    stop(
      "remaining-risk data has non-finite follow-up values.",
      call. = FALSE
    )
  }

  data.table::setorderv(dt, c(path_cols, follow_col))
  dt$.dr_sorted_row <- seq_len(nrow(dt))

  dt[
    ,
    .dr_next_sorted_row := data.table::shift(.dr_sorted_row, type = "lead"),
    by = path_cols
  ]
  dt[
    ,
    .dr_next_follow := data.table::shift(.dr_follow_numeric, type = "lead"),
    by = path_cols
  ]

  n_rows <- nrow(dt)
  v_k <- rep(NA_real_, n_rows)
  v_next <- rep(NA_real_, n_rows)
  has_contiguous_next <- rep(FALSE, n_rows)
  missing_nonterminal_next <- rep(FALSE, n_rows)

  follow_values <- sort(unique(dt$.dr_follow_numeric), decreasing = TRUE)
  terminal_follow <- params$tau - 1L

  for (follow_value in follow_values) {
    cur <- which(dt$.dr_follow_numeric == follow_value)
    next_idx <- dt$.dr_next_sorted_row[cur]
    next_follow <- dt$.dr_next_follow[cur]

    contiguous <- !is.na(next_idx) &
      !is.na(next_follow) &
      next_follow == follow_value + 1

    future_v <- rep(0, length(cur))
    if (any(contiguous)) {
      future_v[contiguous] <- v_k[next_idx[contiguous]]
    }

    is_final_follow <- follow_value >= terminal_follow
    if (has_event_col) {
      is_terminal_event <- as.integer(dt[[event_col]][cur]) == 1L
      is_terminal_event[is.na(is_terminal_event)] <- FALSE
    } else {
      is_terminal_event <- rep(FALSE, length(cur))
    }

    nonterminal_gap <- !contiguous &
      !isTRUE(is_final_follow) &
      !is_terminal_event

    if (has_event_col && any(nonterminal_gap)) {
      future_v[nonterminal_gap] <- NA_real_
    }

    h <- dt$dr_hazard[cur]
    v_next[cur] <- future_v
    v_k[cur] <- h + (1 - h) * future_v
    has_contiguous_next[cur] <- contiguous
    missing_nonterminal_next[cur] <- has_event_col & nonterminal_gap
  }

  dt$dr_V_k <- v_k
  dt$dr_V_next <- v_next
  dt$dr_has_contiguous_next <- has_contiguous_next
  dt$dr_missing_nonterminal_next <- missing_nonterminal_next

  data.table::setorder(dt, .dr_original_row)
  drop_cols <- c(
    ".dr_original_row",
    ".dr_follow_numeric",
    ".dr_sorted_row",
    ".dr_next_sorted_row",
    ".dr_next_follow"
  )
  dt[, (drop_cols) := NULL]

  as.data.frame(dt)
}

.dr_duplicate_key_rows <- function(data, key_cols) {
  .dr_required_columns(data, key_cols, "keyed data")

  data %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(key_cols))) %>%
    dplyr::summarise(n_key = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_key > 1L)
}

.dr_duplicate_key_failure_reason <- function(duplicate_keys, data_name) {
  shown <- utils::head(duplicate_keys, 5L)
  key_cols <- setdiff(names(shown), "n_key")

  key_text <- apply(shown[key_cols], 1L, function(row) {
    paste(paste(key_cols, row, sep = "="), collapse = ", ")
  })

  paste0(
    data_name,
    " has duplicate rows for residual join keys: ",
    paste0(key_text, collapse = "; "),
    if (nrow(duplicate_keys) > nrow(shown)) {
      paste0("; ... plus ", nrow(duplicate_keys) - nrow(shown), " more")
    } else {
      ""
    }
  )
}

.dr_mark_residual_transition_valid <- function(data) {
  valid <- rep(TRUE, nrow(data))

  # dY_next precedes natural censoring in the simulator, so the current
  # interval outcome remains usable when dC_next == 1. Protocol adherence is
  # determined by A_k before dY_next and must still hold for the current row.
  if ("artificial_censor" %in% names(data)) {
    valid <- valid & !is.na(data$artificial_censor) & data$artificial_censor == 0L
  }
  if ("D_adherent" %in% names(data)) {
    valid <- valid & !is.na(data$D_adherent) & data$D_adherent == 1L
  }

  data$dr_residual_transition_valid <- valid
  data
}

.dr_complete_pooled_prediction_paths <- function(
    path_data,
    params,
    path_cols,
    follow_col = "follow_k") {

  .dr_required_columns(
    path_data,
    c(path_cols, follow_col),
    "pooled prediction path data"
  )

  path_starts <- path_data %>%
    dplyr::arrange(
      dplyr::across(tidyselect::all_of(path_cols)),
      .data[[follow_col]]
    ) %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(path_cols))) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  if (nrow(path_starts) == 0L) {
    return(path_starts)
  }
  if (any(path_starts[[follow_col]] != 0L)) {
    stop(
      "Every pooled residual path must have a follow_k = 0 baseline row.",
      call. = FALSE
    )
  }

  prediction_paths <- path_starts[
    rep(seq_len(nrow(path_starts)), each = params$tau),
    ,
    drop = FALSE
  ]
  prediction_paths[[follow_col]] <- rep(
    0:(params$tau - 1L),
    times = nrow(path_starts)
  )
  prediction_paths$follow_f <- factor(
    prediction_paths[[follow_col]],
    levels = 0:(params$tau - 1L)
  )

  prediction_paths
}

.dr_add_aipw_residual <- function(
    data,
    outcome_fit,
    params,
    fallback_hazard,
    path_cols,
    path_data = NULL,
    join_cols = NULL,
    event_col = "dY_next") {

  .dr_required_columns(data, event_col, "residual data")

  if (is.null(path_data)) {
    path_data <- data
  }
  if (is.null(join_cols)) {
    join_cols <- unique(c(path_cols, "follow_k"))
  }

  duplicate_path_keys <- .dr_duplicate_key_rows(path_data, join_cols)
  if (nrow(duplicate_path_keys) > 0L) {
    stop(
      .dr_duplicate_key_failure_reason(duplicate_path_keys, "path_data"),
      call. = FALSE
    )
  }

  duplicate_residual_keys <- .dr_duplicate_key_rows(data, join_cols)
  if (nrow(duplicate_residual_keys) > 0L) {
    stop(
      .dr_duplicate_key_failure_reason(duplicate_residual_keys, "residual data"),
      call. = FALSE
    )
  }

  # The pooled outcome regressions depend on baseline covariates and follow-up
  # time, not future L. Complete each path through tau so censoring or protocol
  # deviation does not silently turn all earlier V_{k+1} predictions into NA.
  path_prediction_data <- .dr_complete_pooled_prediction_paths(
    path_data = path_data,
    params = params,
    path_cols = path_cols
  )
  path_predictions <- .dr_add_remaining_risk_predictions(
    data = path_prediction_data,
    outcome_fit = outcome_fit,
    params = params,
    fallback_hazard = fallback_hazard,
    path_cols = path_cols,
    event_col = ".dr_no_observed_event"
  )

  path_values <- path_predictions %>%
    dplyr::select(
      tidyselect::all_of(join_cols),
      dr_hazard,
      dr_V_k,
      dr_V_next,
      dr_has_contiguous_next,
      dr_missing_nonterminal_next
    )

  data <- data %>%
    dplyr::left_join(path_values, by = join_cols)
  data <- .dr_mark_residual_transition_valid(data)

  y <- as.numeric(data[[event_col]])
  data$dr_residual <- NA_real_
  use_residual <- data$dr_residual_transition_valid
  data$dr_residual[use_residual] <-
    y[use_residual] +
    (1 - y[use_residual]) * data$dr_V_next[use_residual] -
    data$dr_V_k[use_residual]
  data
}

.dr_add_hazard_residual <- function(
    data,
    outcome_fit,
    params,
    fallback_hazard,
    event_col = "dY_next") {

  .dr_required_columns(data, c(event_col, "follow_k"), "hazard residual data")

  data <- as.data.frame(data)
  if (!"follow_f" %in% names(data)) {
    data$follow_f <- factor(data$follow_k, levels = 0:(params$tau - 1L))
  }

  data$dr_hazard <- .predict_prob_safe(
    outcome_fit,
    newdata = data,
    fallback = fallback_hazard,
    eps = 1e-8
  )

  data <- .dr_mark_residual_transition_valid(data)

  y <- as.numeric(data[[event_col]])
  data$dr_residual <- NA_real_
  use_residual <- data$dr_residual_transition_valid &
    is.finite(y) &
    is.finite(data$dr_hazard)

  data$dr_residual[use_residual] <- y[use_residual] - data$dr_hazard[use_residual]
  data
}

.dr_bad_vnext_rows <- function(
    residual_data,
    params,
    event_col = "dY_next") {

  .dr_required_columns(
    residual_data,
    c("follow_k", event_col, "dr_missing_nonterminal_next"),
    "residual data"
  )

  if (!"dr_residual_transition_valid" %in% names(residual_data)) {
    residual_data$dr_residual_transition_valid <- TRUE
  }

  residual_data %>%
    dplyr::filter(
      dr_residual_transition_valid,
      follow_k < params$tau - 1L,
      .data[[event_col]] == 0L,
      dr_missing_nonterminal_next
    ) %>%
    dplyr::mutate(
      dr_vnext_problem = "missing_next_state"
    )
}

.dr_bad_vnext_failure_reason <- function(bad_vnext_rows) {
  strategy_label <- if ("strategy" %in% names(bad_vnext_rows)) {
    as.character(bad_vnext_rows$strategy)
  } else if ("clone" %in% names(bad_vnext_rows)) {
    as.character(bad_vnext_rows$clone)
  } else {
    rep("unknown", nrow(bad_vnext_rows))
  }

  bad_vnext_rows$dr_strategy_label <- strategy_label
  diagnostic <- bad_vnext_rows %>%
    dplyr::count(dr_strategy_label, follow_k, dr_vnext_problem, name = "n_bad") %>%
    dplyr::arrange(dr_strategy_label, follow_k, dr_vnext_problem)

  paste0(
    "Remaining-risk recursion found valid nonterminal residual rows with no next state: ",
    paste0(
      "g=", diagnostic$dr_strategy_label,
      ", k=", diagnostic$follow_k,
      ", problem=", diagnostic$dr_vnext_problem,
      ", n=", diagnostic$n_bad,
      collapse = "; "
    )
  )
}

.dr_risk_from_prediction_rows <- function(
    prediction_rows,
    strategy_col = "strategy") {

  prediction_rows$dr_strategy <- as.integer(prediction_rows[[strategy_col]])

  risk_table <- prediction_rows %>%
    dplyr::filter(follow_k == 0L) %>%
    dplyr::group_by(dr_strategy) %>%
    dplyr::summarise(
      risk = mean(dr_V_k, na.rm = TRUE),
      n_paths = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(strategy = dr_strategy)

  hazard_table <- prediction_rows %>%
    dplyr::group_by(dr_strategy, follow_k) %>%
    dplyr::summarise(
      hazard = mean(dr_hazard, na.rm = TRUE),
      remaining_risk = mean(dr_V_k, na.rm = TRUE),
      n_prediction = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(strategy = dr_strategy)

  risk_g0 <- .extract_scalar_or_na(risk_table$risk[risk_table$strategy == 0L])
  risk_g1 <- .extract_scalar_or_na(risk_table$risk[risk_table$strategy == 1L])

  list(
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    psi_hat = risk_g1 - risk_g0,
    risk_table = risk_table,
    hazard_table = hazard_table,
    prediction_rows = prediction_rows
  )
}

.dr_estimate_state_plugin_risk <- function(
    outcome_fit,
    target_states,
    params,
    fallback_hazard) {

  .dr_required_columns(target_states, c("strategy", "follow_k"), "target states")

  pred_data <- target_states %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  pred_data$hazard <- .predict_prob_safe(
    outcome_fit,
    newdata = pred_data,
    fallback = fallback_hazard,
    eps = 1e-8
  )
  pred_data$dr_hazard <- pred_data$hazard

  hazard_table <- pred_data %>%
    dplyr::group_by(strategy, follow_k) %>%
    dplyr::summarise(
      hazard = mean(hazard, na.rm = TRUE),
      n_target = dplyr::n(),
      .groups = "drop"
    )

  missing_strategy_follow <- tidyr::expand_grid(
    strategy = c(0L, 1L),
    follow_k = 0:(params$tau - 1L)
  ) %>%
    dplyr::anti_join(
      hazard_table %>% dplyr::select(strategy, follow_k),
      by = c("strategy", "follow_k")
    )

  if (nrow(missing_strategy_follow) > 0L) {
    warning(
      "DR state plug-in risk prediction is missing strategy-follow cells: ",
      paste0(
        "g=", missing_strategy_follow$strategy,
        ", k=", missing_strategy_follow$follow_k,
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  risk_table <- hazard_table %>%
    dplyr::arrange(strategy, follow_k) %>%
    dplyr::group_by(strategy) %>%
    dplyr::summarise(
      risk = 1 - prod(1 - hazard),
      .groups = "drop"
    )

  risk_g0 <- .extract_scalar_or_na(risk_table$risk[risk_table$strategy == 0L])
  risk_g1 <- .extract_scalar_or_na(risk_table$risk[risk_table$strategy == 1L])

  list(
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    psi_hat = risk_g1 - risk_g0,
    risk_table = risk_table,
    hazard_table = hazard_table,
    prediction_rows = pred_data
  )
}

.dr_make_strategy_baseline_prediction_rows <- function(
    eligible_base,
    params,
    model_strategy_col = "strategy") {

  eligible_base <- eligible_base %>%
    dplyr::distinct(id_trial, trial_m, X1, X2, X3, L_trial)

  prediction_rows <- lapply(c(0L, 1L), function(strategy_value) {
    pred_grid <- eligible_base[
      rep(seq_len(nrow(eligible_base)), each = params$tau),
      ,
      drop = FALSE
    ]
    pred_grid$follow_k <- rep(0:(params$tau - 1L), times = nrow(eligible_base))
    pred_grid$follow_f <- factor(pred_grid$follow_k, levels = 0:(params$tau - 1L))
    pred_grid$strategy <- strategy_value
    pred_grid[[model_strategy_col]] <- strategy_value
    pred_grid
  })

  dplyr::bind_rows(prediction_rows)
}

.dr_estimate_baseline_plugin_risk <- function(
    outcome_fit,
    eligible_base,
    params,
    fallback_hazard,
    model_strategy_col = "strategy") {

  pred_data <- .dr_make_strategy_baseline_prediction_rows(
    eligible_base = eligible_base,
    params = params,
    model_strategy_col = model_strategy_col
  )

  pred_data <- .dr_add_remaining_risk_predictions(
    data = pred_data,
    outcome_fit = outcome_fit,
    params = params,
    fallback_hazard = fallback_hazard,
    path_cols = c("strategy", "id_trial")
  )

  .dr_risk_from_prediction_rows(pred_data, strategy_col = "strategy")
}

.dr_constant_baseline_counts <- function(n_baseline, strategy_values = c(0L, 1L)) {
  data.frame(
    strategy = as.integer(strategy_values),
    n_baseline = rep(as.numeric(n_baseline), length(strategy_values))
  )
}

.dr_target_cell_mass <- function(target_states, target_deplete_events = TRUE) {
  .dr_required_columns(
    target_states,
    c("strategy", "follow_k"),
    "target states"
  )

  if ("target_at_risk_start" %in% names(target_states)) {
    target_states <- target_states %>%
      dplyr::filter(.data$target_at_risk_start == 1L)
  } else if (identical(target_deplete_events, FALSE)) {
    stop(
      "DR-MDR target-cell mass requires target rows to be depleted after events, ",
      "or target_states_all must include target_at_risk_start. ",
      "target_deplete_events = FALSE retains post-event target rows, so raw ",
      "strategy-follow counts are not risk-set masses.",
      call. = FALSE
    )
  }

  target_states %>%
    dplyr::count(strategy, follow_k, name = "n_target_cell") %>%
    dplyr::left_join(
      target_states %>%
        dplyr::filter(follow_k == 0L) %>%
        dplyr::count(strategy, name = "n_target_baseline"),
      by = "strategy"
    ) %>%
    dplyr::mutate(
      strategy = as.integer(strategy),
      target_cell_mass = n_target_cell / n_target_baseline
    ) %>%
    dplyr::select(strategy, follow_k, n_target_cell, n_target_baseline, target_cell_mass)
}

.dr_cell_support_failure_reason <- function(bad_cells, label = "residual") {
  paste0(
    "Target strategy-follow cells have no observed ",
    label,
    " rows: ",
    paste0("g=", bad_cells$strategy, ", k=", bad_cells$follow_k, collapse = "; ")
  )
}

.dr_weighted_residual_correction <- function(
    residual_data,
    weight_col,
    baseline_counts = NULL,
    target_cell_mass = NULL,
    correction_scale = c("baseline_total", "target_cell_mass"),
    strategy_col = "strategy",
    residual_col = "dr_residual") {

  correction_scale <- match.arg(correction_scale)

  .dr_required_columns(
    residual_data,
    c(weight_col, strategy_col, "follow_k", residual_col),
    "residual data"
  )

  residual_data <- residual_data
  residual_data$dr_strategy <- as.integer(residual_data[[strategy_col]])
  residual_data$dr_weight <- as.numeric(residual_data[[weight_col]])

  if (correction_scale == "baseline_total") {
    .dr_required_columns(
      baseline_counts,
      c("strategy", "n_baseline"),
      "baseline counts"
    )
    baseline_counts$strategy <- as.integer(baseline_counts$strategy)
    baseline_counts$n_baseline <- as.numeric(baseline_counts$n_baseline)

    bad_baseline_counts <- baseline_counts %>%
      dplyr::filter(!is.finite(n_baseline) | n_baseline <= 0)
    if (nrow(bad_baseline_counts) > 0L) {
      return(list(
        failed = TRUE,
        failure_reason = paste0(
          "Invalid baseline denominator for strategy/strategies: ",
          paste(bad_baseline_counts$strategy, collapse = ", ")
        ),
        bad_cells = bad_baseline_counts
      ))
    }
  } else {
    .dr_required_columns(
      target_cell_mass,
      c("strategy", "follow_k", "target_cell_mass"),
      "target cell mass"
    )
    target_cell_mass$strategy <- as.integer(target_cell_mass$strategy)
    target_cell_mass$target_cell_mass <- as.numeric(target_cell_mass$target_cell_mass)
    if (!"n_target_cell" %in% names(target_cell_mass)) {
      target_cell_mass$n_target_cell <- NA_integer_
    }

    bad_target_mass <- target_cell_mass %>%
      dplyr::filter(!is.finite(.data$target_cell_mass) | .data$target_cell_mass < 0)
    if (nrow(bad_target_mass) > 0L) {
      return(list(
        failed = TRUE,
        failure_reason = paste0(
          "Invalid target cell mass for strategy-follow cells: ",
          paste0(
            "g=", bad_target_mass$strategy,
            ", k=", bad_target_mass$follow_k,
            collapse = "; "
          )
        ),
        bad_cells = bad_target_mass
      ))
    }
  }

  residual_data <- residual_data[
    is.finite(residual_data$dr_weight) &
      is.finite(residual_data[[residual_col]]),
    ,
    drop = FALSE
  ]

  if (correction_scale == "target_cell_mass") {
    positive_target_cells <- target_cell_mass %>%
      dplyr::filter(.data$target_cell_mass > 0) %>%
      dplyr::distinct(strategy, follow_k)

    available_cells <- residual_data %>%
      dplyr::distinct(strategy = dr_strategy, follow_k)

    bad_cells <- positive_target_cells %>%
      dplyr::anti_join(available_cells, by = c("strategy", "follow_k")) %>%
      dplyr::arrange(strategy, follow_k)

    if (nrow(bad_cells) > 0L) {
      return(list(
        failed = TRUE,
        failure_reason = .dr_cell_support_failure_reason(bad_cells),
        bad_cells = bad_cells
      ))
    }
  }

  correction_by_cell <- residual_data %>%
    dplyr::group_by(dr_strategy, follow_k) %>%
    dplyr::summarise(
      n_residual = dplyr::n(),
      weight_sum = sum(dr_weight, na.rm = TRUE),
      weighted_residual_sum = sum(dr_weight * .data[[residual_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(strategy = dr_strategy)

  if (correction_scale == "baseline_total") {
    correction_by_cell <- correction_by_cell %>%
      dplyr::left_join(baseline_counts, by = "strategy") %>%
      dplyr::mutate(
        weighted_residual_mean = ifelse(
          weight_sum > 0,
          weighted_residual_sum / weight_sum,
          NA_real_
        ),
        residual_correction = weighted_residual_sum / n_baseline
      )
  } else {
    correction_by_cell <- correction_by_cell %>%
      dplyr::left_join(target_cell_mass, by = c("strategy", "follow_k")) %>%
      dplyr::mutate(
        n_target_cell = dplyr::coalesce(n_target_cell, 0L),
        target_cell_mass = dplyr::coalesce(target_cell_mass, 0)
      ) %>%
      dplyr::mutate(
        weighted_residual_mean = ifelse(
          weight_sum > 0,
          weighted_residual_sum / weight_sum,
          NA_real_
        ),
        residual_correction = ifelse(
          target_cell_mass == 0,
          0,
          target_cell_mass * weighted_residual_mean
        )
      )
  }

  bad_scale_cells <- correction_by_cell %>%
    dplyr::filter(!is.finite(residual_correction))
  if (nrow(bad_scale_cells) > 0L) {
    return(list(
      failed = TRUE,
      failure_reason = paste0(
        "Some strategy-follow residual cells have invalid correction scale: ",
        paste0(
          "g=", bad_scale_cells$strategy,
          ", k=", bad_scale_cells$follow_k,
          collapse = "; "
        )
      ),
      bad_cells = bad_scale_cells
    ))
  }

  correction_by_cell <- correction_by_cell %>%
    dplyr::mutate(
      correction_scale = correction_scale
    )

  bad_weight_cells <- correction_by_cell %>%
    dplyr::filter(!is.finite(weight_sum))

  if (nrow(bad_weight_cells) > 0L) {
    return(list(
      failed = TRUE,
      failure_reason = paste0(
        "Some strategy-follow residual cells have invalid weights: ",
        paste0(
          "g=", bad_weight_cells$strategy,
          ", k=", bad_weight_cells$follow_k,
          collapse = "; "
        )
      ),
      bad_cells = bad_weight_cells
    ))
  }

  correction_table <- correction_by_cell %>%
    dplyr::group_by(strategy) %>%
    dplyr::summarise(
      residual_correction = sum(residual_correction),
      .groups = "drop"
    )

  list(
    failed = FALSE,
    correction_g1 = .extract_scalar_or_na(
      correction_table$residual_correction[correction_table$strategy == 1L]
    ),
    correction_g0 = .extract_scalar_or_na(
      correction_table$residual_correction[correction_table$strategy == 0L]
    ),
    correction_table = correction_table,
    correction_by_cell = correction_by_cell,
    residual_data = residual_data
  )
}

.dr_weighted_hazard_residual_correction <- function(
    residual_data,
    weight_col,
    target_hazard_table,
    strategy_col = "strategy",
    residual_col = "dr_residual") {

  .dr_required_columns(
    residual_data,
    c(weight_col, strategy_col, "follow_k", residual_col),
    "hazard residual data"
  )
  .dr_required_columns(
    target_hazard_table,
    c("strategy", "follow_k", "hazard"),
    "target hazard table"
  )

  target_cells <- target_hazard_table %>%
    dplyr::select(
      strategy,
      follow_k,
      tidyselect::any_of(c("n_target", "n_prediction"))
    ) %>%
    dplyr::distinct()
  target_cells$strategy <- as.integer(target_cells$strategy)

  residual_data <- residual_data
  residual_data$dr_strategy <- as.integer(residual_data[[strategy_col]])
  residual_data$dr_weight <- as.numeric(residual_data[[weight_col]])

  negative_weight_rows <- residual_data %>%
    dplyr::filter(is.finite(dr_weight), dr_weight < 0)
  if (nrow(negative_weight_rows) > 0L) {
    return(list(
      failed = TRUE,
      failure_reason = "Some hazard residual rows have negative weights",
      bad_cells = negative_weight_rows
    ))
  }

  residual_data <- residual_data[
    is.finite(residual_data$dr_weight) &
      is.finite(residual_data[[residual_col]]),
    ,
    drop = FALSE
  ]

  available_cells <- residual_data %>%
    dplyr::distinct(strategy = dr_strategy, follow_k)

  bad_cells <- target_cells %>%
    dplyr::select(strategy, follow_k) %>%
    dplyr::anti_join(available_cells, by = c("strategy", "follow_k")) %>%
    dplyr::arrange(strategy, follow_k)

  if (nrow(bad_cells) > 0L) {
    return(list(
      failed = TRUE,
      failure_reason = .dr_cell_support_failure_reason(bad_cells),
      bad_cells = bad_cells
    ))
  }

  cell_summaries <- residual_data %>%
    dplyr::group_by(dr_strategy, follow_k) %>%
    dplyr::summarise(
      n_residual = dplyr::n(),
      weight_sum = sum(dr_weight, na.rm = TRUE),
      weighted_residual_sum = sum(dr_weight * .data[[residual_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(strategy = dr_strategy)

  correction_by_cell <- target_cells %>%
    dplyr::left_join(cell_summaries, by = c("strategy", "follow_k")) %>%
    dplyr::mutate(
      n_residual = dplyr::coalesce(n_residual, 0L),
      weight_sum = dplyr::coalesce(weight_sum, 0),
      weighted_residual_sum = dplyr::coalesce(weighted_residual_sum, 0),
      weighted_residual_mean = ifelse(
        weight_sum > 0,
        weighted_residual_sum / weight_sum,
        NA_real_
      ),
      hazard_residual_correction = weighted_residual_mean,
      residual_correction = hazard_residual_correction,
      correction_scale = "hazard_cell_mean"
    )

  bad_scale_cells <- correction_by_cell %>%
    dplyr::filter(!is.finite(hazard_residual_correction))
  if (nrow(bad_scale_cells) > 0L) {
    return(list(
      failed = TRUE,
      failure_reason = paste0(
        "Some strategy-follow hazard residual cells have invalid correction scale: ",
        paste0(
          "g=", bad_scale_cells$strategy,
          ", k=", bad_scale_cells$follow_k,
          collapse = "; "
        )
      ),
      bad_cells = bad_scale_cells
    ))
  }

  list(
    failed = FALSE,
    correction_by_cell = correction_by_cell,
    residual_data = residual_data
  )
}

.dr_make_hazard_corrected_estimate <- function(
    plugin_obj,
    correction_obj,
    n_rows_stte = NA_integer_,
    n_rows_target = NA_integer_,
    n_rows_adherent = NA_integer_,
    n_rows_outcome = NA_integer_,
    weight_summary = list(weight_mean = NA_real_, weight_max = NA_real_, weight_cv = NA_real_),
    n_target_missing_cells = NA_integer_,
    n_target_present_adherent_missing_cells = NA_integer_,
    n_target_without_observed_outcome_cells = NA_integer_,
    include_tables = FALSE) {

  correction_cells <- correction_obj$correction_by_cell %>%
    dplyr::select(
      strategy,
      follow_k,
      n_residual,
      weight_sum,
      weighted_residual_sum,
      weighted_residual_mean,
      hazard_residual_correction,
      residual_correction,
      correction_scale
    )

  corrected_hazard_table <- plugin_obj$hazard_table %>%
    dplyr::left_join(correction_cells, by = c("strategy", "follow_k")) %>%
    dplyr::mutate(
      corrected_hazard_raw = hazard + hazard_residual_correction,
      corrected_hazard = pmin(pmax(corrected_hazard_raw, 0), 1),
      hazard_was_clipped = corrected_hazard != corrected_hazard_raw
    )

  corrected_risk_table <- corrected_hazard_table %>%
    dplyr::arrange(strategy, follow_k) %>%
    dplyr::group_by(strategy) %>%
    dplyr::summarise(
      risk = 1 - prod(1 - corrected_hazard),
      plugin_risk = 1 - prod(1 - hazard),
      n_hazard_cells = dplyr::n(),
      n_hazard_cells_clipped = sum(hazard_was_clipped, na.rm = TRUE),
      .groups = "drop"
    )

  risk_g0 <- .extract_scalar_or_na(
    corrected_risk_table$risk[corrected_risk_table$strategy == 0L]
  )
  risk_g1 <- .extract_scalar_or_na(
    corrected_risk_table$risk[corrected_risk_table$strategy == 1L]
  )
  correction_g1 <- risk_g1 - plugin_obj$risk_g1
  correction_g0 <- risk_g0 - plugin_obj$risk_g0

  out <- list(
    psi_hat = risk_g1 - risk_g0,
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    plugin_risk_g1 = plugin_obj$risk_g1,
    plugin_risk_g0 = plugin_obj$risk_g0,
    residual_correction_g1 = correction_g1,
    residual_correction_g0 = correction_g0,
    n_rows_stte = n_rows_stte,
    n_rows_target = n_rows_target,
    n_rows_adherent = n_rows_adherent,
    n_rows_outcome = n_rows_outcome,
    weight_mean = weight_summary$weight_mean,
    weight_max = weight_summary$weight_max,
    weight_cv = weight_summary$weight_cv,
    n_target_missing_cells = n_target_missing_cells,
    n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells,
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    n_hazard_cells_clipped = sum(corrected_hazard_table$hazard_was_clipped, na.rm = TRUE),
    failure_reason = NA_character_
  )

  if (isTRUE(include_tables)) {
    out$corrected_risk_table <- corrected_risk_table
    out$corrected_hazard_table <- corrected_hazard_table
  }

  out
}

.dr_make_estimate <- function(
    plugin_obj,
    correction_obj,
    n_rows_stte = NA_integer_,
    n_rows_target = NA_integer_,
    n_rows_adherent = NA_integer_,
    n_rows_outcome = NA_integer_,
    weight_summary = list(weight_mean = NA_real_, weight_max = NA_real_, weight_cv = NA_real_),
    n_target_missing_cells = NA_integer_,
    n_target_present_adherent_missing_cells = NA_integer_,
    n_target_without_observed_outcome_cells = NA_integer_) {

  risk_g1 <- plugin_obj$risk_g1 + correction_obj$correction_g1
  risk_g0 <- plugin_obj$risk_g0 + correction_obj$correction_g0

  list(
    psi_hat = risk_g1 - risk_g0,
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    plugin_risk_g1 = plugin_obj$risk_g1,
    plugin_risk_g0 = plugin_obj$risk_g0,
    residual_correction_g1 = correction_obj$correction_g1,
    residual_correction_g0 = correction_obj$correction_g0,
    n_rows_stte = n_rows_stte,
    n_rows_target = n_rows_target,
    n_rows_adherent = n_rows_adherent,
    n_rows_outcome = n_rows_outcome,
    weight_mean = weight_summary$weight_mean,
    weight_max = weight_summary$weight_max,
    weight_cv = weight_summary$weight_cv,
    n_target_missing_cells = n_target_missing_cells,
    n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells,
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    failure_reason = NA_character_
  )
}

.dr_mdr_cell_support_counts <- function(weight_obj) {
  n_target_missing_cells <- if (is.null(weight_obj$cell_support) ||
      nrow(weight_obj$cell_support) == 0L) {
    NA_integer_
  } else {
    sum(weight_obj$cell_support$target_missing_but_adherent_present, na.rm = TRUE)
  }

  n_target_present_adherent_missing_cells <- if (is.null(weight_obj$cell_support) ||
      nrow(weight_obj$cell_support) == 0L ||
      !"target_present_but_adherent_missing" %in% names(weight_obj$cell_support)) {
    NA_integer_
  } else {
    sum(weight_obj$cell_support$target_present_but_adherent_missing, na.rm = TRUE)
  }

  list(
    n_target_missing_cells = n_target_missing_cells,
    n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells
  )
}

estimate_STTE_per_protocol_dr_ipcw_pool_from_weighted_stte <- function(
    stte_ipcw,
    params,
    ipcw_obj = NULL,
    return_models = FALSE) {

  if (!"ipcw" %in% names(stte_ipcw)) {
    return(.dr_failure(
      "Stabilized IPCW column is missing from weighted STTE data"
    ))
  }

  outcome_data <- stte_ipcw %>%
    dplyr::filter(
      artificial_censor == 0L,
      !is.na(ipcw),
      is.finite(ipcw)
    ) %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  if (nrow(outcome_data) == 0L || length(unique(outcome_data$dY_next)) < 2L) {
    return(.dr_failure("No outcome variation or no usable DR-IPCW pooled outcome rows"))
  }

  outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
    X1 + X2 + X3 + L_trial

  outcome_fit <- .fit_logistic_safe(
    formula = outcome_formula,
    data = outcome_data
  )

  if (is.null(outcome_fit)) {
    return(.dr_failure("DR-IPCW pooled outcome model failed"))
  }

  eligible_base <- stte_ipcw %>%
    dplyr::filter(follow_k == 0L) %>%
    dplyr::distinct(id_trial, trial_m, X1, X2, X3, L_trial)

  plugin_obj <- .dr_estimate_baseline_plugin_risk(
    outcome_fit = outcome_fit,
    eligible_base = eligible_base,
    params = params,
    fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE),
    model_strategy_col = "strategy"
  )

  residual_data <- .dr_add_aipw_residual(
    data = outcome_data,
    outcome_fit = outcome_fit,
    params = params,
    fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE),
    path_cols = c("strategy", "id_trial"),
    path_data = stte_ipcw,
    join_cols = c("id_trial", "strategy", "follow_k")
  )
  bad_vnext_rows <- .dr_bad_vnext_rows(residual_data, params)
  if (nrow(bad_vnext_rows) > 0L) {
    return(.dr_failure(.dr_bad_vnext_failure_reason(bad_vnext_rows)))
  }

  correction_obj <- .dr_weighted_residual_correction(
    residual_data = residual_data,
    weight_col = "ipcw",
    baseline_counts = .dr_constant_baseline_counts(nrow(eligible_base)),
    correction_scale = "baseline_total",
    strategy_col = "strategy"
  )

  if (isTRUE(correction_obj$failed)) {
    return(.dr_failure(correction_obj$failure_reason))
  }

  out <- .dr_make_estimate(
    plugin_obj = plugin_obj,
    correction_obj = correction_obj,
    n_rows_stte = nrow(stte_ipcw),
    n_rows_outcome = nrow(outcome_data),
    weight_summary = .dr_weight_summary(outcome_data, "ipcw")
  )

  if (return_models) {
    out$ipcw_data <- stte_ipcw
    if (!is.null(ipcw_obj)) {
      out$denominator_fit <- ipcw_obj$denominator_fit
      out$numerator_fit <- ipcw_obj$numerator_fit
    }
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$plugin_risk_table <- plugin_obj$risk_table
    out$plugin_hazard_table <- plugin_obj$hazard_table
    out$residual_correction_by_cell <- correction_obj$correction_by_cell
    out$residual_data <- correction_obj$residual_data
  }

  out
}

estimate_STTE_per_protocol_dr_ipcw_state_from_samples <- function(
    state_samples,
    params,
    return_models = FALSE) {

  if (!"ipcw" %in% names(state_samples$adherent_states_all)) {
    return(.dr_failure(
      "Stabilized IPCW column is missing from adherent state rows"
    ))
  }

  outcome_data <- state_samples$adherent_states_all %>%
    dplyr::filter(
      !is.na(ipcw),
      is.finite(ipcw)
    ) %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  bad_target_outcome_cells <- .check_target_outcome_cell_support(
    target_states = state_samples$target_states_all,
    outcome_data = outcome_data
  )
  n_target_without_observed_outcome_cells <- nrow(bad_target_outcome_cells)
  if (n_target_without_observed_outcome_cells > 0L) {
    return(.dr_failure(
      .target_outcome_support_failure_reason(bad_target_outcome_cells),
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells
    ))
  }

  if (nrow(outcome_data) == 0L || length(unique(outcome_data$dY_next)) < 2L) {
    return(.dr_failure(
      "No outcome variation or no usable DR-IPCW state-gcomp outcome rows",
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells
    ))
  }

  outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
    X1 + X2 + X3 + L0 + L_curr

  outcome_fit <- .fit_logistic_safe(
    formula = outcome_formula,
    data = outcome_data
  )

  if (is.null(outcome_fit)) {
    return(.dr_failure(
      "DR-IPCW state-gcomp outcome model failed",
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells
    ))
  }

  fallback_hazard <- mean(outcome_data$dY_next, na.rm = TRUE)
  plugin_obj <- .dr_estimate_state_plugin_risk(
    outcome_fit = outcome_fit,
    target_states = state_samples$target_states_all,
    params = params,
    fallback_hazard = fallback_hazard
  )

  residual_data <- .dr_add_hazard_residual(
    data = outcome_data,
    outcome_fit = outcome_fit,
    params = params,
    fallback_hazard = fallback_hazard
  )

  correction_obj <- .dr_weighted_hazard_residual_correction(
    residual_data = residual_data,
    weight_col = "ipcw",
    target_hazard_table = plugin_obj$hazard_table,
    strategy_col = "strategy"
  )

  if (isTRUE(correction_obj$failed)) {
    return(.dr_failure(
      correction_obj$failure_reason,
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells
    ))
  }

  out <- .dr_make_hazard_corrected_estimate(
    plugin_obj = plugin_obj,
    correction_obj = correction_obj,
    n_rows_stte = nrow(state_samples$stte_state_data),
    n_rows_target = nrow(state_samples$target_states_all),
    n_rows_adherent = nrow(state_samples$adherent_states_all),
    n_rows_outcome = nrow(outcome_data),
    weight_summary = .dr_weight_summary(outcome_data, "ipcw"),
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    include_tables = return_models
  )

  if (return_models) {
    out$state_samples <- state_samples
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$plugin_risk_table <- plugin_obj$risk_table
    out$plugin_hazard_table <- plugin_obj$hazard_table
    out$residual_correction_by_cell <- correction_obj$correction_by_cell
    out$residual_data <- correction_obj$residual_data
    out$target_without_observed_outcome_cells <- bad_target_outcome_cells
  }

  out
}

estimate_STTE_per_protocol_dr_mdr_from_weight_obj <- function(
    state_samples,
    weight_obj,
    params,
    estimator = c("state_gcomp", "pool"),
    return_models = FALSE) {

  estimator <- match.arg(estimator)
  support_counts <- .dr_mdr_cell_support_counts(weight_obj)

  make_mdr_failure <- function(
      failure_reason,
      n_target_without_observed_outcome_cells = NA_integer_) {
    .dr_failure(
      failure_reason = failure_reason,
      n_target_missing_cells = support_counts$n_target_missing_cells,
      n_target_present_adherent_missing_cells =
        support_counts$n_target_present_adherent_missing_cells,
      n_target_without_observed_outcome_cells =
        n_target_without_observed_outcome_cells
    )
  }

  outcome_data <- weight_obj$weighted_rows %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  weight_col <- "w_mdr"

  outcome_data <- outcome_data[
    !is.na(outcome_data[[weight_col]]) &
      is.finite(outcome_data[[weight_col]]),
    ,
    drop = FALSE
  ]

  bad_target_outcome_cells <- if (estimator == "state_gcomp") {
    .check_target_outcome_cell_support(
      target_states = state_samples$target_states_all,
      outcome_data = outcome_data
    )
  } else {
    data.frame(strategy = integer(), follow_k = integer())
  }
  n_target_without_observed_outcome_cells <- if (estimator == "state_gcomp") {
    nrow(bad_target_outcome_cells)
  } else {
    NA_integer_
  }

  if (nrow(bad_target_outcome_cells) > 0L) {
    return(make_mdr_failure(
      .target_outcome_support_failure_reason(bad_target_outcome_cells),
      n_target_without_observed_outcome_cells
    ))
  }

  if (nrow(outcome_data) == 0L || length(unique(outcome_data$dY_next)) < 2L) {
    return(make_mdr_failure(
      paste0("No outcome variation or no usable DR-MDR ", estimator, " outcome rows"),
      n_target_without_observed_outcome_cells
    ))
  }

  if (estimator == "pool") {
    outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
      X1 + X2 + X3 + L_trial
  } else {
    outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
      X1 + X2 + X3 + L0 + L_curr
  }

  outcome_fit <- .fit_logistic_safe(
    formula = outcome_formula,
    data = outcome_data
  )

  if (is.null(outcome_fit)) {
    return(make_mdr_failure(
      paste0("DR-MDR ", estimator, " outcome model failed"),
      n_target_without_observed_outcome_cells
    ))
  }

  fallback_hazard <- mean(outcome_data$dY_next, na.rm = TRUE)
  if (estimator == "pool") {
    plugin_obj <- .dr_estimate_baseline_plugin_risk(
      outcome_fit = outcome_fit,
      eligible_base = state_samples$eligible_base,
      params = params,
      fallback_hazard = fallback_hazard,
      model_strategy_col = "strategy"
    )
  } else {
    plugin_obj <- .dr_estimate_state_plugin_risk(
      outcome_fit = outcome_fit,
      target_states = state_samples$target_states_all,
      params = params,
      fallback_hazard = fallback_hazard
    )
  }

  if (estimator == "state_gcomp") {
    residual_data <- .dr_add_hazard_residual(
      data = outcome_data,
      outcome_fit = outcome_fit,
      params = params,
      fallback_hazard = fallback_hazard
    )

    correction_obj <- .dr_weighted_hazard_residual_correction(
      residual_data = residual_data,
      weight_col = weight_col,
      target_hazard_table = plugin_obj$hazard_table,
      strategy_col = "strategy"
    )

    if (isTRUE(correction_obj$failed)) {
      return(make_mdr_failure(
        correction_obj$failure_reason,
        n_target_without_observed_outcome_cells
      ))
    }

    out <- .dr_make_hazard_corrected_estimate(
      plugin_obj = plugin_obj,
      correction_obj = correction_obj,
      n_rows_stte = nrow(state_samples$stte_state_data),
      n_rows_target = nrow(state_samples$target_states_all),
      n_rows_adherent = nrow(state_samples$adherent_states_all),
      n_rows_outcome = nrow(outcome_data),
      weight_summary = .dr_weight_summary(outcome_data, weight_col),
      n_target_missing_cells = support_counts$n_target_missing_cells,
      n_target_present_adherent_missing_cells =
        support_counts$n_target_present_adherent_missing_cells,
      n_target_without_observed_outcome_cells =
        n_target_without_observed_outcome_cells,
      include_tables = return_models
    )
  } else {
    residual_data <- .dr_add_aipw_residual(
      data = outcome_data,
      outcome_fit = outcome_fit,
      params = params,
      fallback_hazard = fallback_hazard,
      path_cols = c("strategy", "id_trial"),
      path_data = state_samples$stte_state_data,
      join_cols = c("id_trial", "strategy", "follow_k")
    )
    bad_vnext_rows <- .dr_bad_vnext_rows(residual_data, params)
    if (nrow(bad_vnext_rows) > 0L) {
      return(make_mdr_failure(
        .dr_bad_vnext_failure_reason(bad_vnext_rows),
        n_target_without_observed_outcome_cells
      ))
    }

    target_cell_mass <- try(
      .dr_target_cell_mass(
        state_samples$target_states_all,
        target_deplete_events = state_samples$target_deplete_events
      ),
      silent = TRUE
    )
    if (inherits(target_cell_mass, "try-error")) {
      target_cell_mass_failure <- conditionMessage(attr(target_cell_mass, "condition"))
      if (is.null(target_cell_mass_failure) || length(target_cell_mass_failure) == 0L) {
        target_cell_mass_failure <- as.character(target_cell_mass)
      }
      return(make_mdr_failure(
        target_cell_mass_failure,
        n_target_without_observed_outcome_cells
      ))
    }

    correction_obj <- .dr_weighted_residual_correction(
      residual_data = residual_data,
      weight_col = weight_col,
      target_cell_mass = target_cell_mass,
      correction_scale = "target_cell_mass",
      strategy_col = "strategy"
    )

    if (isTRUE(correction_obj$failed)) {
      return(make_mdr_failure(
        correction_obj$failure_reason,
        n_target_without_observed_outcome_cells
      ))
    }

    out <- .dr_make_estimate(
      plugin_obj = plugin_obj,
      correction_obj = correction_obj,
      n_rows_stte = nrow(state_samples$stte_state_data),
      n_rows_target = nrow(state_samples$target_states_all),
      n_rows_adherent = nrow(state_samples$adherent_states_all),
      n_rows_outcome = nrow(outcome_data),
      weight_summary = .dr_weight_summary(outcome_data, weight_col),
      n_target_missing_cells = support_counts$n_target_missing_cells,
      n_target_present_adherent_missing_cells =
        support_counts$n_target_present_adherent_missing_cells,
      n_target_without_observed_outcome_cells =
        n_target_without_observed_outcome_cells
    )
  }

  if (return_models) {
    out$state_samples <- state_samples
    out$weighted_rows <- weight_obj$weighted_rows
    out$ratio_fits <- weight_obj$ratio_fits
    out$cell_support <- weight_obj$cell_support
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$plugin_risk_table <- plugin_obj$risk_table
    out$plugin_hazard_table <- plugin_obj$hazard_table
    out$residual_correction_by_cell <- correction_obj$correction_by_cell
    out$residual_data <- correction_obj$residual_data
    out$target_without_observed_outcome_cells <- bad_target_outcome_cells
  }

  out
}
