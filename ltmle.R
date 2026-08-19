############################################################
# Longitudinal TMLE helper functions for the STTE simulation
#
# Primary comparison method:
#   ltmle_tmle_glm
#
# The analysis unit is one eligible person-trial. Observed histories are
# represented in wide form using the simulator's within-interval ordering:
#
#   W, L_0,
#   A_0, Y_1, C_1, L_1,
#   ...,
#   A_{tau-1}, Y_tau
#
# Y nodes are cumulative event indicators. Natural censoring occurs after the
# current interval's Y node, so there is no C_tau after the final outcome.
# Treatment-protocol adherence is handled by the A nodes and abar; cloned STTE
# rows and artificial-censoring indicators are deliberately not used here.
############################################################

.ltmle_required_observed_columns <- function() {
  c(
    "id", "k", "X1", "X2", "X3", "L_k", "A_prev", "A_k",
    "Y_start", "C_start", "dC_next", "Y_next"
  )
}

.ltmle_require_columns <- function(data, columns, data_name = "data") {
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

.ltmle_is_complete_binary <- function(values) {
  (is.numeric(values) || is.logical(values)) &&
    !anyNA(values) &&
    all(is.finite(as.numeric(values))) &&
    all(as.numeric(values) %in% c(0, 1))
}

.ltmle_validate_observed_data <- function(obs_data) {
  if (anyNA(obs_data$id)) {
    stop("obs_data$id cannot contain missing values.", call. = FALSE)
  }

  numeric_columns <- c("k", "X1", "X2", "X3", "L_k")
  invalid_numeric <- vapply(
    obs_data[numeric_columns],
    function(values) !is.numeric(values) || any(!is.finite(values)),
    logical(1)
  )
  if (any(invalid_numeric)) {
    stop(
      "obs_data must contain finite numeric values in: ",
      paste(numeric_columns[invalid_numeric], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (any(obs_data$k < 0 | obs_data$k != floor(obs_data$k))) {
    stop("obs_data$k must contain nonnegative integers.", call. = FALSE)
  }

  binary_columns <- c(
    "A_prev", "A_k", "Y_start", "C_start", "dC_next", "Y_next"
  )
  invalid_binary <- !vapply(
    obs_data[binary_columns],
    .ltmle_is_complete_binary,
    logical(1)
  )
  if (any(invalid_binary)) {
    stop(
      "obs_data must contain complete binary values in: ",
      paste(binary_columns[invalid_binary], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

validate_ltmle_gbounds <- function(gbounds = c(0.01, 1)) {
  if (!is.numeric(gbounds) || length(gbounds) != 2L ||
      any(!is.finite(gbounds))) {
    stop("ltmle gbounds must contain two finite numeric values.", call. = FALSE)
  }

  gbounds <- as.numeric(gbounds)
  if (gbounds[[1L]] <= 0 || gbounds[[2L]] > 1 ||
      gbounds[[1L]] >= gbounds[[2L]]) {
    stop(
      "ltmle gbounds must satisfy 0 < lower < upper <= 1.",
      call. = FALSE
    )
  }

  gbounds
}

make_ltmle_stte_node_spec <- function(tau) {
  if (length(tau) != 1L || !is.finite(tau) || tau < 1L ||
      tau != floor(tau)) {
    stop("tau must be a positive integer.", call. = FALSE)
  }
  tau <- as.integer(tau)

  baseline_nodes <- c("trial_m", "X1", "X2", "X3", "L_0")
  A_nodes <- paste0("A_", 0:(tau - 1L))
  Y_nodes <- paste0("Y_", seq_len(tau))
  C_nodes <- if (tau > 1L) paste0("C_", seq_len(tau - 1L)) else character()
  L_nodes <- if (tau > 1L) paste0("L_", seq_len(tau - 1L)) else character()

  ordered_nodes <- baseline_nodes
  for (j in 0:(tau - 1L)) {
    ordered_nodes <- c(ordered_nodes, A_nodes[[j + 1L]], Y_nodes[[j + 1L]])
    if (j < tau - 1L) {
      ordered_nodes <- c(
        ordered_nodes,
        C_nodes[[j + 1L]],
        L_nodes[[j + 1L]]
      )
    }
  }

  list(
    tau = tau,
    baseline_nodes = baseline_nodes,
    Anodes = A_nodes,
    Cnodes = C_nodes,
    Lnodes = L_nodes,
    Ynodes = Y_nodes,
    LYnodes = ordered_nodes[ordered_nodes %in% c(L_nodes, Y_nodes)],
    ACnodes = ordered_nodes[ordered_nodes %in% c(A_nodes, C_nodes)],
    ordered_nodes = ordered_nodes
  )
}

.ltmle_history_key <- function(id, k) {
  paste(as.character(id), as.integer(k), sep = "\r")
}

.ltmle_binary_to_censoring <- function(is_censored) {
  if (requireNamespace("ltmle", quietly = TRUE)) {
    return(ltmle::BinaryToCensoring(is.censored = is_censored))
  }

  factor(
    ifelse(
      is.na(is_censored),
      NA_character_,
      ifelse(is_censored == 1L, "censored", "uncensored")
    ),
    levels = c("censored", "uncensored")
  )
}

prepare_ltmle_stte_data <- function(obs_data, params) {
  .ltmle_require_columns(
    obs_data,
    .ltmle_required_observed_columns(),
    data_name = "obs_data"
  )
  .ltmle_validate_observed_data(obs_data)

  if (!exists("get_eligible_person_trials", mode = "function")) {
    stop(
      "get_eligible_person_trials() is unavailable. Source ",
      "stte_per_protocol_simulation.R before ltmle.R.",
      call. = FALSE
    )
  }
  if (is.null(params$tau)) {
    stop("params$tau is required for LTMLE data preparation.", call. = FALSE)
  }
  if (is.null(params$K) || length(params$K) != 1L ||
      !is.finite(params$K) || params$K < 1L || params$K != floor(params$K)) {
    stop("params$K must be a positive integer.", call. = FALSE)
  }
  if (is.null(params$grace) || length(params$grace) != 1L ||
      !is.finite(params$grace) || params$grace != 0) {
    stop(
      "The LTMLE comparison currently implements only grace = 0: ",
      "always treated versus never treated.",
      call. = FALSE
    )
  }
  if (!isTRUE(params$require_untreated_entry)) {
    stop(
      "The current LTMLE comparison requires require_untreated_entry = TRUE ",
      "so A_prev_entry is fixed at zero for every eligible person-trial.",
      call. = FALSE
    )
  }

  nodes <- make_ltmle_stte_node_spec(params$tau)
  if (params$K < nodes$tau) {
    stop("params$K must be at least params$tau.", call. = FALSE)
  }
  if (any(obs_data$k >= params$K)) {
    stop("obs_data$k must be smaller than params$K.", call. = FALSE)
  }
  eligible <- get_eligible_person_trials(obs_data, params)
  if (nrow(eligible) == 0L) {
    stop("No eligible person-trials were available for LTMLE.", call. = FALSE)
  }
  if (!.ltmle_is_complete_binary(eligible$A_prev_entry) ||
      any(eligible$A_prev_entry != 0)) {
    stop(
      "Eligible person-trials must all have A_prev_entry = 0.",
      call. = FALSE
    )
  }

  obs_key <- .ltmle_history_key(obs_data$id, obs_data$k)
  if (anyDuplicated(obs_key)) {
    stop("obs_data contains duplicate (id, k) rows.", call. = FALSE)
  }

  n_trials <- nrow(eligible)
  wide <- data.frame(
    trial_m = as.integer(eligible$trial_m),
    X1 = eligible$X1,
    X2 = eligible$X2,
    X3 = eligible$X3,
    L_0 = eligible$L_trial,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  event_seen <- rep(FALSE, n_trials)
  censor_seen <- rep(FALSE, n_trials)

  for (j in 0:(nodes$tau - 1L)) {
    calendar_k <- eligible$trial_m + j
    row_index <- match(.ltmle_history_key(eligible$id, calendar_k), obs_key)
    should_have_row <- !event_seen & !censor_seen
    missing_active <- should_have_row & is.na(row_index)
    if (any(missing_active)) {
      examples <- utils::head(eligible$id_trial[missing_active], 5L)
      stop(
        "Observed histories end before an event or censoring indicator for ",
        "eligible person-trial(s): ",
        paste(examples, collapse = ", "),
        call. = FALSE
      )
    }

    A_now <- rep(NA_integer_, n_trials)
    Y_now <- rep(NA_integer_, n_trials)
    dC_now <- rep(NA_integer_, n_trials)

    active <- should_have_row & !is.na(row_index)
    if (any(active)) {
      active_index <- row_index[active]
      A_now[active] <- as.integer(obs_data$A_k[active_index])
      Y_now[active] <- as.integer(obs_data$Y_next[active_index])
      dC_now[active] <- as.integer(obs_data$dC_next[active_index])
    }

    # Preserve the absorbing event outcome. Outcomes after censoring remain
    # missing because they were not observed.
    Y_now[event_seen] <- 1L
    Y_now[censor_seen] <- NA_integer_

    wide[[nodes$Anodes[[j + 1L]]]] <- A_now
    wide[[nodes$Ynodes[[j + 1L]]]] <- Y_now

    event_now <- active & Y_now == 1L
    censor_now <- active & !event_now & dC_now == 1L

    if (j < nodes$tau - 1L) {
      C_now <- rep(NA_integer_, n_trials)
      alive_at_censor_node <- active & !event_now
      C_now[alive_at_censor_node] <- dC_now[alive_at_censor_node]
      wide[[nodes$Cnodes[[j + 1L]]]] <- .ltmle_binary_to_censoring(C_now)

      continues <- active & !event_now & !censor_now
      L_next <- rep(NA_real_, n_trials)
      if (any(continues)) {
        next_index <- match(
          .ltmle_history_key(
            eligible$id[continues],
            calendar_k[continues] + 1L
          ),
          obs_key
        )
        if (anyNA(next_index)) {
          bad <- which(continues)[is.na(next_index)]
          examples <- utils::head(eligible$id_trial[bad], 5L)
          stop(
            "Missing next-L rows for uncensored, event-free person-trial(s): ",
            paste(examples, collapse = ", "),
            call. = FALSE
          )
        }
        L_next[continues] <- obs_data$L_k[next_index]
      }
      wide[[nodes$Lnodes[[j + 1L]]]] <- L_next
    }

    event_seen <- event_seen | event_now
    censor_seen <- censor_seen | censor_now
  }

  wide <- wide[, nodes$ordered_nodes, drop = FALSE]

  out <- list(
    data = wide,
    id = eligible$id,
    id_trial = eligible$id_trial,
    eligible = eligible,
    nodes = nodes
  )
  class(out) <- c("stte_ltmle_data", class(out))
  out
}

.ltmle_formula_string <- function(lhs, rhs) {
  rhs <- unique(rhs)
  paste(lhs, "~", paste(rhs, collapse = " + "))
}

make_ltmle_stte_formulas <- function(nodes) {
  if (is.null(nodes$tau) || is.null(nodes$ordered_nodes)) {
    stop("nodes must be created by make_ltmle_stte_node_spec().", call. = FALSE)
  }

  baseline_rhs <- nodes$baseline_nodes

  Qform <- stats::setNames(character(length(nodes$LYnodes)), nodes$LYnodes)
  for (node in nodes$LYnodes) {
    if (grepl("^Y_", node)) {
      outcome_index <- as.integer(sub("^Y_", "", node))
      j <- outcome_index - 1L
      current_L <- if (j == 0L) character() else paste0("L_", j)
      rhs <- c(baseline_rhs, current_L, paste0("A_", j))
    } else {
      L_index <- as.integer(sub("^L_", "", node))
      previous_j <- L_index - 1L
      previous_L <- if (previous_j == 0L) character() else paste0("L_", previous_j)
      rhs <- c(baseline_rhs, previous_L, paste0("A_", previous_j))
    }
    Qform[[node]] <- .ltmle_formula_string("Q.kplus1", rhs)
  }

  gform <- stats::setNames(character(length(nodes$ACnodes)), nodes$ACnodes)
  for (node in nodes$ACnodes) {
    if (grepl("^A_", node)) {
      j <- as.integer(sub("^A_", "", node))
      current_L <- if (j == 0L) character() else paste0("L_", j)
      previous_A <- if (j == 0L) character() else paste0("A_", j - 1L)
      rhs <- c(baseline_rhs, current_L, previous_A)
    } else {
      censor_index <- as.integer(sub("^C_", "", node))
      j <- censor_index - 1L
      current_L <- if (j == 0L) character() else paste0("L_", j)
      rhs <- c(baseline_rhs, current_L, paste0("A_", j))
    }
    gform[[node]] <- .ltmle_formula_string(node, rhs)
  }

  list(Qform = Qform, gform = gform)
}

.ltmle_validate_formula_order <- function(formulas, expected_nodes, data_names, type) {
  if (!is.character(formulas) || anyNA(formulas)) {
    stop(type, " must be a non-missing character vector.", call. = FALSE)
  }
  if (length(formulas) != length(expected_nodes) ||
      !identical(names(formulas), expected_nodes)) {
    stop(
      type,
      " must be named and ordered as: ",
      paste(expected_nodes, collapse = ", "),
      call. = FALSE
    )
  }

  for (node in expected_nodes) {
    if (identical(formulas[[node]], "IDENTITY")) {
      if (!identical(type, "Qform")) {
        stop("IDENTITY is allowed only in Qform.", call. = FALSE)
      }
      next
    }
    formula_obj <- try(stats::as.formula(formulas[[node]]), silent = TRUE)
    if (inherits(formula_obj, "try-error") || length(formula_obj) != 3L) {
      stop("Invalid ", type, " formula for node ", node, ".", call. = FALSE)
    }
    lhs_variables <- all.vars(formula_obj[[2L]])
    expected_lhs <- if (identical(type, "Qform")) "Q.kplus1" else node
    if (!identical(lhs_variables, expected_lhs)) {
      stop(
        "The left-hand side of ", type, " for ", node,
        " must be ", expected_lhs, ".",
        call. = FALSE
      )
    }
    if ("." %in% all.names(formula_obj[[3L]], functions = FALSE)) {
      stop(
        type, " for ", node,
        " cannot use '.' because it could include the current or future nodes.",
        call. = FALSE
      )
    }
    rhs_variables <- all.vars(formula_obj[[3L]])
    missing_variables <- setdiff(rhs_variables, data_names)
    if (length(missing_variables) > 0L) {
      stop(
        type,
        " for ", node, " references missing variable(s): ",
        paste(missing_variables, collapse = ", "),
        call. = FALSE
      )
    }
    node_position <- match(node, data_names)
    rhs_positions <- match(rhs_variables, data_names)
    if (any(rhs_positions >= node_position)) {
      bad_variables <- rhs_variables[rhs_positions >= node_position]
      stop(
        type,
        " for ", node, " uses variable(s) that do not precede the node: ",
        paste(bad_variables, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

validate_ltmle_stte_data <- function(ltmle_data, Qform = NULL, gform = NULL) {
  if (!is.list(ltmle_data) || is.null(ltmle_data$data) ||
      is.null(ltmle_data$nodes) || is.null(ltmle_data$id)) {
    stop("ltmle_data must be created by prepare_ltmle_stte_data().", call. = FALSE)
  }

  data <- ltmle_data$data
  nodes <- ltmle_data$nodes
  if (!identical(names(data), nodes$ordered_nodes)) {
    stop("LTMLE data columns are not in the required causal order.", call. = FALSE)
  }
  if (nrow(data) != length(ltmle_data$id) ||
      nrow(data) != length(ltmle_data$id_trial)) {
    stop("LTMLE row metadata are not aligned with the wide data.", call. = FALSE)
  }
  if (anyNA(data[, nodes$baseline_nodes, drop = FALSE])) {
    stop("LTMLE baseline variables cannot contain missing values.", call. = FALSE)
  }

  for (node in nodes$Anodes) {
    values <- data[[node]]
    if (any(!is.na(values) & !values %in% c(0L, 1L))) {
      stop(node, " must contain only 0, 1, or NA.", call. = FALSE)
    }
  }
  for (node in nodes$Ynodes) {
    values <- data[[node]]
    if (any(!is.na(values) & !values %in% c(0L, 1L))) {
      stop(node, " must contain only 0, 1, or NA.", call. = FALSE)
    }
  }
  for (node in nodes$Cnodes) {
    values <- data[[node]]
    observed_values <- unique(as.character(values[!is.na(values)]))
    if (!is.factor(values) ||
        any(!observed_values %in% c("censored", "uncensored"))) {
      stop(
        node,
        " must contain only censoring labels 'censored', 'uncensored', or NA.",
        call. = FALSE
      )
    }
  }

  if (!is.null(Qform)) {
    .ltmle_validate_formula_order(
      Qform,
      nodes$LYnodes,
      names(data),
      type = "Qform"
    )
  }
  if (!is.null(gform)) {
    .ltmle_validate_formula_order(
      gform,
      nodes$ACnodes,
      names(data),
      type = "gform"
    )
  }

  invisible(TRUE)
}

.ltmle_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(default)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (!is.finite(value)) default else value
}

summarize_ltmle_cum_g <- function(fit, gbounds) {
  gbounds <- validate_ltmle_gbounds(gbounds)
  lower_bound <- gbounds[[1L]]
  upper_bound <- gbounds[[2L]]
  raw <- as.numeric(fit$cum.g.unbounded)
  bounded <- as.numeric(fit$cum.g)
  used <- as.logical(fit$cum.g.used)

  raw_finite <- raw[is.finite(raw) & raw >= 0]
  bounded_finite <- bounded[is.finite(bounded) & bounded > 0]
  used_finite <- rep(FALSE, length(raw))
  if (length(used) == length(raw)) {
    used_finite <- is.finite(raw) & raw >= 0 & !is.na(used) & used
  }
  bounded_used <- bounded[used_finite & is.finite(bounded) & bounded > 0]

  was_bounded <- function(values) {
    values < lower_bound | values > upper_bound
  }

  list(
    min_cum_g_unbounded = if (length(raw_finite) > 0L) min(raw_finite) else NA_real_,
    prop_cum_g_bounded = if (length(raw_finite) > 0L) {
      mean(was_bounded(raw_finite))
    } else {
      NA_real_
    },
    max_inverse_cum_g = if (length(bounded_finite) > 0L) {
      max(1 / bounded_finite)
    } else {
      NA_real_
    },
    min_cum_g_unbounded_used = if (any(used_finite)) {
      min(raw[used_finite])
    } else {
      NA_real_
    },
    prop_cum_g_bounded_used = if (any(used_finite)) {
      mean(was_bounded(raw[used_finite]))
    } else {
      NA_real_
    },
    max_inverse_cum_g_used = if (length(bounded_used) > 0L) {
      max(1 / bounded_used)
    } else {
      NA_real_
    }
  )
}

estimate_STTE_per_protocol_ltmle <- function(
    obs_data,
    params,
    gbounds = c(0.01, 1),
    SL.library = "glm",
    Qform = NULL,
    gform = NULL,
    return_models = FALSE,
    seed = NULL) {

  # gbounds truncates the cumulative treatment-and-censoring mechanism used by
  # ltmle. It is intentionally separate from params$prob_clip, which clips the
  # simulator's and IPCW nuisance-model probabilities node by node.
  if (!requireNamespace("ltmle", quietly = TRUE)) {
    stop(
      "The ltmle package is required for the ltmle_tmle_glm estimator. ",
      "Install it with install.packages('ltmle').",
      call. = FALSE
    )
  }

  gbounds <- validate_ltmle_gbounds(gbounds)
  prepared <- prepare_ltmle_stte_data(obs_data, params)
  formulas <- make_ltmle_stte_formulas(prepared$nodes)
  if (is.null(Qform)) Qform <- formulas$Qform
  if (is.null(gform)) gform <- formulas$gform
  validate_ltmle_stte_data(prepared, Qform = Qform, gform = gform)

  if (!is.null(seed)) set.seed(seed)
  fit <- ltmle::ltmle(
    data = prepared$data,
    Anodes = prepared$nodes$Anodes,
    Cnodes = prepared$nodes$Cnodes,
    Lnodes = prepared$nodes$Lnodes,
    Ynodes = prepared$nodes$Ynodes,
    survivalOutcome = TRUE,
    Qform = Qform,
    gform = gform,
    abar = list(
      rep(1L, prepared$nodes$tau),
      rep(0L, prepared$nodes$tau)
    ),
    gbounds = gbounds,
    stratify = FALSE,
    SL.library = SL.library,
    gcomp = FALSE,
    estimate.time = FALSE,
    variance.method = "ic",
    id = prepared$id
  )

  fit_summary <- summary(fit, estimator = "tmle")
  effect_measures <- fit_summary$effect.measures
  if (is.null(effect_measures$treatment) ||
      is.null(effect_measures$control) ||
      is.null(effect_measures$ATE)) {
    stop("Unable to extract LTMLE treatment-effect summaries.", call. = FALSE)
  }

  risk_g1 <- .ltmle_scalar(effect_measures$treatment$estimate)
  risk_g0 <- .ltmle_scalar(effect_measures$control$estimate)
  psi_hat <- .ltmle_scalar(effect_measures$ATE$estimate)
  if (any(!is.finite(c(risk_g1, risk_g0, psi_hat)))) {
    stop("LTMLE returned a non-finite risk or risk difference.", call. = FALSE)
  }

  ci <- suppressWarnings(as.numeric(effect_measures$ATE$CI))
  if (length(ci) != 2L || any(!is.finite(ci))) ci <- c(NA_real_, NA_real_)
  cum_g_summary <- summarize_ltmle_cum_g(fit, gbounds)

  out <- c(
    list(
      psi_hat = psi_hat,
      risk_g1 = risk_g1,
      risk_g0 = risk_g0,
      estimate_se = .ltmle_scalar(effect_measures$ATE$std.dev),
      ci_lower = ci[[1L]],
      ci_upper = ci[[2L]],
      n_rows_stte = nrow(prepared$data),
      n_rows_outcome = sum(!is.na(prepared$data[[utils::tail(prepared$nodes$Ynodes, 1L)]])),
      n_unique_ids = length(unique(prepared$id)),
      gbounds = gbounds,
      failure_reason = NA_character_
    ),
    cum_g_summary
  )

  if (isTRUE(return_models)) {
    out$fit <- fit
    out$summary <- fit_summary
    out$prepared_data <- prepared
    out$Qform <- Qform
    out$gform <- gform
  }

  out
}
