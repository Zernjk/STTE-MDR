############################################################
# Simulation study for the ten estimators reported in the paper
#
# This driver uses:
#   - stte_per_protocol_simulation.R for the DGP, non-MDR estimators,
#     performance summaries, and existing visualization functions.
#   - mdr_weight.R for the follow-balanced, cell-normalized MDR weight
#     calculation, with the MDR ratio classifier set to random forest by
#     default in this driver.
#
# The public method registry below is intentionally restricted to the exact
# comparator set used in the manuscript. Pass comparing_methods = c(...) to
# run a subset of these estimators.
############################################################

.all_methods_script_dir <- local({
  frame_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
  if (is.character(frame_file) && length(frame_file) == 1L && !is.na(frame_file)) {
    return(dirname(normalizePath(frame_file, winslash = "/", mustWork = TRUE)))
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)))
  }

  getwd()
})

source(file.path(.all_methods_script_dir, "stte_per_protocol_simulation.R"))
source(file.path(.all_methods_script_dir, "mdr_weight.R"))
source(file.path(.all_methods_script_dir, "stte_doubly_robust_helpers.R"))
source(file.path(.all_methods_script_dir, "ltmle.R"))

mdr_weight_methods_for_paper <- c(
  follow_balanced_cell_normalized = "follow_balanced_cell_normalized"
)

paper_compare_methods <- c(
  "ipcw_pool",
  "dr_ipcw_pool",
  "ipcw_state_gcomp",
  "dr_ipcw_state_gcomp",
  "mdr_state_gcomp_follow_balanced_cell_normalized",
  "mdr_pool_follow_balanced_cell_normalized",
  "dr_mdr_state_gcomp_follow_balanced_cell_normalized",
  "dr_mdr_pool_follow_balanced_cell_normalized",
  "trial_emulation",
  "ltmle_tmle_glm"
)

all_methods_mdr_solutions_method_names <- function(
    include_trial_emulation = TRUE,
    estimate_trial_emulation = include_trial_emulation) {

  include_trial_emulation <- isTRUE(estimate_trial_emulation)
  if (include_trial_emulation) {
    return(paper_compare_methods)
  }
  setdiff(paper_compare_methods, "trial_emulation")
}

.validate_all_methods_comparing_methods <- function(
    comparing_methods = NULL,
    estimate_trial_emulation = TRUE) {

  all_valid_methods <- all_methods_mdr_solutions_method_names(
    estimate_trial_emulation = TRUE
  )

  if (is.null(comparing_methods)) {
    return(all_methods_mdr_solutions_method_names(
      estimate_trial_emulation = estimate_trial_emulation
    ))
  }

  if (!is.character(comparing_methods)) {
    stop(
      "comparing_methods must be NULL or a character vector of method names.",
      call. = FALSE
    )
  }

  if (length(comparing_methods) == 0L) {
    stop("comparing_methods must include at least one method.", call. = FALSE)
  }

  if (any(is.na(comparing_methods)) || any(!nzchar(comparing_methods))) {
    stop("comparing_methods cannot contain missing or empty names.", call. = FALSE)
  }

  comparing_methods <- unique(comparing_methods)
  unknown_methods <- setdiff(comparing_methods, all_valid_methods)
  if (length(unknown_methods) > 0L) {
    stop(
      "Unknown comparing_methods: ",
      paste(unknown_methods, collapse = ", "),
      ". Available methods: ",
      paste(all_valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  comparing_methods
}

.estimate_STTE_per_protocol_ipcw_pool_from_weighted_stte <- function(
    stte_ipcw,
    params,
    ipcw_obj = NULL,
    return_models = FALSE) {

  if (!"ipcw" %in% names(stte_ipcw)) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      failure_reason = "IPCW column is missing from weighted STTE data"
    ))
  }

  # This is the pooled IPCW estimator from
  # estimate_STTE_per_protocol_ipcw_pool(), but it starts from precomputed IPCW
  # weights so a replication fits censoring models only once.
  outcome_data <- stte_ipcw %>%
    dplyr::filter(
      artificial_censor == 0L,
      !is.na(ipcw),
      is.finite(ipcw)
    ) %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  if (nrow(outcome_data) == 0 || length(unique(outcome_data$dY_next)) < 2) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      failure_reason = "No outcome variation or no usable outcome rows"
    ))
  }

  outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
    X1 + X2 + X3 + L_trial

  outcome_fit <- .fit_logistic_safe(
    formula = outcome_formula,
    data = outcome_data,
    weights = outcome_data$ipcw
  )

  if (is.null(outcome_fit)) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      failure_reason = "Outcome model failed"
    ))
  }

  risk_obj <- .estimate_clone_standardized_risks(
    outcome_fit = outcome_fit,
    stte_data = stte_ipcw,
    params = params,
    fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE)
  )

  weight_mean <- mean(outcome_data$ipcw, na.rm = TRUE)

  out <- list(
    psi_hat = risk_obj$psi_hat,
    risk_g1 = risk_obj$risk_g1,
    risk_g0 = risk_obj$risk_g0,
    n_rows_stte = nrow(stte_ipcw),
    n_rows_outcome = nrow(outcome_data),
    weight_mean = weight_mean,
    weight_max = max(outcome_data$ipcw, na.rm = TRUE),
    weight_cv = ifelse(
      abs(weight_mean) < 1e-12,
      NA_real_,
      stats::sd(outcome_data$ipcw, na.rm = TRUE) / weight_mean
    ),
    failure_reason = NA_character_
  )

  if (return_models) {
    out$ipcw_data <- stte_ipcw
    if (!is.null(ipcw_obj)) {
      out$denominator_fit <- ipcw_obj$denominator_fit
      out$numerator_fit <- ipcw_obj$numerator_fit
    }
    out$outcome_fit <- outcome_fit
  }

  out
}

.estimate_STTE_per_protocol_ipcw_state_from_samples <- function(
    state_samples,
    params,
    return_models = FALSE) {

  if (!"ipcw" %in% names(state_samples$adherent_states_all)) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      failure_reason = "IPCW column is missing from adherent state rows"
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
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
      failure_reason = .target_outcome_support_failure_reason(bad_target_outcome_cells)
    ))
  }

  if (nrow(outcome_data) == 0 || length(unique(outcome_data$dY_next)) < 2) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
      failure_reason = "No outcome variation or no usable IPCW state-gcomp outcome rows"
    ))
  }

  outcome_formula <- dY_next ~ strategy + follow_f + trial_m +
    X1 + X2 + X3 + L0 + L_curr
  outcome_fit <- .fit_logistic_safe(
    formula = outcome_formula,
    data = outcome_data,
    weights = outcome_data$ipcw
  )

  if (is.null(outcome_fit)) {
    return(list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
      failure_reason = "IPCW state-gcomp outcome model failed"
    ))
  }

  risk_obj <- .estimate_state_standardized_risks(
    outcome_fit = outcome_fit,
    params = params,
    fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE),
    target_states = state_samples$target_states_all
  )

  weight_mean <- mean(outcome_data$ipcw, na.rm = TRUE)
  weight_max <- max(outcome_data$ipcw, na.rm = TRUE)
  weight_cv <- ifelse(
    abs(weight_mean) < 1e-12,
    NA_real_,
    stats::sd(outcome_data$ipcw, na.rm = TRUE) / weight_mean
  )

  out <- list(
    psi_hat = risk_obj$psi_hat,
    risk_g1 = risk_obj$risk_g1,
    risk_g0 = risk_obj$risk_g0,
    n_rows_stte = nrow(state_samples$stte_state_data),
    n_rows_target = nrow(state_samples$target_states_all),
    n_rows_adherent = nrow(state_samples$adherent_states_all),
    n_rows_outcome = nrow(outcome_data),
    weight_mean = weight_mean,
    weight_max = weight_max,
    weight_cv = weight_cv,
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    failure_reason = NA_character_
  )

  if (return_models) {
    out$state_samples <- state_samples
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$risk_table <- risk_obj$risk_table
    out$hazard_table <- risk_obj$hazard_table
    out$target_without_observed_outcome_cells <- bad_target_outcome_cells
  }

  out
}

.estimate_STTE_per_protocol_mdr_from_weight_obj <- function(
    state_samples,
    weight_obj,
    params,
    estimator = c("state_gcomp", "pool"),
    return_models = FALSE) {

  estimator <- match.arg(estimator)

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

  make_mdr_failure <- function(
      failure_reason,
      n_target_without_observed_outcome_cells = NA_integer_) {
    list(
      psi_hat = NA_real_,
      risk_g1 = NA_real_,
      risk_g0 = NA_real_,
      n_target_missing_cells = n_target_missing_cells,
      n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells,
      n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
      failure_reason = failure_reason
    )
  }

  outcome_data <- weight_obj$weighted_rows %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

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

  if (nrow(outcome_data) == 0 || length(unique(outcome_data$dY_next)) < 2) {
    return(make_mdr_failure(
      paste0("No outcome variation or no usable MDR ", estimator, " outcome rows"),
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
    data = outcome_data,
    weights = outcome_data$w_mdr
  )

  if (is.null(outcome_fit)) {
    return(make_mdr_failure(
      paste0("MDR ", estimator, " outcome model failed"),
      n_target_without_observed_outcome_cells
    ))
  }

  if (estimator == "pool") {
    risk_obj <- .estimate_strategy_baseline_standardized_risks(
      outcome_fit = outcome_fit,
      eligible_base = state_samples$eligible_base,
      params = params,
      fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE)
    )
  } else {
    risk_obj <- .estimate_state_standardized_risks(
      outcome_fit = outcome_fit,
      params = params,
      fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE),
      target_states = state_samples$target_states_all
    )
  }

  weight_mean <- mean(outcome_data$w_mdr, na.rm = TRUE)

  out <- list(
    psi_hat = risk_obj$psi_hat,
    risk_g1 = risk_obj$risk_g1,
    risk_g0 = risk_obj$risk_g0,
    n_rows_stte = nrow(state_samples$stte_state_data),
    n_rows_target = nrow(state_samples$target_states_all),
    n_rows_adherent = nrow(state_samples$adherent_states_all),
    n_rows_outcome = nrow(outcome_data),
    weight_mean = weight_mean,
    weight_max = max(outcome_data$w_mdr, na.rm = TRUE),
    weight_cv = ifelse(
      abs(weight_mean) < 1e-12,
      NA_real_,
      stats::sd(outcome_data$w_mdr, na.rm = TRUE) / weight_mean
    ),
    n_target_missing_cells = n_target_missing_cells,
    n_target_present_adherent_missing_cells = n_target_present_adherent_missing_cells,
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    failure_reason = NA_character_
  )

  if (return_models) {
    out$state_samples <- state_samples
    out$weighted_rows <- weight_obj$weighted_rows
    out$ratio_fits <- weight_obj$ratio_fits
    out$cell_support <- weight_obj$cell_support
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$risk_table <- risk_obj$risk_table
    out$hazard_table <- risk_obj$hazard_table
    out$target_without_observed_outcome_cells <- bad_target_outcome_cells
  }

  out
}

.failed_estimate <- function(err) {
  list(
    psi_hat = NA_real_,
    risk_g1 = NA_real_,
    risk_g0 = NA_real_,
    failure_reason = as.character(err)
  )
}

.try_estimate <- function(expr) {
  est <- try(expr, silent = TRUE)
  if (inherits(est, "try-error")) {
    return(.failed_estimate(est))
  }
  est
}

run_one_replication_all_methods_mdr_solutions <- function(
    b,
    n,
    params,
    seed_base = 10000,
    direct_mdr_target_sim_multiplier = 5,
    direct_mdr_target_sim_size = NULL,
    direct_mdr_eps = 0.001,
    direct_mdr_weight_cap = NULL,
    direct_mdr_target_deplete_events = TRUE,
    direct_mdr_classifier_method = "random_forest",
    direct_mdr_random_forest_ntree = 500,
    direct_mdr_random_forest_mtry = NULL,
    direct_mdr_random_forest_min_node_size = NULL,
    direct_mdr_random_forest_engine = "randomForest",
    comparing_methods = NULL,
    ltmle_gbounds = c(0.01, 1),
    estimate_trial_emulation = TRUE,
    trial_emulation_chunk_size = 500,
    trial_emulation_conf_int = FALSE,
    trial_emulation_samples = 100,
    trial_emulation_model_dir = file.path(tempdir(), "trial_emulation_simulation_models")) {

  direct_mdr_classifier_method <- validate_mdr_weight_classifier_method(
    direct_mdr_classifier_method
  )
  direct_mdr_random_forest_engine <- match.arg(
    direct_mdr_random_forest_engine,
    c("randomForest", "ranger")
  )
  comparing_methods <- .validate_all_methods_comparing_methods(
    comparing_methods = comparing_methods,
    estimate_trial_emulation = estimate_trial_emulation
  )
  if ("ltmle_tmle_glm" %in% comparing_methods) {
    ltmle_gbounds <- validate_ltmle_gbounds(ltmle_gbounds)
  }
  if (identical(direct_mdr_classifier_method, "random_forest") &&
      any(grepl("^(mdr_|dr_mdr_)", comparing_methods))) {
    direct_mdr_random_forest_engine <- .mdr_weight_random_forest_engine(
      direct_mdr_random_forest_engine
    )
  }

  want_method <- function(method) {
    method %in% comparing_methods
  }

  obs_data <- generate_observed_data(
    n = n,
    params = params,
    seed = seed_base + b,
    simulate_censoring = TRUE
  )

  results <- list()

  needs_stte <- any(!comparing_methods %in% c(
    "ltmle_tmle_glm",
    "trial_emulation"
  ))
  needs_ipcw <- any(comparing_methods %in% c(
    "ipcw_pool",
    "ipcw_state_gcomp",
    "dr_ipcw_pool",
    "dr_ipcw_state_gcomp"
  ))
  needs_state_samples <- any(comparing_methods %in% c(
    "ipcw_state_gcomp",
    "dr_ipcw_state_gcomp"
  )) || any(grepl("^(mdr_|dr_mdr_)", comparing_methods))

  stte_data <- NULL
  if (isTRUE(needs_stte)) {
    stte_data <- construct_STTE_dataset(obs_data, params)
  }

  ipcw_obj <- NULL
  ipcw_available <- FALSE
  stte_ipcw <- stte_data
  if (isTRUE(needs_ipcw)) {
    # Fit IPCW once for all selected methods that need censoring weights.
    ipcw_obj <- try(estimate_ipcw(stte_data, params), silent = TRUE)
    ipcw_available <- !inherits(ipcw_obj, "try-error")
    stte_ipcw <- if (ipcw_available) ipcw_obj$data else stte_data
  }

  if (want_method("ipcw_pool")) {
    ipcw_pool_est <- if (ipcw_available) {
      .try_estimate(
        .estimate_STTE_per_protocol_ipcw_pool_from_weighted_stte(
          stte_ipcw = stte_ipcw,
          params = params,
          ipcw_obj = ipcw_obj,
          return_models = FALSE
        )
      )
    } else {
      .failed_estimate(ipcw_obj)
    }
    results[["ipcw_pool"]] <- .make_replication_result_row(
      replication = b,
      method = "ipcw_pool",
      est = ipcw_pool_est,
      n_rows_input = nrow(obs_data)
    )
  }

  if (want_method("dr_ipcw_pool")) {
    dr_ipcw_pool_est <- if (ipcw_available) {
      .try_estimate(
        estimate_STTE_per_protocol_dr_ipcw_pool_from_weighted_stte(
          stte_ipcw = stte_ipcw,
          params = params,
          ipcw_obj = ipcw_obj,
          return_models = FALSE
        )
      )
    } else {
      .failed_estimate(ipcw_obj)
    }
    results[["dr_ipcw_pool"]] <- .make_replication_result_row(
      replication = b,
      method = "dr_ipcw_pool",
      est = dr_ipcw_pool_est,
      n_rows_input = nrow(obs_data)
    )
  }

  state_samples <- NULL
  if (isTRUE(needs_state_samples)) {
    # Build the target and observed-adherent state samples once for all selected
    # state-based and MDR estimators.
    state_samples <- try(
      .prepare_direct_mdr_state_samples(
        obs_data = obs_data,
        params = params,
        stte_data = stte_ipcw,
        target_sim_multiplier = direct_mdr_target_sim_multiplier,
        target_sim_size = direct_mdr_target_sim_size,
        target_deplete_events = direct_mdr_target_deplete_events,
        seed = seed_base + 200000L + b
      ),
      silent = TRUE
    )
  }

  if (isTRUE(needs_state_samples) && inherits(state_samples, "try-error")) {
    failed_state_est <- .failed_estimate(state_samples)

    if (want_method("ipcw_state_gcomp")) {
      results[["ipcw_state_gcomp"]] <- .make_replication_result_row(
        replication = b,
        method = "ipcw_state_gcomp",
        est = failed_state_est,
        n_rows_input = nrow(obs_data)
      )
    }

    if (want_method("dr_ipcw_state_gcomp")) {
      results[["dr_ipcw_state_gcomp"]] <- .make_replication_result_row(
        replication = b,
        method = "dr_ipcw_state_gcomp",
        est = failed_state_est,
        n_rows_input = nrow(obs_data)
      )
    }

    for (weight_method in names(mdr_weight_methods_for_paper)) {
      suffix <- unname(mdr_weight_methods_for_paper[[weight_method]])

      state_method_name <- paste0("mdr_state_gcomp_", suffix)
      dr_state_method_name <- paste0("dr_mdr_state_gcomp_", suffix)
      if (want_method(state_method_name)) {
        results[[state_method_name]] <- .make_replication_result_row(
          replication = b,
          method = state_method_name,
          est = failed_state_est,
          n_rows_input = nrow(obs_data)
        )
      }
      if (want_method(dr_state_method_name)) {
        results[[dr_state_method_name]] <- .make_replication_result_row(
          replication = b,
          method = dr_state_method_name,
          est = failed_state_est,
          n_rows_input = nrow(obs_data)
        )
      }

      pool_method_name <- paste0("mdr_pool_", suffix)
      dr_pool_method_name <- paste0("dr_mdr_pool_", suffix)
      if (want_method(pool_method_name)) {
        results[[pool_method_name]] <- .make_replication_result_row(
          replication = b,
          method = pool_method_name,
          est = failed_state_est,
          n_rows_input = nrow(obs_data)
        )
      }
      if (want_method(dr_pool_method_name)) {
        results[[dr_pool_method_name]] <- .make_replication_result_row(
          replication = b,
          method = dr_pool_method_name,
          est = failed_state_est,
          n_rows_input = nrow(obs_data)
        )
      }

    }
  } else if (isTRUE(needs_state_samples)) {
    if (want_method("ipcw_state_gcomp")) {
      ipcw_state_gcomp_est <- if (ipcw_available) {
        .try_estimate(
          .estimate_STTE_per_protocol_ipcw_state_from_samples(
            state_samples = state_samples,
            params = params,
            return_models = FALSE
          )
        )
      } else {
        .failed_estimate(ipcw_obj)
      }
      results[["ipcw_state_gcomp"]] <- .make_replication_result_row(
        replication = b,
        method = "ipcw_state_gcomp",
        est = ipcw_state_gcomp_est,
        n_rows_input = nrow(obs_data)
      )
    }

    if (want_method("dr_ipcw_state_gcomp")) {
      dr_ipcw_state_gcomp_est <- if (ipcw_available) {
        .try_estimate(
          estimate_STTE_per_protocol_dr_ipcw_state_from_samples(
            state_samples = state_samples,
            params = params,
            return_models = FALSE
          )
        )
      } else {
        .failed_estimate(ipcw_obj)
      }
      results[["dr_ipcw_state_gcomp"]] <- .make_replication_result_row(
        replication = b,
        method = "dr_ipcw_state_gcomp",
        est = dr_ipcw_state_gcomp_est,
        n_rows_input = nrow(obs_data)
      )
    }

    # For each requested MDR weighting variant, compute weights once, then
    # reuse the same weighted rows for selected outcome regressions.
    for (weight_method in names(mdr_weight_methods_for_paper)) {
      suffix <- unname(mdr_weight_methods_for_paper[[weight_method]])
      state_method_name <- paste0("mdr_state_gcomp_", suffix)
      pool_method_name <- paste0("mdr_pool_", suffix)
      dr_state_method_name <- paste0("dr_mdr_state_gcomp_", suffix)
      dr_pool_method_name <- paste0("dr_mdr_pool_", suffix)

      needs_weight_method <- any(c(
        state_method_name,
        pool_method_name,
        dr_state_method_name,
        dr_pool_method_name
      ) %in%
        comparing_methods)
      if (!isTRUE(needs_weight_method)) {
        next
      }

      weight_obj <- try(
        compute_mdr_weights_by_method(
          state_samples = state_samples,
          method = weight_method,
          balance_classes = TRUE,
          eps = direct_mdr_eps,
          weight_cap = direct_mdr_weight_cap,
          classifier_method = direct_mdr_classifier_method,
          random_forest_ntree = direct_mdr_random_forest_ntree,
          random_forest_mtry = direct_mdr_random_forest_mtry,
          random_forest_min_node_size = direct_mdr_random_forest_min_node_size,
          random_forest_engine = direct_mdr_random_forest_engine,
          seed = seed_base + 300000L + b
        ),
        silent = TRUE
      )

      if (want_method(state_method_name)) {
        mdr_state_est <- if (inherits(weight_obj, "try-error")) {
          .failed_estimate(weight_obj)
        } else {
          .try_estimate(
            .estimate_STTE_per_protocol_mdr_from_weight_obj(
              state_samples = state_samples,
              weight_obj = weight_obj,
              params = params,
              estimator = "state_gcomp",
              return_models = FALSE
            )
          )
        }

        results[[state_method_name]] <- .make_replication_result_row(
          replication = b,
          method = state_method_name,
          est = mdr_state_est,
          n_rows_input = nrow(obs_data)
        )
      }

      if (want_method(dr_state_method_name)) {
        dr_mdr_state_est <- if (inherits(weight_obj, "try-error")) {
          .failed_estimate(weight_obj)
        } else {
          .try_estimate(
            estimate_STTE_per_protocol_dr_mdr_from_weight_obj(
              state_samples = state_samples,
              weight_obj = weight_obj,
              params = params,
              estimator = "state_gcomp",
              return_models = FALSE
            )
          )
        }

        results[[dr_state_method_name]] <- .make_replication_result_row(
          replication = b,
          method = dr_state_method_name,
          est = dr_mdr_state_est,
          n_rows_input = nrow(obs_data)
        )
      }

      if (want_method(pool_method_name)) {
        mdr_pool_est <- if (inherits(weight_obj, "try-error")) {
          .failed_estimate(weight_obj)
        } else {
          .try_estimate(
            .estimate_STTE_per_protocol_mdr_from_weight_obj(
              state_samples = state_samples,
              weight_obj = weight_obj,
              params = params,
              estimator = "pool",
              return_models = FALSE
            )
          )
        }

        results[[pool_method_name]] <- .make_replication_result_row(
          replication = b,
          method = pool_method_name,
          est = mdr_pool_est,
          n_rows_input = nrow(obs_data)
        )
      }

      if (want_method(dr_pool_method_name)) {
        dr_mdr_pool_est <- if (inherits(weight_obj, "try-error")) {
          .failed_estimate(weight_obj)
        } else {
          .try_estimate(
            estimate_STTE_per_protocol_dr_mdr_from_weight_obj(
              state_samples = state_samples,
              weight_obj = weight_obj,
              params = params,
              estimator = "pool",
              return_models = FALSE
            )
          )
        }

        results[[dr_pool_method_name]] <- .make_replication_result_row(
          replication = b,
          method = dr_pool_method_name,
          est = dr_mdr_pool_est,
          n_rows_input = nrow(obs_data)
        )
      }

    }
  }

  if (want_method("ltmle_tmle_glm")) {
    if (!requireNamespace("ltmle", quietly = TRUE)) {
      ltmle_est <- list(
        psi_hat = NA_real_,
        risk_g1 = NA_real_,
        risk_g0 = NA_real_,
        failure_reason = "ltmle package is not installed"
      )
    } else {
      ltmle_est <- .try_estimate(
        estimate_STTE_per_protocol_ltmle(
          obs_data = obs_data,
          params = params,
          gbounds = ltmle_gbounds,
          SL.library = "glm",
          return_models = FALSE
        )
      )
    }

    results[["ltmle_tmle_glm"]] <- .make_replication_result_row(
      replication = b,
      method = "ltmle_tmle_glm",
      est = ltmle_est,
      n_rows_input = nrow(obs_data)
    )
  }

  if (want_method("trial_emulation")) {
    if (!requireNamespace("TrialEmulation", quietly = TRUE)) {
      trial_emulation_est <- list(
        psi_hat = NA_real_,
        risk_g1 = NA_real_,
        risk_g0 = NA_real_,
        failure_reason = "TrialEmulation package is not installed"
      )
    } else {
      trial_emulation_est <- .try_estimate(
        estimate_STTE_per_protocol_trial_emulation(
          obs_data = obs_data,
          params = params,
          chunk_size = trial_emulation_chunk_size,
          conf_int = trial_emulation_conf_int,
          samples = trial_emulation_samples,
          model_dir = file.path(trial_emulation_model_dir, paste0("rep_", b)),
          return_trial_sequence = FALSE
        )
      )
    }

    results[["trial_emulation"]] <- .make_replication_result_row(
      replication = b,
      method = "trial_emulation",
      est = trial_emulation_est,
      n_rows_input = nrow(obs_data)
    )
  }

  missing_results <- setdiff(comparing_methods, names(results))
  if (length(missing_results) > 0L) {
    stop(
      "No replication results were produced for selected methods: ",
      paste(missing_results, collapse = ", "),
      call. = FALSE
    )
  }

  dplyr::bind_rows(results[comparing_methods])
}

run_all_methods_mdr_solutions_simulation <- function(
    scenario_names = c("base", "strong_confounding", "poor_positivity"),
    n = 2000,
    n_truth = 100000,
    B = 100,
    truth_target = c("sequential", "baseline"),
    seed = 2026,
    direct_mdr_target_sim_multiplier = 5,
    direct_mdr_target_sim_size = NULL,
    direct_mdr_eps = 0.001,
    direct_mdr_weight_cap = NULL,
    direct_mdr_target_deplete_events = TRUE,
    direct_mdr_classifier_method = "random_forest",
    direct_mdr_random_forest_ntree = 500,
    direct_mdr_random_forest_mtry = NULL,
    direct_mdr_random_forest_min_node_size = NULL,
    direct_mdr_random_forest_engine = "randomForest",
    comparing_methods = NULL,
    ltmle_gbounds = c(0.01, 1),
    estimate_trial_emulation = TRUE,
    trial_emulation_chunk_size = 500,
    trial_emulation_conf_int = FALSE,
    trial_emulation_samples = 100,
    trial_emulation_model_dir = file.path(tempdir(), "trial_emulation_all_methods_mdr_solutions")) {

  truth_target <- match.arg(truth_target)
  direct_mdr_classifier_method <- validate_mdr_weight_classifier_method(
    direct_mdr_classifier_method
  )
  direct_mdr_random_forest_engine <- match.arg(
    direct_mdr_random_forest_engine,
    c("randomForest", "ranger")
  )
  comparing_methods <- .validate_all_methods_comparing_methods(
    comparing_methods = comparing_methods,
    estimate_trial_emulation = estimate_trial_emulation
  )
  if ("ltmle_tmle_glm" %in% comparing_methods) {
    ltmle_gbounds <- validate_ltmle_gbounds(ltmle_gbounds)
  }
  estimate_trial_emulation <- "trial_emulation" %in% comparing_methods
  if (identical(direct_mdr_classifier_method, "random_forest") &&
      any(grepl("^(mdr_|dr_mdr_)", comparing_methods))) {
    direct_mdr_random_forest_engine <- .mdr_weight_random_forest_engine(
      direct_mdr_random_forest_engine
    )
  }

  B <- .as_positive_integer(B, "B")

  all_outputs <- vector("list", length(scenario_names))
  names(all_outputs) <- scenario_names

  for (s in scenario_names) {
    message(
      "Running STTE method comparison scenario: ",
      s,
      " (methods: ",
      paste(comparing_methods, collapse = ", "),
      ")"
    )
    params <- define_parameters(s)

    truth <- calculate_true_effect(
      params = params,
      n_truth = n_truth,
      target = truth_target,
      seed = seed + 1000L
    )

    rep_results <- dplyr::bind_rows(lapply(seq_len(B), function(b) {
      run_one_replication_all_methods_mdr_solutions(
        b = b,
        n = n,
        params = params,
        seed_base = seed + 100000L,
        direct_mdr_target_sim_multiplier = direct_mdr_target_sim_multiplier,
        direct_mdr_target_sim_size = direct_mdr_target_sim_size,
        direct_mdr_eps = direct_mdr_eps,
        direct_mdr_weight_cap = direct_mdr_weight_cap,
        direct_mdr_target_deplete_events = direct_mdr_target_deplete_events,
        direct_mdr_classifier_method = direct_mdr_classifier_method,
        direct_mdr_random_forest_ntree = direct_mdr_random_forest_ntree,
        direct_mdr_random_forest_mtry = direct_mdr_random_forest_mtry,
        direct_mdr_random_forest_min_node_size = direct_mdr_random_forest_min_node_size,
        direct_mdr_random_forest_engine = direct_mdr_random_forest_engine,
        comparing_methods = comparing_methods,
        ltmle_gbounds = ltmle_gbounds,
        estimate_trial_emulation = estimate_trial_emulation,
        trial_emulation_chunk_size = trial_emulation_chunk_size,
        trial_emulation_conf_int = trial_emulation_conf_int,
        trial_emulation_samples = trial_emulation_samples,
        trial_emulation_model_dir = file.path(
          trial_emulation_model_dir,
          paste0("scenario_", s)
        )
      )
    }))
    rep_results$scenario <- s

    performance <- summarize_performance(
      rep_results = rep_results,
      psi_true = truth$psi_true
    )

    performance$scenario <- s
    performance$risk_g1_true <- truth$risk_g1
    performance$risk_g0_true <- truth$risk_g0
    performance$n_truth <- truth$n_truth
    performance$n_eligible_person_trials_truth <- truth$n_eligible_person_trials

    all_outputs[[s]] <- list(
      scenario = s,
      params = params,
      truth = truth,
      rep_results = rep_results,
      performance = performance
    )
  }

  list(
    settings = list(
      scenario_names = scenario_names,
      n = n,
      n_truth = n_truth,
      B = B,
      truth_target = truth_target,
      seed = seed,
      comparing_methods = comparing_methods,
      direct_mdr_target_sim_multiplier = direct_mdr_target_sim_multiplier,
      direct_mdr_target_sim_size = direct_mdr_target_sim_size,
      direct_mdr_eps = direct_mdr_eps,
      direct_mdr_weight_cap = direct_mdr_weight_cap,
      direct_mdr_target_deplete_events = direct_mdr_target_deplete_events,
      direct_mdr_classifier_method = direct_mdr_classifier_method,
      direct_mdr_random_forest_ntree = direct_mdr_random_forest_ntree,
      direct_mdr_random_forest_mtry = direct_mdr_random_forest_mtry,
      direct_mdr_random_forest_min_node_size = direct_mdr_random_forest_min_node_size,
      direct_mdr_random_forest_engine = direct_mdr_random_forest_engine,
      ltmle_gbounds = ltmle_gbounds,
      estimate_trial_emulation = estimate_trial_emulation,
      trial_emulation_chunk_size = trial_emulation_chunk_size,
      trial_emulation_conf_int = trial_emulation_conf_int,
      trial_emulation_samples = trial_emulation_samples
    ),
    outputs = all_outputs,
    performance_table = dplyr::bind_rows(lapply(all_outputs, `[[`, "performance"))
  )
}

if (sys.nframe() == 0L) {
  # Defaults are intentionally modest. Increase n, n_truth, and B for the
  # final paper-quality simulation.
  sim_all_methods_mdr_solutions <- run_all_methods_mdr_solutions_simulation(
    scenario_names = c("base", "strong_confounding", "poor_positivity"),
    n = 300,
    n_truth = 50000,
    B = 2,
    truth_target = "sequential",
    seed = 20260630,
    direct_mdr_target_sim_multiplier = 3,
    direct_mdr_target_sim_size = NULL,
    direct_mdr_eps = 0.001,
    direct_mdr_weight_cap = NULL,
    direct_mdr_target_deplete_events = TRUE,
    direct_mdr_classifier_method = "random_forest",
    direct_mdr_random_forest_ntree = 500,
    direct_mdr_random_forest_mtry = NULL,
    direct_mdr_random_forest_min_node_size = NULL,
    direct_mdr_random_forest_engine = "randomForest",
    comparing_methods = paper_compare_methods
  )

  print(sim_all_methods_mdr_solutions$performance_table)

  estimates_plot <- plot_estimates_by_scenario(sim_all_methods_mdr_solutions)
  performance_plot <- plot_performance_by_scenario(sim_all_methods_mdr_solutions)
  weight_plot <- plot_weight_behavior_by_scenario(sim_all_methods_mdr_solutions)

  print(estimates_plot)
  print(performance_plot)
  print(weight_plot)

#   ggplot2::ggsave(
#     filename = "all_methods_mdr_solutions_estimates.png",
#     plot = estimates_plot,
#     width = 11,
#     height = 6,
#     dpi = 300
#   )
#   ggplot2::ggsave(
#     filename = "all_methods_mdr_solutions_performance.png",
#     plot = performance_plot,
#     width = 11,
#     height = 6,
#     dpi = 300
#   )
#   ggplot2::ggsave(
#     filename = "all_methods_mdr_solutions_weights.png",
#     plot = weight_plot,
#     width = 11,
#     height = 7,
#     dpi = 300
#   )

#   utils::write.csv(
#     sim_all_methods_mdr_solutions$performance_table,
#     file = "all_methods_mdr_solutions_performance.csv",
#     row.names = FALSE
#   )

#   saveRDS(
#     sim_all_methods_mdr_solutions,
#     file = "all_methods_mdr_solutions_simulation.rds"
#   )
}
