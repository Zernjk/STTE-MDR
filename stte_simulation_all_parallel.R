############################################################
# Parallel simulation study comparing all STTE methods with MDR variants
#
# This file sources stte_simulation_all.R and runs the same all-methods
# simulation, including the paper's follow-balanced, cell-normalized MDR and
# doubly robust MDR variants, with
# replications distributed across a PSOCK cluster.
# The returned object mirrors the serial runner's main structure:
#   settings, outputs, performance_table
# and also includes a top-level rep_results convenience table.
#
# Pass comparing_methods = c(...) to run only selected methods. Use
# all_methods_mdr_solutions_method_names() to list the valid names. Include
# "ltmle_tmle_glm" or "trial_emulation" in comparing_methods to run those
# estimators directly from the observed longitudinal data.
############################################################

.parallel_all_methods_script_dir <- local({
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

.parallel_all_methods_simulation_file <- file.path(
  .parallel_all_methods_script_dir,
  "stte_simulation_all.R"
)

source(.parallel_all_methods_simulation_file)

.parallel_all_methods_progress <- function(verbose = TRUE, ...) {
  if (isTRUE(verbose)) {
    cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")
    try(flush.console(), silent = TRUE)
    try(flush(stdout()), silent = TRUE)
  }
  invisible(NULL)
}

.parallel_all_methods_default_workers <- function(max_workers = 4L) {
  n_cores <- parallel::detectCores(logical = TRUE)
  if (is.na(n_cores) || n_cores < 2L) {
    return(1L)
  }

  n_available <- max(1L, n_cores - 1L)
  if (is.null(max_workers) || identical(max_workers, Inf)) {
    return(n_available)
  }

  max_workers <- .as_positive_integer(max_workers, "max_workers")

  min(n_available, max_workers)
}

.parallel_all_methods_worker_outfile <- function(verbose_workers = FALSE) {
  if (isTRUE(verbose_workers)) {
    return("")
  }
  tempfile("stte_all_methods_parallel_workers_", fileext = ".log")
}

.parallel_all_methods_make_cluster <- function(n_workers, worker_outfile) {
  if (n_workers <= 1L) {
    return(NULL)
  }

  parallel::makeCluster(
    n_workers,
    type = "PSOCK",
    outfile = worker_outfile
  )
}

.parallel_all_methods_log_tail <- function(path, n = 30L) {
  if (!nzchar(path) || !file.exists(path)) {
    return(character())
  }

  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (length(lines) == 0L) {
    return(character())
  }
  utils::tail(lines, n)
}

.parallel_all_methods_run_one_replication <- function(task) {
  make_failure_rows <- function(failure_reason) {
    rows <- lapply(task$method_names, function(method) {
      data.frame(
        replication = task$b,
        method = method,
        psi_hat = NA_real_,
        risk_g1_hat = NA_real_,
        risk_g0_hat = NA_real_,
        estimate_se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n_rows_input = NA_integer_,
        n_rows_stte = NA_integer_,
        n_rows_outcome = NA_integer_,
        n_unique_ids = NA_integer_,
        weight_mean = NA_real_,
        weight_max = NA_real_,
        weight_cv = NA_real_,
        n_target_missing_cells = NA_integer_,
        n_target_present_adherent_missing_cells = NA_integer_,
        n_target_without_observed_outcome_cells = NA_integer_,
        min_cum_g_unbounded = NA_real_,
        prop_cum_g_bounded = NA_real_,
        max_inverse_cum_g = NA_real_,
        min_cum_g_unbounded_used = NA_real_,
        prop_cum_g_bounded_used = NA_real_,
        max_inverse_cum_g_used = NA_real_,
        failure_reason = failure_reason,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  }

  tryCatch({
    source(task$simulation_file)

    run_one_replication_all_methods_mdr_solutions(
      b = task$b,
      n = task$n,
      params = task$params,
      seed_base = task$seed_base,
      direct_mdr_target_sim_multiplier = task$direct_mdr_target_sim_multiplier,
      direct_mdr_target_sim_size = task$direct_mdr_target_sim_size,
      direct_mdr_eps = task$direct_mdr_eps,
      direct_mdr_weight_cap = task$direct_mdr_weight_cap,
      direct_mdr_target_deplete_events = task$direct_mdr_target_deplete_events,
      direct_mdr_classifier_method = task$direct_mdr_classifier_method,
      direct_mdr_random_forest_ntree = task$direct_mdr_random_forest_ntree,
      direct_mdr_random_forest_mtry = task$direct_mdr_random_forest_mtry,
      direct_mdr_random_forest_min_node_size = task$direct_mdr_random_forest_min_node_size,
      direct_mdr_random_forest_engine = task$direct_mdr_random_forest_engine,
      comparing_methods = task$comparing_methods,
      ltmle_gbounds = task$ltmle_gbounds,
      estimate_trial_emulation = task$estimate_trial_emulation,
      trial_emulation_chunk_size = task$trial_emulation_chunk_size,
      trial_emulation_conf_int = task$trial_emulation_conf_int,
      trial_emulation_samples = task$trial_emulation_samples,
      trial_emulation_model_dir = task$trial_emulation_model_dir
    )
  }, error = function(e) {
    make_failure_rows(
      paste("parallel worker replication failure:", conditionMessage(e))
    )
  })
}

run_all_methods_mdr_solutions_simulation_parallel <- function(
    scenario_names = c(
      "base",
      "small_size",
      "strong_confounding",
      "poor_positivity",
      "poor_adherence",
      "rare_event",
      "nonlinear"
    ),
    n = 2000,
    n_truth = 50000,
    B = 500,
    truth_target = c("sequential", "baseline"),
    seed = 20260630,
    direct_mdr_target_sim_multiplier = 10,
    direct_mdr_target_sim_size = NULL,
    direct_mdr_eps = 0.001,
    direct_mdr_weight_cap = NULL,
    direct_mdr_target_deplete_events = TRUE,
    direct_mdr_classifier_method = "random_forest",
    direct_mdr_random_forest_ntree = 500,
    direct_mdr_random_forest_mtry = NULL,
    direct_mdr_random_forest_min_node_size = NULL,
    direct_mdr_random_forest_engine = "ranger",
    comparing_methods = NULL,
    ltmle_gbounds = c(0.01, 1),
    estimate_trial_emulation = TRUE,
    trial_emulation_chunk_size = 500,
    trial_emulation_conf_int = FALSE,
    trial_emulation_samples = 100,
    trial_emulation_model_dir = file.path(tempdir(), "trial_emulation_all_methods_mdr_solutions_parallel"),
    n_workers = .parallel_all_methods_default_workers(max_workers = 10L),
    verbose = TRUE,
    verbose_workers = FALSE,
    simulation_file = .parallel_all_methods_simulation_file) {

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

  n_workers <- .as_positive_integer(n_workers, "n_workers")
  n_workers <- min(n_workers, B)

  simulation_file <- normalizePath(simulation_file, winslash = "/", mustWork = TRUE)
  worker_outfile <- .parallel_all_methods_worker_outfile(verbose_workers)
  method_names <- comparing_methods

  .parallel_all_methods_progress(
    verbose,
    "Starting parallel all-methods STTE simulation",
    " (workers=", n_workers,
    "; scenarios=", paste(scenario_names, collapse = ", "),
    "; B=", B,
    "; methods=", paste(comparing_methods, collapse = ", "),
    "; MDR classifier=", direct_mdr_classifier_method, ")"
  )
  if (n_workers > 1L && nzchar(worker_outfile)) {
    .parallel_all_methods_progress(
      verbose,
      "Worker output log: ", worker_outfile
    )
  }

  cl <- .parallel_all_methods_make_cluster(
    n_workers = n_workers,
    worker_outfile = worker_outfile
  )
  if (!is.null(cl)) {
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  all_outputs <- vector("list", length(scenario_names))
  names(all_outputs) <- scenario_names

  for (s in scenario_names) {
    .parallel_all_methods_progress(verbose, "Scenario ", s, ": defining parameters")
    params <- define_parameters(s)

    .parallel_all_methods_progress(
      verbose,
      "Scenario ", s, ": calculating truth",
      " (n_truth=", n_truth, ", target=", truth_target, ")"
    )
    truth <- calculate_true_effect(
      params = params,
      n_truth = n_truth,
      target = truth_target,
      seed = seed + 1000L
    )
    .parallel_all_methods_progress(
      verbose,
      "Scenario ", s, ": truth ready",
      " (psi_true=", signif(truth$psi_true, 4), ")"
    )

    scenario_model_dir <- file.path(
      trial_emulation_model_dir,
      paste0("scenario_", s)
    )

    tasks <- lapply(seq_len(B), function(b) {
      list(
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
        trial_emulation_model_dir = scenario_model_dir,
        method_names = method_names,
        simulation_file = simulation_file
      )
    })

    .parallel_all_methods_progress(
      verbose,
      "Scenario ", s, ": running ", B,
      " replications on ", n_workers, " worker(s)"
    )

    rep_list <- tryCatch({
      if (is.null(cl)) {
        lapply(tasks, .parallel_all_methods_run_one_replication)
      } else {
        parallel::parLapplyLB(cl, tasks, .parallel_all_methods_run_one_replication)
      }
    }, error = function(e) {
      log_tail <- .parallel_all_methods_log_tail(worker_outfile)
      log_message <- if (length(log_tail) > 0L) {
        paste0(
          "\n\nLast lines from worker log (", worker_outfile, "):\n",
          paste(log_tail, collapse = "\n")
        )
      } else if (nzchar(worker_outfile)) {
        paste0("\n\nWorker log path: ", worker_outfile)
      } else {
        ""
      }

      stop(
        "Parallel workers stopped while running scenario '", s, "'. ",
        "Original error: ", conditionMessage(e), "\n",
        "This often means at least one worker was killed by memory pressure. ",
        "Try n_workers = 2 or 4, smaller n/B/n_truth, smaller target_sim_size, ",
        "or lighter random-forest tuning.",
        log_message,
        call. = FALSE
      )
    })

    rep_results <- dplyr::bind_rows(rep_list)
    rep_results$scenario <- s

    .parallel_all_methods_progress(verbose, "Scenario ", s, ": summarizing performance")
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

    .parallel_all_methods_progress(verbose, "Scenario ", s, ": complete")
  }

  out <- list(
    settings = list(
      scenario_names = scenario_names,
      n = n,
      n_truth = n_truth,
      B = B,
      truth_target = truth_target,
      seed = seed,
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
      n_workers = n_workers,
      verbose = verbose,
      verbose_workers = verbose_workers,
      worker_outfile = worker_outfile,
      simulation_file = simulation_file
    ),
    outputs = all_outputs,
    rep_results = dplyr::bind_rows(lapply(all_outputs, `[[`, "rep_results")),
    performance_table = dplyr::bind_rows(lapply(all_outputs, `[[`, "performance"))
  )

  .parallel_all_methods_progress(verbose, "Parallel all-methods STTE simulation complete")
  out
}

############################################################
# Implementation code: run the parallel all-methods simulation
############################################################

if (sys.nframe() == 0L) {
  n_workers <- .parallel_all_methods_default_workers(max_workers = 10L)

  sim_all_methods_mdr_solutions_parallel <- run_all_methods_mdr_solutions_simulation_parallel(
    scenario_names = c(
      "base",
      "small_size",
      "strong_confounding",
      "poor_positivity",
      "poor_adherence",
      "rare_event",
      "nonlinear"
    ),
    n = 2000,
    n_truth = 50000,
    B = 500,
    truth_target = "sequential",
    seed = 20260630,
    direct_mdr_target_sim_multiplier = 10,
    direct_mdr_target_sim_size = NULL,
    direct_mdr_eps = 0.001,
    direct_mdr_weight_cap = NULL,
    direct_mdr_target_deplete_events = TRUE,
    direct_mdr_classifier_method = "random_forest",
    direct_mdr_random_forest_ntree = 500,
    direct_mdr_random_forest_mtry = NULL,
    direct_mdr_random_forest_min_node_size = NULL,
    direct_mdr_random_forest_engine = "ranger",
    comparing_methods = paper_compare_methods,
    n_workers = n_workers,
    verbose = TRUE,
    verbose_workers = TRUE
  )

  print(sim_all_methods_mdr_solutions_parallel$performance_table)

  estimates_plot <- plot_estimates_by_scenario(sim_all_methods_mdr_solutions_parallel)
  performance_plot <- plot_performance_by_scenario(sim_all_methods_mdr_solutions_parallel)
  weight_plot <- plot_weight_behavior_by_scenario(sim_all_methods_mdr_solutions_parallel)

  print(estimates_plot)
  print(performance_plot)
  print(weight_plot)

  # output_prefix <- "all_methods_parallel"

  # ggplot2::ggsave(
  #   filename = paste0(output_prefix, "_estimates.png"),
  #   plot = estimates_plot,
  #   width = 11,
  #   height = 6,
  #   dpi = 300
  # )
  # ggplot2::ggsave(
  #   filename = paste0(output_prefix, "_performance.png"),
  #   plot = performance_plot,
  #   width = 11,
  #   height = 6,
  #   dpi = 300
  # )
  # ggplot2::ggsave(
  #   filename = paste0(output_prefix, "_weights.png"),
  #   plot = weight_plot,
  #   width = 11,
  #   height = 7,
  #   dpi = 300
  # )

  # utils::write.csv(
  #   sim_all_methods_mdr_solutions_parallel$performance_table,
  #   file = paste0(output_prefix, "_performance.csv"),
  #   row.names = FALSE
  # )

  # saveRDS(
  #   sim_all_methods_mdr_solutions_parallel,
  #   file = paste0(output_prefix, "_simulation.rds")
  # )
}
