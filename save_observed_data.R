############################################################
# Save observed data generated for the all-methods MDR simulation
#
# This script mirrors the observed-data generation step inside
# run_one_replication_all_methods_mdr_solutions():
#
#   obs_data <- generate_observed_data(
#     n = n,
#     params = params,
#     seed = seed_base + b,
#     simulate_censoring = TRUE
#   )
#
# The full simulation driver sets seed_base = seed + 100000L, so this script
# uses the same rule. n_truth and truth_target are kept for alignment with the
# full simulation interface; they do not affect observed data generation.
############################################################

.save_obs_script_dir <- local({
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

# Observed-data export needs only the DGP, not every estimator and optional
# comparison backend.
source(file.path(.save_obs_script_dir, "stte_per_protocol_simulation.R"))
load_required_packages()

.safe_observed_data_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "scenario")
}

save_generated_observed_data_all_methods_mdr_solutions <- function(
    scenario_names = c("base", "strong_confounding", "poor_positivity"),
    n = 1500,
    n_truth = 50000,
    B = 10,
    truth_target = "sequential",
    seed = 20260630,
    output_dir = file.path(.save_obs_script_dir, "all_methods_mdr_solutions_observed_data_3"),
    overwrite = TRUE) {

  truth_target <- match.arg(truth_target, c("sequential", "baseline"))

  n <- .as_positive_integer(n, "n")
  n_truth <- .as_positive_integer(n_truth, "n_truth")
  B <- .as_positive_integer(B, "B")
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be a single finite number.", call. = FALSE)
  }
  if (!is.character(scenario_names) || length(scenario_names) < 1L) {
    stop("scenario_names must be a non-empty character vector.", call. = FALSE)
  }

  seed <- as.integer(seed)
  seed_base <- seed + 100000L

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  # Retain the cumulative event/censoring states as well as their increments so
  # saved replications can be analyzed by the wide-history LTMLE helper.
  obs_data_columns <- c(
    "id",
    "k",
    "X1",
    "X2",
    "X3",
    "L_k",
    "A_prev",
    "A_k",
    "Y_start",
    "C_start",
    "dY_next",
    "dC_next",
    "Y_next",
    "C_next"
  )

  saved_files <- list()
  row_index <- 1L

  for (scenario_name in scenario_names) {
    message("Generating observed data for scenario: ", scenario_name)
    params <- define_parameters(scenario_name)
    scenario_n <- .resolve_observed_sample_size(n, params)
    scenario_file_name <- .safe_observed_data_name(scenario_name)
    scenario_dir <- file.path(output_dir, scenario_file_name)
    dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    scenario_dir <- normalizePath(scenario_dir, winslash = "/", mustWork = TRUE)

    for (b in seq_len(B)) {
      observed_seed <- seed_base + b
      message(
        "  replication ", b, "/", B,
        " (n = ", scenario_n, ", seed = ", observed_seed, ")"
      )

      obs_data <- generate_observed_data(
        n = scenario_n,
        params = params,
        seed = observed_seed,
        simulate_censoring = TRUE
      )
      missing_columns <- setdiff(obs_data_columns, names(obs_data))
      if (length(missing_columns) > 0L) {
        stop(
          "Generated observed data is missing expected column(s): ",
          paste(missing_columns, collapse = ", "),
          call. = FALSE
        )
      }
      obs_data <- obs_data[, obs_data_columns, drop = FALSE]

      rds_path <- file.path(scenario_dir, paste0("obs_data_", b, ".rds"))
      if (file.exists(rds_path) && !isTRUE(overwrite)) {
        stop(
          "Refusing to overwrite existing file: ",
          rds_path,
          call. = FALSE
        )
      }

      saveRDS(obs_data, file = rds_path)

      saved_files[[row_index]] <- data.frame(
        scenario = scenario_name,
        replication = b,
        n = scenario_n,
        seed_base = seed_base,
        observed_seed = observed_seed,
        n_rows = nrow(obs_data),
        n_cols = ncol(obs_data),
        rds_file = normalizePath(rds_path, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  message("Saved observed data to: ", output_dir)

  invisible(do.call(rbind, saved_files))
}

if (sys.nframe() == 0L) {
  saved_files <- save_generated_observed_data_all_methods_mdr_solutions(
    scenario_names = c("base", "small_size",
                        "strong_confounding", "poor_positivity", "poor_adherence", 
                        "rare_event", 
                        "nonlinear"),
    n = 2000,
    n_truth = 50000,
    B = 100,
    truth_target = "sequential",
    seed = 688202,
    output_dir = file.path(.save_obs_script_dir, "all_methods_mdr_solutions_observed_data_6")
  )

  # print(saved_files)
}
