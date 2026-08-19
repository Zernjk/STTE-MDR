############################################################
# Core data-generation and estimator helpers for the paper simulation
# Main estimand: 12-month per-protocol cumulative risk difference
# Strategies implemented by default:
#   g1: initiate treatment immediately at trial start and remain treated
#   g0: remain untreated throughout follow-up
############################################################

load_required_packages <- function() {
  needed <- c("dplyr", "tidyr", "data.table")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Please install the following packages before running this script: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
  })
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

clip_prob <- function(p, eps = 1e-4) {
  pmin(pmax(p, eps), 1 - eps)
}

.dgp_prob_clip <- function(params, component) {
  component_name <- paste0(component, "_prob_clip")
  component_clip <- params[[component_name]]
  if (is.null(component_clip)) component_clip <- params$prob_clip

  if (length(component_clip) != 1L || !is.finite(component_clip) ||
      component_clip < 0 || component_clip >= 0.5) {
    stop(
      component_name,
      " (or the fallback prob_clip) must be one finite value in [0, 0.5).",
      call. = FALSE
    )
  }

  component_clip
}

.as_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }

  as.integer(x)
}

#-----------------------------------------------------------
# 1. Define DGP parameters
#-----------------------------------------------------------

define_parameters <- function(scenario = "base") {
  params <- list(
    # Time structure
    K = 24,                    # total observed cohort follow-up months
    tau = 12,                  # target-trial follow-up months
    grace = 0,                 # default exact protocol: no grace period
    require_untreated_entry = TRUE,
    observed_sample_size = NULL, # scenario-specific override; NULL uses caller-supplied n

    # Baseline covariate distribution
    x1_prob = 0.50,
    x2_mean = 0,
    x2_sd = 1,
    x3_intercept = -0.20,
    x3_x1 = 0.30,
    x3_x2 = 0.50,

    # Baseline L0 model
    L0_intercept = 0,
    L0_x1 = 0.50,
    L0_x2 = 0.70,
    L0_x3 = 0.40,
    L0_sd = 1.00,

    # Time-varying confounder L_{k+1}
    L_intercept = 0,
    L_rho = 0.60,
    L_x1 = 0.30,
    L_x2 = 0.40,
    L_x3 = 0.20,
    L_Aprev = -0.50,           # treatment reduces future severity
    L_time = 0.00,
    L_sd = 1.00,

    # Natural treatment assignment model
    alpha0 = -1.20,
    alpha_Aprev = 2.00,
    alpha_L = 0.80,
    alpha_x1 = 0.30,
    alpha_x2 = 0.20,
    alpha_x3 = 0.10,
    alpha_time = 0.02,

    # Outcome discrete-time hazard model
    beta0 = -6.50,
    beta_A = -0.40,            # protective treatment effect
    beta_L = 0.80,
    beta_x1 = 0.30,
    beta_x2 = 0.20,
    beta_x3 = 0.20,
    beta_time = 0.03,

    # Natural censoring model
    delta0 = -5.00,
    delta_A = 0.20,
    delta_L = 0.40,
    delta_x1 = 0.20,
    delta_x2 = 0.20,
    delta_x3 = 0.10,
    delta_time = 0.00,

    # Nonlinear DGP option
    nonlinear = FALSE,
    alpha_L2 = 0.00,
    beta_L2 = 0.00,
    beta_A_L = 0.00,

    # Weight handling
    weight_truncation_quantile = NULL, # e.g., 0.99 if desired
    prob_clip = 0.01,

    # Optional DGP-specific overrides. Keeping these separate prevents a
    # scenario that changes the event rate from also changing treatment and
    # censoring probabilities. NULL retains the legacy prob_clip fallback.
    treatment_prob_clip = NULL,
    outcome_prob_clip = NULL,
    censor_prob_clip = NULL
  )

  if (is.list(scenario)) {
    params <- utils::modifyList(params, scenario)
    return(params)
  }

  if (scenario == "base") {
    return(params)
  }

  if (scenario == "small_size") {
    params$observed_sample_size <- 100L
    return(params)
  }

  if (scenario == "null_effect") {
    params$beta_A <- 0
    # Treatment must not affect either the outcome directly or its mediator L.
    params$L_Aprev <- 0
    return(params)
  }

  if (scenario == "weak_confounding") {
    params$alpha_L <- 0.30
    params$beta_L <- 0.30
    return(params)
  }

  if (scenario == "strong_confounding") {
    params$alpha_L <- 1.20
    params$beta_L <- 1.20
    return(params)
  }

  if (scenario == "poor_positivity") {
    params$alpha_L <- 2.00
    params$beta_L <- 0.80
    return(params)
  }

  if (scenario == "poor_adherence") {
    params$alpha_Aprev <- 0.60
    return(params)
  }

  if (scenario == "rare_event") {
    params$beta0 <- -8.00
    params$outcome_prob_clip <- 1e-4
    return(params)
  }

  if (scenario == "nonlinear") {
    params$nonlinear <- TRUE
    params$alpha_L2 <- 0.35
    params$beta_L2 <- 0.20
    params$beta_A_L <- -0.25
    return(params)
  }

  stop("Unknown scenario: ", scenario, call. = FALSE)
}

#-----------------------------------------------------------
# 2. Core DGP helper functions
#-----------------------------------------------------------

generate_large_baseline_population <- function(n, params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X1 <- rbinom(n, 1, params$x1_prob)
  X2 <- rnorm(n, mean = params$x2_mean, sd = params$x2_sd)
  pX3 <- expit(params$x3_intercept + params$x3_x1 * X1 + params$x3_x2 * X2)
  X3 <- rbinom(n, 1, pX3)

  L0 <- params$L0_intercept +
    params$L0_x1 * X1 +
    params$L0_x2 * X2 +
    params$L0_x3 * X3 +
    rnorm(n, mean = 0, sd = params$L0_sd)

  data.frame(
    id = seq_len(n),
    X1 = X1,
    X2 = X2,
    X3 = X3,
    L0 = L0
  )
}

compute_pA <- function(A_prev, L_k, X1, X2, X3, k, params) {
  lp <- params$alpha0 +
    params$alpha_Aprev * A_prev +
    params$alpha_L * L_k +
    params$alpha_x1 * X1 +
    params$alpha_x2 * X2 +
    params$alpha_x3 * X3 +
    params$alpha_time * k

  if (isTRUE(params$nonlinear)) {
    lp <- lp + params$alpha_L2 * L_k^2
  }

  clip_prob(expit(lp), .dgp_prob_clip(params, "treatment"))
}

compute_pY <- function(A_k, L_k, X1, X2, X3, k, params) {
  lp <- params$beta0 +
    params$beta_A * A_k +
    params$beta_L * L_k +
    params$beta_x1 * X1 +
    params$beta_x2 * X2 +
    params$beta_x3 * X3 +
    params$beta_time * k

  if (isTRUE(params$nonlinear)) {
    lp <- lp + params$beta_L2 * L_k^2 + params$beta_A_L * A_k * L_k
  }

  clip_prob(expit(lp), .dgp_prob_clip(params, "outcome"))
}

compute_pC <- function(A_k, L_k, X1, X2, X3, k, params) {
  lp <- params$delta0 +
    params$delta_A * A_k +
    params$delta_L * L_k +
    params$delta_x1 * X1 +
    params$delta_x2 * X2 +
    params$delta_x3 * X3 +
    params$delta_time * k

  clip_prob(expit(lp), .dgp_prob_clip(params, "censor"))
}

simulate_next_L <- function(L_k, A_k, X1, X2, X3, k, params) {
  params$L_intercept +
    params$L_rho * L_k +
    params$L_x1 * X1 +
    params$L_x2 * X2 +
    params$L_x3 * X3 +
    params$L_Aprev * A_k +
    params$L_time * k +
    rnorm(length(L_k), mean = 0, sd = params$L_sd)
}

#-----------------------------------------------------------
# 3. Generate observed longitudinal data
#-----------------------------------------------------------

simulate_observed_from_baseline <- function(base, params, simulate_censoring = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(base)
  K <- params$K

  X1 <- base$X1
  X2 <- base$X2
  X3 <- base$X3

  L <- matrix(NA_real_, nrow = n, ncol = K + 1)
  A <- matrix(NA_integer_, nrow = n, ncol = K)
  dY <- matrix(0L, nrow = n, ncol = K)
  dC <- matrix(0L, nrow = n, ncol = K)
  Y <- matrix(0L, nrow = n, ncol = K + 1)
  C <- matrix(0L, nrow = n, ncol = K + 1)

  L[, 1] <- base$L0

  for (k in 0:(K - 1)) {
    at_risk <- which(Y[, k + 1] == 0 & C[, k + 1] == 0)

    # Carry forward cumulative statuses by default
    Y[, k + 2] <- Y[, k + 1]
    C[, k + 2] <- C[, k + 1]

    if (length(at_risk) == 0) next

    A_prev <- if (k == 0) rep(0, length(at_risk)) else A[at_risk, k]
    A_prev[is.na(A_prev)] <- 0

    pA <- compute_pA(
      A_prev = A_prev,
      L_k = L[at_risk, k + 1],
      X1 = X1[at_risk],
      X2 = X2[at_risk],
      X3 = X3[at_risk],
      k = k,
      params = params
    )
    A[at_risk, k + 1] <- rbinom(length(at_risk), 1, pA)

    pY <- compute_pY(
      A_k = A[at_risk, k + 1],
      L_k = L[at_risk, k + 1],
      X1 = X1[at_risk],
      X2 = X2[at_risk],
      X3 = X3[at_risk],
      k = k,
      params = params
    )
    dY[at_risk, k + 1] <- rbinom(length(at_risk), 1, pY)
    Y[at_risk[dY[at_risk, k + 1] == 1], k + 2] <- 1L

    survived_interval <- at_risk[dY[at_risk, k + 1] == 0]

    if (simulate_censoring && length(survived_interval) > 0) {
      pC <- compute_pC(
        A_k = A[survived_interval, k + 1],
        L_k = L[survived_interval, k + 1],
        X1 = X1[survived_interval],
        X2 = X2[survived_interval],
        X3 = X3[survived_interval],
        k = k,
        params = params
      )
      dC[survived_interval, k + 1] <- rbinom(length(survived_interval), 1, pC)
      C[survived_interval[dC[survived_interval, k + 1] == 1], k + 2] <- 1L
    }

    # Generate next L for people who were at risk at start of interval.
    # It is only used later if the person remains event-free and uncensored.
    # Drawing for the full risk set preserves the established seeded DGP stream.
    if (k < K - 1) {
      L[at_risk, k + 2] <- simulate_next_L(
        L_k = L[at_risk, k + 1],
        A_k = A[at_risk, k + 1],
        X1 = X1[at_risk],
        X2 = X2[at_risk],
        X3 = X3[at_risk],
        k = k,
        params = params
      )
    }
  }

  rows <- vector("list", K)
  for (k in 0:(K - 1)) {
    at_start <- which(Y[, k + 1] == 0 & C[, k + 1] == 0)
    if (length(at_start) == 0) {
      rows[[k + 1]] <- NULL
      next
    }

    A_prev <- if (k == 0) rep(0, length(at_start)) else A[at_start, k]
    A_prev[is.na(A_prev)] <- 0

    rows[[k + 1]] <- data.frame(
      id = base$id[at_start],
      k = k,
      X1 = X1[at_start],
      X2 = X2[at_start],
      X3 = X3[at_start],
      L_k = L[at_start, k + 1],
      A_prev = A_prev,
      A_k = A[at_start, k + 1],
      Y_start = Y[at_start, k + 1],
      C_start = C[at_start, k + 1],
      dY_next = dY[at_start, k + 1],
      dC_next = dC[at_start, k + 1],
      Y_next = Y[at_start, k + 2],
      C_next = C[at_start, k + 2]
    )
  }

  dplyr::bind_rows(rows)
}

.resolve_observed_sample_size <- function(n, params) {
  scenario_n <- params$observed_sample_size
  resolved_n <- if (is.null(scenario_n)) n else scenario_n

  .as_positive_integer(resolved_n, "The observed sample size")
}

generate_observed_data <- function(
    n,
    params,
    seed = NULL,
    simulate_censoring = TRUE,
    use_scenario_sample_size = TRUE) {

  if (isTRUE(use_scenario_sample_size)) {
    n <- .resolve_observed_sample_size(n, params)
  }

  base <- generate_large_baseline_population(n = n, params = params, seed = seed)
  simulate_observed_from_baseline(
    base = base,
    params = params,
    simulate_censoring = simulate_censoring,
    seed = NULL
  )
}

#-----------------------------------------------------------
# 4. Eligibility and intervention truth
#-----------------------------------------------------------

get_eligible_person_trials <- function(obs_data, params) {
  max_trial <- params$K - params$tau

  elig <- obs_data %>%
    dplyr::filter(k <= max_trial, Y_start == 0, C_start == 0)

  if (isTRUE(params$require_untreated_entry)) {
    elig <- elig %>% dplyr::filter(k == 0 | A_prev == 0)
  }

  elig %>%
    dplyr::transmute(
      id = id,
      trial_m = k,
      id_trial = paste(id, k, sep = "_"),
      X1 = X1,
      X2 = X2,
      X3 = X3,
      L_trial = L_k,
      A_prev_entry = A_prev
    )
}

simulate_future_from_states <- function(states, strategy = c("g1", "g0"), params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  strategy <- match.arg(strategy)

  n <- nrow(states)
  tau <- params$tau

  X1 <- states$X1
  X2 <- states$X2
  X3 <- states$X3
  L_current <- states$L_trial
  trial_m <- states$trial_m

  Y_cum <- rep(0L, n)

  for (j in 0:(tau - 1)) {
    active <- which(Y_cum == 0)
    if (length(active) == 0) break

    calendar_k <- trial_m[active] + j

    # Deterministic per-protocol interventions.
    # g1: immediately treated and remain treated.
    # g0: untreated throughout follow-up.
    A_forced <- if (strategy == "g1") rep(1L, length(active)) else rep(0L, length(active))

    pY <- compute_pY(
      A_k = A_forced,
      L_k = L_current[active],
      X1 = X1[active],
      X2 = X2[active],
      X3 = X3[active],
      k = calendar_k,
      params = params
    )

    dY <- rbinom(length(active), 1, pY)
    event_ids <- active[dY == 1]
    Y_cum[event_ids] <- 1L

    still_active <- active[dY == 0]
    if (length(still_active) > 0 && j < tau - 1) {
      A_still <- if (strategy == "g1") rep(1L, length(still_active)) else rep(0L, length(still_active))
      L_current[still_active] <- simulate_next_L(
        L_k = L_current[still_active],
        A_k = A_still,
        X1 = X1[still_active],
        X2 = X2[still_active],
        X3 = X3[still_active],
        k = trial_m[still_active] + j,
        params = params
      )
    }
  }

  Y_cum
}

simulate_under_intervention <- function(data, strategy = c("g1", "g0"), params, seed = NULL) {
  # This wrapper matches the pseudo-code. It treats input `data` as a baseline
  # population and computes baseline-trial counterfactual outcomes from m = 0.
  # For full sequential-trial truth, use calculate_true_effect(..., target = "sequential").
  strategy <- match.arg(strategy)

  states <- data.frame(
    id = data$id,
    trial_m = 0L,
    id_trial = paste(data$id, 0L, sep = "_"),
    X1 = data$X1,
    X2 = data$X2,
    X3 = data$X3,
    L_trial = data$L0,
    A_prev_entry = 0L
  )

  simulate_future_from_states(
    states = states,
    strategy = strategy,
    params = params,
    seed = seed
  )
}

calculate_true_effect <- function(params, n_truth = 100000, target = c("sequential", "baseline"), seed = 1) {
  target <- match.arg(target)
  n_truth <- .as_positive_integer(n_truth, "n_truth")

  if (params$grace != 0) {
    warning(
      "The truth function currently implements deterministic strategies: ",
      "g1 = immediately treated and g0 = never treated. ",
      "For a grace-period estimand, define a stochastic or deterministic initiation rule explicitly."
    )
  }

  base <- generate_large_baseline_population(n = n_truth, params = params, seed = seed)

  if (target == "baseline") {
    Y_g1 <- simulate_under_intervention(base, strategy = "g1", params = params, seed = seed + 10)
    # Common random numbers reduce Monte Carlo noise in the truth contrast and
    # make a genuinely null DGP return an exactly null paired contrast.
    Y_g0 <- simulate_under_intervention(base, strategy = "g0", params = params, seed = seed + 10)

    return(list(
      psi_true = mean(Y_g1 - Y_g0),
      risk_g1 = mean(Y_g1),
      risk_g0 = mean(Y_g0),
      n_truth = n_truth,
      n_eligible_person_trials = n_truth,
      target = target
    ))
  }

  # Sequential-trial truth: first generate a large natural history without
  # natural censoring, identify eligible person-trials, then simulate each
  # eligible state forward under g1 and g0.
  natural_long <- simulate_observed_from_baseline(
    base = base,
    params = params,
    simulate_censoring = FALSE,
    seed = seed + 5
  )

  states <- get_eligible_person_trials(natural_long, params)

  if (nrow(states) == 0) {
    stop("No eligible person-trials in truth simulation.", call. = FALSE)
  }

  Y_g1 <- simulate_future_from_states(states, strategy = "g1", params = params, seed = seed + 10)
  Y_g0 <- simulate_future_from_states(states, strategy = "g0", params = params, seed = seed + 10)

  list(
    psi_true = mean(Y_g1 - Y_g0),
    risk_g1 = mean(Y_g1),
    risk_g0 = mean(Y_g0),
    n_truth = n_truth,
    n_eligible_person_trials = nrow(states),
    target = target
  )
}

#-----------------------------------------------------------
# 5. Construct STTE strategy/adherence dataset
#-----------------------------------------------------------

.assign_clone_status_one <- function(df, clone_value = NULL, grace = 0) {
  df <- df[order(df$follow_k), ]
  if (is.null(clone_value)) {
    clone_value <- unique(df$clone)
  }
  if (length(clone_value) != 1 || is.na(clone_value)) {
    stop("Each group must contain one clone value.")
  }

  initiated <- FALSE
  keep <- rep(TRUE, nrow(df))
  art_censor <- rep(0L, nrow(df))
  initiated_before <- rep(0L, nrow(df))

  for (r in seq_len(nrow(df))) {
    A_now <- df$A_k[r]
    follow_now <- df$follow_k[r]
    initiated_before[r] <- as.integer(initiated)

    if (clone_value == 0) {
      # g0: remain untreated. Any treatment initiation is a protocol deviation.
      if (!is.na(A_now) && A_now == 1) {
        art_censor[r] <- 1L
      }
    } else {
      # g1: initiate treatment by grace month and remain treated afterward.
      # With grace = 0, this means A must equal 1 at every follow-up month.
      if (!initiated && !is.na(A_now) && A_now == 1) {
        initiated <- TRUE
      } else if (!initiated && !is.na(A_now) && A_now == 0 && follow_now >= grace) {
        art_censor[r] <- 1L
      } else if (initiated && !is.na(A_now) && A_now == 0) {
        art_censor[r] <- 1L
      }
    }

    # Stop clone after first artificial censoring, natural censoring, or outcome.
    if (art_censor[r] == 1L || df$dC_next[r] == 1L || df$dY_next[r] == 1L) {
      if (r < nrow(df)) keep[(r + 1):nrow(df)] <- FALSE
      break
    }
  }

  df$artificial_censor <- art_censor
  df$initiated_before <- initiated_before
  df[keep, , drop = FALSE]
}

construct_STTE_dataset <- function(obs_data, params) {
  elig <- get_eligible_person_trials(obs_data, params)

  if (nrow(elig) == 0) {
    stop("No eligible person-trials found.", call. = FALSE)
  }

  # With an exact protocol (grace = 0), current treatment at trial entry
  # uniquely determines the compatible strategy. Keep one observed copy of
  # each eligible trial and follow it until its first treatment deviation.
  if (isTRUE(params$grace == 0)) {
    follow_grid <- data.frame(
      follow_k = 0:(params$tau - 1L)
    )

    elig_rep <- elig[rep(seq_len(nrow(elig)), each = nrow(follow_grid)), , drop = FALSE]
    grid_rep <- follow_grid[rep(seq_len(nrow(follow_grid)), times = nrow(elig)), , drop = FALSE]

    stte_raw <- dplyr::bind_cols(elig_rep, grid_rep) %>%
      dplyr::mutate(k = trial_m + follow_k) %>%
      dplyr::left_join(
        obs_data %>%
          dplyr::select(
            id, k, L_k, A_prev, A_k,
            dY_next, dC_next, Y_start, C_start,
            Y_next, C_next
          ),
        by = c("id", "k")
      ) %>%
      dplyr::filter(!is.na(A_k)) %>%
      dplyr::arrange(id_trial, follow_k)

    entry_strategy <- stte_raw %>%
      dplyr::filter(follow_k == 0L) %>%
      dplyr::transmute(id_trial, strategy = as.integer(A_k))

    return(
      stte_raw %>%
        dplyr::left_join(entry_strategy, by = "id_trial") %>%
        dplyr::group_by(id_trial, strategy) %>%
        dplyr::group_modify(~ .assign_clone_status_one(
          .x,
          clone_value = .y$strategy,
          grace = 0
        )) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          censor_event = as.integer(artificial_censor == 1L | dC_next == 1L)
        )
    )
  }

  clone_grid <- tidyr::expand_grid(
    clone = c(0L, 1L),
    follow_k = 0:(params$tau - 1)
  )

  # Cross join eligible person-trials with clone/follow-up grid.
  # This base-R construction is used for compatibility across dplyr/tidyr versions.
  elig_rep <- elig[rep(seq_len(nrow(elig)), each = nrow(clone_grid)), , drop = FALSE]
  grid_rep <- clone_grid[rep(seq_len(nrow(clone_grid)), times = nrow(elig)), , drop = FALSE]

  stte_raw <- dplyr::bind_cols(elig_rep, grid_rep) %>%
    dplyr::mutate(k = trial_m + follow_k) %>%
    dplyr::left_join(
      obs_data %>%
        dplyr::select(
          id, k, L_k, A_prev, A_k,
          dY_next, dC_next, Y_start, C_start,
          Y_next, C_next
        ),
      by = c("id", "k")
    ) %>%
    dplyr::filter(!is.na(A_k)) %>%
    dplyr::arrange(id_trial, clone, follow_k)

  stte <- stte_raw %>%
    dplyr::group_by(id_trial, clone) %>%
    dplyr::group_modify(~ .assign_clone_status_one(
      .x,
      clone_value = .y$clone,
      grace = params$grace
    )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      strategy = clone,
      censor_event = as.integer(artificial_censor == 1L | dC_next == 1L)
    )

  stte
}

#-----------------------------------------------------------
# 6a. Estimate IPCW and per-protocol risk difference
#-----------------------------------------------------------

.fit_logistic_safe <- function(formula, data, weights = NULL) {
  model_frame <- try(
    model.frame(formula, data = data),
    silent = TRUE
  )
  if (inherits(model_frame, "try-error")) return(NULL)

  y <- model.response(model_frame)
  if (length(unique(y[!is.na(y)])) < 2) {
    return(NULL)
  }

  if (!is.null(weights)) {
    if (length(weights) != nrow(data)) return(NULL)
    weights <- as.numeric(weights)
    if (any(!is.finite(weights)) || any(weights < 0)) return(NULL)
  }

  glm_args <- list(
    formula = formula,
    family = binomial(),
    data = data,
    control = glm.control(maxit = 50)
  )
  if (!is.null(weights)) {
    glm_args$weights <- weights
  }

  fit <- try(
    suppressWarnings(do.call(glm, glm_args)),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(NULL)
  fit
}

.predict_prob_safe <- function(fit, newdata, fallback, eps = 0.01) {
  fallback <- as.numeric(fallback)[1]
  if (!is.finite(fallback)) fallback <- 0.5
  fallback <- clip_prob(fallback, eps)

  if (is.null(fit)) {
    return(rep(fallback, nrow(newdata)))
  }

  p <- try(suppressWarnings(predict(fit, newdata = newdata, type = "response")), silent = TRUE)
  if (inherits(p, "try-error")) {
    return(rep(fallback, nrow(newdata)))
  }

  p <- as.numeric(p)
  p[!is.finite(p)] <- fallback
  clip_prob(p, eps)
}

.estimate_clone_standardized_risks <- function(
    outcome_fit,
    stte_data,
    params,
    fallback_hazard) {

  eligible_base <- stte_data %>%
    dplyr::filter(follow_k == 0L) %>%
    dplyr::distinct(id_trial, trial_m, X1, X2, X3, L_trial)

  estimate_standardized_risk <- function(g_value) {
    pred_grid <- eligible_base[rep(seq_len(nrow(eligible_base)), each = params$tau), , drop = FALSE]
    pred_grid$follow_k <- rep(0:(params$tau - 1), times = nrow(eligible_base))
    pred_grid$follow_f <- factor(pred_grid$follow_k, levels = 0:(params$tau - 1L))
    pred_grid$strategy <- g_value

    pred_grid$hazard <- .predict_prob_safe(
      outcome_fit,
      newdata = pred_grid,
      fallback = fallback_hazard,
      eps = 1e-8
    )

    risks <- pred_grid %>%
      dplyr::arrange(id_trial, follow_k) %>%
      dplyr::group_by(id_trial) %>%
      dplyr::summarise(
        risk = 1 - prod(1 - hazard),
        .groups = "drop"
      )

    mean(risks$risk)
  }

  risk_g1 <- estimate_standardized_risk(1L)
  risk_g0 <- estimate_standardized_risk(0L)

  list(
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    psi_hat = risk_g1 - risk_g0
  )
}

estimate_ipcw <- function(stte_data, params) {
  stte_model <- stte_data %>%
    dplyr::mutate(
      not_artificial_censored = as.integer(artificial_censor == 0L),
      not_natural_censored = as.integer(dC_next == 0L)
    )

  # Artificial censoring is treatment-protocol deviation. With grace = 0,
  # observed trials are not cloned: strategy is assigned from treatment at
  # entry. Fit the natural treatment probability once and convert it to the
  # probability of following the assigned strategy. This includes the entry
  # treatment probability that was previously represented by immediately
  # censoring the incompatible baseline clone.
  if (isTRUE(params$grace == 0)) {
    art_den_formula <- A_k ~ factor(follow_k) + trial_m +
      X1 + X2 + X3 + L_k + A_prev

    art_num_formula <- A_k ~ factor(follow_k) + trial_m +
      X1 + X2 + X3 + A_prev

    art_den_fit <- .fit_logistic_safe(art_den_formula, data = stte_model)
    art_num_fit <- .fit_logistic_safe(art_num_formula, data = stte_model)

    treatment_fallback <- mean(stte_model$A_k, na.rm = TRUE)
    p_treatment_den <- .predict_prob_safe(
      art_den_fit,
      stte_model,
      fallback = treatment_fallback,
      eps = params$prob_clip
    )
    p_treatment_num <- .predict_prob_safe(
      art_num_fit,
      stte_model,
      fallback = treatment_fallback,
      eps = params$prob_clip
    )

    p_art_den <- ifelse(
      stte_model$strategy == 1L,
      p_treatment_den,
      1 - p_treatment_den
    )
    p_art_num <- ifelse(
      stte_model$strategy == 1L,
      p_treatment_num,
      1 - p_treatment_num
    )
  } else {
    # With a grace period, histories can initially be compatible with both
    # strategies, so retain the original clone-specific censoring models.
    art_den_formula <- not_artificial_censored ~ clone * (
      factor(follow_k) + trial_m + X1 + X2 + X3 + L_k + A_prev
    )

    art_num_formula <- not_artificial_censored ~ clone * (
      factor(follow_k) + trial_m + X1 + X2 + X3
    )

    art_den_fit <- .fit_logistic_safe(art_den_formula, data = stte_model)
    art_num_fit <- .fit_logistic_safe(art_num_formula, data = stte_model)

    art_fallback <- mean(stte_model$not_artificial_censored, na.rm = TRUE)
    p_art_den <- .predict_prob_safe(
      art_den_fit,
      stte_model,
      fallback = art_fallback,
      eps = params$prob_clip
    )
    p_art_num <- .predict_prob_safe(
      art_num_fit,
      stte_model,
      fallback = art_fallback,
      eps = params$prob_clip
    )
  }

  # Natural censoring is generated only after surviving the interval. Fit its
  # model among outcome-free rows that have not already deviated from the
  # assigned treatment protocol. Conditional on actual A_k and history, clone
  # or strategy is not part of the natural-censoring DGP.
  natural_model_data <- stte_model %>%
    dplyr::filter(
      artificial_censor == 0L,
      dY_next == 0L
    )

  nat_den_formula <- not_natural_censored ~ factor(follow_k) + trial_m +
    X1 + X2 + X3 + L_k + A_k

  nat_num_formula <- not_natural_censored ~ strategy + factor(follow_k) + trial_m +
    X1 + X2 + X3

  nat_den_fit <- .fit_logistic_safe(nat_den_formula, data = natural_model_data)
  nat_num_fit <- .fit_logistic_safe(nat_num_formula, data = natural_model_data)

  nat_fallback <- mean(natural_model_data$not_natural_censored, na.rm = TRUE)

  p_nat_den <- .predict_prob_safe(
    nat_den_fit,
    stte_model,
    fallback = nat_fallback,
    eps = params$prob_clip
  )
  p_nat_num <- .predict_prob_safe(
    nat_num_fit,
    stte_model,
    fallback = nat_fallback,
    eps = params$prob_clip
  )

  stte_w <- stte_model %>%
    dplyr::mutate(
      p_art_den = p_art_den,
      p_art_num = p_art_num,
      p_nat_den = p_nat_den,
      p_nat_num = p_nat_num,
      ipcw_art_step = p_art_num / p_art_den,
      ipcw_nat_step = p_nat_num / p_nat_den
    ) %>%
    dplyr::arrange(id_trial, strategy, follow_k) %>%
    dplyr::group_by(id_trial, strategy) %>%
    dplyr::mutate(
      ipcw_artificial = cumprod(ipcw_art_step),
      ipcw_artificial_unstabilized = cumprod(1 / p_art_den),
      # dY_next is generated before dC_next. The outcome at follow_k is
      # therefore weighted by natural censoring only through follow_k - 1.
      ipcw_natural_end = cumprod(ipcw_nat_step),
      ipcw_natural = dplyr::lag(ipcw_natural_end, default = 1),
      ipcw_natural_unstabilized_end = cumprod(1 / p_nat_den),
      ipcw_natural_unstabilized = dplyr::lag(
        ipcw_natural_unstabilized_end,
        default = 1
      ),
      ipcw = ipcw_artificial * ipcw_natural,
      # AIPW corrections require inverse mechanism probabilities rather than
      # stabilized MSM weights.
      ipcw_unstabilized = ipcw_artificial_unstabilized *
        ipcw_natural_unstabilized
    ) %>%
    dplyr::ungroup()

  if (!is.null(params$weight_truncation_quantile)) {
    cap <- stats::quantile(stte_w$ipcw, probs = params$weight_truncation_quantile, na.rm = TRUE)
    stte_w$ipcw <- pmin(stte_w$ipcw, as.numeric(cap))

    unstabilized_cap <- stats::quantile(
      stte_w$ipcw_unstabilized,
      probs = params$weight_truncation_quantile,
      na.rm = TRUE
    )
    stte_w$ipcw_unstabilized <- pmin(
      stte_w$ipcw_unstabilized,
      as.numeric(unstabilized_cap)
    )
  }

  list(
    data = stte_w,
    denominator_fit = list(
      artificial = art_den_fit,
      natural = nat_den_fit
    ),
    numerator_fit = list(
      artificial = art_num_fit,
      natural = nat_num_fit
    ),
    artificial_denominator_fit = art_den_fit,
    artificial_numerator_fit = art_num_fit,
    natural_denominator_fit = nat_den_fit,
    natural_numerator_fit = nat_num_fit
  )
}

estimate_STTE_per_protocol_ipcw_pool <- function(
    stte_data,
    params,
    return_models = FALSE) {
  ipcw_obj <- estimate_ipcw(stte_data, params)
  stte_w <- ipcw_obj$data

  # Current treatment is assigned before dY_next, so protocol deviations make
  # the current outcome incompatible with the strategy. Natural censoring is
  # generated after dY_next, so that interval's outcome remains observed.
  outcome_data <- stte_w %>%
    dplyr::filter(artificial_censor == 0L) %>%
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
    stte_data = stte_data,
    params = params,
    fallback_hazard = mean(outcome_data$dY_next, na.rm = TRUE)
  )

  out <- list(
    psi_hat = risk_obj$psi_hat,
    risk_g1 = risk_obj$risk_g1,
    risk_g0 = risk_obj$risk_g0,
    n_rows_stte = nrow(stte_data),
    n_rows_outcome = nrow(outcome_data),
    weight_mean = mean(outcome_data$ipcw, na.rm = TRUE),
    weight_max = max(outcome_data$ipcw, na.rm = TRUE),
    weight_cv = stats::sd(outcome_data$ipcw, na.rm = TRUE) / mean(outcome_data$ipcw, na.rm = TRUE),
    failure_reason = NA_character_
  )

  if (return_models) {
    out$ipcw_data <- stte_w
    out$denominator_fit <- ipcw_obj$denominator_fit
    out$numerator_fit <- ipcw_obj$numerator_fit
    out$outcome_fit <- outcome_fit
  }

  out
}

#-----------------------------------------------------------
# 6b. Optional estimator using the TrialEmulation package
#-----------------------------------------------------------

prepare_trial_emulation_input <- function(obs_data, params) {
  obs_data %>%
    dplyr::mutate(
      period = k,
      treatment = A_k,
      outcome = dY_next,
      censored = dC_next,
      eligible = as.integer(
        k <= params$K - params$tau &
          Y_start == 0L &
          C_start == 0L &
          (
            !isTRUE(params$require_untreated_entry) |
              k == 0L |
              A_prev == 0L
          )
      )
    ) %>%
    dplyr::select(
      id, period, treatment, outcome, censored, eligible,
      X1, X2, X3, L_k, A_prev
    )
}

.extract_trial_emulation_cum_inc <- function(pred_df, followup_time) {
  if (!"followup_time" %in% names(pred_df) || !"cum_inc" %in% names(pred_df)) {
    return(NA_real_)
  }

  row <- pred_df[pred_df$followup_time == followup_time, , drop = FALSE]
  if (nrow(row) == 0) return(NA_real_)

  as.numeric(row$cum_inc[1])
}

.predict_trial_emulation <- function(object, ...) {
  predict_fun <- get(
    "predict",
    envir = asNamespace("TrialEmulation"),
    inherits = FALSE
  )

  predict_fun(object, ...)
}

.summarize_trial_emulation_weights <- function(outcome_df) {
  weight_col <- if ("weight" %in% names(outcome_df)) {
    "weight"
  } else if ("w" %in% names(outcome_df)) {
    "w"
  } else {
    NA_character_
  }

  if (is.na(weight_col)) {
    return(list(
      weight_mean = NA_real_,
      weight_max = NA_real_,
      weight_cv = NA_real_,
      weight_source = NA_character_
    ))
  }

  weights <- as.numeric(outcome_df[[weight_col]])
  weights <- weights[is.finite(weights)]

  if (length(weights) == 0) {
    return(list(
      weight_mean = NA_real_,
      weight_max = NA_real_,
      weight_cv = NA_real_,
      weight_source = weight_col
    ))
  }

  weight_mean <- mean(weights)

  list(
    weight_mean = weight_mean,
    weight_max = max(weights),
    weight_cv = ifelse(abs(weight_mean) < 1e-12, NA_real_, stats::sd(weights) / weight_mean),
    weight_source = weight_col
  )
}

estimate_STTE_per_protocol_trial_emulation <- function(
    obs_data,
    params,
    chunk_size = 500,
    conf_int = FALSE,
    samples = 100,
    model_dir = file.path(tempdir(), "trial_emulation_models"),
    return_trial_sequence = FALSE) {

  if (!requireNamespace("TrialEmulation", quietly = TRUE)) {
    stop(
      "The TrialEmulation package is required for this estimator. ",
      "Install it with install.packages('TrialEmulation').",
      call. = FALSE
    )
  }

  if (params$grace != 0) {
    warning(
      "This TrialEmulation wrapper is written for the default immediate-initiation ",
      "per-protocol strategies. Check the package setup before using a nonzero grace period."
    )
  }

  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }

  te_input <- prepare_trial_emulation_input(obs_data, params)

  switch_model_dir <- file.path(model_dir, "switch_models")
  censor_model_dir <- file.path(model_dir, "censor_models")
  dir.create(switch_model_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(censor_model_dir, recursive = TRUE, showWarnings = FALSE)

  # TrialEmulation's PP analysis censors people when they deviate from their
  # baseline assigned treatment and uses switching weights for that censoring.
  trial_pp <- TrialEmulation::trial_sequence("PP")

  trial_pp <- TrialEmulation::set_data(
    trial_pp,
    data = te_input,
    id = "id",
    period = "period",
    treatment = "treatment",
    outcome = "outcome",
    eligible = "eligible"
  )

  trial_pp <- TrialEmulation::set_switch_weight_model(
    trial_pp,
    numerator = ~ X1 + X2 + X3,
    denominator = ~ X1 + X2 + X3 + L_k + A_prev + period,
    model_fitter = TrialEmulation::stats_glm_logit(save_path = switch_model_dir)
  )

  trial_pp <- TrialEmulation::set_censor_weight_model(
    trial_pp,
    censor_event = "censored",
    numerator = ~ X1 + X2 + X3,
    denominator = ~ X1 + X2 + X3 + L_k + A_prev + period,
    pool_models = "none",
    model_fitter = TrialEmulation::stats_glm_logit(save_path = censor_model_dir)
  )

  trial_pp <- TrialEmulation::calculate_weights(trial_pp)

  trial_pp <- TrialEmulation::set_outcome_model(
    trial_pp,
    adjustment_terms = ~ X1 + X2 + X3 + L_k,
    followup_time_terms = ~ factor(followup_time),
    trial_period_terms = ~ trial_period,
    model_fitter = TrialEmulation::stats_glm_logit(save_path = NA)
  )

  trial_pp <- TrialEmulation::set_expansion_options(
    trial_pp,
    output = TrialEmulation::save_to_datatable(),
    chunk_size = chunk_size,
    first_period = 0L,
    last_period = params$K - params$tau
  )

  trial_pp <- TrialEmulation::expand_trials(trial_pp)
  trial_pp <- TrialEmulation::load_expanded_data(trial_pp)
  trial_pp <- TrialEmulation::fit_msm(trial_pp)

  outcome_df <- TrialEmulation::outcome_data(trial_pp)
  weight_summary <- .summarize_trial_emulation_weights(outcome_df)
  baseline_newdata <- outcome_df[outcome_df$followup_time == 0L, , drop = FALSE]
  predict_times <- 0:(params$tau - 1L)

  predictions <- .predict_trial_emulation(
    trial_pp,
    newdata = baseline_newdata,
    predict_times = predict_times,
    conf_int = conf_int,
    samples = samples,
    type = "cum_inc"
  )

  tau_followup <- params$tau - 1L

  # TrialEmulation returns cumulative incidence for assigned treatment 0,
  # assigned treatment 1, and then the difference.
  risk_g0 <- .extract_trial_emulation_cum_inc(predictions[[1]], tau_followup)
  risk_g1 <- .extract_trial_emulation_cum_inc(predictions[[2]], tau_followup)

  out <- list(
    psi_hat = risk_g1 - risk_g0,
    risk_g1 = risk_g1,
    risk_g0 = risk_g0,
    predict_followup_time = tau_followup,
    predictions = predictions,
    n_rows_input = nrow(te_input),
    n_rows_outcome = nrow(outcome_df),
    weight_mean = weight_summary$weight_mean,
    weight_max = weight_summary$weight_max,
    weight_cv = weight_summary$weight_cv,
    weight_source = weight_summary$weight_source
  )

  if (return_trial_sequence) {
    out$trial_sequence <- trial_pp
    out$te_input <- te_input
  }

  out
}

#-----------------------------------------------------------
# 6c. Target-state construction for the paper MDR estimators
#-----------------------------------------------------------

.direct_mdr_seed <- function(seed, offset) {
  if (is.null(seed)) return(NULL)
  seed + offset
}

.fit_lm_safe <- function(formula, data) {
  fit <- try(
    suppressWarnings(stats::lm(formula = formula, data = data)),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(NULL)
  fit
}

.direct_mdr_transition_resid_sd <- function(fit) {
  resid_sd <- try(stats::sigma(fit), silent = TRUE)
  if (inherits(resid_sd, "try-error") || !is.finite(resid_sd) || resid_sd <= 0) {
    resid_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)
  }
  if (!is.finite(resid_sd) || resid_sd <= 0) {
    resid_sd <- 1e-8
  }

  resid_sd
}

.predict_lm_safe <- function(fit, newdata, fallback) {
  if (is.null(fit)) {
    return(rep(fallback, nrow(newdata)))
  }

  pred <- try(suppressWarnings(stats::predict(fit, newdata = newdata)), silent = TRUE)
  if (inherits(pred, "try-error")) {
    return(rep(fallback, nrow(newdata)))
  }

  pred <- as.numeric(pred)
  pred[!is.finite(pred)] <- fallback
  pred
}

prepare_direct_mdr_observed_state_data <- function(obs_data) {
  L0_data <- obs_data %>%
    dplyr::arrange(id, k) %>%
    dplyr::group_by(id) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(id = id, L0 = L_k)

  obs_data %>%
    dplyr::left_join(L0_data, by = "id") %>%
    dplyr::arrange(id, k) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      L_prev_state = dplyr::lag(L_k),
      L_prev_state = ifelse(is.na(L_prev_state), L_k, L_prev_state),
      L_curr = L_k,
      A_curr = A_k,
      L_next = dplyr::lead(L_k),
      next_k = dplyr::lead(k)
    ) %>%
    dplyr::ungroup()
}

prepare_direct_mdr_transition_data <- function(obs_data) {
  prepare_direct_mdr_observed_state_data(obs_data) %>%
    dplyr::filter(
      !is.na(L_next),
      next_k == k + 1L,
      dY_next == 0L,
      dC_next == 0L
    ) %>%
    dplyr::select(
      id, k, X1, X2, X3, L0,
      L_prev_state, L_curr, A_curr, L_next
    )
}

fit_direct_mdr_transition_model <- function(obs_data) {
  transition_data <- prepare_direct_mdr_transition_data(obs_data)

  if (nrow(transition_data) == 0) {
    stop("No rows available for direct MDR transition model fitting.", call. = FALSE)
  }

  transition_formula <- L_next ~ X1 + X2 + X3 + L0 +
    L_prev_state + L_curr + A_curr + k

  fit <- .fit_lm_safe(transition_formula, data = transition_data)
  if (is.null(fit)) {
    stop("Direct MDR transition model failed.", call. = FALSE)
  }

  list(
    fit = fit,
    formula = transition_formula,
    resid_sd = .direct_mdr_transition_resid_sd(fit),
    fallback_mean = mean(transition_data$L_next, na.rm = TRUE),
    data = transition_data
  )
}

prepare_direct_mdr_event_data <- function(obs_data) {
  prepare_direct_mdr_observed_state_data(obs_data) %>%
    dplyr::filter(
      Y_start == 0L,
      C_start == 0L,
      !is.na(dY_next)
    ) %>%
    dplyr::select(
      id, k, X1, X2, X3, L0,
      L_prev_state, L_curr, A_curr, dY_next
    )
}

fit_direct_mdr_event_model <- function(obs_data) {
  event_data <- prepare_direct_mdr_event_data(obs_data)

  if (nrow(event_data) == 0) {
    stop("No rows available for direct MDR event model fitting.", call. = FALSE)
  }

  event_formula <- dY_next ~ X1 + X2 + X3 + L0 +
    L_prev_state + L_curr + A_curr + k

  fit <- .fit_logistic_safe(event_formula, data = event_data)
  fallback_hazard <- mean(event_data$dY_next, na.rm = TRUE)
  if (!is.finite(fallback_hazard)) fallback_hazard <- 0

  list(
    fit = fit,
    formula = event_formula,
    fallback_hazard = pmin(pmax(fallback_hazard, 0), 1),
    degenerate = is.null(fit),
    data = event_data
  )
}

.predict_direct_mdr_event_hazard <- function(event_model, newdata, eps = 1e-8) {
  fallback <- event_model$fallback_hazard
  if (!is.finite(fallback)) fallback <- 0
  fallback <- pmin(pmax(fallback, 0), 1)

  if (is.null(event_model$fit)) {
    return(rep(fallback, nrow(newdata)))
  }

  .predict_prob_safe(
    event_model$fit,
    newdata = newdata,
    fallback = fallback,
    eps = eps
  )
}

construct_direct_mdr_state_data <- function(obs_data, params, stte_data = NULL) {
  if (is.null(stte_data)) {
    stte_data <- construct_STTE_dataset(obs_data, params)
  }

  L0_data <- obs_data %>%
    dplyr::arrange(id, k) %>%
    dplyr::group_by(id) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(id = id, L0 = L_k)

  stte_data %>%
    dplyr::left_join(L0_data, by = "id") %>%
    dplyr::mutate(
      L0 = ifelse(is.na(L0), L_trial, L0),
      strategy = as.integer(strategy)
    ) %>%
    dplyr::arrange(id_trial, strategy, follow_k) %>%
    dplyr::group_by(id_trial, strategy) %>%
    dplyr::mutate(
      L_prev_state = dplyr::lag(L_k),
      L_prev_state = ifelse(follow_k == 0L | is.na(L_prev_state), L_trial, L_prev_state),
      L_curr = L_k,
      D_adherent = as.integer(artificial_censor == 0L),
      at_risk_start = as.integer(Y_start == 0L & C_start == 0L)
    ) %>%
    dplyr::ungroup()
}

get_direct_mdr_eligible_base <- function(stte_state_data) {
  stte_state_data %>%
    dplyr::filter(follow_k == 0L) %>%
    dplyr::distinct(
      id_trial, id, trial_m, X1, X2, X3,
      L0, L_trial, A_prev_entry
    )
}

simulate_direct_mdr_target_states <- function(
    eligible_base,
    transition_model,
    strategy_z,
    params,
    B = NULL,
    target_sim_multiplier = 5,
    deterministic_transition = FALSE,
    event_model = NULL,
    target_deplete_events = TRUE,
    target_event_eps = 1e-8,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n_base <- nrow(eligible_base)
  if (n_base == 0) {
    stop("No eligible rows available for direct MDR target simulation.", call. = FALSE)
  }

  if (is.null(B)) {
    if (!is.numeric(target_sim_multiplier) || length(target_sim_multiplier) != 1L ||
        !is.finite(target_sim_multiplier) || target_sim_multiplier <= 0) {
      stop("target_sim_multiplier must be one positive finite number.", call. = FALSE)
    }
    B <- max(n_base, as.integer(ceiling(target_sim_multiplier * n_base)))
  }
  B <- .as_positive_integer(B, "Direct MDR target simulation size")

  if (B == n_base) {
    sims <- eligible_base
  } else {
    sims <- eligible_base[sample(seq_len(n_base), size = B, replace = TRUE), , drop = FALSE]
  }

  sim_id <- seq_len(nrow(sims))
  L_prev_state <- sims$L_trial
  L_curr <- sims$L_trial
  active <- rep(TRUE, nrow(sims))
  target_rows <- vector("list", params$tau)

  if (isTRUE(target_deplete_events) && is.null(event_model)) {
    warning(
      "target_deplete_events = TRUE, but no event_model was supplied. ",
      "Direct MDR target states will not be depleted by simulated events.",
      call. = FALSE
    )
    target_deplete_events <- FALSE
  }

  for (follow_j in 0:(params$tau - 1L)) {
    active_index <- which(active)
    if (length(active_index) == 0) break

    calendar_k <- sims$trial_m[active_index] + follow_j

    target_rows[[follow_j + 1L]] <- data.frame(
      sim_id = sim_id[active_index],
      source_id_trial = sims$id_trial[active_index],
      strategy = as.integer(strategy_z),
      follow_k = follow_j,
      trial_m = sims$trial_m[active_index],
      k = calendar_k,
      X1 = sims$X1[active_index],
      X2 = sims$X2[active_index],
      X3 = sims$X3[active_index],
      L0 = sims$L0[active_index],
      L_trial = sims$L_trial[active_index],
      L_prev_state = L_prev_state[active_index],
      L_curr = L_curr[active_index],
      A_prev = if (follow_j == 0L) {
        sims$A_prev_entry[active_index]
      } else {
        rep(as.integer(strategy_z), length(active_index))
      },
      G = 1L
    )

    transition_index <- active_index

    if (isTRUE(target_deplete_events)) {
      event_newdata <- data.frame(
        X1 = sims$X1[active_index],
        X2 = sims$X2[active_index],
        X3 = sims$X3[active_index],
        L0 = sims$L0[active_index],
        L_prev_state = L_prev_state[active_index],
        L_curr = L_curr[active_index],
        A_curr = as.integer(strategy_z),
        k = calendar_k
      )

      event_hazard <- .predict_direct_mdr_event_hazard(
        event_model = event_model,
        newdata = event_newdata,
        eps = target_event_eps
      )
      dY_next <- stats::rbinom(length(active_index), 1, event_hazard)
      event_index <- active_index[dY_next == 1L]
      if (length(event_index) > 0) {
        active[event_index] <- FALSE
      }
      transition_index <- active_index[dY_next == 0L]
    }

    if (follow_j < params$tau - 1L && length(transition_index) > 0) {
      transition_calendar_k <- sims$trial_m[transition_index] + follow_j
      transition_newdata <- data.frame(
        X1 = sims$X1[transition_index],
        X2 = sims$X2[transition_index],
        X3 = sims$X3[transition_index],
        L0 = sims$L0[transition_index],
        L_prev_state = L_prev_state[transition_index],
        L_curr = L_curr[transition_index],
        A_curr = as.integer(strategy_z),
        k = transition_calendar_k
      )

      L_next_mean <- .predict_lm_safe(
        transition_model$fit,
        newdata = transition_newdata,
        fallback = transition_model$fallback_mean
      )

      if (isTRUE(deterministic_transition)) {
        L_next <- L_next_mean
      } else {
        L_next <- stats::rnorm(
          length(L_next_mean),
          mean = L_next_mean,
          sd = transition_model$resid_sd
        )
      }

      L_prev_state[transition_index] <- L_curr[transition_index]
      L_curr[transition_index] <- L_next
    }
  }

  dplyr::bind_rows(target_rows)
}

extract_direct_mdr_adherent_states <- function(stte_state_data, strategy_z) {
  stte_state_data %>%
    dplyr::filter(
      strategy == as.integer(strategy_z),
      D_adherent == 1L,
      at_risk_start == 1L
    ) %>%
    dplyr::mutate(G = 0L)
}

check_target_adherent_cell_support <- function(target_states_z, adherent_states_z, strategy_z) {
  # The most important estimator-support problem is a target cell with no
  # observed adherent rows. Risk standardization can otherwise ask the outcome
  # model to predict a strategy-month cell with no empirical outcome support.
  cell_counts <- dplyr::full_join(
    target_states_z %>%
      dplyr::count(follow_k, name = "n_target"),
    adherent_states_z %>%
      dplyr::count(follow_k, name = "n_adherent"),
    by = "follow_k"
  ) %>%
    dplyr::mutate(
      n_target = dplyr::coalesce(n_target, 0L),
      n_adherent = dplyr::coalesce(n_adherent, 0L),
      strategy = strategy_z,
      target_missing_but_adherent_present = n_target == 0L & n_adherent > 0L,
      target_present_but_adherent_missing = n_target > 0L & n_adherent == 0L,
      sparse_cell = n_target > 0L & n_adherent > 0L & pmin(n_target, n_adherent) < 50L
    ) %>%
    dplyr::select(
      strategy, follow_k, n_target, n_adherent,
      target_missing_but_adherent_present,
      target_present_but_adherent_missing,
      sparse_cell
    ) %>%
    dplyr::arrange(strategy, follow_k)

  bad <- cell_counts %>%
    dplyr::filter(target_present_but_adherent_missing)

  if (nrow(bad) > 0) {
    warning(
      "Some target strategy-follow cells have no adherent observed rows: ",
      paste0("g=", strategy_z, ", k=", bad$follow_k, collapse = "; "),
      call. = FALSE
    )
  }

  cell_counts
}

.check_target_outcome_cell_support <- function(target_states, outcome_data) {
  target_cells <- target_states %>%
    dplyr::distinct(strategy, follow_k)

  observed_cells <- outcome_data %>%
    dplyr::distinct(strategy, follow_k)

  target_cells %>%
    dplyr::anti_join(observed_cells, by = c("strategy", "follow_k")) %>%
    dplyr::arrange(strategy, follow_k)
}

.target_outcome_support_failure_reason <- function(bad_cells) {
  paste0(
    "Target strategy-follow cells have no observed adherent outcome rows: ",
    paste0("g=", bad_cells$strategy, ", k=", bad_cells$follow_k, collapse = "; ")
  )
}

.extract_scalar_or_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(x[1L])
}

.estimate_state_standardized_risks <- function(
    outcome_fit,
    params,
    fallback_hazard,
    target_states = NULL) {

  if (is.null(target_states)) {
    pred_data <- tidyr::expand_grid(
      strategy = c(0L, 1L),
      follow_k = 0:(params$tau - 1L)
    )
  } else {
    pred_data <- target_states
  }

  pred_data <- pred_data %>%
    dplyr::mutate(
      follow_f = factor(follow_k, levels = 0:(params$tau - 1L))
    )

  pred_data$hazard <- .predict_prob_safe(
    outcome_fit,
    pred_data,
    fallback = fallback_hazard,
    eps = 1e-8
  )

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
      "State-standardized risk prediction is missing strategy-follow cells: ",
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

.estimate_strategy_baseline_standardized_risks <- function(
    outcome_fit,
    eligible_base,
    params,
    fallback_hazard) {

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

    pred_grid$hazard <- .predict_prob_safe(
      outcome_fit,
      pred_grid,
      fallback = fallback_hazard,
      eps = 1e-8
    )

    pred_grid
  }) %>%
    dplyr::bind_rows()

  hazard_table <- prediction_rows %>%
    dplyr::group_by(strategy, follow_k) %>%
    dplyr::summarise(
      hazard = mean(hazard, na.rm = TRUE),
      n_baseline = dplyr::n(),
      .groups = "drop"
    )

  risk_table <- prediction_rows %>%
    dplyr::arrange(strategy, id_trial, follow_k) %>%
    dplyr::group_by(strategy, id_trial) %>%
    dplyr::summarise(
      risk = 1 - prod(1 - hazard),
      .groups = "drop"
    ) %>%
    dplyr::group_by(strategy) %>%
    dplyr::summarise(
      risk = mean(risk, na.rm = TRUE),
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
    prediction_rows = prediction_rows
  )
}

.prepare_direct_mdr_state_samples <- function(
    obs_data,
    params,
    stte_data = NULL,
    target_sim_multiplier = 5,
    target_sim_size = NULL,
    deterministic_transition = FALSE,
    target_deplete_events = TRUE,
    target_event_eps = 1e-8,
    seed = NULL) {

  transition_model <- fit_direct_mdr_transition_model(obs_data)
  event_model <- if (isTRUE(target_deplete_events)) {
    fit_direct_mdr_event_model(obs_data)
  } else {
    NULL
  }

  stte_state_data <- construct_direct_mdr_state_data(obs_data, params, stte_data)
  
  # because stte_state_data came from observed data with censoring, 
  # this eligible base is the observed eligible base, not the no-natural-censoring oracle eligible base.
  eligible_base <- get_direct_mdr_eligible_base(stte_state_data)

  target_states <- list()
  adherent_states <- list()

  for (z in c(0L, 1L)) {
    target_states[[as.character(z)]] <- simulate_direct_mdr_target_states(
      eligible_base = eligible_base,
      transition_model = transition_model,
      strategy_z = z,
      params = params,
      B = target_sim_size,
      target_sim_multiplier = target_sim_multiplier,
      deterministic_transition = deterministic_transition,
      event_model = event_model,
      target_deplete_events = target_deplete_events,
      target_event_eps = target_event_eps,
      seed = .direct_mdr_seed(seed, 100L + z)
    )

    adherent_states[[as.character(z)]] <- extract_direct_mdr_adherent_states(
      stte_state_data = stte_state_data,
      strategy_z = z
    )
  }

  list(
    transition_model = transition_model,
    event_model = event_model,
    stte_state_data = stte_state_data,
    eligible_base = eligible_base,
    target_states = target_states,
    adherent_states = adherent_states,
    target_states_all = dplyr::bind_rows(target_states),
    adherent_states_all = dplyr::bind_rows(adherent_states),
    target_deplete_events = target_deplete_events
  )
}

estimate_STTE_per_protocol_ipcw_state_gcomp <- function(
    obs_data,
    params,
    stte_data = NULL,
    target_sim_multiplier = 5,
    target_sim_size = NULL,
    deterministic_transition = FALSE,
    target_deplete_events = TRUE,
    target_event_eps = 1e-8,
    seed = NULL,
    return_models = FALSE) {

  if (params$grace != 0) {
    warning(
      "The IPCW state-gcomp implementation currently simulates ",
      "deterministic always-treated and never-treated target states. ",
      "Check this before using a nonzero grace period."
    )
  }

  if (is.null(stte_data)) {
    stte_data <- construct_STTE_dataset(obs_data, params)
  }

  ipcw_obj <- NULL
  stte_ipcw <- stte_data
  if (!"ipcw" %in% names(stte_ipcw)) {
    ipcw_obj <- estimate_ipcw(stte_data, params)
    stte_ipcw <- ipcw_obj$data
  }

  state_samples <- .prepare_direct_mdr_state_samples(
    obs_data = obs_data,
    params = params,
    stte_data = stte_ipcw,
    target_sim_multiplier = target_sim_multiplier,
    target_sim_size = target_sim_size,
    deterministic_transition = deterministic_transition,
    target_deplete_events = target_deplete_events,
    target_event_eps = target_event_eps,
    seed = seed
  )

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

  out <- list(
    psi_hat = risk_obj$psi_hat,
    risk_g1 = risk_obj$risk_g1,
    risk_g0 = risk_obj$risk_g0,
    n_rows_stte = nrow(state_samples$stte_state_data),
    n_rows_target = nrow(state_samples$target_states_all),
    n_rows_adherent = nrow(state_samples$adherent_states_all),
    n_rows_outcome = nrow(outcome_data),
    weight_mean = weight_mean,
    weight_max = max(outcome_data$ipcw, na.rm = TRUE),
    weight_cv = ifelse(
      abs(weight_mean) < 1e-12,
      NA_real_,
      stats::sd(outcome_data$ipcw, na.rm = TRUE) / weight_mean
    ),
    n_target_without_observed_outcome_cells = n_target_without_observed_outcome_cells,
    failure_reason = NA_character_
  )

  if (return_models) {
    out$ipcw_data <- stte_ipcw
    out$ipcw_object <- ipcw_obj
    out$transition_model <- state_samples$transition_model
    out$event_model <- state_samples$event_model
    out$stte_state_data <- state_samples$stte_state_data
    out$target_states <- state_samples$target_states_all
    out$adherent_states <- state_samples$adherent_states_all
    out$outcome_data <- outcome_data
    out$outcome_fit <- outcome_fit
    out$risk_table <- risk_obj$risk_table
    out$hazard_table <- risk_obj$hazard_table
    out$target_without_observed_outcome_cells <- bad_target_outcome_cells
  }

  out
}

#-----------------------------------------------------------
# 7. Simulation result rows and performance summaries
#-----------------------------------------------------------

.make_replication_result_row <- function(
    replication,
    method,
    est,
    n_rows_input = NA_integer_) {

  data.frame(
    replication = replication,
    method = method,
    psi_hat = ifelse(is.null(est$psi_hat), NA_real_, est$psi_hat),
    risk_g1_hat = ifelse(is.null(est$risk_g1), NA_real_, est$risk_g1),
    risk_g0_hat = ifelse(is.null(est$risk_g0), NA_real_, est$risk_g0),
    estimate_se = ifelse(is.null(est$estimate_se), NA_real_, est$estimate_se),
    ci_lower = ifelse(is.null(est$ci_lower), NA_real_, est$ci_lower),
    ci_upper = ifelse(is.null(est$ci_upper), NA_real_, est$ci_upper),
    n_rows_input = n_rows_input,
    n_rows_stte = ifelse(is.null(est$n_rows_stte), NA_integer_, est$n_rows_stte),
    n_rows_outcome = ifelse(is.null(est$n_rows_outcome), NA_integer_, est$n_rows_outcome),
    n_unique_ids = ifelse(is.null(est$n_unique_ids), NA_integer_, est$n_unique_ids),
    weight_mean = ifelse(is.null(est$weight_mean), NA_real_, est$weight_mean),
    weight_max = ifelse(is.null(est$weight_max), NA_real_, est$weight_max),
    weight_cv = ifelse(is.null(est$weight_cv), NA_real_, est$weight_cv),
    n_target_missing_cells = ifelse(
      is.null(est$n_target_missing_cells),
      NA_integer_,
      est$n_target_missing_cells
    ),
    n_target_present_adherent_missing_cells = ifelse(
      is.null(est$n_target_present_adherent_missing_cells),
      NA_integer_,
      est$n_target_present_adherent_missing_cells
    ),
    n_target_without_observed_outcome_cells = ifelse(
      is.null(est$n_target_without_observed_outcome_cells),
      NA_integer_,
      est$n_target_without_observed_outcome_cells
    ),
    min_cum_g_unbounded = ifelse(
      is.null(est$min_cum_g_unbounded),
      NA_real_,
      est$min_cum_g_unbounded
    ),
    prop_cum_g_bounded = ifelse(
      is.null(est$prop_cum_g_bounded),
      NA_real_,
      est$prop_cum_g_bounded
    ),
    max_inverse_cum_g = ifelse(
      is.null(est$max_inverse_cum_g),
      NA_real_,
      est$max_inverse_cum_g
    ),
    min_cum_g_unbounded_used = ifelse(
      is.null(est$min_cum_g_unbounded_used),
      NA_real_,
      est$min_cum_g_unbounded_used
    ),
    prop_cum_g_bounded_used = ifelse(
      is.null(est$prop_cum_g_bounded_used),
      NA_real_,
      est$prop_cum_g_bounded_used
    ),
    max_inverse_cum_g_used = ifelse(
      is.null(est$max_inverse_cum_g_used),
      NA_real_,
      est$max_inverse_cum_g_used
    ),
    failure_reason = ifelse(is.null(est$failure_reason), NA_character_, est$failure_reason),
    stringsAsFactors = FALSE
  )
}

summarize_performance <- function(rep_results, psi_true) {
  if (!"method" %in% names(rep_results)) {
    rep_results$method <- "ipcw_pool"
  }
  if (!"n_target_missing_cells" %in% names(rep_results)) {
    rep_results$n_target_missing_cells <- NA_integer_
  }
  if (!"n_target_present_adherent_missing_cells" %in% names(rep_results)) {
    rep_results$n_target_present_adherent_missing_cells <- NA_integer_
  }
  if (!"n_target_without_observed_outcome_cells" %in% names(rep_results)) {
    rep_results$n_target_without_observed_outcome_cells <- NA_integer_
  }

  mean_or_na <- function(x) {
    out <- mean(x, na.rm = TRUE)
    ifelse(is.nan(out), NA_real_, out)
  }

  dplyr::bind_rows(lapply(split(rep_results, rep_results$method), function(method_results) {
    valid <- !is.na(method_results$psi_hat)
    psi_hat <- method_results$psi_hat[valid]

    if (length(psi_hat) == 0) {
      mean_psi_hat <- NA_real_
      bias <- NA_real_
      relative_bias <- NA_real_
      empirical_se <- NA_real_
      rmse <- NA_real_
    } else {
      mean_psi_hat <- mean(psi_hat)
      bias <- mean(psi_hat - psi_true)
      relative_bias <- ifelse(abs(psi_true) < 1e-12, NA_real_, bias / psi_true)
      empirical_se <- stats::sd(psi_hat)
      rmse <- sqrt(mean((psi_hat - psi_true)^2))
    }

    data.frame(
      method = method_results$method[1],
      n_replications = nrow(method_results),
      n_success = sum(valid),
      failure_rate = mean(!valid),
      psi_true = psi_true,
      mean_psi_hat = mean_psi_hat,
      bias = bias,
      relative_bias = relative_bias,
      empirical_se = empirical_se,
      rmse = rmse,
      mean_weight_max = mean_or_na(method_results$weight_max),
      mean_weight_cv = mean_or_na(method_results$weight_cv),
      mean_target_missing_cells = mean_or_na(method_results$n_target_missing_cells),
      mean_target_present_adherent_missing_cells = mean_or_na(
        method_results$n_target_present_adherent_missing_cells
      ),
      mean_target_without_observed_outcome_cells = mean_or_na(
        method_results$n_target_without_observed_outcome_cells
      ),
      stringsAsFactors = FALSE
    )
  }))
}

#-----------------------------------------------------------
# 8. Visualize scenario characteristics and simulation results
#-----------------------------------------------------------

available_scenarios <- function() {
  c(
    "base",
    "small_size",
    "null_effect",
    "weak_confounding",
    "strong_confounding",
    "poor_positivity",
    "poor_adherence",
    "rare_event",
    "nonlinear"
  )
}

load_visualization_packages <- function() {
  load_required_packages()

  needed <- c("ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Please install the following package before using plotting helpers: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(ggplot2)
  })
}

summarize_scenario_parameters <- function(
    scenario_names = available_scenarios(),
    parameter_names = c(
      "K", "tau", "grace", "require_untreated_entry",
      "observed_sample_size",
      "alpha0", "alpha_Aprev", "alpha_L", "alpha_time",
      "beta0", "beta_A", "beta_L", "beta_time",
      "delta0", "delta_A", "delta_L",
      "L_Aprev", "nonlinear", "alpha_L2", "beta_L2", "beta_A_L",
      "prob_clip", "treatment_prob_clip", "outcome_prob_clip",
      "censor_prob_clip", "weight_truncation_quantile"
    )) {

  rows <- lapply(scenario_names, function(scenario_name) {
    params <- define_parameters(scenario_name)

    values <- lapply(parameter_names, function(parameter_name) {
      value <- params[[parameter_name]]
      if (is.null(value)) {
        NA
      } else {
        value
      }
    })
    names(values) <- parameter_names

    data.frame(
      scenario = scenario_name,
      as.data.frame(values, check.names = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

summarize_scenario_observed_characteristics <- function(
    scenario_names = available_scenarios(),
    n = 5000,
    seed = 2026,
    simulate_censoring = TRUE) {

  load_required_packages()

  rows <- lapply(seq_along(scenario_names), function(scenario_index) {
    scenario_name <- scenario_names[[scenario_index]]
    params <- define_parameters(scenario_name)

    obs_data <- generate_observed_data(
      n = n,
      params = params,
      seed = seed + scenario_index,
      simulate_censoring = simulate_censoring
    )

    elig <- get_eligible_person_trials(obs_data, params)

    first_rows <- obs_data %>%
      dplyr::arrange(id, k) %>%
      dplyr::group_by(id) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup()

    person_summary <- obs_data %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        ever_treated = any(A_k == 1L, na.rm = TRUE),
        observed_event = any(dY_next == 1L, na.rm = TRUE),
        observed_censor = any(dC_next == 1L, na.rm = TRUE),
        months_observed = dplyr::n(),
        .groups = "drop"
      )

    data.frame(
      scenario = scenario_name,
      n_people = dplyr::n_distinct(obs_data$id),
      n_rows_observed = nrow(obs_data),
      eligible_person_trials = nrow(elig),
      eligible_trials_per_person = nrow(elig) / dplyr::n_distinct(obs_data$id),
      mean_L0 = mean(first_rows$L_k, na.rm = TRUE),
      sd_L0 = stats::sd(first_rows$L_k, na.rm = TRUE),
      mean_L_observed = mean(obs_data$L_k, na.rm = TRUE),
      treatment_interval_rate = mean(obs_data$A_k, na.rm = TRUE),
      outcome_interval_rate = mean(obs_data$dY_next, na.rm = TRUE),
      censor_interval_rate = mean(obs_data$dC_next, na.rm = TRUE),
      ever_treated_rate = mean(person_summary$ever_treated, na.rm = TRUE),
      observed_event_rate = mean(person_summary$observed_event, na.rm = TRUE),
      observed_censor_rate = mean(person_summary$observed_censor, na.rm = TRUE),
      mean_months_observed = mean(person_summary$months_observed, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

.mean_or_na <- function(x) {
  out <- mean(x, na.rm = TRUE)
  ifelse(is.nan(out), NA_real_, out)
}

.sd_or_na <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  ifelse(is.nan(out), NA_real_, out)
}

.median_or_na <- function(x) {
  out <- stats::median(x, na.rm = TRUE)
  ifelse(is.nan(out), NA_real_, out)
}

.quantile_or_na <- function(x, probs) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

.first_observed_k <- function(k, indicator) {
  index <- which(indicator == 1L)
  if (length(index) == 0) {
    return(NA_integer_)
  }

  as.integer(min(k[index], na.rm = TRUE))
}

.ensure_scenario_column <- function(obs_data, scenario_name = "scenario") {
  if (!"scenario" %in% names(obs_data)) {
    obs_data$scenario <- scenario_name
  }

  obs_data
}

generate_obs_data_by_scenario <- function(
    scenario_names = available_scenarios(),
    n = 5000,
    seed = 2026,
    simulate_censoring = TRUE) {

  load_required_packages()

  rows <- lapply(seq_along(scenario_names), function(scenario_index) {
    scenario_name <- scenario_names[[scenario_index]]
    params <- define_parameters(scenario_name)

    generate_observed_data(
      n = n,
      params = params,
      seed = seed + scenario_index,
      simulate_censoring = simulate_censoring
    ) %>%
      dplyr::mutate(scenario = scenario_name, .before = 1)
  })

  dplyr::bind_rows(rows)
}

summarize_obs_data_person_level <- function(obs_data) {
  load_required_packages()
  obs_data <- .ensure_scenario_column(obs_data)

  baseline_rows <- obs_data %>%
    dplyr::arrange(scenario, id, k) %>%
    dplyr::group_by(scenario, id) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      scenario = scenario,
      id = id,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      L0_observed = L_k
    )

  obs_data %>%
    dplyr::group_by(scenario, id) %>%
    dplyr::summarise(
      months_observed = dplyr::n(),
      ever_treated = any(A_k == 1L, na.rm = TRUE),
      first_treatment_k = .first_observed_k(k, A_k),
      observed_event = any(dY_next == 1L, na.rm = TRUE),
      event_k = .first_observed_k(k, dY_next),
      observed_censor = any(dC_next == 1L, na.rm = TRUE),
      censor_k = .first_observed_k(k, dC_next),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      baseline_rows,
      by = c("scenario", "id")
    )
}

summarize_obs_data_scenario_table <- function(obs_data) {
  load_required_packages()
  obs_data <- .ensure_scenario_column(obs_data)

  person_summary <- summarize_obs_data_person_level(obs_data)

  interval_summary <- obs_data %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      person_months_observed = dplyr::n(),
      mean_L_observed = .mean_or_na(L_k),
      p10_L_observed = .quantile_or_na(L_k, 0.10),
      p90_L_observed = .quantile_or_na(L_k, 0.90),
      treatment_interval_rate = .mean_or_na(A_k),
      outcome_interval_rate = .mean_or_na(dY_next),
      censor_interval_rate = .mean_or_na(dC_next),
      .groups = "drop"
    )

  person_table <- person_summary %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_people = dplyr::n(),
      mean_months_observed = .mean_or_na(months_observed),
      x1_rate = .mean_or_na(X1),
      x2_mean = .mean_or_na(X2),
      x2_sd = .sd_or_na(X2),
      x3_rate = .mean_or_na(X3),
      mean_L0_observed = .mean_or_na(L0_observed),
      sd_L0_observed = .sd_or_na(L0_observed),
      ever_treated_rate = .mean_or_na(ever_treated),
      median_first_treatment_k = .median_or_na(first_treatment_k),
      observed_event_rate = .mean_or_na(observed_event),
      median_event_k = .median_or_na(event_k),
      observed_censor_rate = .mean_or_na(observed_censor),
      median_censor_k = .median_or_na(censor_k),
      .groups = "drop"
    )

  dplyr::left_join(
    person_table,
    interval_summary,
    by = "scenario"
  ) %>%
    dplyr::select(
      scenario,
      n_people,
      person_months_observed,
      mean_months_observed,
      x1_rate,
      x2_mean,
      x2_sd,
      x3_rate,
      mean_L0_observed,
      sd_L0_observed,
      mean_L_observed,
      p10_L_observed,
      p90_L_observed,
      treatment_interval_rate,
      ever_treated_rate,
      median_first_treatment_k,
      outcome_interval_rate,
      observed_event_rate,
      median_event_k,
      censor_interval_rate,
      observed_censor_rate,
      median_censor_k
    )
}

summarize_obs_data_time_profile <- function(obs_data) {
  load_required_packages()
  obs_data <- .ensure_scenario_column(obs_data)

  n_people <- obs_data %>%
    dplyr::distinct(scenario, id) %>%
    dplyr::count(scenario, name = "n_people")

  obs_data %>%
    dplyr::group_by(scenario, k) %>%
    dplyr::summarise(
      n_at_risk = dplyr::n(),
      mean_L_k = .mean_or_na(L_k),
      p10_L_k = .quantile_or_na(L_k, 0.10),
      p90_L_k = .quantile_or_na(L_k, 0.90),
      treatment_rate = .mean_or_na(A_k),
      outcome_hazard = .mean_or_na(dY_next),
      censor_hazard = .mean_or_na(dC_next),
      event_count = sum(dY_next == 1L, na.rm = TRUE),
      censor_count = sum(dC_next == 1L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario, k) %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(
      cumulative_events = cumsum(event_count),
      cumulative_censors = cumsum(censor_count)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(n_people, by = "scenario") %>%
    dplyr::mutate(
      risk_set_fraction = n_at_risk / n_people,
      cumulative_event_rate = cumulative_events / n_people,
      cumulative_censor_rate = cumulative_censors / n_people
    )
}

make_obs_data_comparison_tables <- function(
    obs_data,
    time_points = c(0, 3, 6, 9, 12, 18, 23),
    time_metrics = c(
      "n_at_risk",
      "risk_set_fraction",
      "mean_L_k",
      "treatment_rate",
      "outcome_hazard",
      "censor_hazard",
      "cumulative_event_rate",
      "cumulative_censor_rate"
    )) {

  scenario_table <- summarize_obs_data_scenario_table(obs_data)
  time_profile <- summarize_obs_data_time_profile(obs_data)

  time_table <- time_profile %>%
    dplyr::filter(k %in% time_points) %>%
    dplyr::select(dplyr::all_of(c("scenario", "k", time_metrics))) %>%
    dplyr::arrange(scenario, k)

  list(
    overall = scenario_table,
    by_time = time_table,
    time_profile = time_profile
  )
}

make_scenario_probability_grid <- function(
    scenario_names = available_scenarios(),
    L_values = seq(-3, 3, length.out = 121),
    k_values = c(0, 6, 12),
    x1 = 0,
    x2 = 0,
    x3 = 0,
    A_prev_values = c(0L, 1L),
    A_values = c(0L, 1L)) {

  rows <- list()
  row_index <- 1L

  for (scenario_name in scenario_names) {
    params <- define_parameters(scenario_name)

    for (k_value in k_values) {
      for (A_prev_value in A_prev_values) {
        rows[[row_index]] <- data.frame(
          scenario = scenario_name,
          component = "Treatment probability",
          conditioning = paste0("A_prev = ", A_prev_value),
          k = k_value,
          L_k = L_values,
          probability = compute_pA(
            A_prev = rep(A_prev_value, length(L_values)),
            L_k = L_values,
            X1 = rep(x1, length(L_values)),
            X2 = rep(x2, length(L_values)),
            X3 = rep(x3, length(L_values)),
            k = k_value,
            params = params
          ),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }

      for (A_value in A_values) {
        rows[[row_index]] <- data.frame(
          scenario = scenario_name,
          component = "Outcome hazard",
          conditioning = paste0("A_k = ", A_value),
          k = k_value,
          L_k = L_values,
          probability = compute_pY(
            A_k = rep(A_value, length(L_values)),
            L_k = L_values,
            X1 = rep(x1, length(L_values)),
            X2 = rep(x2, length(L_values)),
            X3 = rep(x3, length(L_values)),
            k = k_value,
            params = params
          ),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L

        rows[[row_index]] <- data.frame(
          scenario = scenario_name,
          component = "Censoring hazard",
          conditioning = paste0("A_k = ", A_value),
          k = k_value,
          L_k = L_values,
          probability = compute_pC(
            A_k = rep(A_value, length(L_values)),
            L_k = L_values,
            X1 = rep(x1, length(L_values)),
            X2 = rep(x2, length(L_values)),
            X3 = rep(x3, length(L_values)),
            k = k_value,
            params = params
          ),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }

  dplyr::bind_rows(rows)
}

plot_scenario_probability_curves <- function(
    probability_grid,
    components = c("Treatment probability", "Outcome hazard", "Censoring hazard")) {

  load_visualization_packages()

  plot_data <- probability_grid %>%
    dplyr::filter(component %in% components)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = L_k,
      y = probability,
      color = scenario,
      linetype = conditioning
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(component ~ k, scales = "free_y", labeller = ggplot2::label_both) +
    ggplot2::labs(
      title = "Scenario probability curves",
      x = "Time-varying confounder L_k",
      y = "Probability",
      color = "Scenario",
      linetype = "Conditioning"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_scenario_parameter_heatmap <- function(
    parameter_table = NULL,
    scenario_names = available_scenarios(),
    parameters = c(
      "alpha_Aprev", "alpha_L", "beta_A", "beta_L", "beta0",
      "delta_L", "L_Aprev", "nonlinear", "alpha_L2", "beta_L2", "beta_A_L"
    ),
    standardize = TRUE) {

  load_visualization_packages()

  if (is.null(parameter_table)) {
    parameter_table <- summarize_scenario_parameters(scenario_names = scenario_names)
  }

  plot_data <- parameter_table %>%
    dplyr::select(dplyr::all_of(c("scenario", parameters))) %>%
    tidyr::pivot_longer(
      cols = -scenario,
      names_to = "parameter",
      values_to = "value"
    ) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::group_by(parameter) %>%
    dplyr::mutate(
      fill_value = if (isTRUE(standardize)) {
        value_sd <- stats::sd(value, na.rm = TRUE)
        ifelse(
          is.na(value_sd) | value_sd < 1e-12,
          0,
          (value - mean(value, na.rm = TRUE)) / value_sd
        )
      } else {
        value
      },
      label = ifelse(is.na(value), "NA", sprintf("%.2f", value))
    ) %>%
    dplyr::ungroup()

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = scenario, y = parameter, fill = fill_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      na.value = "grey85"
    ) +
    ggplot2::labs(
      title = "Scenario-defining parameters",
      x = NULL,
      y = NULL,
      fill = ifelse(isTRUE(standardize), "Scaled value", "Value")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

plot_scenario_observed_characteristics <- function(
    characteristics,
    metrics = c(
      "eligible_trials_per_person",
      "treatment_interval_rate",
      "ever_treated_rate",
      "observed_event_rate",
      "observed_censor_rate",
      "mean_months_observed"
    )) {

  load_visualization_packages()

  plot_data <- characteristics %>%
    dplyr::select(dplyr::all_of(c("scenario", metrics))) %>%
    tidyr::pivot_longer(
      cols = -scenario,
      names_to = "metric",
      values_to = "value"
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = scenario, y = value, fill = scenario)
  ) +
    ggplot2::geom_col(show.legend = FALSE, width = 0.75) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(
      title = "Observed-data characteristics by scenario",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_obs_data_time_profiles <- function(
    time_profile,
    metrics = c(
      "risk_set_fraction",
      "mean_L_k",
      "treatment_rate",
      "outcome_hazard",
      "censor_hazard",
      "cumulative_event_rate",
      "cumulative_censor_rate"
    )) {

  load_visualization_packages()

  if (!all(c("scenario", "k", metrics) %in% names(time_profile))) {
    time_profile <- summarize_obs_data_time_profile(time_profile)
  }

  plot_data <- time_profile %>%
    dplyr::select(dplyr::all_of(c("scenario", "k", metrics))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metrics),
      names_to = "metric",
      values_to = "value"
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = k, y = value, color = scenario)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.3, alpha = 0.75) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(
      title = "Observed data time profiles by scenario",
      x = "Calendar month k",
      y = NULL,
      color = "Scenario"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_obs_data_L_distribution <- function(
    obs_data,
    k_values = c(0, 6, 12, 18),
    plot_type = c("density", "boxplot")) {

  load_visualization_packages()
  plot_type <- match.arg(plot_type)

  obs_data <- .ensure_scenario_column(obs_data)
  plot_data <- obs_data %>%
    dplyr::filter(k %in% k_values) %>%
    dplyr::mutate(
      k_label = factor(
        paste0("k = ", k),
        levels = paste0("k = ", k_values)
      )
    )

  if (plot_type == "density") {
    return(
      ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = L_k, color = scenario, fill = scenario)
      ) +
        ggplot2::geom_density(alpha = 0.16, linewidth = 0.8) +
        ggplot2::facet_wrap(~ k_label, scales = "free_y") +
        ggplot2::labs(
          title = "Observed L_k distributions by scenario",
          x = "Time-varying confounder L_k",
          y = "Density",
          color = "Scenario",
          fill = "Scenario"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank()
        )
    )
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = scenario, y = L_k, fill = scenario)
  ) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.25, show.legend = FALSE) +
    ggplot2::facet_wrap(~ k_label, scales = "free_y") +
    ggplot2::labs(
      title = "Observed L_k distributions by scenario",
      x = NULL,
      y = "Time-varying confounder L_k"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.extract_simulation_replication_results <- function(sim) {
  if (!is.null(sim$outputs)) {
    return(dplyr::bind_rows(lapply(sim$outputs, function(output) {
      rep_results <- output$rep_results
      rep_results$scenario <- output$scenario
      rep_results$psi_true <- output$truth$psi_true
      rep_results$risk_g1_true <- output$truth$risk_g1
      rep_results$risk_g0_true <- output$truth$risk_g0
      rep_results
    })))
  }

  sim
}

.extract_simulation_performance_table <- function(sim) {
  if (!is.null(sim$performance_table)) {
    return(sim$performance_table)
  }

  sim
}

plot_estimates_by_scenario <- function(
    sim,
    methods = NULL,
    include_truth = TRUE) {

  load_visualization_packages()

  rep_results <- .extract_simulation_replication_results(sim)
  if (!is.null(methods)) {
    rep_results <- rep_results %>% dplyr::filter(method %in% methods)
  }

  rep_results <- rep_results %>%
    dplyr::filter(!is.na(psi_hat))

  truth_rows <- rep_results %>%
    dplyr::distinct(scenario, psi_true)

  p <- ggplot2::ggplot(
    rep_results,
    ggplot2::aes(x = method, y = psi_hat, fill = method)
  ) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.35) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.12, height = 0),
      alpha = 0.35,
      size = 1.4,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ scenario, scales = "free_y") +
    ggplot2::labs(
      title = "Estimated per-protocol risk difference by scenario",
      x = NULL,
      y = "Estimated risk difference",
      fill = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(include_truth) && "psi_true" %in% names(rep_results)) {
    p <- p +
      ggplot2::geom_hline(
        data = truth_rows,
        ggplot2::aes(yintercept = psi_true),
        linetype = "dashed",
        linewidth = 0.5
      )
  }

  p
}

plot_performance_by_scenario <- function(
    sim,
    metrics = c(
      "bias",
      "rmse",
      "failure_rate",
      "mean_weight_cv",
      "mean_weight_max",
      "mean_target_present_adherent_missing_cells",
      "mean_target_without_observed_outcome_cells"
    ),
    methods = NULL) {

  load_visualization_packages()

  performance_table <- .extract_simulation_performance_table(sim)
  if (!is.null(methods)) {
    performance_table <- performance_table %>% dplyr::filter(method %in% methods)
  }
  metrics <- intersect(metrics, names(performance_table))
  if (length(metrics) == 0L) {
    stop("No requested performance metrics are available in the simulation object.", call. = FALSE)
  }

  plot_data <- performance_table %>%
    dplyr::select(dplyr::all_of(c("scenario", "method", metrics))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metrics),
      names_to = "metric",
      values_to = "value"
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = scenario, y = value, fill = method)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(
      title = "Simulation performance by scenario",
      x = NULL,
      y = NULL,
      fill = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_weight_behavior_by_scenario <- function(
    sim,
    metrics = c("weight_cv", "weight_max"),
    methods = NULL) {

  load_visualization_packages()

  rep_results <- .extract_simulation_replication_results(sim)
  if (!is.null(methods)) {
    rep_results <- rep_results %>% dplyr::filter(method %in% methods)
  }

  plot_data <- rep_results %>%
    dplyr::select(dplyr::all_of(c("scenario", "method", metrics))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value))

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = method, y = value, fill = method)
  ) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.35) +
    ggplot2::facet_grid(metric ~ scenario, scales = "free_y") +
    ggplot2::labs(
      title = "Weight behavior by scenario",
      x = NULL,
      y = NULL,
      fill = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}
