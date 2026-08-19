############################################################
# MDR weights for the paper implementation
#
# This module implements the follow-balanced, cell-normalized MDR weights
# shared by the four MDR estimators reported in the paper:
#   - mdr_state_gcomp_follow_balanced_cell_normalized
#   - mdr_pool_follow_balanced_cell_normalized
#   - dr_mdr_state_gcomp_follow_balanced_cell_normalized
#   - dr_mdr_pool_follow_balanced_cell_normalized
#
# State construction and outcome estimation remain in
# stte_per_protocol_simulation.R and stte_doubly_robust_helpers.R.
############################################################

.mdr_weight_script_dir <- local({
  frame_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
  if (is.character(frame_file) && length(frame_file) == 1L && !is.na(frame_file)) {
    return(dirname(normalizePath(frame_file, winslash = "/", mustWork = TRUE)))
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )))
  }

  getwd()
})

# Make the focused weight module safe to source directly as well as through
# stte_simulation_all.R.
if (!exists("load_required_packages", mode = "function") ||
    !exists("check_target_adherent_cell_support", mode = "function") ||
    !exists(".fit_logistic_safe", mode = "function")) {
  source(file.path(.mdr_weight_script_dir, "stte_per_protocol_simulation.R"))
}
load_required_packages()

mdr_weight_methods <- "follow_balanced_cell_normalized"

mdr_weight_method_aliases <- c(
  solution2_follow_balanced_normalized = "follow_balanced_cell_normalized"
)

mdr_weight_classifier_methods <- c(
  "logistic_regression",
  "random_forest"
)

mdr_weight_classifier_method_aliases <- c(
  glm = "logistic_regression",
  logistic = "logistic_regression",
  logit = "logistic_regression",
  current_logistic = "logistic_regression",
  rf = "random_forest",
  randomForest = "random_forest"
)

.mdr_weight_classifier_method_seed_offsets <- c(
  logistic_regression = 0L,
  random_forest = 10000L
)

.mdr_variant_seed <- function(seed, offset) {
  if (is.null(seed)) return(NULL)
  as.integer(seed + offset)
}

.mdr_weight_classifier_seed <- function(seed, classifier_method) {
  if (is.null(seed)) return(NULL)
  method_offset <- .mdr_weight_classifier_method_seed_offsets[[classifier_method]]
  if (is.null(method_offset) || is.na(method_offset)) method_offset <- 0L
  as.integer(seed + method_offset)
}

.check_mdr_weight_cap <- function(weight_cap) {
  if (!is.null(weight_cap) &&
      (!is.numeric(weight_cap) || length(weight_cap) != 1L ||
       !is.finite(weight_cap) || weight_cap <= 0)) {
    stop("weight_cap must be NULL or one positive finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

.normalize_mdr_weight_vector <- function(weights, cap = NULL) {
  weights <- as.numeric(weights)
  if (length(weights) == 0L) return(weights)
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("MDR weights must be finite and nonnegative before normalization.", call. = FALSE)
  }
  if (all(weights == 0)) {
    stop("MDR weights cannot all be zero within a strategy-follow cell.", call. = FALSE)
  }

  if (is.null(cap)) {
    return(weights / mean(weights))
  }

  .check_mdr_weight_cap(cap)
  if (cap < 1) {
    stop(
      "A cell-normalized weight cap must be at least 1 because cell weights have mean 1.",
      call. = FALSE
    )
  }

  max_attainable_mean <- cap * mean(weights > 0)
  if (max_attainable_mean < 1 - 1e-12) {
    stop(
      "The requested weight cap cannot produce mean-one weights in this cell.",
      call. = FALSE
    )
  }
  if (cap == 1) return(rep(1, length(weights)))
  if (abs(max_attainable_mean - 1) <= 1e-12) {
    return(ifelse(weights > 0, cap, 0))
  }

  capped_mean_minus_one <- function(scale) {
    mean(pmin(scale * weights, cap)) - 1
  }
  upper <- max(1, 1 / mean(weights))
  while (capped_mean_minus_one(upper) < 0) {
    upper <- upper * 2
    if (!is.finite(upper)) {
      stop("Unable to normalize MDR weights under the requested cap.", call. = FALSE)
    }
  }

  scale <- stats::uniroot(
    capped_mean_minus_one,
    interval = c(0, upper),
    tol = 1e-12
  )$root
  pmin(scale * weights, cap)
}

validate_mdr_weight_methods <- function(methods) {
  if (length(methods) == 0L) {
    stop("At least one MDR weight method must be supplied.", call. = FALSE)
  }

  methods <- as.character(methods)
  alias_match <- mdr_weight_method_aliases[methods]
  methods <- ifelse(!is.na(alias_match), unname(alias_match), methods)

  bad_methods <- setdiff(methods, mdr_weight_methods)
  if (length(bad_methods) > 0L) {
    stop(
      "Unknown MDR weight method(s): ",
      paste(bad_methods, collapse = ", "),
      ". This paper implementation supports only: ",
      mdr_weight_methods,
      call. = FALSE
    )
  }

  unique(methods)
}

validate_mdr_weight_classifier_method <- function(classifier_method) {
  if (length(classifier_method) == 0L) {
    stop("At least one MDR classifier method must be supplied.", call. = FALSE)
  }

  classifier_method <- as.character(classifier_method)[[1L]]
  alias_match <- mdr_weight_classifier_method_aliases[classifier_method]
  if (!is.na(alias_match)) {
    classifier_method <- unname(alias_match)
  }

  if (!(classifier_method %in% mdr_weight_classifier_methods)) {
    stop(
      "Unknown MDR classifier method: ",
      classifier_method,
      ". Valid methods are: ",
      paste(mdr_weight_classifier_methods, collapse = ", "),
      call. = FALSE
    )
  }

  classifier_method
}

sample_target_by_follow <- function(target_states_z, adherent_states_z, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  follow_levels <- sort(unique(adherent_states_z$follow_k))

  sampled <- lapply(follow_levels, function(k) {
    target_k <- target_states_z[target_states_z$follow_k == k, , drop = FALSE]
    n_k <- sum(adherent_states_z$follow_k == k)

    if (nrow(target_k) == 0L || n_k == 0L) {
      return(target_k[0, , drop = FALSE])
    }

    target_k[
      sample(
        seq_len(nrow(target_k)),
        size = n_k,
        replace = nrow(target_k) < n_k
      ),
      ,
      drop = FALSE
    ]
  })

  dplyr::bind_rows(sampled)
}

normalize_mdr_weights_by_cell <- function(weighted_rows, cap = NULL) {
  .check_mdr_weight_cap(cap)

  weighted_rows %>%
    dplyr::group_by(strategy, follow_k) %>%
    dplyr::mutate(
      w_mdr_joint = w_mdr_raw,
      cell_mean_raw = mean(w_mdr_raw, na.rm = TRUE),
      w_mdr = .normalize_mdr_weight_vector(w_mdr_raw, cap = cap)
    ) %>%
    dplyr::ungroup()
}

.mdr_weight_random_forest_engine <- function(preferred = c("randomForest", "ranger")) {
  preferred <- match.arg(preferred)
  if (preferred == "randomForest" && requireNamespace("randomForest", quietly = TRUE)) {
    return("randomForest")
  }
  if (preferred == "ranger" && requireNamespace("ranger", quietly = TRUE)) {
    return("ranger")
  }
  if (requireNamespace("randomForest", quietly = TRUE)) {
    return("randomForest")
  }
  if (requireNamespace("ranger", quietly = TRUE)) {
    return("ranger")
  }

  stop(
    "Random forest MDR weights require package 'randomForest' or 'ranger'.",
    call. = FALSE
  )
}

.mdr_weight_classifier_feature_formula <- function(include_follow = TRUE) {
  if (isTRUE(include_follow)) {
    return(~ factor(follow_k) + trial_m +
             X1 + X2 + X3 +
             L0 + L_prev_state + L_curr +
             I(L_curr - L_prev_state) +
             I(L_curr^2) + I(L_curr^3) - 1)
  }

  ~ trial_m + X1 + X2 + X3 + L0 + L_prev_state + L_curr - 1
}

.mdr_weight_classifier_model_matrix <- function(data, feature_formula, feature_cols = NULL) {
  x <- stats::model.matrix(feature_formula, data = data)
  x <- as.matrix(x)
  colnames(x) <- make.names(colnames(x), unique = TRUE)

  if (is.null(feature_cols)) {
    return(x)
  }

  missing_cols <- setdiff(feature_cols, colnames(x))
  if (length(missing_cols) > 0L) {
    x <- cbind(
      x,
      matrix(
        0,
        nrow = nrow(x),
        ncol = length(missing_cols),
        dimnames = list(NULL, missing_cols)
      )
    )
  }

  extra_cols <- setdiff(colnames(x), feature_cols)
  if (length(extra_cols) > 0L) {
    x <- x[, setdiff(colnames(x), extra_cols), drop = FALSE]
  }

  x[, feature_cols, drop = FALSE]
}

.fit_mdr_weight_random_forest_ratio_classifier <- function(
    ratio_data,
    include_follow = TRUE,
    seed = NULL,
    ntree = 500,
    mtry = NULL,
    min_node_size = NULL,
    engine = c("randomForest", "ranger")) {

  engine <- .mdr_weight_random_forest_engine(match.arg(engine))
  feature_formula <- .mdr_weight_classifier_feature_formula(include_follow = include_follow)
  x <- .mdr_weight_classifier_model_matrix(ratio_data, feature_formula)
  y <- factor(ratio_data$G, levels = c(0L, 1L))

  if (!is.null(seed)) set.seed(seed)

  if (engine == "randomForest") {
    rf_args <- list(x = x, y = y, ntree = ntree)
    if (!is.null(mtry)) rf_args$mtry <- mtry
    if (!is.null(min_node_size)) rf_args$nodesize <- min_node_size

    fit <- try(
      suppressWarnings(do.call(randomForest::randomForest, rf_args)),
      silent = TRUE
    )
  } else {
    train_df <- as.data.frame(x)
    train_df$G <- y
    rf_args <- list(
      formula = G ~ .,
      data = train_df,
      probability = TRUE,
      num.trees = ntree,
      seed = seed
    )
    if (!is.null(mtry)) rf_args$mtry <- mtry
    if (!is.null(min_node_size)) rf_args$min.node.size <- min_node_size

    fit <- try(
      suppressWarnings(do.call(ranger::ranger, rf_args)),
      silent = TRUE
    )
  }

  if (inherits(fit, "try-error") || is.null(fit)) {
    stop("MDR random forest ratio classifier failed.", call. = FALSE)
  }

  list(
    classifier_method = "random_forest",
    fit = fit,
    rho = mean(ratio_data$G),
    formula = feature_formula,
    feature_cols = colnames(x),
    engine = engine,
    ntree = ntree,
    mtry = mtry,
    min_node_size = min_node_size
  )
}

.predict_mdr_weight_random_forest_prob <- function(ratio_fit_z, newdata, eps = 0.001) {
  x <- .mdr_weight_classifier_model_matrix(
    data = newdata,
    feature_formula = ratio_fit_z$formula,
    feature_cols = ratio_fit_z$feature_cols
  )

  if (identical(ratio_fit_z$engine, "randomForest")) {
    p <- try(
      suppressWarnings(predict(ratio_fit_z$fit, newdata = x, type = "prob")),
      silent = TRUE
    )
    if (inherits(p, "try-error")) {
      return(rep(clip_prob(ratio_fit_z$rho, eps), nrow(newdata)))
    }
    p <- p[, "1"]
  } else {
    pred <- try(
      suppressWarnings(predict(ratio_fit_z$fit, data = as.data.frame(x))),
      silent = TRUE
    )
    if (inherits(pred, "try-error") || is.null(pred$predictions)) {
      return(rep(clip_prob(ratio_fit_z$rho, eps), nrow(newdata)))
    }
    p <- pred$predictions[, "1"]
  }

  p <- as.numeric(p)
  fallback <- clip_prob(ratio_fit_z$rho, eps)
  p[!is.finite(p)] <- fallback
  clip_prob(p, eps)
}

.fit_mdr_ratio_classifier_pooled <- function(
    target_states_z,
    adherent_states_z,
    balance_classes = TRUE,
    seed = NULL,
    classifier_method = mdr_weight_classifier_methods,
    random_forest_ntree = 500,
    random_forest_mtry = NULL,
    random_forest_min_node_size = NULL,
    random_forest_engine = c("randomForest", "ranger")) {

  classifier_method <- validate_mdr_weight_classifier_method(classifier_method)
  if (!is.null(seed)) set.seed(seed)

  if (nrow(target_states_z) == 0L || nrow(adherent_states_z) == 0L) {
    stop("Target and adherent state samples must both be nonempty.", call. = FALSE)
  }

  target_states_z <- target_states_z %>% dplyr::mutate(G = 1L)
  adherent_states_z <- adherent_states_z %>% dplyr::mutate(G = 0L)

  target_train <- if (isTRUE(balance_classes)) {
    sample_target_by_follow(
      target_states_z = target_states_z,
      adherent_states_z = adherent_states_z,
      seed = seed
    )
  } else {
    target_states_z
  }
  ratio_data <- dplyr::bind_rows(target_train, adherent_states_z)

  balance_counts <- ratio_data %>%
    dplyr::count(follow_k, G, name = "n") %>%
    tidyr::pivot_wider(
      names_from = G,
      values_from = n,
      names_prefix = "n_G",
      values_fill = 0
    ) %>%
    dplyr::rename(
      n_adherent_train_k = n_G0,
      n_target_train_k = n_G1
    ) %>%
    dplyr::full_join(
      target_states_z %>% dplyr::count(follow_k, name = "n_target_available_k"),
      by = "follow_k"
    ) %>%
    dplyr::full_join(
      adherent_states_z %>% dplyr::count(follow_k, name = "n_adherent_available_k"),
      by = "follow_k"
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("n_"), ~ dplyr::coalesce(.x, 0L)),
      balance_within_follow = TRUE
    )

  ratio_formula <- G ~ factor(follow_k) + trial_m +
    X1 + X2 + X3 + L0 + L_prev_state + L_curr

  if (classifier_method == "logistic_regression") {
    fit <- .fit_logistic_safe(ratio_formula, data = ratio_data)
    if (is.null(fit)) {
      stop("Direct MDR pooled ratio classifier failed.", call. = FALSE)
    }

    fit_obj <- list(
      classifier_method = "logistic_regression",
      fit = fit,
      rho = mean(ratio_data$G),
      formula = ratio_formula
    )
  } else {
    fit_obj <- .fit_mdr_weight_random_forest_ratio_classifier(
      ratio_data = ratio_data,
      include_follow = TRUE,
      seed = seed,
      ntree = random_forest_ntree,
      mtry = random_forest_mtry,
      min_node_size = random_forest_min_node_size,
      engine = random_forest_engine
    )
  }

  c(
    fit_obj,
    list(
      balance_within_follow = TRUE,
      n_target_train = sum(ratio_data$G == 1L),
      n_adherent_train = sum(ratio_data$G == 0L),
      n_target_available = nrow(target_states_z),
      n_adherent_available = nrow(adherent_states_z),
      balance_counts = balance_counts,
      ratio_data = ratio_data
    )
  )
}

.compute_mdr_weights_from_fit <- function(
    adherent_rows_z,
    ratio_fit_z,
    eps = 0.001,
    cap = NULL) {

  .check_mdr_weight_cap(cap)

  classifier_method <- ratio_fit_z$classifier_method
  if (is.null(classifier_method)) classifier_method <- "logistic_regression"

  if (identical(classifier_method, "random_forest")) {
    h_hat <- .predict_mdr_weight_random_forest_prob(
      ratio_fit_z = ratio_fit_z,
      newdata = adherent_rows_z,
      eps = eps
    )
  } else {
    h_hat <- .predict_prob_safe(
      ratio_fit_z$fit,
      adherent_rows_z,
      fallback = ratio_fit_z$rho,
      eps = eps
    )
  }

  rho <- ratio_fit_z$rho
  w_raw <- ((1 - rho) / rho) * h_hat / (1 - h_hat)
  w <- if (is.null(cap)) w_raw else pmin(w_raw, cap)

  adherent_rows_z$h_hat <- h_hat
  adherent_rows_z$rho <- rho
  adherent_rows_z$w_mdr_raw <- w_raw
  adherent_rows_z$w_mdr <- w
  adherent_rows_z
}

.compute_mdr_weighted_rows_pooled_raw <- function(
    state_samples,
    balance_classes = TRUE,
    eps = 0.001,
    seed = NULL,
    classifier_method = mdr_weight_classifier_methods,
    random_forest_ntree = 500,
    random_forest_mtry = NULL,
    random_forest_min_node_size = NULL,
    random_forest_engine = c("randomForest", "ranger")) {

  classifier_method <- validate_mdr_weight_classifier_method(classifier_method)
  weighted_rows <- list()
  ratio_fits <- list()
  cell_support <- list()

  for (z in c(0L, 1L)) {
    target_z <- state_samples$target_states[[as.character(z)]]
    adherent_z <- state_samples$adherent_states[[as.character(z)]]

    cell_support[[as.character(z)]] <- check_target_adherent_cell_support(
      target_states_z = target_z,
      adherent_states_z = adherent_z,
      strategy_z = z
    )

    ratio_fit_z <- .fit_mdr_ratio_classifier_pooled(
      target_states_z = target_z,
      adherent_states_z = adherent_z,
      balance_classes = balance_classes,
      seed = .mdr_weight_classifier_seed(
        .mdr_variant_seed(seed, 200L + z),
        classifier_method
      ),
      classifier_method = classifier_method,
      random_forest_ntree = random_forest_ntree,
      random_forest_mtry = random_forest_mtry,
      random_forest_min_node_size = random_forest_min_node_size,
      random_forest_engine = random_forest_engine
    )

    weighted_z <- .compute_mdr_weights_from_fit(
      adherent_rows_z = adherent_z,
      ratio_fit_z = ratio_fit_z,
      eps = eps,
      cap = NULL
    )
    ratio_source <- "pooled_follow_balanced"
    if (!identical(classifier_method, "logistic_regression")) {
      ratio_source <- paste0(ratio_source, "_", classifier_method)
    }
    weighted_z$ratio_source <- ratio_source
    weighted_z$classifier_method <- classifier_method

    weighted_rows[[as.character(z)]] <- weighted_z
    ratio_fits[[as.character(z)]] <- ratio_fit_z
  }

  list(
    weighted_rows = dplyr::bind_rows(weighted_rows),
    ratio_fits = ratio_fits,
    cell_support = dplyr::bind_rows(cell_support)
  )
}

compute_mdr_weights_follow_balanced_cell_normalized <- function(
    state_samples,
    balance_classes = TRUE,
    eps = 0.001,
    weight_cap = NULL,
    seed = NULL,
    classifier_method = mdr_weight_classifier_methods,
    random_forest_ntree = 500,
    random_forest_mtry = NULL,
    random_forest_min_node_size = NULL,
    random_forest_engine = c("randomForest", "ranger")) {

  raw_obj <- .compute_mdr_weighted_rows_pooled_raw(
    state_samples = state_samples,
    balance_classes = balance_classes,
    eps = eps,
    seed = seed,
    classifier_method = classifier_method,
    random_forest_ntree = random_forest_ntree,
    random_forest_mtry = random_forest_mtry,
    random_forest_min_node_size = random_forest_min_node_size,
    random_forest_engine = random_forest_engine
  )

  list(
    weighted_rows = normalize_mdr_weights_by_cell(
      raw_obj$weighted_rows,
      cap = weight_cap
    ),
    ratio_fits = raw_obj$ratio_fits,
    cell_support = raw_obj$cell_support,
    weight_method = "follow_balanced_cell_normalized",
    classifier_method = validate_mdr_weight_classifier_method(classifier_method)
  )
}

compute_mdr_weights_by_method <- function(
    state_samples,
    method = mdr_weight_methods,
    balance_classes = TRUE,
    eps = 0.001,
    weight_cap = NULL,
    seed = NULL,
    classifier_method = mdr_weight_classifier_methods,
    random_forest_ntree = 500,
    random_forest_mtry = NULL,
    random_forest_min_node_size = NULL,
    random_forest_engine = c("randomForest", "ranger")) {

  method <- validate_mdr_weight_methods(method)[[1L]]
  classifier_method <- validate_mdr_weight_classifier_method(classifier_method)

  out <- compute_mdr_weights_follow_balanced_cell_normalized(
    state_samples = state_samples,
    balance_classes = balance_classes,
    eps = eps,
    weight_cap = weight_cap,
    seed = seed,
    classifier_method = classifier_method,
    random_forest_ntree = random_forest_ntree,
    random_forest_mtry = random_forest_mtry,
    random_forest_min_node_size = random_forest_min_node_size,
    random_forest_engine = random_forest_engine
  )
  out$weight_method <- method
  out$classifier_method <- classifier_method
  out
}
