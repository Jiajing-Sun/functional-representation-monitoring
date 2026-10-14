#!/usr/bin/env Rscript

# Paired comparison of the restored LDL componentwise adjusted-range monitor
# with Shao's quadratic self-normalization benchmark. Every pair uses
# the same simulated curves, FPCA fit, retained dimension, LDL transformation,
# normalization block, monitoring block, boundary, and nominal alpha.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath("simulations/componentwise_range")
project_dir <- normalizePath(file.path(script_dir, "../.."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "simulations", "formal", "R", "null_dgp.R"))
source(file.path(project_dir, "simulations", "formal", "R", "null_worker.R"))

dc_parse_shao_args <- function(args) {
  defaults <- list(
    nrep_primary = 5000L,
    nrep_stress = 5000L,
    nrep_power = 5000L,
    workers = max(1L, min(10L, parallel::detectCores(logical = TRUE) - 1L)),
    chunk_size = 100L,
    cv_nrep = 250000L,
    cv_train_grid = 1000L,
    stage = "all",
    panel = "all",
    scenario = "all",
    seed = 20260830L
  )
  character_keys <- c("stage", "panel", "scenario")
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[1L] %in% names(defaults)) {
      stop("Unknown argument: ", arg)
    }
    key <- pieces[1L]
    defaults[[key]] <- if (key %in% character_keys) pieces[2L] else as.integer(pieces[2L])
  }
  if (!defaults$stage %in% c("all", "null", "power")) {
    stop("stage must be one of all, null, or power.")
  }
  if (!defaults$panel %in% c("all", "primary", "stress")) {
    stop("panel must be one of all, primary, or stress.")
  }
  defaults
}

dc_comparison_scenario <- function(id, design, panel, dgp_parameters,
                                   eigen_profile = "fixed", K_mode = "fixed",
                                   K_fixed = 2L, nrep = 5000L) {
  list(
    id = id,
    panel = panel,
    assumption_alignment = if (panel == "primary") {
      "componentwise score condition"
    } else "deliberate misspecification stress",
    dgp_parameters = dgp_parameters,
    design = design,
    eigen_profile = eigen_profile,
    K_mode = K_mode,
    K_fixed = as.integer(K_fixed),
    fve = 0.90,
    d_grid = 8L,
    m_fit = 500L,
    m_calibration = 1000L,
    m_monitor = 500L,
    alpha = 0.05,
    gamma = 0,
    prewhiten_mode = "demix",
    feature_var_order = 0L,
    ldl_ridge = 1e-8,
    sn_ridge = 0,
    lrv_diagnostic_bandwidth = 10L,
    nrep = as.integer(nrep)
  )
}

dc_comparison_scenarios <- function(nrep_primary, nrep_stress) {
  list(
    dc_comparison_scenario(
      "primary_gaussian_iid", "gaussian_iid", "primary",
      "iid Gaussian functional scores", nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_hong_diag_var", "hong_var_diag", "primary",
      "VAR(1): diagonal 0.50; iid N(0,I) errors", nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_component_lrv", "component_ar_lrv", "primary",
      "independent heterogeneous AR scores; non-diagonal LRV operator on the observed function grid",
      nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_component_garch", "component_garch", "primary",
      "independent score GARCH alpha=0.10, beta=0.79 with own AR coefficient 0.35",
      nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_component_sv", "component_sv", "primary",
      "independent score SV rho=0.90, sigma_eta=0.15 with own AR coefficient 0.35",
      nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_component_t5", "component_t5", "primary",
      "independent standardized t5 score errors with heterogeneous own AR coefficients",
      nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_randomK_stable", "component_ar_lrv", "primary",
      "component AR scores; 90% FVE cutoff separated from population boundary",
      eigen_profile = "fve_stable", K_mode = "random_fve", nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "primary_randomK_boundary", "component_ar_lrv", "primary",
      "component AR scores; population FVE approximately at 90% boundary",
      eigen_profile = "fve_boundary", K_mode = "random_fve", nrep = nrep_primary
    ),
    dc_comparison_scenario(
      "stress_hong_cross_var", "hong_var_cross", "stress",
      "banded VAR(1): diagonal 0.50, adjacent 0.10; error correlation rho=0.10",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_hong_cross_garch", "hong_garch", "stress",
      "banded VAR(1); GARCH alpha=0.10, beta=0.79; shock correlation rho=0.10",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_general_nondiag_lrv", "nondiag_lrv", "stress",
      "general VAR rotation and correlated errors; non-diagonal score LRV after lag-zero LDL",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_common_persistent_garch", "garch", "stress",
      "common-volatility GARCH alpha=0.08, beta=0.87 plus general VAR rotation",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_common_persistent_sv", "sv", "stress",
      "common-factor SV rho=0.96, sigma_eta=0.18 plus general VAR rotation",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_general_t5", "heavy_t5", "stress",
      "standardized t5 with correlated innovations and general VAR rotation",
      nrep = nrep_stress
    ),
    dc_comparison_scenario(
      "stress_general_t3", "heavy_t3_stress", "stress",
      "standardized t3 stress outside the preferred moment regime",
      nrep = nrep_stress
    )
  )
}

dc_comparison_critical_cache <- function(cache_file, K_values, scenario,
                                         nrep, train_grid, seed) {
  requested <- list(
    K_values = as.integer(K_values),
    horizon = scenario$m_monitor / scenario$m_calibration,
    gamma = scenario$gamma,
    alpha = scenario$alpha,
    nrep = as.integer(nrep),
    train_grid = as.integer(train_grid)
  )
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$settings, requested)) return(cached$values)
  }
  audit_file <- file.path(
    project_dir, "simulations", "componentwise_range", "outputs",
    "critical_value_audit", "critical_value_grid_and_monte_carlo_audit.csv"
  )
  if (file.exists(audit_file) && requested$horizon == 0.5 &&
      requested$gamma == 0 && requested$alpha == 0.05 &&
      requested$nrep == 250000L && requested$train_grid == 1000L) {
    audit <- utils::read.csv(audit_file, stringsAsFactors = FALSE)
    audit <- audit[
      audit$configuration == "exact1000_500_T0.5" & audit$K %in% K_values,
      , drop = FALSE
    ]
    range_rows <- audit[audit$method_id == "adjusted_range", , drop = FALSE]
    full_rows <- audit[audit$method_id == "full_matrix", , drop = FALSE]
    range_rows <- range_rows[match(K_values, range_rows$K), , drop = FALSE]
    full_rows <- full_rows[match(K_values, full_rows$K), , drop = FALSE]
    if (nrow(range_rows) == length(K_values) && nrow(full_rows) == length(K_values) &&
        !anyNA(range_rows$critical_value) && !anyNA(full_rows$critical_value)) {
      adjusted_range <- as.list(range_rows$critical_value)
      shao <- as.list(full_rows$critical_value)
      names(adjusted_range) <- names(shao) <- paste0("K", K_values)
      values <- list(adjusted_range = adjusted_range, shao = shao)
      saveRDS(list(settings = requested, values = values, source = audit_file), cache_file)
      return(values)
    }
  }
  adjusted_range <- lapply(K_values, function(K) {
    dc_simulate_adjusted_range_critical_value(
      K,
      T = requested$horizon,
      gamma = requested$gamma,
      alpha = requested$alpha,
      nrep = requested$nrep,
      train_grid = requested$train_grid,
      seed = seed + 100L * K
    )$critical_value
  })
  shao <- lapply(K_values, function(K) {
    dc_simulate_sn_critical_value(
      K,
      T = requested$horizon,
      gamma = requested$gamma,
      alpha = requested$alpha,
      nrep = requested$nrep,
      train_grid = requested$train_grid,
      seed = seed + 100L * K
    )$critical_value
  })
  names(adjusted_range) <- names(shao) <- paste0("K", K_values)
  values <- list(adjusted_range = adjusted_range, shao = shao)
  saveRDS(list(settings = requested, values = values), cache_file)
  message(paste(
    sprintf(
      "K=%d: range %.4f; Shao %.4f",
      K_values,
      unlist(adjusted_range, use.names = FALSE),
      unlist(shao, use.names = FALSE)
    ),
    collapse = "; "
  ))
  values
}

dc_fit_and_monitor_pair <- function(x, scenario, critical_values) {
  fit_index <- seq_len(scenario$m_fit)
  cal_index <- scenario$m_fit + seq_len(scenario$m_calibration)
  monitor_index <- scenario$m_fit + scenario$m_calibration + seq_len(scenario$m_monitor)
  fit <- dc_fit_score_comparison(
    x_fit = x[fit_index, , drop = FALSE],
    x_calibration = x[cal_index, , drop = FALSE],
    K = if (scenario$K_mode == "fixed") scenario$K_fixed else NULL,
    fve = scenario$fve,
    K_max = 5L,
    prewhiten_mode = scenario$prewhiten_mode,
    ldl_ridge = scenario$ldl_ridge,
    sn_ridge = scenario$sn_ridge,
    feature_var_order = scenario$feature_var_order
  )
  K_hat <- fit$settings$K
  key <- paste0("K", K_hat)
  critical <- c(
    adjusted_range = critical_values$adjusted_range[[key]],
    shao = critical_values$shao[[key]]
  )
  if (any(!is.finite(critical))) stop("No critical values cached for K = ", K_hat, ".")
  monitor <- dc_monitor_score_comparison(
    fit,
    x[monitor_index, , drop = FALSE],
    alpha = scenario$alpha,
    gamma = scenario$gamma,
    critical_values = critical
  )
  list(fit = fit, monitor = monitor)
}

dc_paired_null_replicate <- function(replicate_id, scenario, critical_values, seed_base) {
  set.seed(seed_base + as.integer(replicate_id))
  total_n <- scenario$m_fit + scenario$m_calibration + scenario$m_monitor
  x <- dc_simulate_null_scores(
    n = total_n,
    d = scenario$d_grid,
    design = scenario$design,
    eigen_profile = scenario$eigen_profile
  )
  tryCatch({
    result <- dc_fit_and_monitor_pair(x, scenario, critical_values)
    fit <- result$fit
    monitor <- result$monitor
    demixed <- dc_adjusted_range_demixed_features(
      fit$adjusted_range$sn,
      fit$adjusted_range$calibration_sn_scores
    )
    data.frame(
      replicate = replicate_id,
      K_hat = fit$settings$K,
      adjusted_range_reject = !is.na(monitor$alarm[["adjusted_range"]]),
      shao_reject = !is.na(monitor$alarm[["shao"]]),
      adjusted_range_first_alarm = monitor$alarm[["adjusted_range"]],
      shao_first_alarm = monitor$alarm[["shao"]],
      adjusted_range_max = max(monitor$path$adjusted_range),
      shao_max = max(monitor$path$shao),
      adjusted_range_critical = monitor$critical_values[["adjusted_range"]],
      shao_critical = monitor$critical_values[["shao"]],
      adjusted_range_condition = fit$adjusted_range$sn$condition_number,
      shao_condition = fit$shao$sn$condition_number,
      offdiag_lrv = dc_bartlett_lrv_offdiag_ratio(
        demixed, bandwidth = scenario$lrv_diagnostic_bandwidth
      ),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = replicate_id,
      K_hat = NA_integer_,
      adjusted_range_reject = NA,
      shao_reject = NA,
      adjusted_range_first_alarm = NA_integer_,
      shao_first_alarm = NA_integer_,
      adjusted_range_max = NA_real_,
      shao_max = NA_real_,
      adjusted_range_critical = NA_real_,
      shao_critical = NA_real_,
      adjusted_range_condition = NA_real_,
      shao_condition = NA_real_,
      offdiag_lrv = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

dc_run_chunks <- function(replicate_function, scenario, critical_values,
                          checkpoint_dir, workers, chunk_size, seed_base, ...) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  starts <- seq.int(1L, scenario$nrep, by = chunk_size)
  chunks <- vector("list", length(starts))
  for (chunk_index in seq_along(starts)) {
    first <- starts[chunk_index]
    last <- min(scenario$nrep, first + chunk_size - 1L)
    checkpoint <- file.path(checkpoint_dir, sprintf("rep_%05d_%05d.rds", first, last))
    if (file.exists(checkpoint)) {
      chunks[[chunk_index]] <- readRDS(checkpoint)
      next
    }
    pieces <- parallel::mclapply(
      seq.int(first, last),
      replicate_function,
      scenario = scenario,
      critical_values = critical_values,
      seed_base = seed_base,
      ...,
      mc.cores = workers,
      mc.preschedule = TRUE
    )
    chunks[[chunk_index]] <- do.call(rbind, pieces)
    saveRDS(chunks[[chunk_index]], checkpoint)
    message(sprintf("%s: completed %d/%d", scenario$id, last, scenario$nrep))
  }
  do.call(rbind, chunks)
}

dc_method_null_summary <- function(raw, scenario, method) {
  valid <- raw[is.na(raw$error), , drop = FALSE]
  reject <- valid[[paste0(method, "_reject")]]
  size <- mean(reject)
  mcse <- sqrt(size * (1 - size) / nrow(valid))
  k_table <- table(valid$K_hat)
  data.frame(
    panel = scenario$panel,
    scenario_id = scenario$id,
    design = scenario$design,
    dgp_parameters = scenario$dgp_parameters,
    method = if (method == "adjusted_range") {
      "LDL componentwise adjusted-range"
    } else "Shao (2010) quadratic self-normalization",
    method_id = method,
    n_requested = scenario$nrep,
    n_valid = nrow(valid),
    n_errors = sum(!is.na(raw$error)),
    empirical_size = size,
    monte_carlo_se = mcse,
    ci95_low = max(0, size - 1.96 * mcse),
    ci95_high = min(1, size + 1.96 * mcse),
    mean_K = mean(valid$K_hat),
    K_distribution = paste0(
      names(k_table), ":", sprintf("%.4f", as.numeric(k_table) / nrow(valid)),
      collapse = ";"
    ),
    median_condition_number = stats::median(valid[[paste0(method, "_condition")]]),
    median_offdiag_lrv = stats::median(valid$offdiag_lrv),
    critical_value = stats::median(valid[[paste0(method, "_critical")]]),
    stringsAsFactors = FALSE
  )
}

dc_paired_null_summary <- function(raw, scenario) {
  valid <- raw[is.na(raw$error), , drop = FALSE]
  range_reject <- as.integer(valid$adjusted_range_reject)
  shao_reject <- as.integer(valid$shao_reject)
  difference <- range_reject - shao_reject
  se <- stats::sd(difference) / sqrt(length(difference))
  data.frame(
    panel = scenario$panel,
    scenario_id = scenario$id,
    n_valid = nrow(valid),
    adjusted_range_size = mean(range_reject),
    shao_size = mean(shao_reject),
    paired_difference = mean(difference),
    paired_difference_se = se,
    paired_ci95_low = mean(difference) - 1.96 * se,
    paired_ci95_high = mean(difference) + 1.96 * se,
    both_reject = mean(range_reject == 1L & shao_reject == 1L),
    range_only = mean(range_reject == 1L & shao_reject == 0L),
    shao_only = mean(range_reject == 0L & shao_reject == 1L),
    neither = mean(range_reject == 0L & shao_reject == 0L),
    stringsAsFactors = FALSE
  )
}

dc_power_specs <- function(nrep) {
  abrupt <- expand.grid(
    design = c(
      "gaussian_iid", "hong_var_diag", "component_garch",
      "component_sv", "component_t5"
    ),
    direction = c("sparse_first", "dense_first_two"),
    delta = c(0.05, 0.15, 0.30),
    change_shape = "abrupt",
    stringsAsFactors = FALSE
  )
  ramp <- expand.grid(
    design = c("hong_var_diag", "component_garch", "component_t5"),
    direction = "sparse_first",
    delta = c(0.20, 0.40, 0.60),
    change_shape = "ramp",
    stringsAsFactors = FALSE
  )
  base <- rbind(abrupt, ramp)
  null_ids <- c(
    gaussian_iid = "primary_gaussian_iid",
    hong_var_diag = "primary_hong_diag_var",
    component_garch = "primary_component_garch",
    component_sv = "primary_component_sv",
    component_t5 = "primary_component_t5"
  )
  lapply(seq_len(nrow(base)), function(i) {
    row <- base[i, ]
    scenario <- dc_comparison_scenario(
      id = sprintf(
        "power_%s_%s_%s_delta_%0.2f",
        row$change_shape, row$design, row$direction, row$delta
      ),
      design = row$design,
      panel = "power",
      dgp_parameters = "standardized retained-score mean shift",
      K_mode = "fixed",
      K_fixed = 2L,
      nrep = nrep
    )
    scenario$direction <- row$direction
    scenario$delta <- as.numeric(row$delta)
    scenario$change_shape <- row$change_shape
    scenario$break_after <- 100L
    scenario$null_scenario_id <- unname(null_ids[[row$design]])
    scenario
  })
}

dc_shift_vector <- function(design, d, eigen_profile, K, direction, delta) {
  lambda <- dc_score_eigenvalues(d, profile = eigen_profile)
  rotation <- if (design == "gaussian_iid") diag(d) else dc_functional_basis_rotation(d)
  standardized_direction <- if (direction == "sparse_first") {
    c(1, rep(0, K - 1L))
  } else if (direction == "dense_first_two") {
    rep(1 / sqrt(K), K)
  } else stop("Unknown shift direction: ", direction)
  score_shift <- delta * standardized_direction * sqrt(lambda[seq_len(K)])
  as.numeric(score_shift %*% rotation[seq_len(K), , drop = FALSE])
}

dc_first_crossing <- function(path, critical) {
  crossed <- which(path > critical)
  if (length(crossed)) crossed[1L] else NA_integer_
}

dc_paired_power_replicate <- function(replicate_id, scenario, critical_values,
                                      seed_base, empirical_critical) {
  set.seed(seed_base + as.integer(replicate_id))
  total_n <- scenario$m_fit + scenario$m_calibration + scenario$m_monitor
  x <- dc_simulate_null_scores(
    n = total_n,
    d = scenario$d_grid,
    design = scenario$design,
    eigen_profile = scenario$eigen_profile
  )
  change_start <- scenario$m_fit + scenario$m_calibration + scenario$break_after + 1L
  shift <- dc_shift_vector(
    scenario$design,
    scenario$d_grid,
    scenario$eigen_profile,
    scenario$K_fixed,
    scenario$direction,
    scenario$delta
  )
  change_rows <- change_start:total_n
  change_weight <- if (scenario$change_shape == "abrupt") {
    rep(1, length(change_rows))
  } else if (scenario$change_shape == "ramp") {
    seq_len(length(change_rows)) / length(change_rows)
  } else stop("Unknown change shape: ", scenario$change_shape)
  x[change_rows, ] <- x[change_rows, , drop = FALSE] +
    outer(change_weight, shift)
  tryCatch({
    result <- dc_fit_and_monitor_pair(x, scenario, critical_values)
    path <- result$monitor$path
    data.frame(
      replicate = replicate_id,
      adjusted_range_alarm_nominal = dc_first_crossing(
        path$adjusted_range, result$monitor$critical_values[["adjusted_range"]]
      ),
      shao_alarm_nominal = dc_first_crossing(
        path$shao, result$monitor$critical_values[["shao"]]
      ),
      adjusted_range_alarm_size_adjusted = dc_first_crossing(
        path$adjusted_range, empirical_critical[["adjusted_range"]]
      ),
      shao_alarm_size_adjusted = dc_first_crossing(
        path$shao, empirical_critical[["shao"]]
      ),
      adjusted_range_max = max(path$adjusted_range),
      shao_max = max(path$shao),
      K_hat = result$fit$settings$K,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = replicate_id,
      adjusted_range_alarm_nominal = NA_integer_,
      shao_alarm_nominal = NA_integer_,
      adjusted_range_alarm_size_adjusted = NA_integer_,
      shao_alarm_size_adjusted = NA_integer_,
      adjusted_range_max = NA_real_,
      shao_max = NA_real_,
      K_hat = NA_integer_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

dc_power_summary <- function(raw, scenario, method, calibration) {
  valid <- raw[is.na(raw$error), , drop = FALSE]
  alarm <- valid[[paste0(method, "_alarm_", calibration)]]
  prechange <- !is.na(alarm) & alarm <= scenario$break_after
  postchange <- !is.na(alarm) & alarm > scenario$break_after
  at_risk <- !prechange
  delay <- alarm[postchange] - scenario$break_after
  post_probability <- mean(postchange)
  pre_probability <- mean(prechange)
  conditional_probability <- sum(postchange) / sum(at_risk)
  data.frame(
    scenario_id = scenario$id,
    design = scenario$design,
    change_shape = scenario$change_shape,
    direction = scenario$direction,
    delta_standardized = scenario$delta,
    break_after = scenario$break_after,
    method = if (method == "adjusted_range") {
      "LDL componentwise adjusted-range"
    } else "Shao (2010) quadratic self-normalization",
    method_id = method,
    calibration = calibration,
    n_requested = scenario$nrep,
    n_valid = nrow(valid),
    n_errors = sum(!is.na(raw$error)),
    rejection_by_horizon = mean(!is.na(alarm)),
    prechange_false_alarm = pre_probability,
    prechange_false_alarm_mcse = sqrt(
      pre_probability * (1 - pre_probability) / nrow(valid)
    ),
    postchange_detection_probability = post_probability,
    postchange_detection_mcse = sqrt(
      post_probability * (1 - post_probability) / nrow(valid)
    ),
    conditional_postchange_power = conditional_probability,
    conditional_postchange_power_mcse = sqrt(
      conditional_probability * (1 - conditional_probability) / sum(at_risk)
    ),
    restricted_mean_delay = mean(ifelse(
      postchange[at_risk],
      alarm[at_risk] - scenario$break_after,
      scenario$m_monitor - scenario$break_after + 1L
    )),
    mean_delay_conditional = if (length(delay)) mean(delay) else NA_real_,
    median_delay_conditional = if (length(delay)) stats::median(delay) else NA_real_,
    stringsAsFactors = FALSE
  )
}

dc_paired_power_summary <- function(raw, scenario, calibration) {
  valid <- raw[is.na(raw$error), , drop = FALSE]
  range_alarm <- valid[[paste0("adjusted_range_alarm_", calibration)]]
  shao_alarm <- valid[[paste0("shao_alarm_", calibration)]]
  range_pre <- !is.na(range_alarm) & range_alarm <= scenario$break_after
  shao_pre <- !is.na(shao_alarm) & shao_alarm <= scenario$break_after
  range_post <- !is.na(range_alarm) & range_alarm > scenario$break_after
  shao_post <- !is.na(shao_alarm) & shao_alarm > scenario$break_after
  post_difference <- as.integer(range_post) - as.integer(shao_post)
  post_difference_se <- stats::sd(post_difference) / sqrt(nrow(valid))
  common_risk <- !range_pre & !shao_pre
  cap <- scenario$m_monitor - scenario$break_after + 1L
  range_restricted <- ifelse(
    range_post[common_risk], range_alarm[common_risk] - scenario$break_after, cap
  )
  shao_restricted <- ifelse(
    shao_post[common_risk], shao_alarm[common_risk] - scenario$break_after, cap
  )
  delay_difference <- range_restricted - shao_restricted
  delay_difference_se <- stats::sd(delay_difference) / sqrt(length(delay_difference))
  data.frame(
    scenario_id = scenario$id,
    design = scenario$design,
    change_shape = scenario$change_shape,
    direction = scenario$direction,
    delta_standardized = scenario$delta,
    calibration = calibration,
    n_valid = nrow(valid),
    range_postchange_detection = mean(range_post),
    shao_postchange_detection = mean(shao_post),
    paired_power_difference = mean(post_difference),
    paired_power_difference_se = post_difference_se,
    paired_power_ci95_low = mean(post_difference) - 1.96 * post_difference_se,
    paired_power_ci95_high = mean(post_difference) + 1.96 * post_difference_se,
    common_risk_n = sum(common_risk),
    range_restricted_mean_delay = mean(range_restricted),
    shao_restricted_mean_delay = mean(shao_restricted),
    paired_restricted_delay_difference = mean(delay_difference),
    paired_restricted_delay_difference_se = delay_difference_se,
    paired_restricted_delay_ci95_low = mean(delay_difference) - 1.96 * delay_difference_se,
    paired_restricted_delay_ci95_high = mean(delay_difference) + 1.96 * delay_difference_se,
    stringsAsFactors = FALSE
  )
}

args <- dc_parse_shao_args(commandArgs(trailingOnly = TRUE))
output_dir <- file.path(script_dir, "outputs", "shao_comparison")
checkpoint_root <- file.path(
  script_dir,
  "checkpoints",
  sprintf(
    "shao_comparison_formula_v2_m500_1000_500_cv%d_grid%d_seed%d",
    args$cv_nrep, args$cv_train_grid, args$seed
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)

all_null_scenarios <- dc_comparison_scenarios(args$nrep_primary, args$nrep_stress)
reference_scenario <- all_null_scenarios[[1L]]
critical_file <- file.path(
  output_dir,
  sprintf("paired_critical_T0.5_cv%d_grid%d.rds", args$cv_nrep, args$cv_train_grid)
)
critical_values <- dc_comparison_critical_cache(
  critical_file,
  K_values = 1:5,
  scenario = reference_scenario,
  nrep = args$cv_nrep,
  train_grid = args$cv_train_grid,
  seed = args$seed + 900000L
)
critical_table <- rbind(
  data.frame(
    method_id = "adjusted_range",
    method_display = "LDL componentwise adjusted-range",
    K = 1:5,
    T = reference_scenario$m_monitor / reference_scenario$m_calibration,
    gamma = reference_scenario$gamma,
    alpha = reference_scenario$alpha,
    cv_nrep = args$cv_nrep,
    train_grid = args$cv_train_grid,
    seed = args$seed + 900000L + 100L * (1:5),
    critical_value = unlist(critical_values$adjusted_range, use.names = FALSE)
  ),
  data.frame(
    method_id = "shao",
    method_display = "Shao (2010) quadratic self-normalization",
    K = 1:5,
    T = reference_scenario$m_monitor / reference_scenario$m_calibration,
    gamma = reference_scenario$gamma,
    alpha = reference_scenario$alpha,
    cv_nrep = args$cv_nrep,
    train_grid = args$cv_train_grid,
    seed = args$seed + 900000L + 100L * (1:5),
    critical_value = unlist(critical_values$shao, use.names = FALSE)
  )
)
utils::write.csv(
  critical_table,
  file.path(
    output_dir,
    sprintf("paired_critical_values_cv%d_grid%d.csv", args$cv_nrep, args$cv_train_grid)
  ),
  row.names = FALSE
)

null_scenarios <- all_null_scenarios
if (args$panel != "all") {
  null_scenarios <- Filter(function(x) identical(x$panel, args$panel), null_scenarios)
}
if (args$scenario != "all" && args$stage %in% c("all", "null")) {
  null_scenarios <- Filter(function(x) identical(x$id, args$scenario), null_scenarios)
}

null_raw_map <- list()
null_method_summaries <- list()
null_paired_summaries <- list()
if (args$stage %in% c("all", "null")) {
  scenario_ids <- vapply(all_null_scenarios, `[[`, character(1L), "id")
  for (i in seq_along(null_scenarios)) {
    scenario <- null_scenarios[[i]]
    message("Starting paired null comparison: ", scenario$id)
    raw <- dc_run_chunks(
      dc_paired_null_replicate,
      scenario,
      critical_values,
      checkpoint_dir = file.path(checkpoint_root, "null", scenario$id),
      workers = args$workers,
      chunk_size = args$chunk_size,
      seed_base = args$seed + match(scenario$id, scenario_ids) * 1000000L
    )
    null_raw_map[[scenario$id]] <- raw
    saveRDS(raw, file.path(output_dir, sprintf("%s_paired_null_raw_%drep.rds", scenario$id, scenario$nrep)))
    null_method_summaries[[length(null_method_summaries) + 1L]] <-
      dc_method_null_summary(raw, scenario, "adjusted_range")
    null_method_summaries[[length(null_method_summaries) + 1L]] <-
      dc_method_null_summary(raw, scenario, "shao")
    null_paired_summaries[[length(null_paired_summaries) + 1L]] <-
      dc_paired_null_summary(raw, scenario)
  }
  if (length(null_method_summaries)) {
    null_long <- do.call(rbind, null_method_summaries)
    null_paired <- do.call(rbind, null_paired_summaries)
    utils::write.csv(
      null_long,
      file.path(
        output_dir,
        sprintf("paired_null_size_primary%d_stress%d.csv", args$nrep_primary, args$nrep_stress)
      ),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      null_paired,
      file.path(
        output_dir,
        sprintf("paired_null_differences_primary%d_stress%d.csv", args$nrep_primary, args$nrep_stress)
      ),
      row.names = FALSE,
      na = ""
    )
    print(null_long)
    print(null_paired)
  }
}

power_summaries <- list()
paired_power_summaries <- list()
if (args$stage %in% c("all", "power")) {
  all_power_specs <- dc_power_specs(args$nrep_power)
  power_specs <- all_power_specs
  if (args$scenario != "all") {
    power_specs <- Filter(function(x) identical(x$id, args$scenario), power_specs)
  }
  if (!length(power_specs)) stop("No power scenarios selected.")
  all_null_lookup <- stats::setNames(all_null_scenarios, vapply(all_null_scenarios, `[[`, character(1L), "id"))
  power_design_order <- unique(vapply(all_power_specs, `[[`, character(1L), "design"))
  for (i in seq_along(power_specs)) {
    scenario <- power_specs[[i]]
    null_id <- scenario$null_scenario_id
    null_raw <- null_raw_map[[null_id]]
    if (is.null(null_raw)) {
      null_scenario <- all_null_lookup[[null_id]]
      null_file <- file.path(
        output_dir,
        sprintf("%s_paired_null_raw_%drep.rds", null_id, null_scenario$nrep)
      )
      if (!file.exists(null_file)) {
        stop("Power stage requires paired null output: ", null_file)
      }
      null_raw <- readRDS(null_file)
    }
    null_valid <- null_raw[is.na(null_raw$error) & null_raw$K_hat == 2L, , drop = FALSE]
    empirical_critical <- c(
      adjusted_range = as.numeric(stats::quantile(
        null_valid$adjusted_range_max, 0.95, type = 8, names = FALSE
      )),
      shao = as.numeric(stats::quantile(
        null_valid$shao_max, 0.95, type = 8, names = FALSE
      ))
    )
    message("Starting paired power comparison: ", scenario$id)
    raw <- dc_run_chunks(
      dc_paired_power_replicate,
      scenario,
      critical_values,
      checkpoint_dir = file.path(checkpoint_root, "power", scenario$id),
      workers = args$workers,
      chunk_size = args$chunk_size,
      seed_base = args$seed + 100000000L +
        match(scenario$design, power_design_order) * 1000000L,
      empirical_critical = empirical_critical
    )
    saveRDS(raw, file.path(output_dir, sprintf("%s_paired_power_raw_%drep.rds", scenario$id, scenario$nrep)))
    for (method in c("adjusted_range", "shao")) {
      for (calibration in c("nominal", "size_adjusted")) {
        power_summaries[[length(power_summaries) + 1L]] <-
          dc_power_summary(raw, scenario, method, calibration)
      }
    }
    for (calibration in c("nominal", "size_adjusted")) {
      paired_power_summaries[[length(paired_power_summaries) + 1L]] <-
        dc_paired_power_summary(raw, scenario, calibration)
    }
  }
  power_table <- do.call(rbind, power_summaries)
  utils::write.csv(
    power_table,
    file.path(output_dir, sprintf("paired_power_summary_%drep.csv", args$nrep_power)),
    row.names = FALSE,
    na = ""
  )
  paired_power_table <- do.call(rbind, paired_power_summaries)
  utils::write.csv(
    paired_power_table,
    file.path(output_dir, sprintf("paired_power_differences_%drep.csv", args$nrep_power)),
    row.names = FALSE,
    na = ""
  )
  print(power_table)
  print(paired_power_table)
}

manifest <- list(
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  arguments = args,
  comparison = list(
    target = "retained FPCA score mean",
    methods = c(
      "LDL componentwise adjusted-range",
      "Shao (2010) quadratic self-normalization, implemented as in Chan et al. (2021)"
    ),
    common_design = list(
      m_fit = 500L,
      m_calibration = 1000L,
      m_monitor = 500L,
      alpha = 0.05,
      gamma = 0,
      prewhiten_mode = "demix",
      feature_var_order = 0L
    ),
    power_design = list(
      break_after = 100L,
      abrupt_delta_standardized = c(0.05, 0.15, 0.30),
      ramp_terminal_delta_standardized = c(0.20, 0.40, 0.60),
      directions = c("sparse_first", "dense_first_two"),
      change_shapes = c("abrupt", "ramp")
    )
  ),
  critical_value_file = normalizePath(critical_file),
  session_info = utils::capture.output(sessionInfo())
)
saveRDS(
  manifest,
  file.path(
    output_dir,
    sprintf(
      "shao_comparison_manifest_primary%d_stress%d_power%d.rds",
      args$nrep_primary, args$nrep_stress, args$nrep_power
    )
  )
)
