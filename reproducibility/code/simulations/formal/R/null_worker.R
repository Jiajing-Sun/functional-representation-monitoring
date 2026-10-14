dc_bartlett_lrv_offdiag_ratio <- function(features, bandwidth = 10L) {
  features <- as.matrix(features)
  n <- nrow(features)
  d <- ncol(features)
  bandwidth <- min(as.integer(bandwidth), n - 1L)
  if (d == 1L || bandwidth < 1L) return(0)
  centered <- sweep(features, 2L, colMeans(features), FUN = "-")
  lrv <- crossprod(centered) / n
  for (lag in seq_len(bandwidth)) {
    weight <- 1 - lag / (bandwidth + 1)
    gamma_lag <- crossprod(
      centered[seq_len(n - lag), , drop = FALSE],
      centered[lag + seq_len(n - lag), , drop = FALSE]
    ) / n
    lrv <- lrv + weight * (gamma_lag + t(gamma_lag))
  }
  offdiag <- lrv - diag(diag(lrv), nrow = d)
  diagonal <- diag(diag(lrv), nrow = d)
  norm(offdiag, type = "F") / max(norm(diagonal, type = "F"), .Machine$double.eps)
}

dc_adjusted_range_demixed_features <- function(sn_fit, features) {
  centered <- sweep(features, 2L, sn_fit$feature_mean, FUN = "-")
  if (!is.null(sn_fit$ldl_C)) centered %*% t(solve(sn_fit$ldl_C)) else centered
}

dc_null_replicate <- function(replicate_id, scenario, critical_values, seed_base) {
  set.seed(seed_base + as.integer(replicate_id))
  total_n <- scenario$m_fit + scenario$m_calibration + scenario$m_monitor
  x <- dc_simulate_null_scores(
    n = total_n,
    d = scenario$d_grid,
    design = scenario$design,
    eigen_profile = scenario$eigen_profile
  )
  fit_index <- seq_len(scenario$m_fit)
  cal_index <- scenario$m_fit + seq_len(scenario$m_calibration)
  monitor_index <- scenario$m_fit + scenario$m_calibration + seq_len(scenario$m_monitor)

  result <- tryCatch({
    fit <- dc_fit_dual_channel(
      x_fit = x[fit_index, , drop = FALSE],
      x_calibration = x[cal_index, , drop = FALSE],
      K = if (scenario$K_mode == "fixed") scenario$K_fixed else NULL,
      L = scenario$L,
      fve = scenario$fve,
      prewhiten_mode = scenario$prewhiten_mode,
      ldl_ridge = scenario$ldl_ridge,
      sn_ridge = scenario$sn_ridge,
      sn_type = scenario$sn_type,
      feature_var_order = scenario$feature_var_order,
      r_huber_threshold = scenario$r_huber_threshold,
      ae_hull_grid = scenario$ae_hull_grid
    )
    K_hat <- fit$fpca$K
    diagnostic_bandwidth <- if (is.null(scenario$lrv_diagnostic_bandwidth)) {
      10L
    } else as.integer(scenario$lrv_diagnostic_bandwidth)
    diagnostic_S <- dc_adjusted_range_demixed_features(
      fit$sn_S, fit$calibration_sn_features$channel_S
    )
    diagnostic_R <- dc_adjusted_range_demixed_features(
      fit$sn_R, fit$calibration_sn_features$channel_R
    )
    offdiag_lrv_S <- dc_bartlett_lrv_offdiag_ratio(
      diagnostic_S, bandwidth = diagnostic_bandwidth
    )
    offdiag_lrv_R <- dc_bartlett_lrv_offdiag_ratio(
      diagnostic_R, bandwidth = diagnostic_bandwidth
    )
    if (identical(scenario$critical_method, "stationary_bootstrap")) {
      monitored <- dc_monitor(
        fit,
        x[monitor_index, , drop = FALSE],
        alpha = scenario$alpha,
        alpha_weights = scenario$alpha_weights,
        gamma = scenario$gamma,
        critical_nrep = scenario$bootstrap_nrep,
        seed = seed_base + 500000000L + as.integer(replicate_id),
        critical_method = "stationary_bootstrap",
        bootstrap_block_length = scenario$bootstrap_block_length
      )
    } else {
      cv_key <- paste0("K", K_hat)
      cv <- critical_values[[cv_key]]
      if (is.null(cv)) stop("No critical value cached for selected K = ", K_hat, ".")
      monitored <- dc_monitor(
        fit,
        x[monitor_index, , drop = FALSE],
        alpha = scenario$alpha,
        alpha_weights = scenario$alpha_weights,
        gamma = scenario$gamma,
        critical_values = cv
      )
    }
    selected_block_length <- if (identical(scenario$critical_method, "stationary_bootstrap")) {
      monitored$critical_metadata$block_length
    } else NA_integer_
    data.frame(
      replicate = replicate_id,
      reject_global = !is.na(monitored$alarm$stop),
      reject_S = !is.na(monitored$alarm$first_S),
      reject_R = !is.na(monitored$alarm$first_R),
      first_alarm = ifelse(is.na(monitored$alarm$stop), NA_integer_, monitored$alarm$stop),
      label = monitored$alarm$label,
      K_hat = K_hat,
      d_S = fit$sn_S$dimension,
      d_R = fit$sn_R$dimension,
      cond_S = fit$sn_S$condition_number,
      cond_R = fit$sn_R$condition_number,
      ridge_S = fit$sn_S$ridge_added,
      ridge_R = fit$sn_R$ridge_added,
      offdiag_lrv_S = offdiag_lrv_S,
      offdiag_lrv_R = offdiag_lrv_R,
      critical_S = unname(monitored$critical_values["S"]),
      critical_R = unname(monitored$critical_values["R"]),
      block_length = selected_block_length,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = replicate_id,
      reject_global = NA,
      reject_S = NA,
      reject_R = NA,
      first_alarm = NA_integer_,
      label = "error",
      K_hat = NA_integer_,
      d_S = NA_integer_,
      d_R = NA_integer_,
      cond_S = NA_real_,
      cond_R = NA_real_,
      ridge_S = NA_real_,
      ridge_R = NA_real_,
      offdiag_lrv_S = NA_real_,
      offdiag_lrv_R = NA_real_,
      critical_S = NA_real_,
      critical_R = NA_real_,
      block_length = NA_integer_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  result
}

dc_summarize_null_results <- function(raw, scenario) {
  valid <- raw[is.na(raw$error), , drop = FALSE]
  n_valid <- nrow(valid)
  size <- if (n_valid) mean(valid$reject_global) else NA_real_
  mcse <- if (n_valid) sqrt(size * (1 - size) / n_valid) else NA_real_
  target_mcse <- sqrt(scenario$alpha * (1 - scenario$alpha) / max(n_valid, 1L))
  k_table <- if (n_valid) table(valid$K_hat) else integer()
  k_distribution <- if (length(k_table)) {
    paste0(names(k_table), ":", sprintf("%.4f", as.numeric(k_table) / n_valid), collapse = ";")
  } else NA_character_

  data.frame(
    scenario_id = scenario$id,
    design = scenario$design,
    eigen_profile = scenario$eigen_profile,
    K_mode = scenario$K_mode,
    K_fixed = if (scenario$K_mode == "fixed") scenario$K_fixed else NA_integer_,
    fve = scenario$fve,
    m_fit = scenario$m_fit,
    m_calibration = scenario$m_calibration,
    m_monitor = scenario$m_monitor,
    n_requested = scenario$nrep,
    n_valid = n_valid,
    n_errors = sum(!is.na(raw$error)),
    empirical_global_size = size,
    monte_carlo_se = mcse,
    nominal_mcse = target_mcse,
    ci95_low = max(0, size - 1.96 * mcse),
    ci95_high = min(1, size + 1.96 * mcse),
    channel_S_crossing = if (n_valid) mean(valid$reject_S) else NA_real_,
    channel_R_crossing = if (n_valid) mean(valid$reject_R) else NA_real_,
    label_S = if (n_valid) mean(valid$label == "S") else NA_real_,
    label_R = if (n_valid) mean(valid$label == "R") else NA_real_,
    label_both = if (n_valid) mean(valid$label == "both") else NA_real_,
    mean_K = if (n_valid) mean(valid$K_hat) else NA_real_,
    sd_K = if (n_valid) stats::sd(valid$K_hat) else NA_real_,
    modal_K = if (length(k_table)) as.integer(names(k_table)[which.max(k_table)]) else NA_integer_,
    modal_K_probability = if (length(k_table)) max(k_table) / n_valid else NA_real_,
    K_distribution = k_distribution,
    median_cond_S = if (n_valid) stats::median(valid$cond_S) else NA_real_,
    median_cond_R = if (n_valid) stats::median(valid$cond_R) else NA_real_,
    median_offdiag_lrv_S = if (n_valid && "offdiag_lrv_S" %in% names(valid)) {
      stats::median(valid$offdiag_lrv_S)
    } else NA_real_,
    median_offdiag_lrv_R = if (n_valid && "offdiag_lrv_R" %in% names(valid)) {
      stats::median(valid$offdiag_lrv_R)
    } else NA_real_,
    median_critical_S = if (n_valid) stats::median(valid$critical_S) else NA_real_,
    median_critical_R = if (n_valid) stats::median(valid$critical_R) else NA_real_,
    median_block_length = if (n_valid && any(is.finite(valid$block_length))) {
      stats::median(valid$block_length, na.rm = TRUE)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}
