#!/usr/bin/env Rscript

# Controlled post-LDL long-run-correlation experiment.  The bivariate process
# x1_t=e1_t, x2_t=rho*e1_{t-1}+sqrt(1-rho^2)*e2_t has identity lag-zero
# covariance and LRV correlation rho.  Hence population LDL is the identity and
# rho directly controls the misspecification relevant to componentwise ranges.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else normalizePath("simulations/componentwise_range")
project_dir <- normalizePath(file.path(script_dir, "../.."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "simulations", "formal", "R", "null_worker.R"))

parse_args <- function(args) {
  out <- list(
    nrep = 5000L,
    workers = max(1L, min(10L, parallel::detectCores(logical = TRUE) - 1L)),
    seed = 20260831L
  )
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out)) stop("Unknown argument: ", arg)
    out[[pair[1L]]] <- as.integer(pair[2L])
  }
  out
}

simulate_scores <- function(n, rho) {
  e1 <- rnorm(n + 1L)
  e2 <- rnorm(n)
  cbind(e1[-1L], rho * e1[-length(e1)] + sqrt(1 - rho^2) * e2)
}

first_crossing <- function(path, critical) {
  hit <- which(path > critical)
  if (length(hit)) hit[1L] else NA_integer_
}

one_replicate <- function(rep_id, rho, delta, critical, seed, m = 1000L,
                          h = 500L, break_after = 100L) {
  set.seed(seed + rep_id)
  scores <- simulate_scores(m + h, rho)
  if (delta > 0) {
    scores[m + break_after + seq_len(h - break_after), 1L] <-
      scores[m + break_after + seq_len(h - break_after), 1L] + delta
  }
  calibration <- scores[seq_len(m), , drop = FALSE]
  monitoring <- scores[m + seq_len(h), , drop = FALSE]
  outer_ldl <- dc_fit_ldl_prewhitener(calibration, ridge = 1e-8, mode = "demix")
  cal_demixed <- dc_apply_ldl_prewhitener(outer_ldl, calibration)
  mon_demixed <- dc_apply_ldl_prewhitener(outer_ldl, monitoring)
  range_fit <- dc_fit_adjusted_range_sn(cal_demixed, ridge = 0)
  full_fit <- dc_fit_full_matrix_sn(cal_demixed, ridge = 0)
  range_path <- dc_sn_monitor_path(range_fit, mon_demixed, gamma = 0)
  full_path <- dc_sn_monitor_path(full_fit, mon_demixed, gamma = 0)
  range_alarm <- first_crossing(range_path, critical[["adjusted_range"]])
  full_alarm <- first_crossing(full_path, critical[["full_matrix"]])
  data.frame(
    replicate = rep_id,
    rho = rho,
    delta = delta,
    adjusted_range_alarm = range_alarm,
    full_matrix_alarm = full_alarm,
    adjusted_range_max = max(range_path),
    full_matrix_max = max(full_path),
    estimated_finite_lag_offdiag = dc_bartlett_lrv_offdiag_ratio(
      cal_demixed, bandwidth = 10L
    ),
    stringsAsFactors = FALSE
  )
}

summarize_cell <- function(raw, break_after = 100L) {
  null_cell <- raw$delta[1L] == 0
  range_reject <- !is.na(raw$adjusted_range_alarm)
  full_reject <- !is.na(raw$full_matrix_alarm)
  if (null_cell) {
    range_out <- range_reject
    full_out <- full_reject
  } else {
    range_out <- !is.na(raw$adjusted_range_alarm) & raw$adjusted_range_alarm > break_after
    full_out <- !is.na(raw$full_matrix_alarm) & raw$full_matrix_alarm > break_after
  }
  difference <- as.integer(range_out) - as.integer(full_out)
  se <- stats::sd(difference) / sqrt(nrow(raw))
  data.frame(
    rho = raw$rho[1L],
    delta = raw$delta[1L],
    estimand = if (null_cell) "null_rejection" else "post_change_alarm_probability",
    nrep = nrow(raw),
    adjusted_range = mean(range_out),
    full_matrix = mean(full_out),
    paired_difference = mean(difference),
    paired_ci95_low = mean(difference) - 1.96 * se,
    paired_ci95_high = mean(difference) + 1.96 * se,
    median_estimated_finite_lag_offdiag = stats::median(raw$estimated_finite_lag_offdiag),
    q10_estimated_finite_lag_offdiag = stats::quantile(raw$estimated_finite_lag_offdiag, 0.10),
    q90_estimated_finite_lag_offdiag = stats::quantile(raw$estimated_finite_lag_offdiag, 0.90),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
critical_file <- file.path(
  script_dir, "outputs", "critical_value_audit",
  "critical_value_grid_and_monte_carlo_audit.csv"
)
critical_table <- utils::read.csv(critical_file, stringsAsFactors = FALSE)
rows <- critical_table[
  critical_table$configuration == "exact1000_500_T0.5" &
    critical_table$K == 2L,
  , drop = FALSE
]
critical <- c(
  adjusted_range = rows$critical_value[rows$method_id == "adjusted_range"],
  full_matrix = rows$critical_value[rows$method_id == "full_matrix"]
)

cells <- rbind(
  expand.grid(
    rho = c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.70),
    delta = 0,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    rho = c(0.20, 0.40),
    delta = c(0.05, 0.15, 0.30),
    stringsAsFactors = FALSE
  )
)

output_dir <- file.path(script_dir, "outputs", "post_ldl_correlation_sweep")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- vector("list", nrow(cells))
for (cell_id in seq_len(nrow(cells))) {
  rho <- cells$rho[cell_id]
  delta <- cells$delta[cell_id]
  cache <- file.path(
    output_dir,
    sprintf("rho_%0.2f_delta_%0.2f_%drep.rds", rho, delta, args$nrep)
  )
  if (file.exists(cache)) {
    raw <- readRDS(cache)
  } else {
    message(sprintf("Starting rho=%0.2f, delta=%0.2f", rho, delta))
    pieces <- parallel::mclapply(
      seq_len(args$nrep),
      one_replicate,
      rho = rho,
      delta = delta,
      critical = critical,
      seed = args$seed + cell_id * 1000000L,
      mc.cores = args$workers,
      mc.preschedule = TRUE
    )
    raw <- do.call(rbind, pieces)
    saveRDS(raw, cache, compress = "xz")
  }
  summaries[[cell_id]] <- summarize_cell(raw)
}

summary_table <- do.call(rbind, summaries)
utils::write.csv(
  summary_table,
  file.path(output_dir, "post_ldl_correlation_sweep_5000rep.csv"),
  row.names = FALSE
)
print(summary_table)
