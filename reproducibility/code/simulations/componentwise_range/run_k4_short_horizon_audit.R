#!/usr/bin/env Rscript

# Five-thousand-replication audit at the empirical dimension and horizon.
# Each two-coordinate block has identity lag-zero covariance and post-LDL LRV
# correlation rho, so rho=0 is assumption-aligned and rho>0 is controlled
# misspecification.

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
    seed = 20260901L
  )
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out)) stop("Unknown argument: ", arg)
    out[[pair[1L]]] <- as.integer(pair[2L])
  }
  out
}

simulate_pair <- function(n, rho) {
  e1 <- stats::rnorm(n + 1L)
  e2 <- stats::rnorm(n)
  cbind(e1[-1L], rho * e1[-length(e1)] + sqrt(1 - rho^2) * e2)
}

first_crossing <- function(path, critical) {
  hit <- which(path > critical)
  if (length(hit)) hit[1L] else NA_integer_
}

one_replicate <- function(rep_id, rho, critical, seed, m = 997L, h = 251L) {
  set.seed(seed + rep_id)
  scores <- cbind(simulate_pair(m + h, rho), simulate_pair(m + h, rho))
  calibration <- scores[seq_len(m), , drop = FALSE]
  monitoring <- scores[m + seq_len(h), , drop = FALSE]
  ldl <- dc_fit_ldl_prewhitener(calibration, ridge = 1e-8, mode = "demix")
  cal_demixed <- dc_apply_ldl_prewhitener(ldl, calibration)
  mon_demixed <- dc_apply_ldl_prewhitener(ldl, monitoring)
  range_fit <- dc_fit_adjusted_range_sn(cal_demixed, ridge = 0)
  full_fit <- dc_fit_full_matrix_sn(cal_demixed, ridge = 0)
  range_path <- dc_sn_monitor_path(range_fit, mon_demixed, gamma = 0)
  full_path <- dc_sn_monitor_path(full_fit, mon_demixed, gamma = 0)
  data.frame(
    replicate = rep_id,
    rho = rho,
    adjusted_range_alarm = first_crossing(range_path, critical[["adjusted_range"]]),
    full_matrix_alarm = first_crossing(full_path, critical[["full_matrix"]]),
    adjusted_range_max = max(range_path),
    full_matrix_max = max(full_path),
    estimated_finite_lag_offdiag = dc_bartlett_lrv_offdiag_ratio(cal_demixed, bandwidth = 10L)
  )
}

summarize_cell <- function(raw) {
  range_reject <- !is.na(raw$adjusted_range_alarm)
  full_reject <- !is.na(raw$full_matrix_alarm)
  difference <- as.integer(range_reject) - as.integer(full_reject)
  se <- stats::sd(difference) / sqrt(nrow(raw))
  data.frame(
    K = 4L,
    m = 997L,
    h = 251L,
    T = 251 / 997,
    rho = raw$rho[1L],
    nrep = nrow(raw),
    adjusted_range_size = mean(range_reject),
    full_matrix_size = mean(full_reject),
    paired_difference = mean(difference),
    paired_ci95_low = mean(difference) - 1.96 * se,
    paired_ci95_high = mean(difference) + 1.96 * se,
    median_estimated_finite_lag_offdiag = stats::median(raw$estimated_finite_lag_offdiag),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
critical_table <- utils::read.csv(file.path(
  script_dir, "outputs", "critical_value_audit",
  "critical_value_grid_and_monte_carlo_audit.csv"
), stringsAsFactors = FALSE)
critical_rows <- critical_table[
  critical_table$configuration == "exact997_251_empirical" &
    critical_table$K == 4L,
  , drop = FALSE
]
critical <- c(
  adjusted_range = critical_rows$critical_value[critical_rows$method_id == "adjusted_range"],
  full_matrix = critical_rows$critical_value[critical_rows$method_id == "full_matrix"]
)

output_dir <- file.path(script_dir, "outputs", "k4_short_horizon_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()
for (cell_id in seq_along(c(0, 0.20))) {
  rho <- c(0, 0.20)[cell_id]
  cache <- file.path(output_dir, sprintf("rho_%0.2f_%drep.rds", rho, args$nrep))
  raw <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.data.frame(raw)) {
    pieces <- parallel::mclapply(
      seq_len(args$nrep), one_replicate, rho = rho, critical = critical,
      seed = args$seed + cell_id * 1000000L,
      mc.cores = args$workers, mc.preschedule = TRUE
    )
    raw <- do.call(rbind, pieces)
    saveRDS(raw, cache, compress = "xz")
  }
  summaries[[cell_id]] <- summarize_cell(raw)
}
summary_table <- do.call(rbind, summaries)
utils::write.csv(
  summary_table,
  file.path(output_dir, "k4_short_horizon_size_5000rep.csv"),
  row.names = FALSE
)
print(summary_table)
