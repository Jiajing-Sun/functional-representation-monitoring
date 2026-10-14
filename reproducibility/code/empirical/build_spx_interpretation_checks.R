#!/usr/bin/env Rscript

# Reproducible interpretation checks for the S&P 500 retained-score analysis.
# This script uses the saved primary-analysis object and does not require the
# licensed minute-level source data.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else {
  normalizePath("reproducibility/code/empirical")
}
project_dir <- normalizePath(file.path(script_dir, ".."))
paper_dir <- normalizePath(file.path(project_dir, ".."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))

analysis_file <- file.path(
  script_dir, "outputs", "spx_score_only_range_vs_shao_analysis.rds"
)
if (!file.exists(analysis_file)) stop("Primary SPX analysis object not found.")
analysis <- readRDS(analysis_file)

data_dir <- file.path(paper_dir, "data")
figure_dir <- file.path(dirname(paper_dir), "figures", "revised")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fit <- analysis$fit
monitor <- analysis$monitor
fpca <- fit$adjusted_range$fpca
K_primary <- fit$settings$K
monitor_dates <- analysis$split_dates$monitor_dates

# Retained eigenfunctions and fitted-spectrum diagnostics.
shares <- fpca$eigenvalues / sum(fpca$eigenvalues)
eigenfunction_table <- data.frame(
  bin = seq_len(fpca$n_grid),
  minutes_from_session_start = seq(2.5, by = 5, length.out = fpca$n_grid),
  fpca$eigenfunctions[, seq_len(K_primary), drop = FALSE],
  check.names = FALSE
)
utils::write.csv(
  eigenfunction_table,
  file.path(data_dir, "spx_retained_eigenfunctions.csv"),
  row.names = FALSE
)

spectrum_table <- data.frame(
  component = seq_len(min(10L, length(shares))),
  eigenvalue = fpca$eigenvalues[seq_len(min(10L, length(shares)))],
  variance_share = shares[seq_len(min(10L, length(shares)))],
  cumulative_fve = fpca$fve_path[seq_len(min(10L, length(shares)))],
  next_eigengap = c(
    head(fpca$eigenvalues, min(10L, length(shares))) -
      fpca$eigenvalues[2:(min(10L, length(shares)) + 1L)]
  ),
  relative_next_eigengap = c(
    (head(fpca$eigenvalues, min(10L, length(shares))) -
       fpca$eigenvalues[2:(min(10L, length(shares)) + 1L)]) /
      head(fpca$eigenvalues, min(10L, length(shares)))
  )
)
utils::write.csv(
  spectrum_table,
  file.path(data_dir, "spx_fpca_dimension_diagnostics.csv"),
  row.names = FALSE
)

# Fixed-K monitoring checks use the leading subvector of the same FPCA score
# stream. The leading LDL factors are unchanged by dropping trailing scores.
critical_file <- file.path(data_dir, "critical_value_grid_and_monte_carlo_audit.csv")
critical_audit <- utils::read.csv(critical_file, stringsAsFactors = FALSE)
critical_audit <- critical_audit[
  critical_audit$configuration == "exact997_251_empirical", , drop = FALSE
]

calibration_scores <- fit$adjusted_range$calibration_sn_scores
monitor_scores <- monitor$adjusted_range$sn_scores
m <- nrow(calibration_scores)
h <- nrow(monitor_scores)

fixed_k_rows <- lapply(2:K_primary, function(K) {
  calibration_K <- calibration_scores[, seq_len(K), drop = FALSE]
  monitor_K <- monitor_scores[, seq_len(K), drop = FALSE]
  range_sn <- dc_fit_adjusted_range_sn(calibration_K, ridge = 0)
  quadratic_sn <- dc_fit_full_matrix_sn(calibration_K, ridge = 0)
  range_path <- dc_sn_monitor_path(range_sn, monitor_K, gamma = 0)
  quadratic_path <- dc_sn_monitor_path(quadratic_sn, monitor_K, gamma = 0)
  range_cv <- critical_audit$critical_value[
    critical_audit$K == K & critical_audit$method_id == "adjusted_range"
  ]
  quadratic_cv <- critical_audit$critical_value[
    critical_audit$K == K & critical_audit$method_id == "full_matrix"
  ]
  alarm_summary <- function(method, path, critical) {
    crossing <- which(path > critical)
    alarm_index <- if (length(crossing)) crossing[1L] else NA_integer_
    data.frame(
      K = K,
      method = method,
      critical_value = critical,
      alarm_index = alarm_index,
      alarm_date = if (is.na(alarm_index)) NA_character_ else
        as.character(monitor_dates[alarm_index]),
      maximum_threshold_ratio = max(path / critical),
      maximum_date = as.character(monitor_dates[which.max(path)]),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    alarm_summary("Adjusted-range SN", range_path, range_cv),
    alarm_summary("Quadratic SN", quadratic_path, quadratic_cv)
  )
})
fixed_k_table <- do.call(rbind, fixed_k_rows)
utils::write.csv(
  fixed_k_table,
  file.path(data_dir, "spx_fixed_k_sensitivity.csv"),
  row.names = FALSE, na = ""
)

# Exact additive coordinate decomposition at each method's first alarm.
coordinate_decomposition <- function(method, sn_fit, scores, alarm_index) {
  cumulative <- colSums(scores[seq_len(alarm_index), , drop = FALSE])
  scale <- m * dc_boundary_sq(alarm_index / m, gamma = 0)
  contribution <- as.numeric(cumulative * (sn_fit$inverse %*% cumulative) / scale)
  statistic <- sum(contribution)
  data.frame(
    method = method,
    alarm_index = alarm_index,
    alarm_date = as.character(monitor_dates[alarm_index]),
    coordinate = seq_along(contribution),
    contribution = contribution,
    signed_share = contribution / statistic,
    statistic = statistic,
    stringsAsFactors = FALSE
  )
}
contribution_table <- rbind(
  coordinate_decomposition(
    "Adjusted-range SN", fit$adjusted_range$sn,
    monitor$adjusted_range$sn_scores, monitor$alarm[["adjusted_range"]]
  ),
  coordinate_decomposition(
    "Quadratic SN", fit$shao$sn,
    monitor$shao$sn_scores, monitor$alarm[["shao"]]
  )
)
utils::write.csv(
  contribution_table,
  file.path(data_dir, "spx_alarm_coordinate_contributions.csv"),
  row.names = FALSE
)

# Mean change in the original FPCA coordinates at the adjusted-range alarm.
range_alarm <- monitor$alarm[["adjusted_range"]]
demixed_mean <- colMeans(
  monitor$adjusted_range$sn_scores[seq_len(range_alarm), , drop = FALSE]
)
raw_score_mean <- as.numeric(demixed_mean %*% t(fit$adjusted_range$prewhitener$C))
reconstructed_mean <- as.numeric(
  fpca$eigenfunctions[, seq_len(K_primary), drop = FALSE] %*% raw_score_mean
)
mean_shift_table <- data.frame(
  coordinate = seq_len(K_primary),
  mean_score_shift = raw_score_mean,
  eigenfunction_variance_share = shares[seq_len(K_primary)],
  stringsAsFactors = FALSE
)
attr(mean_shift_table, "reconstructed_min") <- min(reconstructed_mean)
attr(mean_shift_table, "reconstructed_max") <- max(reconstructed_mean)
utils::write.csv(
  mean_shift_table,
  file.path(data_dir, "spx_retained_mean_shift_at_range_alarm.csv"),
  row.names = FALSE
)

# Reader-facing 2-by-2 plot of the retained eigenfunctions.
plot_file <- file.path(figure_dir, "spx_retained_eigenfunctions.png")
grDevices::png(plot_file, width = 1800, height = 1250, res = 180)
old_par <- graphics::par(no.readonly = TRUE)
on.exit({
  graphics::par(old_par)
  grDevices::dev.off()
}, add = TRUE)
graphics::par(
  mfrow = c(2, 2), mar = c(4.0, 4.2, 2.4, 1.0),
  oma = c(0.7, 0.8, 0.3, 0.2), las = 1
)
x_minutes <- eigenfunction_table$minutes_from_session_start
tick_minutes <- c(0, 90, 180, 270, 390)
tick_labels <- c("09:30", "11:00", "12:30", "14:00", "16:00")
for (j in seq_len(K_primary)) {
  y <- fpca$eigenfunctions[, j]
  graphics::plot(
    x_minutes, y, type = "l", lwd = 2.4, col = "black",
    xaxt = "n", xlab = "Intraday time", ylab = "Eigenfunction value",
    main = sprintf("PC%d (%.2f%%)", j, 100 * shares[j])
  )
  graphics::axis(1, at = tick_minutes, labels = tick_labels)
  graphics::abline(h = 0, col = "grey55", lty = 2, lwd = 1)
  graphics::box()
}
grDevices::dev.off()
on.exit(NULL, add = FALSE)

message("Wrote interpretation outputs to ", data_dir)
message("Wrote eigenfunction figure to ", plot_file)
