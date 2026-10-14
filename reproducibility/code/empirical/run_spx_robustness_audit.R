#!/usr/bin/env Rscript

# Prespecified empirical checks for JTSA preparation: reference-period
# summaries, a historical 2019 monitoring exercise, preprocessing comparisons,
# and exact reconstruction metadata for the TickData minute interface.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else normalizePath("empirical")
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "R", "intraday_interface.R"))

data_file <- Sys.getenv(
  "SPX_DATA_FILE",
  unset = file.path(
    dirname(dirname(project_dir)), "GLMY-20260405", "paper_r_code_bundle_current",
    "empirical", "legacy_source", "data", "SPX.csv"
  )
)
data_file <- normalizePath(data_file, mustWork = TRUE)
output_dir <- file.path(script_dir, "outputs", "spx_robustness_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", data_file)
raw <- dc_read_intraday_ohlc(data_file)
primary_curves <- dc_build_intraday_volatility_curves(
  raw, start_date = "2014-01-01", end_date = "2020-12-31",
  bin_minutes = 5L, min_day_coverage = 0.95, spline_df = 12L
)
primary_split <- dc_split_intraday_curves(
  primary_curves,
  fit_dates = c("2014-01-01", "2015-12-31"),
  calibration_dates = c("2016-01-01", "2019-12-31"),
  monitor_dates = c("2020-01-01", "2020-12-31")
)

# Normalization-sample score audit.
fpca <- dc_fit_fpca(
  primary_split$x_fit, K = NULL, L = 1L, fve = 0.90,
  K_max = 5L,
  grid = primary_split$grid
)
projected_cal <- dc_project_fpca(fpca, primary_split$x_calibration)$retained
centered_cal <- sweep(projected_cal, 2L, colMeans(projected_cal), FUN = "-")
score_table <- data.frame(
  date = primary_split$calibration_dates,
  centered_cal,
  stringsAsFactors = FALSE
)
names(score_table)[-1L] <- paste0("score", seq_len(ncol(centered_cal)))
utils::write.csv(score_table, file.path(output_dir, "spx_calibration_scores_2016_2019.csv"), row.names = FALSE)

annual_rows <- list()
years <- format(primary_split$calibration_dates, "%Y")
for (year in sort(unique(years))) {
  idx <- which(years == year)
  for (j in seq_len(ncol(centered_cal))) {
    annual_rows[[length(annual_rows) + 1L]] <- data.frame(
      year = year, score = j, n = length(idx),
      mean = mean(centered_cal[idx, j]), sd = stats::sd(centered_cal[idx, j]),
      stringsAsFactors = FALSE
    )
  }
}
annual_summary <- do.call(rbind, annual_rows)
utils::write.csv(annual_summary, file.path(output_dir, "spx_calibration_annual_summary.csv"), row.names = FALSE)

acf_rows <- list()
for (j in seq_len(ncol(centered_cal))) {
  acf_score <- stats::acf(centered_cal[, j], lag.max = 20L, plot = FALSE)$acf[-1L]
  acf_square <- stats::acf(centered_cal[, j]^2, lag.max = 20L, plot = FALSE)$acf[-1L]
  acf_rows[[j]] <- data.frame(
    score = j, lag = seq_len(20L), acf_score = as.numeric(acf_score),
    acf_squared_score = as.numeric(acf_square)
  )
}
utils::write.csv(do.call(rbind, acf_rows), file.path(output_dir, "spx_calibration_acf.csv"), row.names = FALSE)

bridge <- apply(centered_cal, 2L, cumsum) / sqrt(nrow(centered_cal))
grDevices::png(file.path(output_dir, "spx_calibration_scores_and_bridges.png"),
               width = 1800, height = 1400, res = 180)
graphics::par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1))
graphics::matplot(primary_split$calibration_dates, centered_cal, type = "l", lty = 1,
                  col = c("#00008B", "#2E5AAC", "#5A7FC0", "#88A4D0"),
                  xlab = "Date", ylab = "Centered retained score",
                  main = "Normalization-sample retained scores")
graphics::matplot(primary_split$calibration_dates, bridge, type = "l", lty = 1,
                  col = c("#00008B", "#2E5AAC", "#5A7FC0", "#88A4D0"),
                  xlab = "Date", ylab = "Scaled cumulative bridge",
                  main = "Normalization-block bridge paths")
grDevices::dev.off()

half <- floor(nrow(centered_cal) / 2L)
normalizer_audit <- lapply(
  list(first_half = centered_cal[seq_len(half), , drop = FALSE],
       second_half = centered_cal[(half + 1L):nrow(centered_cal), , drop = FALSE]),
  function(z) {
    ldl <- dc_fit_ldl_prewhitener(z, ridge = 1e-8, mode = "demix")
    demixed <- dc_apply_ldl_prewhitener(ldl, z)
    list(range = dc_fit_adjusted_range_sn(demixed, ridge = 0),
         full = dc_fit_full_matrix_sn(demixed, ridge = 0))
  }
)
normalizer_summary <- data.frame(
  half = c("first_half", "second_half"),
  adjusted_range_condition = vapply(normalizer_audit, function(x) x$range$condition_number, numeric(1L)),
  full_matrix_condition = vapply(normalizer_audit, function(x) x$full$condition_number, numeric(1L)),
  range1 = vapply(normalizer_audit, function(x) x$range$adjusted_ranges[1L], numeric(1L)),
  range2 = vapply(normalizer_audit, function(x) x$range$adjusted_ranges[2L], numeric(1L)),
  range3 = vapply(normalizer_audit, function(x) x$range$adjusted_ranges[3L], numeric(1L)),
  range4 = vapply(normalizer_audit, function(x) x$range$adjusted_ranges[4L], numeric(1L))
)
utils::write.csv(normalizer_summary, file.path(output_dir, "spx_calibration_half_normalizers.csv"), row.names = FALSE)

# Historical 2019 monitoring exercise. Legacy variable and output names retain
# "placebo" only for compatibility with the frozen critical-value identifiers.
placebo_split <- dc_split_intraday_curves(
  primary_curves,
  fit_dates = c("2014-01-01", "2015-12-31"),
  calibration_dates = c("2016-01-01", "2018-12-31"),
  monitor_dates = c("2019-01-01", "2019-12-31")
)
placebo_fit <- dc_fit_score_comparison(
  placebo_split$x_fit, placebo_split$x_calibration, K = NULL, fve = 0.90,
  K_max = 5L,
  grid = placebo_split$grid, prewhiten_mode = "demix", feature_var_order = 0L
)
critical_table <- utils::read.csv(file.path(
  project_dir, "simulations", "componentwise_range", "outputs",
  "critical_value_audit", "critical_value_grid_and_monte_carlo_audit.csv"
), stringsAsFactors = FALSE)
placebo_cv_rows <- critical_table[
  critical_table$configuration == "exact748_249_placebo" &
    critical_table$K == placebo_fit$settings$K,
  , drop = FALSE
]
placebo_cv <- c(
  adjusted_range = placebo_cv_rows$critical_value[placebo_cv_rows$method_id == "adjusted_range"],
  shao = placebo_cv_rows$critical_value[placebo_cv_rows$method_id == "full_matrix"]
)
placebo_monitor <- dc_monitor_score_comparison(
  placebo_fit, placebo_split$x_monitor, alpha = 0.05, gamma = 0,
  critical_values = placebo_cv
)
placebo_path <- data.frame(
  date = placebo_split$monitor_dates,
  k = placebo_monitor$path$k,
  adjusted_range_ratio = placebo_monitor$path$adjusted_range / placebo_cv[["adjusted_range"]],
  full_matrix_ratio = placebo_monitor$path$shao / placebo_cv[["shao"]]
)
utils::write.csv(placebo_path, file.path(output_dir, "spx_2019_placebo_path.csv"), row.names = FALSE)
placebo_alarm <- data.frame(
  method = c("adjusted_range", "full_matrix"),
  alarm_index = unname(placebo_monitor$alarm[c("adjusted_range", "shao")]),
  alarm_date = c(
    if (is.na(placebo_monitor$alarm[["adjusted_range"]])) NA else as.character(placebo_split$monitor_dates[placebo_monitor$alarm[["adjusted_range"]]]),
    if (is.na(placebo_monitor$alarm[["shao"]])) NA else as.character(placebo_split$monitor_dates[placebo_monitor$alarm[["shao"]]])
  ),
  maximum_ratio = c(max(placebo_path$adjusted_range_ratio), max(placebo_path$full_matrix_ratio)),
  K = placebo_fit$settings$K,
  n_calibration = nrow(placebo_split$x_calibration),
  n_monitor = nrow(placebo_split$x_monitor)
)
utils::write.csv(placebo_alarm, file.path(output_dir, "spx_2019_placebo_alarm.csv"), row.names = FALSE, na = "")

# Prespecified preprocessing sensitivity grid.  The main 18-cell grid crosses
# FVE, spline dimension and bin width at 95% coverage; two additional rows use
# complete-bin days for the primary FVE/spline choices.
specs <- expand.grid(
  fve = c(0.85, 0.90, 0.95), spline_df = c(8L, 12L, 16L),
  bin_minutes = c(5L, 10L), min_day_coverage = 0.95,
  stringsAsFactors = FALSE
)
specs <- rbind(specs, data.frame(
  fve = 0.90, spline_df = 12L, bin_minutes = c(5L, 10L),
  min_day_coverage = 1.00
))
critical_cache <- new.env(parent = emptyenv())
get_sensitivity_cv <- function(K, m, h, seed) {
  key <- paste(K, m, h, sep = "_")
  if (exists(key, envir = critical_cache, inherits = FALSE)) return(get(key, envir = critical_cache))
  if (K == 4L && m == 997L && h == 251L) {
    exact_rows <- critical_table[
      critical_table$configuration == "exact997_251_empirical" &
        critical_table$K == K,
      , drop = FALSE
    ]
    cv <- c(
      adjusted_range = exact_rows$critical_value[exact_rows$method_id == "adjusted_range"],
      shao = exact_rows$critical_value[exact_rows$method_id == "full_matrix"]
    )
    assign(key, cv, envir = critical_cache)
    return(cv)
  }
  cv <- c(
    adjusted_range = dc_simulate_adjusted_range_critical_value(
      K, T = h / m, alpha = 0.05, nrep = 20000L, train_grid = m, seed = seed
    )$critical_value,
    shao = dc_simulate_sn_critical_value(
      K, T = h / m, alpha = 0.05, nrep = 20000L, train_grid = m, seed = seed + 500000L,
      ridge = 0
    )$critical_value
  )
  assign(key, cv, envir = critical_cache)
  cv
}

curve_cache <- new.env(parent = emptyenv())
sensitivity_rows <- vector("list", nrow(specs))
for (i in seq_len(nrow(specs))) {
  spec <- specs[i, ]
  curve_key <- paste(spec$bin_minutes, spec$min_day_coverage, spec$spline_df, sep = "_")
  if (exists(curve_key, envir = curve_cache, inherits = FALSE)) {
    curves_i <- get(curve_key, envir = curve_cache)
  } else {
    message("Building sensitivity curves ", curve_key)
    curves_i <- dc_build_intraday_volatility_curves(
      raw, start_date = "2014-01-01", end_date = "2020-12-31",
      bin_minutes = spec$bin_minutes,
      min_day_coverage = spec$min_day_coverage,
      spline_df = spec$spline_df
    )
    assign(curve_key, curves_i, envir = curve_cache)
  }
  split_i <- dc_split_intraday_curves(
    curves_i,
    fit_dates = c("2014-01-01", "2015-12-31"),
    calibration_dates = c("2016-01-01", "2019-12-31"),
    monitor_dates = c("2020-01-01", "2020-12-31")
  )
  fit_i <- dc_fit_score_comparison(
    split_i$x_fit, split_i$x_calibration, K = NULL, fve = spec$fve,
    K_max = as.integer(spec$spline_df) - 1L,
    grid = split_i$grid, prewhiten_mode = "demix", feature_var_order = 0L
  )
  cv_i <- get_sensitivity_cv(
    fit_i$settings$K, nrow(split_i$x_calibration), nrow(split_i$x_monitor),
    seed = 20260831L + i * 10000L
  )
  monitor_i <- dc_monitor_score_comparison(
    fit_i, split_i$x_monitor, alpha = 0.05, gamma = 0, critical_values = cv_i
  )
  alarm_date <- function(index) if (is.na(index)) NA_character_ else as.character(split_i$monitor_dates[index])
  sensitivity_rows[[i]] <- data.frame(
    fve = spec$fve, spline_df = spec$spline_df,
    bin_minutes = spec$bin_minutes, min_day_coverage = spec$min_day_coverage,
    n_fit = nrow(split_i$x_fit), n_calibration = nrow(split_i$x_calibration),
    n_monitor = nrow(split_i$x_monitor), K = fit_i$settings$K,
    adjusted_range_alarm = alarm_date(monitor_i$alarm[["adjusted_range"]]),
    full_matrix_alarm = alarm_date(monitor_i$alarm[["shao"]]),
    adjusted_range_max_ratio = max(monitor_i$path$adjusted_range / cv_i[["adjusted_range"]]),
    full_matrix_max_ratio = max(monitor_i$path$shao / cv_i[["shao"]]),
    stringsAsFactors = FALSE
  )
}
utils::write.csv(do.call(rbind, sensitivity_rows),
                 file.path(output_dir, "spx_preprocessing_sensitivity.csv"), row.names = FALSE, na = "")

# Reconstruction metadata and quality-control details.
removed <- primary_curves$quality[primary_curves$quality$coverage < 0.95, , drop = FALSE]
utils::write.csv(removed, file.path(output_dir, "spx_removed_dates_primary.csv"), row.names = FALSE)
utils::write.csv(primary_curves$kept_quality,
                 file.path(output_dir, "spx_retained_day_interpolation_audit.csv"), row.names = FALSE)
basis <- splines::bs(primary_curves$grid, df = 12L, degree = 3L, intercept = TRUE)
metadata <- data.frame(
  field = c(
    "vendor", "instrument", "timestamp_convention", "session",
    "bin_endpoints", "interpolation_scale", "spline_fit", "spline_df",
    "boundary_knots", "interior_knots", "fpca_inner_product"
  ),
  value = c(
    "TickData", "S&P 500 cash index (SPX)",
    "exchange-local minute stamps as recorded; 09:31 is treated as the first regular-session minute",
    "09:31--16:00 New York exchange-local time; trading calendar carries DST",
    "[09:31,09:35], [09:36,09:40], ..., [15:56,16:00] for five-minute bins",
    "realized variance before logarithm", "unpenalized least squares", "12",
    paste(attr(basis, "Boundary.knots"), collapse = ";"),
    paste(attr(basis, "knots"), collapse = ";"),
    "trapezoidal weights on normalized bin-midpoint grid"
  )
)
utils::write.csv(metadata, file.path(output_dir, "spx_reconstruction_metadata.csv"), row.names = FALSE)

print(placebo_alarm)
print(do.call(rbind, sensitivity_rows))
