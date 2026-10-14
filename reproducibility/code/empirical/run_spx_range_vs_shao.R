#!/usr/bin/env Rscript

# Fair actual-data comparison on one retained-score stream. The two methods
# share the FPCA fit, FVE-selected K, LDL demixing, reference period,
# monitoring period, alpha, gamma, and finite horizon. Each method receives
# its own Brownian-limit critical value at the actual horizon 251/997.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else normalizePath("empirical")
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "R", "intraday_interface.R"))

dc_parse_spx_comparison_args <- function(args) {
  defaults <- list(
    data_file = "",
    fit_start = "2014-01-01",
    fit_end = "2015-12-31",
    calibration_start = "2016-01-01",
    calibration_end = "2019-12-31",
    monitor_start = "2020-01-01",
    monitor_end = "2020-12-31",
    fve = 0.90,
    cv_nrep = 250000L,
    cv_train_grid = 997L,
    critical_value_file = "",
    seed = 20260830L
  )
  integer_keys <- c("cv_nrep", "cv_train_grid", "seed")
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[1L] %in% names(defaults)) {
      stop("Unknown argument: ", arg)
    }
    key <- pieces[1L]
    defaults[[key]] <- if (key %in% integer_keys) {
      as.integer(pieces[2L])
    } else if (key == "fve") as.numeric(pieces[2L]) else pieces[2L]
  }
  defaults
}

dc_locate_comparison_spx_data <- function(explicit_path, project_dir) {
  candidates <- c(
    explicit_path,
    Sys.getenv("SPX_DATA_FILE", unset = ""),
    file.path(
      dirname(dirname(project_dir)), "GLMY-20260405", "paper_r_code_bundle_current",
      "empirical", "legacy_source", "data", "SPX.csv"
    )
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop("Could not locate SPX.csv. Supply --data_file=/absolute/path/SPX.csv.")
  }
  normalizePath(existing[1L])
}

args <- dc_parse_spx_comparison_args(commandArgs(trailingOnly = TRUE))
data_file <- dc_locate_comparison_spx_data(args$data_file, project_dir)
message("Reading intraday data from ", data_file)
raw <- dc_read_intraday_ohlc(data_file)
curve_data <- dc_build_intraday_volatility_curves(
  raw,
  start_date = args$fit_start,
  end_date = args$monitor_end,
  bin_minutes = 5L,
  min_day_coverage = 0.95
)
split <- dc_split_intraday_curves(
  curve_data,
  fit_dates = c(args$fit_start, args$fit_end),
  calibration_dates = c(args$calibration_start, args$calibration_end),
  monitor_dates = c(args$monitor_start, args$monitor_end)
)

fit <- dc_fit_score_comparison(
  split$x_fit,
  split$x_calibration,
  K = NULL,
  fve = args$fve,
  K_max = 5L,
  grid = split$grid,
  prewhiten_mode = "demix",
  feature_var_order = 0L
)
message("Selected shared K = ", fit$settings$K)
default_critical_candidates <- c(
  file.path(
    project_dir, "simulations", "componentwise_range", "outputs",
    "critical_value_audit", "critical_value_grid_and_monte_carlo_audit.csv"
  ),
  file.path(
    dirname(project_dir), "data",
    "critical_value_grid_and_monte_carlo_audit.csv"
  )
)
default_critical_file <- default_critical_candidates[
  file.exists(default_critical_candidates)
][1L]
critical_file <- if (nzchar(args$critical_value_file)) {
  normalizePath(args$critical_value_file, mustWork = TRUE)
} else {
  if (is.na(default_critical_file)) {
    stop("High-precision critical-value audit file not found.")
  }
  normalizePath(default_critical_file, mustWork = TRUE)
}
critical_audit <- utils::read.csv(critical_file, stringsAsFactors = FALSE)
critical_rows <- critical_audit[
  critical_audit$configuration == "exact997_251_empirical" &
    critical_audit$K == fit$settings$K,
  , drop = FALSE
]
if (!all(c("adjusted_range", "full_matrix") %in% critical_rows$method_id)) {
  stop("High-precision critical-value audit lacks the selected K.")
}
critical_values <- c(
  adjusted_range = critical_rows$critical_value[
    critical_rows$method_id == "adjusted_range"
  ],
  shao = critical_rows$critical_value[
    critical_rows$method_id == "full_matrix"
  ]
)
monitor <- dc_monitor_score_comparison(
  fit,
  split$x_monitor,
  alpha = 0.05,
  gamma = 0,
  critical_values = critical_values,
  seed = args$seed
)

output_dir <- file.path(script_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- "spx_score_only_range_vs_shao"
normalized_range <- monitor$path$adjusted_range /
  monitor$critical_values[["adjusted_range"]]
normalized_shao <- monitor$path$shao / monitor$critical_values[["shao"]]
path_table <- data.frame(
  date = as.character(split$monitor_dates),
  k = monitor$path$k,
  adjusted_range_statistic = monitor$path$adjusted_range,
  adjusted_range_critical = monitor$critical_values[["adjusted_range"]],
  adjusted_range_threshold_ratio = normalized_range,
  shao_statistic = monitor$path$shao,
  shao_critical = monitor$critical_values[["shao"]],
  shao_threshold_ratio = normalized_shao,
  stringsAsFactors = FALSE
)
utils::write.csv(
  path_table, file.path(output_dir, paste0(prefix, "_path.csv")), row.names = FALSE
)

alarm_row <- function(method, alarm_index, statistic, critical) {
  alarm_date <- if (is.na(alarm_index)) NA_character_ else {
    as.character(split$monitor_dates[alarm_index])
  }
  max_index <- which.max(statistic)
  data.frame(
    method = method,
    alarm_index = alarm_index,
    alarm_date = alarm_date,
    critical_value = critical,
    maximum_statistic = statistic[max_index],
    maximum_threshold_ratio = statistic[max_index] / critical,
    maximum_date = as.character(split$monitor_dates[max_index]),
    stringsAsFactors = FALSE
  )
}
alarm_table <- rbind(
  alarm_row(
    "Adjusted-range SN",
    monitor$alarm[["adjusted_range"]],
    monitor$path$adjusted_range,
    monitor$critical_values[["adjusted_range"]]
  ),
  alarm_row(
    "Shao (2010) quadratic SN",
    monitor$alarm[["shao"]],
    monitor$path$shao,
    monitor$critical_values[["shao"]]
  )
)
utils::write.csv(
  alarm_table, file.path(output_dir, paste0(prefix, "_alarms.csv")),
  row.names = FALSE, na = ""
)

fit_summary <- data.frame(
  data_file = data_file,
  curve_definition = "daily 5-minute log realized-variance profile",
  fit_period = paste(args$fit_start, args$fit_end, sep = " to "),
  calibration_period = paste(args$calibration_start, args$calibration_end, sep = " to "),
  monitoring_period = paste(args$monitor_start, args$monitor_end, sep = " to "),
  n_fit = nrow(split$x_fit),
  n_calibration = nrow(split$x_calibration),
  n_monitor = nrow(split$x_monitor),
  horizon_T = nrow(split$x_monitor) / nrow(split$x_calibration),
  n_grid = ncol(split$x_fit),
  fve_target = args$fve,
  K_selected = fit$settings$K,
  K_max = fit$settings$K_max,
  ldl_mode = fit$settings$prewhiten_mode,
  feature_var_order = fit$settings$feature_var_order,
  adjusted_range_condition = fit$adjusted_range$sn$condition_number,
  shao_condition = fit$shao$sn$condition_number,
  adjusted_range_internal_ldl_applied = FALSE,
  adjusted_range_normalizer_ridge = fit$adjusted_range$sn$ridge_added,
  full_matrix_normalizer_ridge = fit$shao$sn$ridge_added,
  critical_value_file = critical_file,
  critical_value_nrep = unique(critical_rows$nrep),
  critical_value_train_grid = unique(critical_rows$train_grid),
  critical_value_monitor_grid = unique(critical_rows$monitor_grid),
  stringsAsFactors = FALSE
)
utils::write.csv(
  fit_summary, file.path(output_dir, paste0(prefix, "_fit_summary.csv")), row.names = FALSE
)
saveRDS(
  list(
    fit = fit,
    monitor = monitor,
    split_dates = split[c("fit_dates", "calibration_dates", "monitor_dates")],
    curve_settings = curve_data$settings,
    arguments = args
  ),
  file.path(output_dir, paste0(prefix, "_analysis.rds"))
)

source(file.path(script_dir, "plot_spx_audit_figures.R"))
plot_spx_monitor(split$monitor_dates, normalized_range, normalized_shao,
  file.path(output_dir, paste0(prefix, "_monitor.png")))

print(fit_summary)
print(alarm_table)
