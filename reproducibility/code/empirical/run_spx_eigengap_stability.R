#!/usr/bin/env Rscript

# Stability checks for the empirical FPCA directions used in the paper.
#
# The first check rotates the fitted PC3--PC4 plane while holding its span
# fixed.  The second applies a stationary block bootstrap to the 2014--2015
# fit block, re-estimates FPCA, re-applies the 90% FVE rule, refits the LDL map
# and both self-normalizers, and then monitors the unchanged 2020 curves.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1L])))
} else {
  normalizePath("reproducibility/code/empirical")
}
project_dir <- normalizePath(file.path(script_dir, ".."))
paper_dir <- normalizePath(file.path(project_dir, ".."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "R", "intraday_interface.R"))

parse_args <- function(args) {
  defaults <- list(
    data_file = "",
    nrep = 1000L,
    workers = max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L)),
    seed = 20260902L,
    rotation_step = 1
  )
  integer_keys <- c("nrep", "workers", "seed")
  numeric_keys <- "rotation_step"
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[1L] %in% names(defaults)) {
      stop("Unknown argument: ", arg)
    }
    key <- pieces[1L]
    defaults[[key]] <- if (key %in% integer_keys) {
      as.integer(pieces[2L])
    } else if (key %in% numeric_keys) {
      as.numeric(pieces[2L])
    } else pieces[2L]
  }
  if (defaults$nrep < 100L) stop("nrep must be at least 100.")
  if (defaults$workers < 1L) stop("workers must be positive.")
  if (!is.finite(defaults$rotation_step) || defaults$rotation_step <= 0 ||
      defaults$rotation_step > 15) stop("rotation_step must lie in (0,15].")
  defaults
}

locate_data <- function(explicit_path, project_dir) {
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

first_crossing <- function(path, critical) {
  crossing <- which(path > critical)
  if (length(crossing)) crossing[1L] else NA_integer_
}

fit_score_paths <- function(reference_scores, monitoring_scores, critical_values) {
  ldl <- dc_fit_ldl_prewhitener(reference_scores, ridge = 1e-8, mode = "demix")
  reference_demixed <- dc_apply_ldl_prewhitener(ldl, reference_scores)
  monitoring_demixed <- dc_apply_ldl_prewhitener(ldl, monitoring_scores)
  range_fit <- dc_fit_adjusted_range_sn(reference_demixed, ridge = 0)
  quadratic_fit <- dc_fit_full_matrix_sn(reference_demixed, ridge = 0)
  range_path <- dc_sn_monitor_path(range_fit, monitoring_demixed, gamma = 0)
  quadratic_path <- dc_sn_monitor_path(quadratic_fit, monitoring_demixed, gamma = 0)
  list(
    range_path = range_path,
    quadratic_path = quadratic_path,
    range_alarm = first_crossing(range_path, critical_values[["adjusted_range"]]),
    quadratic_alarm = first_crossing(quadratic_path, critical_values[["quadratic"]]),
    range_max_ratio = max(range_path / critical_values[["adjusted_range"]]),
    quadratic_max_ratio = max(quadratic_path / critical_values[["quadratic"]])
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
analysis_file <- file.path(
  script_dir, "outputs", "spx_score_only_range_vs_shao_analysis.rds"
)
if (!file.exists(analysis_file)) stop("Primary SPX analysis object not found.")
analysis <- readRDS(analysis_file)

requested_data_file <- if (nzchar(args$data_file)) {
  args$data_file
} else analysis$arguments$data_file
data_file <- locate_data(requested_data_file, project_dir)
message("Reading intraday data from ", data_file)
raw <- dc_read_intraday_ohlc(data_file)
curve_data <- dc_build_intraday_volatility_curves(
  raw,
  start_date = "2014-01-01",
  end_date = "2020-12-31",
  bin_minutes = 5L,
  min_day_coverage = 0.95,
  spline_df = 12L
)
split <- dc_split_intraday_curves(
  curve_data,
  fit_dates = c("2014-01-01", "2015-12-31"),
  calibration_dates = c("2016-01-01", "2019-12-31"),
  monitor_dates = c("2020-01-01", "2020-12-31")
)

primary_fpca <- analysis$fit$adjusted_range$fpca
if (primary_fpca$K != 4L) stop("The stored primary analysis does not use K=4.")
critical_file <- file.path(
  paper_dir, "data", "critical_value_grid_and_monte_carlo_audit.csv"
)
critical_table <- utils::read.csv(critical_file, stringsAsFactors = FALSE)
critical_table <- critical_table[
  critical_table$configuration == "exact997_251_empirical", , drop = FALSE
]
critical_for_K <- function(K) {
  rows <- critical_table[critical_table$K == K, , drop = FALSE]
  values <- c(
    adjusted_range = rows$critical_value[rows$method_id == "adjusted_range"],
    quadratic = rows$critical_value[rows$method_id == "full_matrix"]
  )
  if (length(values) != 2L || any(!is.finite(values))) {
    stop("Missing matched-grid critical values for K=", K, ".")
  }
  values
}

# Rotation of the empirically close PC3--PC4 directions.
reference_projection <- dc_project_fpca(primary_fpca, split$x_calibration)
monitor_projection <- dc_project_fpca(primary_fpca, split$x_monitor)
reference_scores <- reference_projection$retained
monitor_scores <- monitor_projection$retained
angles <- seq(0, 90, by = args$rotation_step)
rotation_rows <- lapply(angles, function(angle) {
  theta <- angle * pi / 180
  rotation <- diag(4L)
  rotation[3:4, 3:4] <- matrix(
    c(cos(theta), sin(theta), -sin(theta), cos(theta)),
    nrow = 2L, byrow = TRUE
  )
  result <- fit_score_paths(
    reference_scores %*% rotation,
    monitor_scores %*% rotation,
    critical_for_K(4L)
  )
  data.frame(
    angle_degrees = angle,
    adjusted_range_alarm_index = result$range_alarm,
    adjusted_range_alarm_date = if (is.na(result$range_alarm)) NA_character_ else
      as.character(split$monitor_dates[result$range_alarm]),
    adjusted_range_maximum_ratio = result$range_max_ratio,
    quadratic_alarm_index = result$quadratic_alarm,
    quadratic_alarm_date = if (is.na(result$quadratic_alarm)) NA_character_ else
      as.character(split$monitor_dates[result$quadratic_alarm]),
    quadratic_maximum_ratio = result$quadratic_max_ratio,
    stringsAsFactors = FALSE
  )
})
rotation_table <- do.call(rbind, rotation_rows)

# Stationary block bootstrap of the fit block.  The block length is selected
# once from the fitted score stream and then held fixed across replications.
fit_projection <- dc_project_fpca(primary_fpca, split$x_fit)$retained
block_length <- dc_select_stationary_block_length(fit_projection)
set.seed(args$seed)
bootstrap_indices <- replicate(
  args$nrep,
  dc_stationary_bootstrap_indices(nrow(split$x_fit), nrow(split$x_fit), block_length),
  simplify = FALSE
)
primary_basis <- primary_fpca$basis_weighted[, seq_len(4L), drop = FALSE]

bootstrap_one <- function(rep_index) {
  tryCatch({
    fitted <- dc_fit_score_comparison(
      split$x_fit[bootstrap_indices[[rep_index]], , drop = FALSE],
      split$x_calibration,
      K = NULL,
      fve = 0.90,
      K_max = 5L,
      grid = split$grid,
      prewhiten_mode = "demix",
      feature_var_order = 0L
    )
    K <- fitted$settings$K
    monitored <- dc_monitor_score_comparison(
      fitted,
      split$x_monitor,
      alpha = 0.05,
      gamma = 0,
      critical_values = c(
        adjusted_range = critical_for_K(K)[["adjusted_range"]],
        shao = critical_for_K(K)[["quadratic"]]
      ),
      seed = args$seed + rep_index
    )
    fpca <- fitted$adjusted_range$fpca
    shares <- fpca$eigenvalues / sum(fpca$eigenvalues)
    comparison_K <- min(4L, K)
    overlap <- crossprod(
      primary_basis[, seq_len(comparison_K), drop = FALSE],
      fpca$basis_weighted[, seq_len(comparison_K), drop = FALSE]
    )
    smallest_cosine <- min(svd(overlap, nu = 0L, nv = 0L)$d)
    smallest_cosine <- min(1, max(0, smallest_cosine))
    data.frame(
      replication = rep_index,
      success = TRUE,
      K_selected = K,
      fve_at_3 = fpca$fve_path[3L],
      fve_at_4 = fpca$fve_path[4L],
      pc3_share = shares[3L],
      pc4_share = shares[4L],
      pc5_share = shares[5L],
      gap_3_4 = shares[3L] - shares[4L],
      gap_4_5 = shares[4L] - shares[5L],
      retained_subspace_max_angle_degrees = acos(smallest_cosine) * 180 / pi,
      adjusted_range_alarm_index = monitored$alarm[["adjusted_range"]],
      adjusted_range_alarm_date = if (is.na(monitored$alarm[["adjusted_range"]]))
        NA_character_ else as.character(split$monitor_dates[monitored$alarm[["adjusted_range"]]]),
      adjusted_range_maximum_ratio = max(
        monitored$path$adjusted_range / monitored$critical_values[["adjusted_range"]]
      ),
      quadratic_alarm_index = monitored$alarm[["shao"]],
      quadratic_alarm_date = if (is.na(monitored$alarm[["shao"]]))
        NA_character_ else as.character(split$monitor_dates[monitored$alarm[["shao"]]]),
      quadratic_maximum_ratio = max(
        monitored$path$shao / monitored$critical_values[["shao"]]
      ),
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replication = rep_index,
      success = FALSE,
      K_selected = NA_integer_,
      fve_at_3 = NA_real_, fve_at_4 = NA_real_,
      pc3_share = NA_real_, pc4_share = NA_real_, pc5_share = NA_real_,
      gap_3_4 = NA_real_, gap_4_5 = NA_real_,
      retained_subspace_max_angle_degrees = NA_real_,
      adjusted_range_alarm_index = NA_integer_,
      adjusted_range_alarm_date = NA_character_,
      adjusted_range_maximum_ratio = NA_real_,
      quadratic_alarm_index = NA_integer_,
      quadratic_alarm_date = NA_character_,
      quadratic_maximum_ratio = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

message(
  "Running ", args$nrep, " fit-block bootstrap replications with block length ",
  block_length, " on ", args$workers, " worker(s)."
)
bootstrap_rows <- if (args$workers > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(
    seq_len(args$nrep), bootstrap_one,
    mc.cores = args$workers, mc.preschedule = TRUE
  )
} else lapply(seq_len(args$nrep), bootstrap_one)
bootstrap_table <- do.call(rbind, bootstrap_rows)

data_dir <- file.path(paper_dir, "data")
utils::write.csv(
  rotation_table,
  file.path(data_dir, "spx_pc34_rotation_sensitivity.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  bootstrap_table,
  file.path(data_dir, "spx_fit_block_bootstrap_stability.csv"),
  row.names = FALSE, na = ""
)

successful <- bootstrap_table[bootstrap_table$success, , drop = FALSE]
alarm_summary <- function(index) {
  observed <- index[!is.na(index)]
  c(
    alarm_rate = length(observed) / length(index),
    alarm_index_median = if (length(observed)) stats::median(observed) else NA_real_,
    alarm_index_q10 = if (length(observed)) stats::quantile(observed, 0.10, names = FALSE) else NA_real_,
    alarm_index_q90 = if (length(observed)) stats::quantile(observed, 0.90, names = FALSE) else NA_real_
  )
}
range_boot <- alarm_summary(successful$adjusted_range_alarm_index)
quadratic_boot <- alarm_summary(successful$quadratic_alarm_index)
dimension_counts <- table(factor(successful$K_selected, levels = 1:5))
rotation_range_alarm <- rotation_table$adjusted_range_alarm_index
rotation_quadratic_alarm <- rotation_table$quadratic_alarm_index
summary_table <- data.frame(
  item = c(
    "rotation_grid_size", "rotation_adjusted_range_alarm_fraction",
    "rotation_adjusted_range_earliest_alarm", "rotation_adjusted_range_latest_alarm",
    "rotation_quadratic_alarm_fraction", "rotation_quadratic_path_max_abs_change",
    "bootstrap_replications", "bootstrap_successes", "bootstrap_block_length",
    paste0("bootstrap_K", 1:5, "_count"),
    "bootstrap_adjusted_range_alarm_rate", "bootstrap_adjusted_range_alarm_index_median",
    "bootstrap_adjusted_range_alarm_index_q10", "bootstrap_adjusted_range_alarm_index_q90",
    "bootstrap_quadratic_alarm_rate", "bootstrap_quadratic_alarm_index_median",
    "bootstrap_quadratic_alarm_index_q10", "bootstrap_quadratic_alarm_index_q90",
    "bootstrap_retained_subspace_angle_median", "bootstrap_retained_subspace_angle_q90"
  ),
  value = c(
    nrow(rotation_table), mean(!is.na(rotation_range_alarm)),
    if (all(is.na(rotation_range_alarm))) NA else min(rotation_range_alarm, na.rm = TRUE),
    if (all(is.na(rotation_range_alarm))) NA else max(rotation_range_alarm, na.rm = TRUE),
    mean(!is.na(rotation_quadratic_alarm)),
    max(abs(rotation_table$quadratic_maximum_ratio -
      rotation_table$quadratic_maximum_ratio[1L])),
    args$nrep, nrow(successful), block_length, as.numeric(dimension_counts),
    unname(range_boot), unname(quadratic_boot),
    stats::median(successful$retained_subspace_max_angle_degrees, na.rm = TRUE),
    stats::quantile(
      successful$retained_subspace_max_angle_degrees, 0.90,
      names = FALSE, na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary_table,
  file.path(data_dir, "spx_eigengap_stability_summary.csv"),
  row.names = FALSE
)

message("Wrote eigengap stability outputs to ", data_dir)
