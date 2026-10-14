#!/usr/bin/env Rscript

# Power audit at a nonregular 90% FVE boundary.  The post-change shift is in
# population PC3, so it is visible when K-hat=3 and omitted when K-hat=2.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else normalizePath("simulations/componentwise_range")
project_dir <- normalizePath(file.path(script_dir, "../.."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))
source(file.path(project_dir, "simulations", "formal", "R", "null_dgp.R"))

parse_args <- function(args) {
  out <- list(
    nrep = 5000L,
    workers = max(1L, min(10L, parallel::detectCores(logical = TRUE) - 1L)),
    seed = 20260902L
  )
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out)) stop("Unknown argument: ", arg)
    out[[pair[1L]]] <- as.integer(pair[2L])
  }
  out
}

first_crossing <- function(path, critical) {
  hit <- which(path > critical)
  if (length(hit)) hit[1L] else NA_integer_
}

one_replicate <- function(rep_id, delta, critical, seed) {
  set.seed(seed + rep_id)
  d <- 8L
  n_fit <- 500L
  m <- 1000L
  h <- 500L
  break_after <- 100L
  total_n <- n_fit + m + h
  x <- dc_simulate_null_scores(
    n = total_n, d = d, design = "component_ar_lrv",
    eigen_profile = "fve_boundary"
  )
  lambda <- dc_score_eigenvalues(d, profile = "fve_boundary")
  rotation <- dc_functional_basis_rotation(d)
  shift <- delta * sqrt(lambda[3L]) * rotation[3L, ]
  changed <- n_fit + m + break_after + seq_len(h - break_after)
  x[changed, ] <- x[changed, , drop = FALSE] +
    matrix(shift, nrow = length(changed), ncol = d, byrow = TRUE)
  fit_index <- seq_len(n_fit)
  cal_index <- n_fit + seq_len(m)
  mon_index <- n_fit + m + seq_len(h)
  result <- tryCatch({
    fit <- dc_fit_score_comparison(
      x[fit_index, , drop = FALSE], x[cal_index, , drop = FALSE],
      K = NULL, fve = 0.90, K_max = 5L, prewhiten_mode = "demix",
      ldl_ridge = 1e-8, sn_ridge = 0, feature_var_order = 0L
    )
    key <- paste0("K", fit$settings$K)
    cv <- c(
      adjusted_range = critical$adjusted_range[[key]],
      shao = critical$full_matrix[[key]]
    )
    monitor <- dc_monitor_score_comparison(
      fit, x[mon_index, , drop = FALSE], alpha = 0.05, gamma = 0,
      critical_values = cv
    )
    data.frame(
      replicate = rep_id, delta = delta, K_hat = fit$settings$K,
      adjusted_range_alarm = first_crossing(monitor$path$adjusted_range, cv[["adjusted_range"]]),
      full_matrix_alarm = first_crossing(monitor$path$shao, cv[["shao"]]),
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = rep_id, delta = delta, K_hat = NA_integer_,
      adjusted_range_alarm = NA_integer_, full_matrix_alarm = NA_integer_,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  result
}

summarize_group <- function(raw, K_group) {
  z <- if (K_group == "all") raw else raw[raw$K_hat == as.integer(K_group), , drop = FALSE]
  range_pre <- !is.na(z$adjusted_range_alarm) & z$adjusted_range_alarm <= 100L
  full_pre <- !is.na(z$full_matrix_alarm) & z$full_matrix_alarm <= 100L
  range_post <- !is.na(z$adjusted_range_alarm) & z$adjusted_range_alarm > 100L
  full_post <- !is.na(z$full_matrix_alarm) & z$full_matrix_alarm > 100L
  data.frame(
    delta_standardized_pc3 = z$delta[1L],
    K_group = K_group,
    n = nrow(z),
    adjusted_range_prechange_alarm = mean(range_pre),
    full_matrix_prechange_alarm = mean(full_pre),
    adjusted_range_postchange_detection = mean(range_post),
    full_matrix_postchange_detection = mean(full_post),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
critical_table <- utils::read.csv(file.path(
  script_dir, "outputs", "critical_value_audit",
  "critical_value_grid_and_monte_carlo_audit.csv"
), stringsAsFactors = FALSE)
critical_rows <- critical_table[
  critical_table$configuration == "exact1000_500_T0.5",
  , drop = FALSE
]
critical <- list(
  adjusted_range = stats::setNames(
    as.list(critical_rows$critical_value[critical_rows$method_id == "adjusted_range"]),
    paste0("K", critical_rows$K[critical_rows$method_id == "adjusted_range"])
  ),
  full_matrix = stats::setNames(
    as.list(critical_rows$critical_value[critical_rows$method_id == "full_matrix"]),
    paste0("K", critical_rows$K[critical_rows$method_id == "full_matrix"])
  )
)

output_dir <- file.path(script_dir, "outputs", "boundary_fve_pc3_power")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()
for (cell_id in seq_along(c(0.30, 0.60))) {
  delta <- c(0.30, 0.60)[cell_id]
  cache <- file.path(output_dir, sprintf("pc3_delta_%0.2f_%drep.rds", delta, args$nrep))
  if (file.exists(cache)) {
    raw <- readRDS(cache)
  } else {
    pieces <- parallel::mclapply(
      seq_len(args$nrep), one_replicate, delta = delta, critical = critical,
      seed = args$seed + cell_id * 1000000L,
      mc.cores = args$workers, mc.preschedule = TRUE
    )
    raw <- do.call(rbind, pieces)
    saveRDS(raw, cache, compress = "xz")
  }
  valid <- raw[is.na(raw$error), , drop = FALSE]
  for (group in c("all", "2", "3")) {
    summary_rows[[length(summary_rows) + 1L]] <- summarize_group(valid, group)
  }
}
summary_table <- do.call(rbind, summary_rows)
utils::write.csv(
  summary_table,
  file.path(output_dir, "boundary_fve_pc3_power_5000rep.csv"),
  row.names = FALSE
)
print(summary_table)
