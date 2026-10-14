#!/usr/bin/env Rscript

# Finite-sample audit for the substantive eigengap requirement.  The second
# and third eigenvalues are deliberately close while K=2 is fixed.  Shifts in
# PC2 are retained; shifts in PC3 lie just outside the target eigenspace.

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
    seed = 20260903L
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

one_replicate <- function(rep_id, target_pc, delta, critical, seed) {
  set.seed(seed + rep_id)
  d <- 8L
  n_fit <- 500L
  m <- 1000L
  h <- 500L
  break_after <- 100L
  lambda <- c(4.5, 2.60, 2.50, 0.45, 0.28, 0.18, 0.12, 0.07)
  rotation <- dc_functional_basis_rotation(d)
  latent <- matrix(stats::rnorm((n_fit + m + h) * d), ncol = d)
  latent <- sweep(latent, 2L, sqrt(lambda), FUN = "*")
  x <- latent %*% rotation
  if (target_pc > 0L) {
    changed <- n_fit + m + break_after + seq_len(h - break_after)
    shift <- delta * sqrt(lambda[target_pc]) * rotation[target_pc, ]
    x[changed, ] <- x[changed, , drop = FALSE] +
      matrix(shift, nrow = length(changed), ncol = d, byrow = TRUE)
  }
  fit_index <- seq_len(n_fit)
  cal_index <- n_fit + seq_len(m)
  mon_index <- n_fit + m + seq_len(h)
  result <- tryCatch({
    fit <- dc_fit_score_comparison(
      x[fit_index, , drop = FALSE], x[cal_index, , drop = FALSE],
      K = 2L, fve = 0.90, K_max = 5L, prewhiten_mode = "demix",
      ldl_ridge = 1e-8, sn_ridge = 0, feature_var_order = 0L
    )
    monitor <- dc_monitor_score_comparison(
      fit, x[mon_index, , drop = FALSE], alpha = 0.05, gamma = 0,
      critical_values = critical
    )
    estimated_basis <- fit$adjusted_range$fpca$basis_weighted[, 1:2, drop = FALSE]
    true_pc2 <- rotation[2L, ]
    true_pc3 <- rotation[3L, ]
    data.frame(
      replicate = rep_id, target_pc = target_pc, delta = delta,
      pc2_retained_share = sum(crossprod(estimated_basis, true_pc2)^2),
      pc3_retained_share = sum(crossprod(estimated_basis, true_pc3)^2),
      adjusted_range_alarm = first_crossing(
        monitor$path$adjusted_range, critical[["adjusted_range"]]
      ),
      full_matrix_alarm = first_crossing(
        monitor$path$shao, critical[["shao"]]
      ),
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      replicate = rep_id, target_pc = target_pc, delta = delta,
      pc2_retained_share = NA_real_, pc3_retained_share = NA_real_,
      adjusted_range_alarm = NA_integer_, full_matrix_alarm = NA_integer_,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  result
}

summarize_cell <- function(raw) {
  null_cell <- raw$target_pc[1L] == 0L
  range_pre <- !is.na(raw$adjusted_range_alarm) & raw$adjusted_range_alarm <= 100L
  full_pre <- !is.na(raw$full_matrix_alarm) & raw$full_matrix_alarm <= 100L
  range_post <- !is.na(raw$adjusted_range_alarm) & raw$adjusted_range_alarm > 100L
  full_post <- !is.na(raw$full_matrix_alarm) & raw$full_matrix_alarm > 100L
  data.frame(
    target = if (null_cell) "null" else paste0("PC", raw$target_pc[1L]),
    delta_standardized = raw$delta[1L],
    nrep = nrow(raw),
    adjusted_range_rejection = if (null_cell) mean(!is.na(raw$adjusted_range_alarm)) else mean(range_post),
    full_matrix_rejection = if (null_cell) mean(!is.na(raw$full_matrix_alarm)) else mean(full_post),
    adjusted_range_prechange_alarm = mean(range_pre),
    full_matrix_prechange_alarm = mean(full_pre),
    median_pc2_retained_share = stats::median(raw$pc2_retained_share),
    median_pc3_retained_share = stats::median(raw$pc3_retained_share),
    q10_pc3_retained_share = stats::quantile(raw$pc3_retained_share, 0.10),
    q90_pc3_retained_share = stats::quantile(raw$pc3_retained_share, 0.90),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
critical_table <- utils::read.csv(file.path(
  script_dir, "outputs", "critical_value_audit",
  "critical_value_grid_and_monte_carlo_audit.csv"
), stringsAsFactors = FALSE)
critical_rows <- critical_table[
  critical_table$configuration == "exact1000_500_T0.5" & critical_table$K == 2L,
  , drop = FALSE
]
critical <- c(
  adjusted_range = critical_rows$critical_value[critical_rows$method_id == "adjusted_range"],
  shao = critical_rows$critical_value[critical_rows$method_id == "full_matrix"]
)

cells <- data.frame(target_pc = c(0L, 2L, 3L), delta = c(0, 0.30, 0.30))
output_dir <- file.path(script_dir, "outputs", "near_eigengap_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()
for (cell_id in seq_len(nrow(cells))) {
  target_pc <- cells$target_pc[cell_id]
  delta <- cells$delta[cell_id]
  cache <- file.path(output_dir, sprintf("pc%d_delta_%0.2f_%drep.rds", target_pc, delta, args$nrep))
  if (file.exists(cache)) {
    raw <- readRDS(cache)
  } else {
    pieces <- parallel::mclapply(
      seq_len(args$nrep), one_replicate,
      target_pc = target_pc, delta = delta, critical = critical,
      seed = args$seed + cell_id * 1000000L,
      mc.cores = args$workers, mc.preschedule = TRUE
    )
    raw <- do.call(rbind, pieces)
    saveRDS(raw, cache, compress = "xz")
  }
  valid <- raw[is.na(raw$error), , drop = FALSE]
  summaries[[cell_id]] <- summarize_cell(valid)
}
summary_table <- do.call(rbind, summaries)
utils::write.csv(
  summary_table,
  file.path(output_dir, "near_eigengap_null_power_5000rep.csv"),
  row.names = FALSE
)
print(summary_table)
