#!/usr/bin/env Rscript

# Matched-grid K=5 critical values needed when the empirical fit-block
# bootstrap selects five retained components.  The primary specification uses
# K=4, so this file extends only the empirical 997/251 grid without rerunning
# the other critical-value configurations.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else normalizePath("reproducibility/code/simulations/componentwise_range")
project_dir <- normalizePath(file.path(script_dir, "../.."))
paper_dir <- normalizePath(file.path(project_dir, ".."))

parse_args <- function(args) {
  out <- list(
    nrep = 250000L,
    workers = max(1L, min(10L, parallel::detectCores(logical = TRUE) - 1L)),
    chunk_size = 500L,
    seed = 20270902L
  )
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out)) stop("Unknown argument: ", arg)
    out[[pair[1L]]] <- as.integer(pair[2L])
  }
  out
}

one_path <- function(train_grid = 997L, monitor_grid = 251L, K = 5L) {
  increments <- matrix(
    stats::rnorm((train_grid + monitor_grid) * K),
    ncol = K
  ) / sqrt(train_grid)
  brownian <- apply(increments, 2L, cumsum)
  terminal <- brownian[train_grid, ]
  r <- seq_len(train_grid) / train_grid
  s <- seq_len(monitor_grid) / train_grid
  bridge <- brownian[seq_len(train_grid), , drop = FALSE] - outer(r, terminal)
  monitoring <- brownian[train_grid + seq_len(monitor_grid), , drop = FALSE] -
    outer(1 + s, terminal)
  ranges <- apply(rbind(rep(0, K), bridge), 2L, max) -
    apply(rbind(rep(0, K), bridge), 2L, min)
  denominator <- (1 + s)^2
  range_max <- max(rowSums(sweep(monitoring^2, 2L, ranges^2, FUN = "/")) /
    denominator)
  V_inverse <- chol2inv(chol(crossprod(bridge) / train_grid))
  quadratic_max <- max(rowSums((monitoring %*% V_inverse) * monitoring) /
    denominator)
  c(adjusted_range = range_max, full_matrix = quadratic_max)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
output_dir <- file.path(script_dir, "outputs", "critical_value_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cache_file <- file.path(
  output_dir,
  sprintf("critical_audit_exact997_251_empirical_n%d_m997_h251_K5_only.rds", args$nrep)
)

if (file.exists(cache_file)) {
  maxima <- readRDS(cache_file)
} else {
  starts <- seq.int(1L, args$nrep, by = args$chunk_size)
  chunks <- parallel::mclapply(
    seq_along(starts),
    function(chunk_id) {
      first <- starts[chunk_id]
      last <- min(args$nrep, first + args$chunk_size - 1L)
      set.seed(args$seed + chunk_id)
      ans <- matrix(NA_real_, nrow = last - first + 1L, ncol = 2L)
      for (i in seq_len(nrow(ans))) ans[i, ] <- one_path()
      if (chunk_id %% 20L == 0L) message("completed ", last, "/", args$nrep)
      ans
    },
    mc.cores = args$workers,
    mc.preschedule = FALSE
  )
  maxima <- do.call(rbind, chunks)
  colnames(maxima) <- c("adjusted_range", "full_matrix")
  saveRDS(maxima, cache_file, compress = "xz")
}

lower_rank <- max(1L, stats::qbinom(0.025, args$nrep, 0.95))
upper_rank <- min(args$nrep, stats::qbinom(0.975, args$nrep, 0.95) + 1L)
new_rows <- do.call(rbind, lapply(seq_len(ncol(maxima)), function(j) {
  values <- sort(maxima[, j])
  data.frame(
    configuration = "exact997_251_empirical",
    method_id = colnames(maxima)[j],
    K = 5L,
    train_grid = 997L,
    monitor_grid = 251L,
    horizon_T = 251 / 997,
    nrep = args$nrep,
    critical_value = as.numeric(stats::quantile(values, 0.95, type = 8, names = FALSE)),
    quantile_ci95_low = values[lower_rank],
    quantile_ci95_high = values[upper_rank],
    stringsAsFactors = FALSE
  )
}))

audit_file <- file.path(paper_dir, "data", "critical_value_grid_and_monte_carlo_audit.csv")
checks_file <- file.path(paper_dir, "data", "critical_value_grid_and_monte_carlo_checks.csv")
existing <- utils::read.csv(audit_file, stringsAsFactors = FALSE)
existing <- existing[!(existing$configuration == "exact997_251_empirical" & existing$K == 5L), ]
combined <- rbind(existing, new_rows)
configuration_order <- unique(existing$configuration)
combined$configuration <- factor(combined$configuration, levels = configuration_order)
combined <- combined[order(combined$configuration, combined$method_id, combined$K), ]
combined$configuration <- as.character(combined$configuration)
utils::write.csv(combined, audit_file, row.names = FALSE)
utils::write.csv(combined, checks_file, row.names = FALSE)
print(new_rows)
