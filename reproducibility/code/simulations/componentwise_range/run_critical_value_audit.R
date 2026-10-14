#!/usr/bin/env Rscript

# High-precision Brownian critical-value audit for the adjusted-range and
# Shao-type quadratic self-normalized monitoring statistics. Each Brownian
# path supplies all K=1,...,Kmax values, and the final configurations use the
# exact discrete normalization and monitoring grid lengths used in the paper.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else normalizePath("simulations/componentwise_range")
project_dir <- normalizePath(file.path(script_dir, "../.."))
source(file.path(project_dir, "R", "dual_channel_monitor.R"))

parse_args <- function(args) {
  out <- list(
    nrep_final = 250000L,
    nrep_grid = 50000L,
    workers = max(1L, min(10L, parallel::detectCores(logical = TRUE) - 1L)),
    chunk_size = 250L,
    seed = 20260831L
  )
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L || !pair[1L] %in% names(out)) stop("Unknown argument: ", arg)
    out[[pair[1L]]] <- as.integer(pair[2L])
  }
  out
}

one_path <- function(train_grid, monitor_grid, Kmax) {
  total <- train_grid + monitor_grid
  increments <- matrix(rnorm(total * Kmax), nrow = total, ncol = Kmax) /
    sqrt(train_grid)
  brownian <- apply(increments, 2L, cumsum)
  if (is.null(dim(brownian))) brownian <- matrix(brownian, ncol = 1L)
  terminal <- brownian[train_grid, ]
  r <- seq_len(train_grid) / train_grid
  s <- seq_len(monitor_grid) / train_grid
  bridge <- brownian[seq_len(train_grid), , drop = FALSE] - outer(r, terminal)
  bridge_with_zero <- rbind(rep(0, Kmax), bridge)
  monitoring <- brownian[train_grid + seq_len(monitor_grid), , drop = FALSE] -
    outer(1 + s, terminal)
  denominator <- (1 + s)^2
  out <- numeric(2L * Kmax)
  for (K in seq_len(Kmax)) {
    bridge_K <- bridge[, seq_len(K), drop = FALSE]
    monitor_K <- monitoring[, seq_len(K), drop = FALSE]
    ranges <- apply(bridge_with_zero[, seq_len(K), drop = FALSE], 2L, max) -
      apply(bridge_with_zero[, seq_len(K), drop = FALSE], 2L, min)
    range_q <- rowSums(sweep(monitor_K^2, 2L, ranges^2, FUN = "/"))
    V <- crossprod(bridge_K) / train_grid
    V_inv <- chol2inv(chol(V))
    full_q <- rowSums((monitor_K %*% V_inv) * monitor_K)
    out[K] <- max(range_q / denominator)
    out[Kmax + K] <- max(full_q / denominator)
  }
  out
}

run_configuration <- function(id, train_grid, monitor_grid, nrep, Kmax,
                              seed, workers, chunk_size, output_dir) {
  cache_file <- file.path(
    output_dir,
    sprintf("critical_audit_%s_n%d_m%d_h%d_K%d.rds", id, nrep, train_grid,
            monitor_grid, Kmax)
  )
  if (file.exists(cache_file)) return(readRDS(cache_file))
  starts <- seq.int(1L, nrep, by = chunk_size)
  chunks <- parallel::mclapply(
    seq_along(starts),
    function(chunk_id) {
      first <- starts[chunk_id]
      last <- min(nrep, first + chunk_size - 1L)
      set.seed(seed + chunk_id)
      ans <- matrix(NA_real_, nrow = last - first + 1L, ncol = 2L * Kmax)
      for (i in seq_len(nrow(ans))) {
        ans[i, ] <- one_path(train_grid, monitor_grid, Kmax)
      }
      if (chunk_id %% 20L == 0L) {
        message(id, ": completed ", last, "/", nrep)
      }
      ans
    },
    mc.cores = workers,
    mc.preschedule = FALSE
  )
  maxima <- do.call(rbind, chunks)
  colnames(maxima) <- c(
    paste0("adjusted_range_K", seq_len(Kmax)),
    paste0("full_matrix_K", seq_len(Kmax))
  )
  result <- list(
    id = id,
    settings = list(
      train_grid = train_grid, monitor_grid = monitor_grid, nrep = nrep,
      Kmax = Kmax, seed = seed
    ),
    maxima = maxima
  )
  saveRDS(result, cache_file, compress = "xz")
  result
}

summarize_configuration <- function(result) {
  nrep <- nrow(result$maxima)
  lower_rank <- max(1L, stats::qbinom(0.025, nrep, 0.95))
  upper_rank <- min(nrep, stats::qbinom(0.975, nrep, 0.95) + 1L)
  rows <- lapply(seq_len(ncol(result$maxima)), function(j) {
    values <- sort(result$maxima[, j])
    parts <- strsplit(colnames(result$maxima)[j], "_K", fixed = TRUE)[[1L]]
    data.frame(
      configuration = result$id,
      method_id = parts[1L],
      K = as.integer(parts[2L]),
      train_grid = result$settings$train_grid,
      monitor_grid = result$settings$monitor_grid,
      horizon_T = result$settings$monitor_grid / result$settings$train_grid,
      nrep = nrep,
      critical_value = as.numeric(stats::quantile(values, 0.95, type = 8, names = FALSE)),
      quantile_ci95_low = values[lower_rank],
      quantile_ci95_high = values[upper_rank],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
output_dir <- file.path(script_dir, "outputs", "critical_value_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

configurations <- list(
  list("grid400_T0.5", 400L, 200L, args$nrep_grid, 5L),
  list("grid800_T0.5", 800L, 400L, args$nrep_grid, 5L),
  list("exact1000_500_T0.5", 1000L, 500L, args$nrep_final, 5L),
  list("exact997_251_empirical", 997L, 251L, args$nrep_final, 4L),
  list("exact748_249_placebo", 748L, 249L, args$nrep_grid, 4L)
)

results <- vector("list", length(configurations))
for (i in seq_along(configurations)) {
  cfg <- configurations[[i]]
  message("Starting ", cfg[[1L]])
  results[[i]] <- run_configuration(
    id = cfg[[1L]], train_grid = cfg[[2L]], monitor_grid = cfg[[3L]],
    nrep = cfg[[4L]], Kmax = cfg[[5L]], seed = args$seed + i * 1000000L,
    workers = args$workers, chunk_size = args$chunk_size,
    output_dir = output_dir
  )
}

summary_table <- do.call(rbind, lapply(results, summarize_configuration))
utils::write.csv(
  summary_table,
  file.path(output_dir, "critical_value_grid_and_monte_carlo_audit.csv"),
  row.names = FALSE
)
print(summary_table)
