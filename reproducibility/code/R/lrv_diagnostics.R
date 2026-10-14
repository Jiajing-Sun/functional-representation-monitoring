# Bartlett long-run covariance estimate for empirical correlation diagnostics.
det_bartlett_lrv <- function(features, bandwidth = 10L) {
  features <- dc_assert_matrix(features, "features", min_rows = 20L)
  n <- nrow(features)
  bandwidth <- min(as.integer(bandwidth), n - 1L)
  center <- colMeans(features)
  centered <- sweep(features, 2L, center, FUN = "-")
  lrv <- crossprod(centered) / n
  for (lag in seq_len(bandwidth)) {
    weight <- 1 - lag / (bandwidth + 1)
    gamma_lag <- crossprod(
      centered[seq_len(n - lag), , drop = FALSE],
      centered[lag + seq_len(n - lag), , drop = FALSE]
    ) / n
    lrv <- lrv + weight * (gamma_lag + t(gamma_lag))
  }
  inverse <- dc_inverse_spd(lrv, ridge = 1e-8)
  list(center = center, lrv = lrv, inverse = inverse$inverse,
       ridge_added = inverse$ridge_added,
       condition_number = inverse$condition_number,
       n_reference = n, dimension = ncol(features), bandwidth = bandwidth)
}
