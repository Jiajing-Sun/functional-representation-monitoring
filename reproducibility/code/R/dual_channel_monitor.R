# Paper-facing monitoring of the mean of retained FPCA scores.
#
# The manuscript-facing entry points fit one frozen FPCA representation, apply
# lag-zero LDL triangular demixing, and compare componentwise adjusted-range
# normalization with Shao's quadratic self-normalization benchmark. The
# paper fixes the unweighted boundary (gamma = 0). General boundary arguments
# and the older dual-channel, orthogonal-audit, and affine-range helpers remain
# in this file only for backward compatibility; the paper-facing scripts do not
# invoke those extensions.

dc_assert_matrix <- function(x, name, min_rows = 2L, min_cols = 1L) {
  if (!is.numeric(x) || is.null(dim(x)) || length(dim(x)) != 2L) {
    stop(name, " must be a numeric matrix.")
  }
  x <- as.matrix(x)
  if (nrow(x) < min_rows || ncol(x) < min_cols) {
    stop(name, " must have at least ", min_rows, " rows and ", min_cols, " columns.")
  }
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x
}

dc_symmetrize <- function(a) (a + t(a)) / 2

dc_relative_ridge <- function(a, ridge) {
  scale <- mean(abs(diag(a)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  as.numeric(ridge) * scale
}

dc_inverse_spd <- function(a, ridge = 1e-8, max_tries = 8L) {
  a <- dc_assert_matrix(a, "a", min_rows = 1L)
  if (nrow(a) != ncol(a)) stop("a must be square.")
  a <- dc_symmetrize(a)
  base_ridge <- dc_relative_ridge(a, ridge)

  for (attempt in 0:max_tries) {
    added <- if (attempt == 0L) 0 else base_ridge * (10^(attempt - 1L))
    candidate <- a + added * diag(nrow(a))
    chol_factor <- tryCatch(chol(candidate), error = function(e) NULL)
    if (!is.null(chol_factor)) {
      return(list(
        inverse = chol2inv(chol_factor),
        matrix = candidate,
        ridge_added = added,
        condition_number = kappa(candidate, exact = TRUE)
      ))
    }
  }
  stop("Unable to regularize matrix to positive definiteness.")
}

# For a positive-definite covariance Sigma, R <- chol(Sigma) satisfies
# Sigma = R'R. Setting C = R' diag(1 / diag(R)) and D = diag(diag(R)^2)
# gives the unique unit-lower-triangular LDL representation Sigma = C D C'.
dc_ldl_factor <- function(sigma, ridge = 1e-8, max_tries = 8L) {
  sigma <- dc_assert_matrix(sigma, "sigma", min_rows = 1L)
  if (nrow(sigma) != ncol(sigma)) stop("sigma must be square.")
  sigma <- dc_symmetrize(sigma)
  base_ridge <- dc_relative_ridge(sigma, ridge)

  for (attempt in 0:max_tries) {
    added <- if (attempt == 0L) 0 else base_ridge * (10^(attempt - 1L))
    sigma_reg <- sigma + added * diag(nrow(sigma))
    r <- tryCatch(chol(sigma_reg), error = function(e) NULL)
    if (!is.null(r)) {
      diagonal <- diag(r)
      cmat <- t(r) %*% diag(1 / diagonal, nrow = length(diagonal))
      dmat <- diag(diagonal^2, nrow = length(diagonal))
      return(list(
        C = cmat,
        D = dmat,
        sigma_regularized = sigma_reg,
        ridge_added = added,
        relative_reconstruction_error =
          norm(sigma_reg - cmat %*% dmat %*% t(cmat), type = "F") /
          max(norm(sigma_reg, type = "F"), .Machine$double.eps)
      ))
    }
  }
  stop("Unable to compute an LDL factorization after ridge regularization.")
}

# Hong et al. use C^{-1}X_t: this removes triangular cross-coordinate mixing
# but deliberately retains own-coordinate scale and serial dependence. The
# optional "whiten" mode also applies D^{-1/2}; it is not the paper default.
dc_fit_ldl_prewhitener <- function(z, ridge = 1e-8,
                                   mode = c("demix", "whiten", "none")) {
  z <- dc_assert_matrix(z, "z")
  mode <- match.arg(mode)
  center <- colMeans(z)
  centered <- sweep(z, 2L, center, FUN = "-")

  if (mode == "none") {
    transform <- diag(ncol(z))
    return(list(
      center = center,
      transform = transform,
      C = diag(ncol(z)),
      D = diag(ncol(z)),
      mode = mode,
      covariance = stats::cov(centered),
      ridge_added = 0
    ))
  }

  sigma <- stats::cov(centered)
  factor <- dc_ldl_factor(sigma, ridge = ridge)
  demix <- solve(factor$C)
  transform <- if (mode == "demix") {
    demix
  } else {
    diag(1 / sqrt(diag(factor$D)), nrow = ncol(z)) %*% demix
  }

  list(
    center = center,
    transform = transform,
    C = factor$C,
    D = factor$D,
    mode = mode,
    covariance = sigma,
    ridge_added = factor$ridge_added,
    relative_reconstruction_error = factor$relative_reconstruction_error
  )
}

dc_apply_ldl_prewhitener <- function(prewhitener, z) {
  z <- dc_assert_matrix(z, "z", min_rows = 1L)
  if (ncol(z) != length(prewhitener$center)) {
    stop("z has a different number of columns from the fitted prewhitener.")
  }
  sweep(z, 2L, prewhitener$center, FUN = "-") %*% t(prewhitener$transform)
}

dc_grid_weights <- function(n_grid, grid = NULL, weights = NULL) {
  if (!is.null(weights)) {
    if (length(weights) != n_grid || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("weights must contain one positive finite value per grid point.")
    }
    return(as.numeric(weights) / sum(weights))
  }
  if (is.null(grid)) return(rep(1 / n_grid, n_grid))
  if (length(grid) != n_grid || any(!is.finite(grid)) || is.unsorted(grid, strictly = TRUE)) {
    stop("grid must be a strictly increasing finite vector with one value per column.")
  }
  if (n_grid == 1L) return(1)
  delta <- diff(grid)
  w <- c(delta[1L] / 2, (head(delta, -1L) + tail(delta, -1L)) / 2, tail(delta, 1L) / 2)
  w / sum(w)
}

dc_fit_fpca <- function(x_fit, K = NULL, L = 1L, fve = 0.90, K_max = NULL,
                        grid = NULL, weights = NULL) {
  x_fit <- dc_assert_matrix(x_fit, "x_fit", min_rows = 3L, min_cols = 2L)
  L <- as.integer(L)
  if (L < 1L) stop("L must be at least one.")
  if (!is.numeric(fve) || length(fve) != 1L || fve <= 0 || fve >= 1) {
    stop("fve must lie strictly between zero and one.")
  }

  w <- dc_grid_weights(ncol(x_fit), grid = grid, weights = weights)
  sqrt_w <- sqrt(w)
  mean_function <- colMeans(x_fit)
  centered <- sweep(x_fit, 2L, mean_function, FUN = "-")
  weighted <- sweep(centered, 2L, sqrt_w, FUN = "*")
  rank_max <- min(nrow(weighted) - 1L, ncol(weighted))
  sv <- svd(weighted, nu = 0L, nv = rank_max)
  eigenvalues <- sv$d^2 / max(nrow(weighted) - 1L, 1L)
  fve_path <- cumsum(eigenvalues) / sum(eigenvalues)

  selection_cap <- if (is.null(K_max)) rank_max - L else as.integer(K_max)
  if (length(selection_cap) != 1L || !is.finite(selection_cap) || selection_cap < 1L) {
    stop("K_max must be one positive integer.")
  }
  if (selection_cap + L > rank_max) {
    stop("K_max + L exceeds the estimable FPCA rank.")
  }

  if (is.null(K)) {
    eligible <- which(seq_along(fve_path) <= selection_cap & fve_path >= fve)
    if (!length(eligible)) {
      stop("The target FVE is not reached by the prespecified K_max.")
    }
    K <- eligible[1L]
  } else {
    K <- as.integer(K)
  }
  if (K < 1L) stop("K must be at least one.")
  if (K > selection_cap) stop("K exceeds the prespecified K_max.")
  if (K + L > rank_max) {
    stop("K + L exceeds the estimable FPCA rank; reduce K or L, or enlarge x_fit.")
  }

  basis_weighted <- sv$v[, seq_len(K + L), drop = FALSE]
  eigenfunctions <- sweep(basis_weighted, 1L, sqrt_w, FUN = "/")
  colnames(eigenfunctions) <- paste0("phi", seq_len(K + L))

  list(
    mean_function = mean_function,
    weights = w,
    sqrt_weights = sqrt_w,
    basis_weighted = basis_weighted,
    eigenfunctions = eigenfunctions,
    eigenvalues = eigenvalues,
    fve_path = fve_path,
    K = K,
    K_max = selection_cap,
    L = L,
    fve_target = fve,
    n_fit = nrow(x_fit),
    n_grid = ncol(x_fit)
  )
}

dc_project_fpca <- function(fpca, x) {
  x <- dc_assert_matrix(x, "x", min_rows = 1L, min_cols = fpca$n_grid)
  if (ncol(x) != fpca$n_grid) stop("x does not use the fitted functional grid.")
  centered <- sweep(x, 2L, fpca$mean_function, FUN = "-")
  weighted <- sweep(centered, 2L, fpca$sqrt_weights, FUN = "*")
  scores <- weighted %*% fpca$basis_weighted
  retained <- scores[, seq_len(fpca$K), drop = FALSE]
  audit_index <- fpca$K + seq_len(fpca$L)
  audit <- scores[, audit_index, drop = FALSE]
  retained_audit_reconstruction <- scores %*% t(fpca$basis_weighted)
  far_tail_weighted <- weighted - retained_audit_reconstruction
  far_tail_energy <- rowSums(far_tail_weighted^2)
  total_reconstruction_loss <- rowSums(audit^2) + far_tail_energy

  colnames(retained) <- paste0("s", seq_len(fpca$K))
  colnames(audit) <- paste0("u", seq_len(fpca$L))
  list(
    retained = retained,
    audit = audit,
    far_tail_energy = far_tail_energy,
    total_reconstruction_loss = total_reconstruction_loss
  )
}

dc_svec_outer <- function(u) {
  u <- dc_assert_matrix(u, "u", min_rows = 1L)
  L <- ncol(u)
  out <- vector("list", L * (L + 1L) / 2L)
  names_out <- character(length(out))
  idx <- 1L
  for (j in seq_len(L)) {
    for (i in seq_len(j)) {
      multiplier <- if (i == j) 1 else sqrt(2)
      out[[idx]] <- multiplier * u[, i] * u[, j]
      names_out[idx] <- if (i == j) paste0("u", i, "_sq") else paste0("sqrt2_u", i, "u", j)
      idx <- idx + 1L
    }
  }
  ans <- do.call(cbind, out)
  colnames(ans) <- names_out
  ans
}

dc_cross_features <- function(s, u) {
  s <- dc_assert_matrix(s, "s", min_rows = 1L)
  u <- dc_assert_matrix(u, "u", min_rows = 1L)
  if (nrow(s) != nrow(u)) stop("s and u must have the same number of rows.")
  out <- matrix(NA_real_, nrow(s), ncol(s) * ncol(u))
  names_out <- character(ncol(out))
  idx <- 1L
  for (j in seq_len(ncol(u))) {
    for (i in seq_len(ncol(s))) {
      out[, idx] <- s[, i] * u[, j]
      names_out[idx] <- paste0("s", i, "u", j)
      idx <- idx + 1L
    }
  }
  colnames(out) <- names_out
  out
}

dc_build_features <- function(fpca_projection, prewhitener, K, L) {
  base <- cbind(fpca_projection$retained, fpca_projection$audit)
  base_pw <- dc_apply_ldl_prewhitener(prewhitener, base)
  s <- base_pw[, seq_len(K), drop = FALSE]
  u <- base_pw[, K + seq_len(L), drop = FALSE]
  colnames(s) <- paste0("S_score_", seq_len(K))
  colnames(u) <- paste0("R_audit_location_", seq_len(L))

  audit_covariance <- dc_svec_outer(u)
  colnames(audit_covariance) <- paste0("R_audit_cov_", colnames(audit_covariance))
  far_tail <- matrix(fpca_projection$far_tail_energy, ncol = 1L)
  colnames(far_tail) <- "R_far_tail_energy"
  leakage <- dc_cross_features(s, u)
  colnames(leakage) <- paste0("R_leakage_", colnames(leakage))

  list(
    channel_S = s,
    channel_R = cbind(u, audit_covariance, far_tail, leakage),
    diagnostics = data.frame(
      total_reconstruction_loss = fpca_projection$total_reconstruction_loss,
      far_tail_energy = fpca_projection$far_tail_energy
    )
  )
}

dc_fit_var_ldl_features <- function(features, order = 0L, ridge = 1e-8,
                                    ldl_mode = c("demix", "whiten")) {
  features <- dc_assert_matrix(features, "features", min_rows = 5L)
  order <- as.integer(order)
  ldl_mode <- match.arg(ldl_mode)
  n <- nrow(features)
  d <- ncol(features)
  if (order < 0L || order > 3L) stop("VAR order must be between zero and three.")
  if (order == 0L) {
    return(list(
      transformed = features,
      model = list(order = 0L, dimension = d, feature_names = colnames(features))
    ))
  }
  if (n - order <= 2L + order * d) {
    stop("Normalization block is too short for the requested feature VAR.")
  }

  response <- features[(order + 1L):n, , drop = FALSE]
  design <- matrix(1, nrow = n - order, ncol = 1L + order * d)
  column <- 2L
  for (lag in seq_len(order)) {
    design[, column:(column + d - 1L)] <- features[(order + 1L - lag):(n - lag), , drop = FALSE]
    column <- column + d
  }
  gram <- crossprod(design)
  penalty <- dc_relative_ridge(gram, ridge)
  coefficient <- solve(gram + penalty * diag(ncol(gram)), crossprod(design, response))
  residuals <- response - design %*% coefficient
  ldl <- dc_fit_ldl_prewhitener(residuals, ridge = ridge, mode = ldl_mode)
  transformed <- dc_apply_ldl_prewhitener(ldl, residuals)
  colnames(transformed) <- colnames(features)
  list(
    transformed = transformed,
    model = list(
      order = order,
      dimension = d,
      coefficient = coefficient,
      ldl = ldl,
      history = tail(features, order),
      feature_names = colnames(features),
      regression_ridge = penalty
    )
  )
}

dc_apply_var_ldl_features <- function(model, features) {
  features <- dc_assert_matrix(features, "features", min_rows = 1L)
  if (ncol(features) != model$dimension) stop("Feature dimension does not match the fitted VAR.")
  if (model$order == 0L) return(features)
  combined <- rbind(model$history, features)
  d <- model$dimension
  design <- matrix(1, nrow = nrow(features), ncol = 1L + model$order * d)
  column <- 2L
  for (lag in seq_len(model$order)) {
    row_index <- model$order + seq_len(nrow(features)) - lag
    design[, column:(column + d - 1L)] <- combined[row_index, , drop = FALSE]
    column <- column + d
  }
  residuals <- features - design %*% model$coefficient
  transformed <- dc_apply_ldl_prewhitener(model$ldl, residuals)
  colnames(transformed) <- model$feature_names
  transformed
}

dc_fit_huberizer <- function(features, threshold = Inf) {
  features <- dc_assert_matrix(features, "features", min_rows = 5L)
  threshold <- as.numeric(threshold)
  if (length(threshold) != 1L || (!is.infinite(threshold) && threshold <= 0)) {
    stop("Huber threshold must be positive or Inf.")
  }
  location <- apply(features, 2L, stats::median)
  robust_scale <- apply(features, 2L, stats::mad, constant = 1.4826)
  fallback <- apply(features, 2L, stats::sd)
  bad <- !is.finite(robust_scale) | robust_scale <= sqrt(.Machine$double.eps)
  robust_scale[bad] <- fallback[bad]
  robust_scale[!is.finite(robust_scale) | robust_scale <= sqrt(.Machine$double.eps)] <- 1
  model <- list(
    location = location,
    scale = robust_scale,
    threshold = threshold,
    dimension = ncol(features),
    feature_names = colnames(features)
  )
  list(transformed = dc_apply_huberizer(model, features), model = model)
}

dc_apply_huberizer <- function(model, features) {
  features <- dc_assert_matrix(features, "features", min_rows = 1L)
  if (ncol(features) != model$dimension) stop("Feature dimension does not match the Huberizer.")
  if (is.infinite(model$threshold)) return(features)
  standardized <- sweep(features, 2L, model$location, FUN = "-")
  standardized <- sweep(standardized, 2L, model$scale, FUN = "/")
  standardized[standardized > model$threshold] <- model$threshold
  standardized[standardized < -model$threshold] <- -model$threshold
  transformed <- sweep(standardized, 2L, model$scale, FUN = "*")
  transformed <- sweep(transformed, 2L, model$location, FUN = "+")
  colnames(transformed) <- model$feature_names
  transformed
}

dc_fit_full_matrix_sn <- function(features, ridge = 0) {
  features <- dc_assert_matrix(features, "features", min_rows = 3L)
  n <- nrow(features)
  d <- ncol(features)
  structurally_singular <- d > n - 1L
  if (structurally_singular && isTRUE(as.numeric(ridge) == 0)) {
    stop(paste(
      "Without ridge regularization, quadratic self-normalization requires",
      "feature dimension <= calibration size - 1."
    ))
  }
  feature_mean <- colMeans(features)
  centered <- sweep(features, 2L, feature_mean, FUN = "-")
  partial_sums <- apply(centered, 2L, cumsum)
  if (is.null(dim(partial_sums))) partial_sums <- matrix(partial_sums, ncol = 1L)
  self_normalizer <- crossprod(partial_sums) / (n^2)
  # The quadratic-SN statistic is exactly invariant to a common nonsingular
  # coordinate transformation only when the bridge matrix is inverted without
  # an additive ridge.  Use the exact inverse by default; a positive ridge is
  # an explicit finite-machine fallback rather than part of the paper's
  # reported statistic.
  if (structurally_singular) {
    forced_ridge <- dc_relative_ridge(self_normalizer, ridge)
    inverse_info <- dc_inverse_spd(
      self_normalizer + forced_ridge * diag(d), ridge = 0
    )
    inverse_info$ridge_added <- inverse_info$ridge_added + forced_ridge
  } else {
    inverse_info <- dc_inverse_spd(self_normalizer, ridge = ridge)
  }
  list(
    feature_mean = feature_mean,
    self_normalizer = self_normalizer,
    inverse = inverse_info$inverse,
    regularized_self_normalizer = inverse_info$matrix,
    ridge_added = inverse_info$ridge_added,
    condition_number = inverse_info$condition_number,
    n_calibration = n,
    dimension = d,
    feature_names = colnames(features),
    normalizer_type = "full_matrix_sequential"
  )
}

# Componentwise adjusted-range normalization for an already LDL-demixed feature
# stream.  The outer LDL map is fitted once in dc_fit_score_monitor(); applying
# a second internal LDL factor here would define a different finite-sample
# statistic from Equation (2.7) of the paper.
dc_fit_adjusted_range_sn <- function(features, ridge = 0) {
  features <- dc_assert_matrix(features, "features", min_rows = 3L)
  n <- nrow(features)
  d <- ncol(features)
  feature_mean <- colMeans(features)
  centered <- sweep(features, 2L, feature_mean, FUN = "-")
  partial_sums <- apply(centered, 2L, cumsum) / sqrt(n)
  if (is.null(dim(partial_sums))) partial_sums <- matrix(partial_sums, ncol = 1L)
  partial_sums <- rbind(rep(0, d), partial_sums)
  adjusted_ranges <- apply(partial_sums, 2L, max) - apply(partial_sums, 2L, min)
  adjusted_ranges <- pmax(adjusted_ranges, sqrt(.Machine$double.eps))
  self_normalizer <- diag(adjusted_ranges^2, nrow = d)
  inverse <- diag(1 / adjusted_ranges^2, nrow = d)
  list(
    feature_mean = feature_mean,
    self_normalizer = self_normalizer,
    inverse = inverse,
    regularized_self_normalizer = self_normalizer,
    ridge_added = 0,
    condition_number = max(adjusted_ranges^2) / min(adjusted_ranges^2),
    n_calibration = n,
    dimension = d,
    feature_names = colnames(features),
    normalizer_type = "adjusted_range",
    adjusted_ranges = adjusted_ranges
  )
}

dc_ae_gauge_lp <- function(G, z, eps = .Machine$double.eps) {
  G <- dc_assert_matrix(G, "G", min_rows = 2L)
  z <- as.numeric(z)
  q <- length(z)
  if (ncol(G) != q) stop("ncol(G) must equal length(z).")
  if (q == 1L) {
    denominator <- max(G[, 1L]) - min(G[, 1L])
    return(abs(z[1L]) / max(denominator, eps))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required for affine-range gauges in dimension above one.")
  }

  objective <- c(z, -z, 0, 0)
  upper_rows <- cbind(G, -G, -1, 1)
  lower_rows <- cbind(-G, G, 1, -1)
  constraint_matrix <- rbind(upper_rows, lower_rows)
  constraint_rhs <- c(rep(1, nrow(G)), rep(0, nrow(G)))
  fit <- lpSolve::lp(
    direction = "max",
    objective.in = objective,
    const.mat = constraint_matrix,
    const.dir = rep("<=", length(constraint_rhs)),
    const.rhs = constraint_rhs,
    scale = 196
  )
  if (fit$status != 0L || !is.finite(fit$objval)) return(NA_real_)
  max(fit$objval, 0)
}

dc_fit_affine_range_sn <- function(features, ridge = NULL, rank_tolerance = 1e-10,
                                   hull_grid = NULL) {
  features <- dc_assert_matrix(features, "features", min_rows = 3L)
  n <- nrow(features)
  d <- ncol(features)
  feature_mean <- colMeans(features)
  centered <- sweep(features, 2L, feature_mean, FUN = "-")
  partial_sums <- apply(centered, 2L, cumsum)
  if (is.null(dim(partial_sums))) partial_sums <- matrix(partial_sums, ncol = 1L)
  G <- rbind(rep(0, d), partial_sums / sqrt(n))
  if (!is.null(hull_grid)) {
    hull_grid <- as.integer(hull_grid)
    if (hull_grid < d + 1L) stop("hull_grid must exceed the channel dimension.")
    selected_rows <- unique(as.integer(round(seq(1, nrow(G), length.out = hull_grid + 1L))))
    G <- G[selected_rows, , drop = FALSE]
  }
  singular_values <- svd(G, nu = 0L, nv = 0L)$d
  numerical_rank <- sum(singular_values > max(singular_values) * rank_tolerance)
  if (numerical_rank < d) {
    stop("The centered training path does not span the full channel dimension.")
  }
  list(
    feature_mean = feature_mean,
    training_path = G,
    n_calibration = n,
    dimension = d,
    feature_names = colnames(features),
    normalizer_type = "affine_range",
    hull_grid = if (is.null(hull_grid)) n else hull_grid,
    numerical_rank = numerical_rank,
    smallest_relative_singular_value = min(singular_values) / max(singular_values),
    condition_number = max(singular_values) / min(singular_values),
    ridge_added = 0
  )
}

dc_affine_range_monitor_path <- function(sn_fit, monitor_features, gamma = 0) {
  monitor_features <- dc_assert_matrix(monitor_features, "monitor_features", min_rows = 1L)
  if (ncol(monitor_features) != sn_fit$dimension) {
    stop("monitor_features has the wrong number of columns.")
  }
  centered <- sweep(monitor_features, 2L, sn_fit$feature_mean, FUN = "-")
  cumulative <- apply(centered, 2L, cumsum)
  if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, ncol = 1L)
  Z <- cumulative / sqrt(sn_fit$n_calibration)
  gauge_squared <- apply(Z, 1L, function(z) dc_ae_gauge_lp(sn_fit$training_path, z)^2)
  s <- seq_len(nrow(monitor_features)) / sn_fit$n_calibration
  as.numeric(gauge_squared / dc_boundary_sq(s, gamma = gamma))
}

dc_boundary_sq <- function(s, gamma = 0) {
  if (any(s <= 0)) stop("s must be positive.")
  (1 + s)^2 * (s / (1 + s))^(2 * gamma)
}

dc_sn_monitor_path <- function(sn_fit, monitor_features, gamma = 0) {
  monitor_features <- dc_assert_matrix(monitor_features, "monitor_features", min_rows = 1L)
  if (ncol(monitor_features) != sn_fit$dimension) {
    stop("monitor_features has the wrong number of columns.")
  }
  centered <- sweep(monitor_features, 2L, sn_fit$feature_mean, FUN = "-")
  cumulative <- apply(centered, 2L, cumsum)
  if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, ncol = 1L)
  quadratic <- rowSums((cumulative %*% sn_fit$inverse) * cumulative) / sn_fit$n_calibration
  s <- seq_len(nrow(monitor_features)) / sn_fit$n_calibration
  as.numeric(quadratic / dc_boundary_sq(s, gamma = gamma))
}

dc_simulate_sn_critical_value <- function(dimension, T = 2, gamma = 0,
                                          alpha = 0.025, nrep = 2000L,
                                          train_grid = 400L, seed = 20260830L,
                                          ridge = 0) {
  dimension <- as.integer(dimension)
  nrep <- as.integer(nrep)
  train_grid <- as.integer(train_grid)
  if (dimension < 1L || nrep < 20L || train_grid < 50L) stop("Invalid simulation settings.")
  if (!is.numeric(T) || T <= 0 || !is.finite(T)) stop("T must be positive and finite.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must lie in (0,1).")
  monitor_grid <- max(1L, as.integer(round(T * train_grid)))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  maxima <- numeric(nrep)
  total_steps <- train_grid + monitor_grid
  r <- seq_len(train_grid) / train_grid
  s <- seq_len(monitor_grid) / train_grid
  denominator <- dc_boundary_sq(s, gamma = gamma)

  for (replicate_index in seq_len(nrep)) {
    increments <- matrix(
      stats::rnorm(total_steps * dimension),
      nrow = total_steps,
      ncol = dimension
    ) / sqrt(train_grid)
    brownian <- apply(increments, 2L, cumsum)
    if (is.null(dim(brownian))) brownian <- matrix(brownian, ncol = 1L)
    terminal_training <- brownian[train_grid, ]
    bridge <- brownian[seq_len(train_grid), , drop = FALSE] - outer(r, terminal_training)
    V <- crossprod(bridge) / train_grid
    V_inv <- dc_inverse_spd(V, ridge = ridge)$inverse
    after_training <- brownian[train_grid + seq_len(monitor_grid), , drop = FALSE]
    monitoring_process <- after_training - outer(1 + s, terminal_training)
    quadratic <- rowSums((monitoring_process %*% V_inv) * monitoring_process)
    maxima[replicate_index] <- max(quadratic / denominator)
  }

  value <- as.numeric(stats::quantile(maxima, probs = 1 - alpha, type = 8, names = FALSE))
  list(
    critical_value = value,
    alpha = alpha,
    dimension = dimension,
    T = T,
    gamma = gamma,
    nrep = nrep,
    train_grid = train_grid,
    monte_carlo_se = sqrt(alpha * (1 - alpha) / nrep)
  )
}

dc_simulate_adjusted_range_critical_value <- function(dimension, T = 2, gamma = 0,
                                                      alpha = 0.025, nrep = 2000L,
                                                      train_grid = 400L,
                                                      seed = 20260830L) {
  dimension <- as.integer(dimension)
  nrep <- as.integer(nrep)
  train_grid <- as.integer(train_grid)
  if (dimension < 1L || nrep < 20L || train_grid < 50L) stop("Invalid simulation settings.")
  if (!is.numeric(T) || T <= 0 || !is.finite(T)) stop("T must be positive and finite.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must lie in (0,1).")
  monitor_grid <- max(1L, as.integer(round(T * train_grid)))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  maxima <- numeric(nrep)
  total_steps <- train_grid + monitor_grid
  r <- seq_len(train_grid) / train_grid
  s <- seq_len(monitor_grid) / train_grid
  denominator <- dc_boundary_sq(s, gamma = gamma)
  for (replicate_index in seq_len(nrep)) {
    increments <- matrix(
      stats::rnorm(total_steps * dimension), nrow = total_steps, ncol = dimension
    ) / sqrt(train_grid)
    brownian <- apply(increments, 2L, cumsum)
    if (is.null(dim(brownian))) brownian <- matrix(brownian, ncol = 1L)
    terminal_training <- brownian[train_grid, ]
    bridge <- brownian[seq_len(train_grid), , drop = FALSE] - outer(r, terminal_training)
    bridge <- rbind(rep(0, dimension), bridge)
    ranges <- apply(bridge, 2L, max) - apply(bridge, 2L, min)
    inverse_diagonal <- 1 / pmax(ranges^2, .Machine$double.eps)
    after_training <- brownian[train_grid + seq_len(monitor_grid), , drop = FALSE]
    monitoring_process <- after_training - outer(1 + s, terminal_training)
    quadratic <- rowSums(sweep(monitoring_process^2, 2L, inverse_diagonal, FUN = "*"))
    maxima[replicate_index] <- max(quadratic / denominator)
  }
  value <- as.numeric(stats::quantile(maxima, probs = 1 - alpha, type = 8, names = FALSE))
  list(
    critical_value = value,
    alpha = alpha,
    dimension = dimension,
    T = T,
    gamma = gamma,
    nrep = nrep,
    train_grid = train_grid,
    normalizer_type = "adjusted_range"
  )
}

dc_simulate_affine_range_critical_value <- function(dimension, T = 2, gamma = 0,
                                                    alpha = 0.025, nrep = 1000L,
                                                    train_grid = 100L,
                                                    seed = 20260830L,
                                                    workers = 1L) {
  dimension <- as.integer(dimension)
  nrep <- as.integer(nrep)
  train_grid <- as.integer(train_grid)
  workers <- as.integer(workers)
  if (dimension < 1L || nrep < 20L || train_grid < 30L) stop("Invalid simulation settings.")
  if (!is.numeric(T) || T <= 0 || !is.finite(T)) stop("T must be positive and finite.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must lie in (0,1).")
  monitor_grid <- max(1L, as.integer(round(T * train_grid)))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  replicate_seeds <- sample.int(.Machine$integer.max, nrep)
  total_steps <- train_grid + monitor_grid
  r <- seq_len(train_grid) / train_grid
  s <- seq_len(monitor_grid) / train_grid
  denominator <- dc_boundary_sq(s, gamma = gamma)
  simulate_one <- function(replicate_index) {
    set.seed(replicate_seeds[replicate_index])
    increments <- matrix(
      stats::rnorm(total_steps * dimension), nrow = total_steps, ncol = dimension
    ) / sqrt(train_grid)
    brownian <- apply(increments, 2L, cumsum)
    if (is.null(dim(brownian))) brownian <- matrix(brownian, ncol = 1L)
    terminal_training <- brownian[train_grid, ]
    bridge <- brownian[seq_len(train_grid), , drop = FALSE] - outer(r, terminal_training)
    G <- rbind(rep(0, dimension), bridge)
    after_training <- brownian[train_grid + seq_len(monitor_grid), , drop = FALSE]
    monitoring_process <- after_training - outer(1 + s, terminal_training)
    gauge_squared <- apply(
      monitoring_process, 1L, function(z) dc_ae_gauge_lp(G, z)^2
    )
    max(gauge_squared / denominator)
  }
  maxima <- if (workers > 1L) {
    unlist(parallel::mclapply(
      seq_len(nrep), simulate_one,
      mc.cores = workers, mc.preschedule = TRUE
    ), use.names = FALSE)
  } else vapply(seq_len(nrep), simulate_one, numeric(1L))
  value <- as.numeric(stats::quantile(maxima, probs = 1 - alpha, type = 8, names = FALSE))
  list(
    critical_value = value,
    alpha = alpha,
    dimension = dimension,
    T = T,
    gamma = gamma,
    nrep = nrep,
    train_grid = train_grid,
    normalizer_type = "affine_range"
  )
}

dc_stationary_bootstrap_indices <- function(n_source, n_draw, block_length) {
  n_source <- as.integer(n_source)
  n_draw <- as.integer(n_draw)
  restart_probability <- 1 / as.numeric(block_length)
  index <- integer(n_draw)
  index[1L] <- sample.int(n_source, 1L)
  if (n_draw > 1L) {
    for (tt in 2:n_draw) {
      if (stats::runif(1L) < restart_probability) {
        index[tt] <- sample.int(n_source, 1L)
      } else {
        index[tt] <- if (index[tt - 1L] == n_source) 1L else index[tt - 1L] + 1L
      }
    }
  }
  index
}

dc_select_stationary_block_length <- function(features, max_lag = NULL) {
  features <- dc_assert_matrix(features, "features", min_rows = 20L)
  n <- nrow(features)
  if (is.null(max_lag)) max_lag <- min(30L, max(5L, floor(10 * log10(n))))
  standardized <- scale(features)
  bad_columns <- !is.finite(colSums(standardized))
  if (any(bad_columns)) standardized[, bad_columns] <- 0
  positive_envelope <- numeric(max_lag)
  for (lag in seq_len(max_lag)) {
    left <- standardized[seq_len(n - lag), , drop = FALSE]
    right <- standardized[lag + seq_len(n - lag), , drop = FALSE]
    correlations <- colMeans(left * right)
    positive_envelope[lag] <- max(0, stats::median(correlations, na.rm = TRUE))
  }
  taper <- 1 - seq_len(max_lag) / (max_lag + 1)
  dependence_factor <- max(1, 1 + 2 * sum(taper * positive_envelope))
  selected <- ceiling(n^(1 / 3) * dependence_factor^(2 / 3))
  as.integer(max(5L, min(floor(n / 4), selected)))
}

# Finite-sample critical values obtained exclusively from the normalization block.
# Each stationary-bootstrap draw supplies a pseudo normalization block and a
# pseudo continuation; the full self-normalizer and normalization mean are refit
# on the former. This reproduces the estimation term in sequential monitoring.
# The observed monitoring sample is never used to choose the boundary.
dc_stationary_bootstrap_critical_values <- function(fit, horizon_rows,
                                                    alpha = 0.05,
                                                    alpha_weights = c(S = 0.5, R = 0.5),
                                                    gamma = 0,
                                                    nrep = 999L,
                                                    block_length = NULL,
                                                    seed = 20260830L) {
  if (!inherits(fit, "dc_monitor_fit")) stop("fit must be produced by dc_fit_dual_channel().")
  nrep <- as.integer(nrep)
  horizon_rows <- as.integer(horizon_rows)
  if (nrep < 99L) stop("nrep must be at least 99 for bootstrap critical values.")
  if (horizon_rows < 1L) stop("horizon_rows must be positive.")
  alpha_weights <- alpha_weights / sum(alpha_weights)
  names(alpha_weights) <- c("S", "R")
  source_S <- fit$calibration_sn_features$channel_S
  source_R <- fit$calibration_sn_features$channel_R
  m <- nrow(source_S)
  joint <- cbind(scale(source_S), scale(source_R))
  joint[!is.finite(joint)] <- 0
  if (is.null(block_length)) block_length <- dc_select_stationary_block_length(joint)
  block_length <- as.integer(block_length)
  if (block_length < 2L || block_length > floor(m / 2)) {
    stop("block_length must be between two and half the normalization-sample size.")
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  maxima_S <- numeric(nrep)
  maxima_R <- numeric(nrep)
  for (replicate_index in seq_len(nrep)) {
    index <- dc_stationary_bootstrap_indices(m, m + horizon_rows, block_length)
    bootstrap_S <- source_S[index, , drop = FALSE]
    bootstrap_R <- source_R[index, , drop = FALSE]
    fit_S <- dc_fit_full_matrix_sn(
      bootstrap_S[seq_len(m), , drop = FALSE], ridge = fit$settings$sn_ridge
    )
    fit_R <- dc_fit_full_matrix_sn(
      bootstrap_R[seq_len(m), , drop = FALSE], ridge = fit$settings$sn_ridge
    )
    monitor_index <- m + seq_len(horizon_rows)
    maxima_S[replicate_index] <- max(dc_sn_monitor_path(
      fit_S, bootstrap_S[monitor_index, , drop = FALSE], gamma = gamma
    ))
    maxima_R[replicate_index] <- max(dc_sn_monitor_path(
      fit_R, bootstrap_R[monitor_index, , drop = FALSE], gamma = gamma
    ))
  }

  critical_values <- c(
    S = as.numeric(stats::quantile(
      maxima_S, 1 - alpha * alpha_weights["S"], type = 8, names = FALSE
    )),
    R = as.numeric(stats::quantile(
      maxima_R, 1 - alpha * alpha_weights["R"], type = 8, names = FALSE
    ))
  )
  list(
    critical_values = critical_values,
    block_length = block_length,
    nrep = nrep,
    horizon_rows = horizon_rows,
    alpha = alpha,
    alpha_weights = alpha_weights,
    gamma = gamma,
    maxima = list(S = maxima_S, R = maxima_R)
  )
}

dc_fit_dual_channel <- function(x_fit, x_calibration, K = NULL, L = 1L,
                                fve = 0.90, grid = NULL, weights = NULL,
                                prewhiten_mode = c("demix", "whiten", "none"),
                                ldl_ridge = 1e-8, sn_ridge = 1e-8,
                                sn_type = c("partial_sum", "adjusted_range", "affine_range"),
                                feature_var_order = 0L,
                                r_huber_threshold = Inf,
                                ae_hull_grid = NULL) {
  x_fit <- dc_assert_matrix(x_fit, "x_fit", min_rows = 3L, min_cols = 2L)
  x_calibration <- dc_assert_matrix(x_calibration, "x_calibration", min_rows = 5L, min_cols = ncol(x_fit))
  if (ncol(x_fit) != ncol(x_calibration)) stop("Fit and normalization curves must share a grid.")
  prewhiten_mode <- match.arg(prewhiten_mode)
  sn_type <- match.arg(sn_type)

  fpca <- dc_fit_fpca(x_fit, K = K, L = L, fve = fve, grid = grid, weights = weights)
  projected_calibration <- dc_project_fpca(fpca, x_calibration)
  base_calibration <- cbind(projected_calibration$retained, projected_calibration$audit)
  prewhitener <- dc_fit_ldl_prewhitener(
    base_calibration,
    ridge = ldl_ridge,
    mode = prewhiten_mode
  )
  calibration_features <- dc_build_features(
    projected_calibration,
    prewhitener,
    K = fpca$K,
    L = fpca$L
  )
  robust_R <- dc_fit_huberizer(
    calibration_features$channel_R,
    threshold = r_huber_threshold
  )

  prewhitened_S <- dc_fit_var_ldl_features(
    calibration_features$channel_S,
    order = feature_var_order,
    ridge = ldl_ridge,
    ldl_mode = "demix"
  )
  prewhitened_R <- dc_fit_var_ldl_features(
    robust_R$transformed,
    order = feature_var_order,
    ridge = ldl_ridge,
    ldl_mode = "demix"
  )
  calibration_sn_features <- list(
    channel_S = prewhitened_S$transformed,
    channel_R = prewhitened_R$transformed
  )

  if (sn_type == "affine_range") {
    sn_S <- dc_fit_affine_range_sn(
      calibration_sn_features$channel_S, hull_grid = ae_hull_grid
    )
    sn_R <- dc_fit_affine_range_sn(
      calibration_sn_features$channel_R, hull_grid = ae_hull_grid
    )
  } else {
    fit_normalizer <- if (sn_type == "partial_sum") dc_fit_full_matrix_sn else dc_fit_adjusted_range_sn
    sn_S <- fit_normalizer(calibration_sn_features$channel_S, ridge = sn_ridge)
    sn_R <- fit_normalizer(calibration_sn_features$channel_R, ridge = sn_ridge)
  }

  structure(list(
    fpca = fpca,
    prewhitener = prewhitener,
    feature_prewhitener_S = prewhitened_S$model,
    feature_prewhitener_R = prewhitened_R$model,
    R_huberizer = robust_R$model,
    sn_S = sn_S,
    sn_R = sn_R,
    calibration_features = calibration_features,
    calibration_sn_features = calibration_sn_features,
    settings = list(
      K = fpca$K,
      L = fpca$L,
      fve = fve,
      prewhiten_mode = prewhiten_mode,
      sn_type = sn_type,
      feature_var_order = as.integer(feature_var_order),
      r_huber_threshold = as.numeric(r_huber_threshold),
      ae_hull_grid = ae_hull_grid,
      ldl_ridge = ldl_ridge,
      sn_ridge = sn_ridge
    )
  ), class = "dc_monitor_fit")
}

# Restored score-only monitor used by the paper's original main method. The
# extra FPCA direction is estimated only so the same FPCA routine can be
# reused; it does not enter the statistic or the stopping rule.
dc_fit_score_monitor <- function(x_fit, x_calibration, K = NULL, fve = 0.90,
                                 K_max = NULL,
                                 grid = NULL, weights = NULL,
                                 prewhiten_mode = c("demix", "whiten", "none"),
                                 ldl_ridge = 1e-8, sn_ridge = 1e-8,
                                 sn_type = c("adjusted_range", "partial_sum"),
                                 feature_var_order = 0L) {
  x_fit <- dc_assert_matrix(x_fit, "x_fit", min_rows = 3L, min_cols = 2L)
  x_calibration <- dc_assert_matrix(
    x_calibration, "x_calibration", min_rows = 5L, min_cols = ncol(x_fit)
  )
  if (ncol(x_fit) != ncol(x_calibration)) {
    stop("Fit and normalization curves must share a grid.")
  }
  prewhiten_mode <- match.arg(prewhiten_mode)
  sn_type <- match.arg(sn_type)

  fpca <- dc_fit_fpca(
    x_fit, K = K, L = 1L, fve = fve, K_max = K_max,
    grid = grid, weights = weights
  )
  projected_calibration <- dc_project_fpca(fpca, x_calibration)
  prewhitener <- dc_fit_ldl_prewhitener(
    projected_calibration$retained,
    ridge = ldl_ridge,
    mode = prewhiten_mode
  )
  calibration_scores <- dc_apply_ldl_prewhitener(
    prewhitener, projected_calibration$retained
  )
  colnames(calibration_scores) <- paste0("score_", seq_len(fpca$K))
  feature_prewhitener <- dc_fit_var_ldl_features(
    calibration_scores,
    order = feature_var_order,
    ridge = ldl_ridge,
    ldl_mode = "demix"
  )
  calibration_sn_scores <- feature_prewhitener$transformed
  fit_normalizer <- if (sn_type == "adjusted_range") {
    dc_fit_adjusted_range_sn
  } else dc_fit_full_matrix_sn
  sn <- fit_normalizer(calibration_sn_scores, ridge = sn_ridge)

  structure(list(
    fpca = fpca,
    prewhitener = prewhitener,
    feature_prewhitener = feature_prewhitener$model,
    sn = sn,
    calibration_scores = calibration_scores,
    calibration_sn_scores = calibration_sn_scores,
    settings = list(
      K = fpca$K,
      K_max = fpca$K_max,
      fve = fve,
      prewhiten_mode = prewhiten_mode,
      sn_type = sn_type,
      feature_var_order = as.integer(feature_var_order),
      ldl_ridge = ldl_ridge,
      sn_ridge = sn_ridge
    )
  ), class = "dc_score_monitor_fit")
}

# Fit the paper's adjusted-range method and Shao's quadratic self-normalization
# benchmark from one common FPCA projection and one common LDL-demixed
# retained-score stream. The quadratic-SN statistic is affine
# invariant, so applying the shared nonsingular LDL transformation does not
# substantively change it; the shared map keeps preprocessing exactly identical
# across the two reported methods.
dc_fit_score_comparison <- function(x_fit, x_calibration, K = NULL, fve = 0.90,
                                    K_max = NULL,
                                    grid = NULL, weights = NULL,
                                    prewhiten_mode = c("demix", "whiten", "none"),
                                    ldl_ridge = 1e-8, sn_ridge = 0,
                                    feature_var_order = 0L) {
  prewhiten_mode <- match.arg(prewhiten_mode)
  adjusted_range <- dc_fit_score_monitor(
    x_fit = x_fit,
    x_calibration = x_calibration,
    K = K,
    fve = fve,
    K_max = K_max,
    grid = grid,
    weights = weights,
    prewhiten_mode = prewhiten_mode,
    ldl_ridge = ldl_ridge,
    sn_ridge = sn_ridge,
    sn_type = "adjusted_range",
    feature_var_order = feature_var_order
  )
  shao <- adjusted_range
  shao$sn <- dc_fit_full_matrix_sn(
    adjusted_range$calibration_sn_scores,
    ridge = 0
  )
  shao$settings$sn_type <- "partial_sum"

  structure(list(
    adjusted_range = adjusted_range,
    shao = shao,
    settings = list(
      K = adjusted_range$fpca$K,
      K_max = adjusted_range$fpca$K_max,
      fve = fve,
      prewhiten_mode = prewhiten_mode,
      feature_var_order = as.integer(feature_var_order),
      ldl_ridge = ldl_ridge,
      sn_ridge = sn_ridge
    )
  ), class = "dc_score_comparison_fit")
}

dc_monitor_score <- function(fit, x_monitor, alpha = 0.05, gamma = 0,
                             critical_value = NULL, critical_nrep = 2000L,
                             critical_train_grid = 400L,
                             seed = 20260830L) {
  if (!inherits(fit, "dc_score_monitor_fit")) {
    stop("fit must be produced by dc_fit_score_monitor().")
  }
  x_monitor <- dc_assert_matrix(
    x_monitor, "x_monitor", min_rows = 1L, min_cols = fit$fpca$n_grid
  )
  if (ncol(x_monitor) != fit$fpca$n_grid) {
    stop("Monitoring curves must share the fitted grid.")
  }
  projected <- dc_project_fpca(fit$fpca, x_monitor)
  scores <- dc_apply_ldl_prewhitener(fit$prewhitener, projected$retained)
  colnames(scores) <- paste0("score_", seq_len(fit$fpca$K))
  sn_scores <- dc_apply_var_ldl_features(fit$feature_prewhitener, scores)
  path <- dc_sn_monitor_path(fit$sn, sn_scores, gamma = gamma)
  horizon <- nrow(x_monitor) / fit$sn$n_calibration

  if (is.null(critical_value)) {
    simulator <- if (fit$settings$sn_type == "adjusted_range") {
      dc_simulate_adjusted_range_critical_value
    } else dc_simulate_sn_critical_value
    critical_metadata <- simulator(
      fit$sn$dimension,
      T = horizon,
      gamma = gamma,
      alpha = alpha,
      nrep = critical_nrep,
      train_grid = critical_train_grid,
      seed = seed
    )
    critical_value <- critical_metadata$critical_value
  } else {
    critical_value <- as.numeric(critical_value)
    if (length(critical_value) != 1L || !is.finite(critical_value)) {
      stop("critical_value must be one finite number.")
    }
    critical_metadata <- list(method = "supplied")
  }
  crossing <- which(path > critical_value)
  stop_time <- if (length(crossing)) crossing[1L] else NA_integer_
  list(
    path = data.frame(
      k = seq_along(path), statistic = path, boundary = critical_value
    ),
    alarm = list(stop = stop_time),
    critical_value = critical_value,
    critical_metadata = critical_metadata,
    alpha = alpha,
    gamma = gamma,
    scores = scores,
    sn_scores = sn_scores,
    diagnostics = list(
      condition_number = fit$sn$condition_number,
      prewhitener_mode = fit$prewhitener$mode,
      feature_var_order = fit$settings$feature_var_order,
      K = fit$fpca$K
    )
  )
}

dc_monitor_score_comparison <- function(fit, x_monitor, alpha = 0.05,
                                        gamma = 0, critical_values = NULL,
                                        critical_nrep = 2000L,
                                        critical_train_grid = 400L,
                                        seed = 20260830L) {
  if (!inherits(fit, "dc_score_comparison_fit")) {
    stop("fit must be produced by dc_fit_score_comparison().")
  }
  if (!is.null(critical_values)) {
    if (is.null(names(critical_values)) ||
        !all(c("adjusted_range", "shao") %in% names(critical_values))) {
      stop("critical_values must be named adjusted_range and shao.")
    }
  }
  range_result <- dc_monitor_score(
    fit$adjusted_range,
    x_monitor,
    alpha = alpha,
    gamma = gamma,
    critical_value = if (is.null(critical_values)) NULL else critical_values[["adjusted_range"]],
    critical_nrep = critical_nrep,
    critical_train_grid = critical_train_grid,
    seed = seed
  )
  shao_result <- dc_monitor_score(
    fit$shao,
    x_monitor,
    alpha = alpha,
    gamma = gamma,
    critical_value = if (is.null(critical_values)) NULL else critical_values[["shao"]],
    critical_nrep = critical_nrep,
    critical_train_grid = critical_train_grid,
    seed = seed
  )
  list(
    path = data.frame(
      k = range_result$path$k,
      adjusted_range = range_result$path$statistic,
      shao = shao_result$path$statistic,
      critical_adjusted_range = range_result$critical_value,
      critical_shao = shao_result$critical_value
    ),
    alarm = c(
      adjusted_range = range_result$alarm$stop,
      shao = shao_result$alarm$stop
    ),
    critical_values = c(
      adjusted_range = range_result$critical_value,
      shao = shao_result$critical_value
    ),
    adjusted_range = range_result,
    shao = shao_result,
    alpha = alpha,
    gamma = gamma
  )
}

dc_channel_label <- function(path_S, path_R, critical_S, critical_R) {
  crossed_S <- which(path_S > critical_S)
  crossed_R <- which(path_R > critical_R)
  first_S <- if (length(crossed_S)) crossed_S[1L] else NA_integer_
  first_R <- if (length(crossed_R)) crossed_R[1L] else NA_integer_
  first_global <- suppressWarnings(min(c(first_S, first_R), na.rm = TRUE))
  if (!is.finite(first_global)) {
    return(list(stop = NA_integer_, label = "none", first_S = first_S, first_R = first_R))
  }
  implicated_S <- path_S[first_global] > critical_S
  implicated_R <- path_R[first_global] > critical_R
  label <- if (implicated_S && implicated_R) "both" else if (implicated_S) "S" else "R"
  list(stop = as.integer(first_global), label = label, first_S = first_S, first_R = first_R)
}

dc_monitor <- function(fit, x_monitor, alpha = 0.05,
                       alpha_weights = c(S = 0.5, R = 0.5), gamma = 0,
                       critical_values = NULL, critical_nrep = 2000L,
                       critical_train_grid = 400L, seed = 20260830L,
                       critical_method = c("brownian", "stationary_bootstrap"),
                       bootstrap_block_length = NULL,
                       critical_workers = 1L) {
  if (!inherits(fit, "dc_monitor_fit")) stop("fit must be produced by dc_fit_dual_channel().")
  x_monitor <- dc_assert_matrix(x_monitor, "x_monitor", min_rows = 1L, min_cols = fit$fpca$n_grid)
  if (ncol(x_monitor) != fit$fpca$n_grid) stop("Monitoring curves must share the fitted grid.")
  if (length(alpha_weights) != 2L || any(alpha_weights <= 0)) {
    stop("alpha_weights must contain two positive values.")
  }
  alpha_weights <- alpha_weights / sum(alpha_weights)
  names(alpha_weights) <- c("S", "R")
  critical_method <- match.arg(critical_method)

  projected <- dc_project_fpca(fit$fpca, x_monitor)
  features <- dc_build_features(projected, fit$prewhitener, fit$fpca$K, fit$fpca$L)
  robust_features_R <- dc_apply_huberizer(fit$R_huberizer, features$channel_R)
  sn_features <- list(
    channel_S = dc_apply_var_ldl_features(fit$feature_prewhitener_S, features$channel_S),
    channel_R = dc_apply_var_ldl_features(fit$feature_prewhitener_R, robust_features_R)
  )
  monitor_path_function <- if (identical(fit$settings$sn_type, "affine_range")) {
    dc_affine_range_monitor_path
  } else dc_sn_monitor_path
  path_S <- monitor_path_function(fit$sn_S, sn_features$channel_S, gamma = gamma)
  path_R <- monitor_path_function(fit$sn_R, sn_features$channel_R, gamma = gamma)
  horizon <- nrow(x_monitor) / fit$sn_S$n_calibration

  if (is.null(critical_values)) {
    if (critical_method == "brownian") {
      critical_simulator <- switch(
        fit$settings$sn_type,
        partial_sum = dc_simulate_sn_critical_value,
        adjusted_range = dc_simulate_adjusted_range_critical_value,
        affine_range = dc_simulate_affine_range_critical_value
      )
      simulate_channel_critical <- function(dimension, channel_alpha, channel_seed) {
        arguments <- list(
          dimension = dimension, T = horizon, gamma = gamma,
          alpha = channel_alpha, nrep = critical_nrep,
          train_grid = critical_train_grid, seed = channel_seed
        )
        if (identical(fit$settings$sn_type, "affine_range")) {
          arguments$workers <- critical_workers
        }
        do.call(critical_simulator, arguments)
      }
      cv_S <- simulate_channel_critical(
        fit$sn_S$dimension, alpha * alpha_weights["S"], seed
      )
      cv_R <- simulate_channel_critical(
        fit$sn_R$dimension, alpha * alpha_weights["R"], seed + 1L
      )
      critical_values <- c(S = cv_S$critical_value, R = cv_R$critical_value)
      critical_metadata <- list(method = "brownian", S = cv_S, R = cv_R)
    } else {
      bootstrap <- dc_stationary_bootstrap_critical_values(
        fit,
        horizon_rows = nrow(x_monitor),
        alpha = alpha,
        alpha_weights = alpha_weights,
        gamma = gamma,
        nrep = critical_nrep,
        block_length = bootstrap_block_length,
        seed = seed
      )
      critical_values <- bootstrap$critical_values
      critical_metadata <- c(list(method = "stationary_bootstrap"), bootstrap)
    }
  } else {
    if (is.null(names(critical_values)) || !all(c("S", "R") %in% names(critical_values))) {
      stop("critical_values must be a named vector with entries S and R.")
    }
    critical_values <- critical_values[c("S", "R")]
    critical_metadata <- list(method = "supplied")
  }

  alarm <- dc_channel_label(path_S, path_R, critical_values["S"], critical_values["R"])
  list(
    path = data.frame(
      k = seq_along(path_S),
      S = path_S,
      R = path_R,
      boundary_S = unname(critical_values["S"]),
      boundary_R = unname(critical_values["R"])
    ),
    alarm = alarm,
    critical_values = critical_values,
    critical_metadata = critical_metadata,
    alpha = alpha,
    alpha_weights = alpha_weights,
    gamma = gamma,
    critical_method = if (is.null(critical_metadata$method)) "supplied" else critical_metadata$method,
    features = features,
    robust_features_R = robust_features_R,
    sn_features = sn_features,
    diagnostics = list(
      S_condition_number = fit$sn_S$condition_number,
      R_condition_number = fit$sn_R$condition_number,
      prewhitener_mode = fit$prewhitener$mode,
      feature_var_order = fit$settings$feature_var_order,
      r_huber_threshold = fit$settings$r_huber_threshold,
      K = fit$fpca$K,
      L = fit$fpca$L,
      d_S = fit$sn_S$dimension,
      d_R = fit$sn_R$dimension
    )
  )
}

print.dc_monitor_fit <- function(x, ...) {
  cat("Dual-channel fixed-FPCA monitor\n")
  cat("  K retained scores:", x$fpca$K, "\n")
  cat("  L audit scores:", x$fpca$L, "\n")
  cat("  Channel dimensions: S =", x$sn_S$dimension, ", R =", x$sn_R$dimension, "\n")
  cat("  LDL mode:", x$prewhitener$mode, "\n")
  cat("  SN type:", x$settings$sn_type, "\n")
  cat("  Feature VAR order:", x$settings$feature_var_order, "\n")
  cat("  R Huber threshold:", x$settings$r_huber_threshold, "\n")
  cat("  SN condition numbers: S =", format(x$sn_S$condition_number, digits = 4),
      ", R =", format(x$sn_R$condition_number, digits = 4), "\n")
  invisible(x)
}
