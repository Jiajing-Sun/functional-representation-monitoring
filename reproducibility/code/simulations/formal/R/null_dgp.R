# Null data-generating processes for the retained-score comparison study.

# Locate the archived completion relative to this source, not the working directory.
dc_frozen_rotation_path <- tryCatch({
  source_path <- normalizePath(sys.frame(1L)$ofile, mustWork = TRUE)
  file.path(dirname(source_path), "..", "..", "..", "..", "data",
            "fixed_functional_rotation_8.csv")
}, error = function(e) character())

dc_make_correlation <- function(d, rho = 0.35) {
  index <- seq_len(d)
  rho^abs(outer(index, index, FUN = "-"))
}

dc_score_eigenvalues <- function(d = 8L,
                                 profile = c("fixed", "fve_stable", "fve_boundary")) {
  d <- as.integer(d)
  profile <- match.arg(profile)
  if (d < 5L) stop("d must be at least five.")

  if (profile == "fixed") {
    base <- c(4.5, 2.6, 0.80, 0.45, 0.28, 0.18, 0.12, 0.07)
  } else if (profile == "fve_stable") {
    # The first two population components explain 93.40659%.
    base <- c(5.2, 3.3, 0.20, 0.14, 0.10, 0.07, 0.05, 0.04)
  } else {
    # The first two population components sit essentially on the 90% cutoff.
    base <- c(5.0, 3.1, 0.32, 0.22, 0.14, 0.09, 0.07, 0.06)
    base[2L] <- 9 * sum(base[3:length(base)]) - base[1L]
  }

  if (d <= length(base)) return(base[seq_len(d)])
  c(base, tail(base, 1L) * 0.72^seq_len(d - length(base)))
}

dc_standard_t <- function(n, d, df) {
  raw <- matrix(stats::rt(n * d, df = df), nrow = n, ncol = d)
  raw / sqrt(df / (df - 2))
}

dc_correlated_innovations <- function(n, d, distribution = c("gaussian", "t5", "t3"),
                                      rho = 0.35) {
  distribution <- match.arg(distribution)
  raw <- switch(
    distribution,
    gaussian = matrix(stats::rnorm(n * d), nrow = n, ncol = d),
    t5 = dc_standard_t(n, d, df = 5),
    t3 = dc_standard_t(n, d, df = 3)
  )
  raw %*% chol(dc_make_correlation(d, rho = rho))
}

dc_var_transition <- function(d) {
  A <- diag(0.28, d)
  for (j in seq_len(d - 1L)) {
    A[j, j + 1L] <- 0.14
    A[j + 1L, j] <- -0.08
  }
  A[1L, d] <- 0.06
  A[d, 1L] <- -0.04
  spectral_radius <- max(Mod(eigen(A, only.values = TRUE)$values))
  if (spectral_radius >= 0.80) A <- A * (0.75 / spectral_radius)
  A
}

dc_apply_var <- function(innovations, A, burn = 500L) {
  n <- nrow(innovations)
  d <- ncol(innovations)
  if (n <= burn) stop("innovations must include observations after burn-in.")
  out <- matrix(0, nrow = n, ncol = d)
  for (tt in 2:n) out[tt, ] <- A %*% out[tt - 1L, ] + innovations[tt, ]
  out[(burn + 1L):n, , drop = FALSE]
}

# Hong et al. (2024) use a bivariate VAR with 0.5 on the diagonal and 0.1
# cross-dynamics. This banded extension keeps those local coefficients when
# d > 2 and remains stable for every dimension used in this project.
dc_hong_var_transition <- function(d, cross = 0.10, diagonal = 0.50) {
  d <- as.integer(d)
  A <- diag(diagonal, d)
  if (d > 1L && cross != 0) {
    for (j in seq_len(d - 1L)) {
      A[j, j + 1L] <- cross
      A[j + 1L, j] <- cross
    }
  }
  spectral_radius <- max(Mod(eigen(A, only.values = TRUE)$values))
  if (spectral_radius >= 0.95) {
    stop("Hong-type VAR transition is not sufficiently stable.")
  }
  A
}

dc_equicorrelation <- function(d, rho = 0.10) {
  d <- as.integer(d)
  out <- matrix(rho, nrow = d, ncol = d)
  diag(out) <- 1
  if (min(eigen(out, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Requested equicorrelation matrix is not positive definite.")
  }
  out
}

dc_hong_gaussian_shocks <- function(n, d, rho = 0.10) {
  raw <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  raw %*% chol(dc_equicorrelation(d, rho = rho))
}

# DGP3 in Hong et al. (2024): componentwise GARCH(1,1) with
# (alpha, beta) = (0.10, 0.79) and contemporaneously correlated Gaussian
# standardized shocks. The intercept makes every marginal unconditional
# innovation variance equal to one.
dc_hong_garch_innovations <- function(n, d, alpha = 0.10, beta = 0.79,
                                      rho = 0.10) {
  if (alpha <= 0 || beta < 0 || alpha + beta >= 1) {
    stop("GARCH parameters must satisfy alpha > 0, beta >= 0, alpha + beta < 1.")
  }
  standardized <- dc_hong_gaussian_shocks(n, d, rho = rho)
  variance <- matrix(1, nrow = n, ncol = d)
  innovations <- matrix(0, nrow = n, ncol = d)
  innovations[1L, ] <- standardized[1L, ]
  intercept <- 1 - alpha - beta
  for (tt in 2:n) {
    variance[tt, ] <- intercept + alpha * innovations[tt - 1L, ]^2 +
      beta * variance[tt - 1L, ]
    innovations[tt, ] <- sqrt(variance[tt, ]) * standardized[tt, ]
  }
  innovations
}

# A moderate stochastic-volatility design for the assumption-aligned panel.
# Volatility states are component-specific; the old common-factor rho=.96
# design remains unchanged below as an explicit misspecification stress case.
dc_moderate_sv_innovations <- function(n, d, rho_sv = 0.90,
                                       sigma_eta = 0.15, shock_rho = 0.10) {
  if (rho_sv <= -1 || rho_sv >= 1 || sigma_eta <= 0) {
    stop("SV parameters require |rho_sv| < 1 and sigma_eta > 0.")
  }
  log_variance <- matrix(0, nrow = n, ncol = d)
  stationary_variance <- sigma_eta^2 / (1 - rho_sv^2)
  log_variance[1L, ] <- stats::rnorm(d, sd = sqrt(stationary_variance))
  for (tt in 2:n) {
    log_variance[tt, ] <- rho_sv * log_variance[tt - 1L, ] +
      stats::rnorm(d, sd = sigma_eta)
  }
  # E[scale^2] = 1 under the stationary Gaussian log-volatility law.
  scale <- exp(0.5 * log_variance - 0.25 * stationary_variance)
  dc_hong_gaussian_shocks(n, d, rho = shock_rho) * scale
}

dc_hong_t5_shocks <- function(n, d, rho = 0.10) {
  dc_standard_t(n, d, df = 5) %*% chol(dc_equicorrelation(d, rho = rho))
}

dc_functional_basis_rotation <- function(d) {
  if (d == 8L && length(dc_frozen_rotation_path) == 1L &&
      file.exists(dc_frozen_rotation_path)) {
    q <- as.matrix(utils::read.csv(dc_frozen_rotation_path))
    if (!identical(dim(q), c(8L, 8L)) || any(!is.finite(q)) ||
        max(abs(crossprod(q) - diag(8))) > 1e-12) {
      stop("The archived functional rotation is not an orthogonal 8 by 8 matrix.")
    }
    return(q)
  }
  qr.Q(qr(matrix(sin(seq_len(d * d) * 0.73), nrow = d)))
}

# Component-structured functional DGPs. The latent coordinates are population
# FPCA scores with pairwise zero cross-covariance at every lead and lag. A
# fixed orthogonal basis makes the LRV operator non-diagonal on the observed
# function grid while preserving the componentwise score condition used by
# adjusted-range normalization.
dc_component_ar_transition <- function(d, lower = 0.15, upper = 0.55) {
  diag(seq(lower, upper, length.out = d), nrow = d)
}

dc_component_process <- function(total, d, lambda,
                                 innovation = c("gaussian", "garch", "sv", "t5")) {
  innovation <- match.arg(innovation)
  transition <- if (innovation %in% c("garch", "sv")) {
    diag(0.35, d)
  } else dc_component_ar_transition(d)
  phi <- diag(transition)
  shocks <- switch(
    innovation,
    gaussian = matrix(stats::rnorm(total * d), nrow = total, ncol = d),
    garch = dc_hong_garch_innovations(
      total, d, alpha = 0.10, beta = 0.79, rho = 0
    ),
    sv = dc_moderate_sv_innovations(
      total, d, rho_sv = 0.90, sigma_eta = 0.15, shock_rho = 0
    ),
    t5 = dc_standard_t(total, d, df = 5)
  )
  marginal_scale <- sqrt(lambda * (1 - phi^2))
  scaled <- sweep(shocks, 2L, marginal_scale, FUN = "*")
  list(innovations = scaled, transition = transition)
}

dc_garch_scales <- function(n, d, omega = 0.05, alpha = 0.08, beta = 0.87) {
  common_shock <- stats::rnorm(n)
  idiosyncratic_shock <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  driver <- 0.55 * common_shock + sqrt(1 - 0.55^2) * idiosyncratic_shock
  h <- matrix(1, nrow = n, ncol = d)
  for (tt in 2:n) {
    h[tt, ] <- omega + alpha * driver[tt - 1L, ]^2 * h[tt - 1L, ] + beta * h[tt - 1L, ]
  }
  sqrt(h)
}

dc_sv_scales <- function(n, d, rho = 0.96, sigma_eta = 0.18) {
  common <- numeric(n)
  idiosyncratic <- matrix(0, nrow = n, ncol = d)
  common_sd <- sigma_eta / sqrt(1 - rho^2)
  common[1L] <- stats::rnorm(1L, sd = common_sd)
  idiosyncratic[1L, ] <- stats::rnorm(d, sd = 0.55 * common_sd)
  for (tt in 2:n) {
    common[tt] <- rho * common[tt - 1L] + stats::rnorm(1L, sd = sigma_eta)
    idiosyncratic[tt, ] <- 0.90 * idiosyncratic[tt - 1L, ] +
      stats::rnorm(d, sd = 0.55 * sigma_eta)
  }
  log_variance <- 0.65 * common + sqrt(1 - 0.65^2) * idiosyncratic
  # The centering keeps the unconditional scale near one.
  exp(0.5 * sweep(log_variance, 2L, colMeans(log_variance), FUN = "-"))
}

dc_simulate_null_scores <- function(n, d = 8L,
                                    design = c("gaussian_iid", "nondiag_lrv",
                                               "garch", "sv", "heavy_t5",
                                               "heavy_t3_stress",
                                               "hong_var_diag", "hong_var_cross",
                                               "hong_garch", "sv_moderate",
                                               "heavy_t5_ldl",
                                               "component_ar_lrv", "component_garch",
                                               "component_sv", "component_t5"),
                                    eigen_profile = c("fixed", "fve_stable", "fve_boundary"),
                                    burn = 500L) {
  n <- as.integer(n)
  d <- as.integer(d)
  burn <- as.integer(burn)
  design <- match.arg(design)
  eigen_profile <- match.arg(eigen_profile)
  total <- n + burn
  lambda <- dc_score_eigenvalues(d, profile = eigen_profile)
  scale_lambda <- matrix(sqrt(lambda), nrow = total, ncol = d, byrow = TRUE)

  if (design == "gaussian_iid") {
    innovations <- dc_correlated_innovations(total, d, "gaussian", rho = 0)
    scores <- innovations * scale_lambda
    scores <- scores[(burn + 1L):total, , drop = FALSE]
  } else if (design %in% c("component_ar_lrv", "component_garch",
                           "component_sv", "component_t5")) {
    innovation <- switch(
      design,
      component_ar_lrv = "gaussian",
      component_garch = "garch",
      component_sv = "sv",
      component_t5 = "t5"
    )
    component <- dc_component_process(total, d, lambda, innovation = innovation)
    scores <- dc_apply_var(component$innovations, component$transition, burn = burn)
    scores <- scores %*% dc_functional_basis_rotation(d)
  } else if (design %in% c("hong_var_diag", "hong_var_cross", "hong_garch",
                           "sv_moderate", "heavy_t5_ldl")) {
    cross <- if (design == "hong_var_diag") 0 else 0.10
    transition <- dc_hong_var_transition(d, cross = cross, diagonal = 0.50)
    raw_innovations <- switch(
      design,
      hong_var_diag = dc_hong_gaussian_shocks(total, d, rho = 0),
      hong_var_cross = dc_hong_gaussian_shocks(total, d, rho = 0.10),
      hong_garch = dc_hong_garch_innovations(
        total, d, alpha = 0.10, beta = 0.79, rho = 0.10
      ),
      sv_moderate = dc_moderate_sv_innovations(
        total, d, rho_sv = 0.90, sigma_eta = 0.15, shock_rho = 0.10
      ),
      heavy_t5_ldl = dc_hong_t5_shocks(total, d, rho = 0.10)
    )
    # The common 0.5 own-lag coefficient would otherwise inflate all
    # population eigenvalues by 1/(1-0.5^2). This fixed scale preserves the
    # intended eigenvalue profile exactly in the diagonal VAR design and
    # approximately in the cross-VAR designs.
    innovation_scale <- sqrt(1 - 0.50^2) * scale_lambda
    scores <- dc_apply_var(raw_innovations * innovation_scale, transition, burn = burn)
    # A deterministic orthogonal functional basis avoids making the raw grid
    # axes coincide with the population eigenfunctions. It does not alter the
    # null or the VAR order and is estimated anew by FPCA in every replication.
    scores <- scores %*% dc_functional_basis_rotation(d)
  } else if (design == "nondiag_lrv") {
    innovations <- dc_correlated_innovations(total, d, "gaussian", rho = 0.45) * scale_lambda
    scores <- dc_apply_var(innovations, dc_var_transition(d), burn = burn)
  } else if (design == "garch") {
    innovations <- dc_correlated_innovations(total, d, "gaussian", rho = 0.40)
    innovations <- innovations * dc_garch_scales(total, d) * scale_lambda
    scores <- dc_apply_var(innovations, 0.45 * dc_var_transition(d), burn = burn)
  } else if (design == "sv") {
    innovations <- dc_correlated_innovations(total, d, "gaussian", rho = 0.40)
    innovations <- innovations * dc_sv_scales(total, d) * scale_lambda
    scores <- dc_apply_var(innovations, 0.35 * dc_var_transition(d), burn = burn)
  } else {
    distribution <- if (design == "heavy_t5") "t5" else "t3"
    innovations <- dc_correlated_innovations(total, d, distribution, rho = 0.40) * scale_lambda
    scores <- dc_apply_var(innovations, 0.55 * dc_var_transition(d), burn = burn)
  }

  # A fixed orthogonal rotation makes both lag-zero covariance and the LRV
  # non-diagonal in the observed coordinates without changing the null.
  if (!design %in% c("gaussian_iid", "hong_var_diag", "hong_var_cross",
                     "hong_garch", "sv_moderate", "heavy_t5_ldl",
                     "component_ar_lrv", "component_garch", "component_sv",
                     "component_t5")) {
    rotation <- qr.Q(qr(matrix(sin(seq_len(d * d) * 0.73), nrow = d)))
    scores <- scores %*% rotation
  }
  colnames(scores) <- paste0("grid", seq_len(d))
  scores
}
