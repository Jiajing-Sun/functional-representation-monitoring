file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
} else {
  normalizePath("JTSA_method_rebuild/tests/test_dual_channel_monitor.R", mustWork = TRUE)
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(project_dir, "R", "dual_channel_monitor.R"))

message("1. LDL factorization and Hong-type demixing")
set.seed(101)
mixing <- matrix(c(1, 0, 0, 0.7, 1, 0, -0.4, 0.5, 1), 3, 3, byrow = TRUE)
latent <- matrix(rnorm(1200), ncol = 3) %*% diag(c(1, 2, 0.5))
observed <- latent %*% t(mixing)
sigma <- cov(observed)
ldl <- dc_ldl_factor(sigma, ridge = 1e-12)
stopifnot(ldl$relative_reconstruction_error < 1e-10)
pw <- dc_fit_ldl_prewhitener(observed, ridge = 1e-12, mode = "demix")
demixed <- dc_apply_ldl_prewhitener(pw, observed)
offdiag <- cov(demixed) - diag(diag(cov(demixed)))
stopifnot(max(abs(offdiag)) < 1e-8)

message("2. Quadratic SN affine invariance")
set.seed(102)
calibration <- matrix(rnorm(900), nrow = 300, ncol = 3)
monitoring <- matrix(rnorm(240), nrow = 80, ncol = 3)
affine_map <- matrix(c(1.2, 0.2, -0.1, 0, 0.8, 0.3, 0, 0, 1.5), 3, 3, byrow = TRUE)
sn_original <- dc_fit_full_matrix_sn(calibration, ridge = 0)
sn_mapped <- dc_fit_full_matrix_sn(calibration %*% t(affine_map), ridge = 0)
path_original <- dc_sn_monitor_path(sn_original, monitoring, gamma = 0)
path_mapped <- dc_sn_monitor_path(sn_mapped, monitoring %*% t(affine_map), gamma = 0)
stopifnot(max(abs(path_original - path_mapped)) < 1e-5)

message("2b. Exact quadratic-SN rank boundary and componentwise wide-dimension fit")
rank_boundary_features <- rbind(diag(4), rep(-1, 4))
rank_boundary_fit <- dc_fit_full_matrix_sn(rank_boundary_features, ridge = 0)
stopifnot(
  rank_boundary_fit$dimension == 4L,
  all(is.finite(rank_boundary_fit$inverse))
)
rank_failure <- tryCatch(
  dc_fit_full_matrix_sn(matrix(seq_len(25), nrow = 5, ncol = 5), ridge = 0),
  error = function(e) e
)
stopifnot(
  inherits(rank_failure, "error"),
  grepl("feature dimension <= calibration size - 1", conditionMessage(rank_failure), fixed = TRUE)
)
regularized_rank_fit <- dc_fit_full_matrix_sn(
  matrix(seq_len(25), nrow = 5, ncol = 5), ridge = 1e-8
)
stopifnot(
  regularized_rank_fit$dimension == 5L,
  regularized_rank_fit$ridge_added > 0,
  all(is.finite(regularized_rank_fit$inverse))
)
wide_range_fit <- dc_fit_adjusted_range_sn(
  matrix(seq_len(30), nrow = 5, ncol = 6), ridge = 0
)
stopifnot(
  wide_range_fit$dimension == 6L,
  all(is.finite(wide_range_fit$adjusted_ranges)),
  all(wide_range_fit$adjusted_ranges > 0)
)

message("3. Fixed-FPCA projection and orthogonal-loss identity")
set.seed(103)
grid <- seq(0, 1, length.out = 81)
weights <- dc_grid_weights(length(grid), grid = grid)
raw_basis <- cbind(
  rep(1, length(grid)),
  sqrt(2) * sin(2 * pi * grid),
  sqrt(2) * cos(2 * pi * grid),
  sqrt(2) * sin(4 * pi * grid),
  sqrt(2) * cos(4 * pi * grid),
  sqrt(2) * sin(6 * pi * grid)
)
weighted_basis <- sweep(raw_basis, 1L, sqrt(weights), FUN = "*")
orthogonal_weighted <- qr.Q(qr(weighted_basis))
functional_basis <- sweep(orthogonal_weighted, 1L, sqrt(weights), FUN = "/")
coefficients <- matrix(rnorm(500 * 6), 500, 6) %*% diag(sqrt(c(9, 5, 3, 1.5, 0.8, 0.4)))
curves <- coefficients %*% t(functional_basis)

fpca <- dc_fit_fpca(curves[1:180, ], K = 3, L = 1, grid = grid)
projected <- dc_project_fpca(fpca, curves[181:260, ])
reconstructed_loss <- rowSums(projected$audit^2) + projected$far_tail_energy
stopifnot(max(abs(reconstructed_loss - projected$total_reconstruction_loss)) < 1e-10)
fpca_capped <- dc_fit_fpca(
  curves[1:180, ], K = NULL, L = 1, fve = 0.90, K_max = 5, grid = grid
)
stopifnot(fpca_capped$K <= 5L, fpca_capped$K_max == 5L)
cap_failure <- try(
  dc_fit_fpca(curves[1:180, ], K = NULL, L = 1, fve = 0.999, K_max = 1, grid = grid),
  silent = TRUE
)
stopifnot(inherits(cap_failure, "try-error"))

message("4. End-to-end two-channel fit and monitoring")
fit <- dc_fit_dual_channel(
  curves[1:180, ], curves[181:360, ], K = 3, L = 1,
  grid = grid, prewhiten_mode = "demix", sn_ridge = 1e-8
)
stopifnot(fit$sn_S$dimension == 3L)
stopifnot(fit$sn_R$dimension == 6L)
result <- dc_monitor(
  fit, curves[361:500, ],
  critical_values = c(S = 1e9, R = 1e9), gamma = 0
)
stopifnot(nrow(result$path) == 140L)
stopifnot(all(is.finite(result$path$S)), all(is.finite(result$path$R)))
stopifnot(result$alarm$label == "none")

message("5. Formal channel labels are assigned only at the global stopping time")
stopifnot(dc_channel_label(c(1, 6, 7), c(1, 2, 3), 5, 5)$label == "S")
stopifnot(dc_channel_label(c(1, 2, 3), c(1, 6, 7), 5, 5)$label == "R")
stopifnot(dc_channel_label(c(1, 6, 7), c(1, 8, 9), 5, 5)$label == "both")
stopifnot(dc_channel_label(c(1, 2, 3), c(1, 2, 3), 5, 5)$label == "none")

message("6. VAR--LDL feature residualization and normalization-only Huberization")
set.seed(104)
feature_process <- matrix(0, nrow = 240, ncol = 4)
innovations <- matrix(rt(240 * 4, df = 5), nrow = 240, ncol = 4)
for (tt in 2:240) feature_process[tt, ] <- 0.55 * feature_process[tt - 1L, ] + innovations[tt, ]
feature_process[120, 3] <- 100
huber <- dc_fit_huberizer(feature_process[1:160, ], threshold = 4.5)
stopifnot(max(huber$transformed[, 3]) < max(feature_process[1:160, 3]))
var_fit <- dc_fit_var_ldl_features(huber$transformed, order = 1L)
var_future <- dc_apply_var_ldl_features(
  var_fit$model,
  dc_apply_huberizer(huber$model, feature_process[161:240, , drop = FALSE])
)
stopifnot(nrow(var_fit$transformed) == 159L, nrow(var_future) == 80L)
stopifnot(all(is.finite(var_fit$transformed)), all(is.finite(var_future)))

message("7. Componentwise adjusted-range normalizer matches its exact formula")
range_fit <- dc_fit_adjusted_range_sn(calibration, ridge = 0)
stopifnot(all(eigen(range_fit$regularized_self_normalizer, symmetric = TRUE)$values > 0))
centered_calibration <- sweep(calibration, 2L, colMeans(calibration), FUN = "-")
calibration_bridge <- rbind(0, apply(centered_calibration, 2L, cumsum) / sqrt(nrow(calibration)))
direct_ranges <- apply(calibration_bridge, 2L, max) - apply(calibration_bridge, 2L, min)
stopifnot(
  is.null(range_fit$ldl_C),
  range_fit$ridge_added == 0,
  max(abs(range_fit$adjusted_ranges - direct_ranges)) < 1e-12
)
range_cv <- dc_simulate_adjusted_range_critical_value(
  dimension = 2L, T = 1, alpha = 0.025, nrep = 20L,
  train_grid = 50L, seed = 105L
)
stopifnot(is.finite(range_cv$critical_value), range_cv$critical_value > 0)

message("8. Increment-hull gauge has scalar reduction and affine invariance")
set.seed(106)
G <- matrix(rnorm(303), nrow = 101, ncol = 3)
G <- sweep(G, 2L, colMeans(G), FUN = "-")
z <- rnorm(3)
H <- matrix(c(1.2, 0.2, -0.1, 0.1, 0.9, 0.3, 0, 0.2, 1.1),
            nrow = 3, byrow = TRUE)
gauge_original <- dc_ae_gauge_lp(G, z)
gauge_mapped <- dc_ae_gauge_lp(G %*% t(H), as.numeric(H %*% z))
stopifnot(abs(gauge_original - gauge_mapped) < 1e-9)
scalar_training <- matrix(rnorm(80), ncol = 1L)
scalar_fit <- dc_fit_affine_range_sn(scalar_training)
scalar_monitor <- matrix(rnorm(30), ncol = 1L)
scalar_path <- dc_affine_range_monitor_path(scalar_fit, scalar_monitor)
training_cs <- cumsum(scalar_training[, 1L] - mean(scalar_training[, 1L])) / sqrt(80)
training_range <- max(c(0, training_cs)) - min(c(0, training_cs))
expected_path <- (cumsum(scalar_monitor[, 1L] - mean(scalar_training[, 1L])) /
                    sqrt(80) / training_range)^2 /
  dc_boundary_sq(seq_len(30) / 80)
stopifnot(max(abs(scalar_path - expected_path)) < 1e-10)

message("9. Restored score-only LDL adjusted-range monitor")
score_fit <- dc_fit_score_monitor(
  x_fit = curves[1:180, ],
  x_calibration = curves[181:360, ],
  K = 3L,
  grid = grid,
  prewhiten_mode = "demix",
  sn_type = "adjusted_range",
  feature_var_order = 0L
)
score_result <- dc_monitor_score(
  score_fit,
  curves[361:500, ],
  critical_value = 1e9,
  gamma = 0
)
stopifnot(
  inherits(score_fit, "dc_score_monitor_fit"),
  score_fit$fpca$K == 3L,
  score_fit$sn$dimension == 3L,
  nrow(score_result$path) == 140L,
  all(is.finite(score_result$path$statistic)),
  is.na(score_result$alarm$stop)
)

message("10. Paired adjusted-range and quadratic-SN monitors share one score stream")
comparison_fit <- dc_fit_score_comparison(
  x_fit = curves[1:180, ],
  x_calibration = curves[181:360, ],
  K = 3L,
  grid = grid,
  prewhiten_mode = "demix",
  feature_var_order = 0L
)
comparison_result <- dc_monitor_score_comparison(
  comparison_fit,
  curves[361:500, ],
  critical_values = c(adjusted_range = 1e9, shao = 1e9),
  gamma = 0
)
stopifnot(
  inherits(comparison_fit, "dc_score_comparison_fit"),
  comparison_fit$adjusted_range$fpca$K == comparison_fit$shao$fpca$K,
  identical(
    comparison_fit$adjusted_range$calibration_sn_scores,
    comparison_fit$shao$calibration_sn_scores
  ),
  comparison_fit$shao$sn$normalizer_type == "full_matrix_sequential",
  is.null(comparison_fit$adjusted_range$sn$ldl_C),
  comparison_fit$adjusted_range$sn$ridge_added == 0,
  comparison_fit$shao$sn$ridge_added == 0,
  nrow(comparison_result$path) == 140L,
  all(is.finite(comparison_result$path$adjusted_range)),
  all(is.finite(comparison_result$path$shao)),
  all(is.na(comparison_result$alarm))
)

message("All monitoring checks passed.")
