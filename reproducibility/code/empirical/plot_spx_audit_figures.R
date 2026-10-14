# Presentation-only plotters. Source from the empirical runners or regenerate
# both manuscript figures from frozen CSVs with:
# Rscript reproducibility/code/empirical/plot_spx_audit_figures.R

plot_spx_demixed_loadings <- function(u, phi, g, output_file) {
  stopifnot(ncol(phi) == 4L, identical(dim(phi), dim(g)), length(u) == nrow(phi))
  grDevices::png(output_file, width = 2000, height = 1500, res = 200)
  on.exit(grDevices::dev.off())
  # A dedicated fifth panel keeps the legend outside both lower axis labels.
  graphics::layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
                   heights = c(1, 1, 0.18))
  graphics::par(mar = c(4.5, 4, 2.5, 1), oma = c(0, 0, 0, 0), cex = 1)
  for (j in seq_len(4L)) {
    graphics::matplot(u, cbind(phi[, j], g[, j]), type = "l", lty = c(1, 2),
      lwd = 2, col = c("#202020", "#0072B2"),
      xlab = "Fraction of regular session", ylab = "Loading",
      main = paste("Direction", j))
    graphics::abline(h = 0, col = "grey75", lty = 3)
  }
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend("center", c("Original FPCA eigenfunction", "LDL-transformed loading"),
    col = c("#202020", "#0072B2"), lty = c(1, 2), lwd = 2,
    horiz = TRUE, bty = "n")
  invisible(output_file)
}

plot_spx_monitor <- function(dates, normalized_range, normalized_shao, output_file) {
  stopifnot(length(dates) == length(normalized_range),
            length(dates) == length(normalized_shao))
  grDevices::png(output_file, width = 1800, height = 900, res = 180)
  on.exit(grDevices::dev.off())
  graphics::matplot(dates, cbind(normalized_range, normalized_shao), type = "l",
    lty = c(1, 1), lwd = c(1.8, 1.6), col = c("#1f5a99", "#d47a17"),
    xlab = "Date", ylab = "Statistic / method-specific 5% critical value",
    main = "S&P 500 monitoring of retained FPC scores", col.main = "black")
  graphics::abline(h = 1, col = "#b22222", lty = 2, lwd = 1.5)
  graphics::legend("topleft",
    legend = c("Adjusted-range SN", "Quadratic SN", "alarm threshold"),
    col = c("#1f5a99", "#d47a17", "#b22222"), text.col = "black",
    lty = c(1, 1, 2), lwd = c(1.8, 1.6, 1.5), bty = "n")
  invisible(output_file)
}

if (sys.nframe() == 0L) {
  loadings <- read.csv("reproducibility/data/spx_original_and_demixed_loadings.csv")
  plot_spx_demixed_loadings(loadings$u, as.matrix(loadings[paste0("phi", 1:4)]),
    as.matrix(loadings[paste0("X", 1:4)]),
    "figures/revised/spx_original_and_demixed_loadings.png")
  path <- read.csv("reproducibility/data/spx_score_only_range_vs_shao_path.csv")
  plot_spx_monitor(as.Date(path$date), path$adjusted_range_threshold_ratio,
    path$shao_threshold_ratio, "figures/revised/spx_score_only_range_vs_shao_monitor.png")
}
