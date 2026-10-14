#!/usr/bin/env Rscript

root <- normalizePath(dirname(grep("^--file=", commandArgs(FALSE), value = TRUE)
  |> sub("^--file=", "", x = _)), mustWork = TRUE)
data_dir <- file.path(root, "reproducibility", "data")
table_dir <- file.path(root, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

escape_tex <- function(x) {
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

design_label <- c(
  gaussian_iid = "DGP1",
  hong_var_diag = "DGP2",
  component_garch = "DGP4",
  component_sv = "DGP5",
  component_t5 = "DGP6"
)
direction_label <- c(
  sparse_first = "PC1",
  dense_first_two = "PC1--PC2"
)

power <- read.csv(file.path(data_dir, "paired_power_differences_5000rep.csv"))
power <- subset(power, calibration == "size_adjusted")

write_power_tables <- function(data, shape, output, base_label) {
  z <- subset(data, change_shape == shape)
  deltas <- sort(unique(z$delta_standardized))
  lines <- character()
  for (index in seq_along(deltas)) {
    delta <- deltas[index]
    d <- z[z$delta_standardized == delta, ]
    d <- d[order(match(d$design, names(design_label)),
                 match(d$direction, names(direction_label))), ]
    lines <- c(
      lines,
      "\\begin{table}[!htbp]",
      "\\centering",
      "\\footnotesize",
      "\\color{revisiondeepred}",
      "\\setlength{\\tabcolsep}{2.5pt}",
      sprintf("\\caption{%s changes, standardized signal %.2f (5,000 paired replications per cell).}",
              if (shape == "abrupt") "Abrupt" else "Ramp", delta),
      sprintf("\\label{%s%s}", base_label, if (index == 1L) "" else paste0("-", index)),
      "\\begin{tabular}{llrrrrrr}",
      "\\toprule",
      "Design & Direction & \\multicolumn{3}{c}{Post-change detection} & \\multicolumn{3}{c}{Capped mean delay} \\\\",
      " & & Range & Quadratic & Difference [95\\% interval] & Range & Quadratic & Difference \\\\",
      "\\midrule"
    )
    for (row in seq_len(nrow(d))) {
      ci <- sprintf("%.4f [%.4f, %.4f]",
                    d$paired_power_difference[row],
                    d$paired_power_ci95_low[row],
                    d$paired_power_ci95_high[row])
      lines <- c(
        lines,
        sprintf(
          "%s & %s & %.4f & %.4f & %s & %.2f & %.2f & %.2f \\\\",
          design_label[[d$design[row]]],
          direction_label[[d$direction[row]]],
          d$range_postchange_detection[row],
          d$shao_postchange_detection[row],
          ci,
          d$range_restricted_mean_delay[row],
          d$shao_restricted_mean_delay[row],
          d$paired_restricted_delay_difference[row]
        )
      )
    }
    lines <- c(
      lines,
      "\\bottomrule",
      "\\end{tabular}",
      "\\begin{flushleft}\\small Range and Quadratic denote adjusted-range SN and quadratic SN. Power is the unconditional post-change alarm probability after size adjustment. Capped mean delay uses the common no-pre-alarm risk set, with missed alarms assigned 401. Positive power differences and negative delay differences favor adjusted-range SN.\\end{flushleft}",
      "\\end{table}",
      ""
    )
  }
  writeLines(lines, file.path(table_dir, output))
}

write_power_tables(power, "abrupt", "revised_power_abrupt.tex", "tab:power-abrupt")
write_power_tables(power, "ramp", "revised_power_ramp.tex", "tab:power-ramp")

null_pair <- read.csv(file.path(
  data_dir, "paired_null_differences_primary5000_stress5000.csv"
))
scenario_label <- c(
  primary_gaussian_iid = "DGP1: Gaussian IID",
  primary_hong_diag_var = "DGP2: Diagonal VAR",
  primary_component_lrv = "DGP3: Component AR",
  primary_component_garch = "DGP4: Component GARCH",
  primary_component_sv = "DGP5: Component SV",
  primary_component_t5 = "DGP6: Component $t_5$",
  primary_randomK_stable = "DGP7: Random $K$, separated",
  primary_randomK_boundary = "DGP8: Random $K$, boundary",
  stress_hong_cross_var = "DGP9: Cross VAR",
  stress_hong_cross_garch = "DGP10: Cross VAR + GARCH",
  stress_general_nondiag_lrv = "DGP11: General non-diagonal LRV",
  stress_common_persistent_garch = "DGP12: Common GARCH",
  stress_common_persistent_sv = "DGP13: Common SV",
  stress_general_t5 = "DGP14: Cholesky-mixed $t_5$",
  stress_general_t3 = "DGP15: Cholesky-mixed $t_3$"
)
lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\small",
  "\\color{revisiondeepblue}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\caption{Paired null differences (5,000 replications per design).}",
  "\\label{tab:null-paired}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{llrrrrr}",
  "\\toprule",
  "Panel & Design & \\prev{Adjusted-range SN} & \\grev{Quadratic SN} & Difference [95\\% interval] & \\prev{Adjusted-range SN only} & \\grev{Quadratic SN only} \\\\",
  "\\midrule"
)
for (row in seq_len(nrow(null_pair))) {
  ci <- sprintf("%.4f [%.4f, %.4f]",
                null_pair$paired_difference[row],
                null_pair$paired_ci95_low[row],
                null_pair$paired_ci95_high[row])
  lines <- c(
    lines,
    sprintf(
      "%s & \\prev{%s} & %.4f & %.4f & %s & %.4f & %.4f \\\\",
      if (null_pair$panel[row] == "primary") "Primary" else "Stress",
      scenario_label[[null_pair$scenario_id[row]]],
      null_pair$adjusted_range_size[row],
      null_pair$shao_size[row],
      ci,
      null_pair$range_only[row],
      null_pair$shao_only[row]
    )
  )
}
lines <- c(
  lines,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\begin{flushleft}\\small Difference is \\prev{adjusted-range SN} minus \\grev{quadratic SN}. The last two columns are discordant rejection proportions from the paired experiments.\\end{flushleft}",
  "\\end{table}"
)
writeLines(lines, file.path(table_dir, "revised_null_paired.tex"))

# Rebuild the null-size figure from the frozen method-level audit file.
null_size <- read.csv(file.path(
  data_dir, "paired_null_size_primary5000_stress5000.csv"
))
figure_dir <- file.path(root, "figures", "revised")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
method_ids <- c("adjusted_range", "shao")
method_labels <- c(
  adjusted_range = "Adjusted-range SN",
  shao = "Quadratic SN"
)
method_colors <- c(adjusted_range = "#1f5a99", shao = "#d47a17")
method_pch <- c(adjusted_range = 16, shao = 17)
scenario_order <- unique(null_size$scenario_id)
plot_labels <- c(
  primary_gaussian_iid = "DGP1: Gaussian IID",
  primary_hong_diag_var = "DGP2: Diagonal VAR",
  primary_component_lrv = "DGP3: Component AR",
  primary_component_garch = "DGP4: Component GARCH",
  primary_component_sv = "DGP5: Component SV",
  primary_component_t5 = "DGP6: Component t5",
  primary_randomK_stable = "DGP7: Random K separated",
  primary_randomK_boundary = "DGP8: Random K boundary",
  stress_hong_cross_var = "DGP9: Cross VAR",
  stress_hong_cross_garch = "DGP10: Cross VAR+GARCH",
  stress_general_nondiag_lrv = "DGP11: General LRV",
  stress_common_persistent_garch = "DGP12: Common GARCH",
  stress_common_persistent_sv = "DGP13: Common SV",
  stress_general_t5 = "DGP14: Mixed t5",
  stress_general_t3 = "DGP15: Mixed t3"
)
png(
  file.path(figure_dir, "paired_null_size_5000.png"),
  width = 2300, height = 1200, res = 180
)
par(mar = c(9.5, 4.5, 3.5, 1), xpd = NA)
plot(
  NA,
  xlim = c(0.5, length(scenario_order) + 0.5),
  ylim = c(0.025, 0.085),
  xaxt = "n", xlab = "",
  ylab = "Empirical null rejection probability",
  main = "Paired null-size comparison (5,000 replications)"
)
abline(h = 0.05, col = "#b22222", lty = 2, lwd = 1.5)
abline(v = 8.5, col = "grey75", lty = 3)
for (method in method_ids) {
  rows <- null_size[null_size$method_id == method, ]
  rows <- rows[match(scenario_order, rows$scenario_id), ]
  offset <- if (method == "adjusted_range") -0.11 else 0.11
  x <- seq_along(scenario_order) + offset
  segments(
    x, rows$ci95_low, x, rows$ci95_high,
    col = method_colors[[method]], lwd = 1.2
  )
  points(
    x, rows$empirical_size,
    col = method_colors[[method]], pch = method_pch[[method]], cex = 1.05
  )
}
axis(
  1, at = seq_along(scenario_order),
  labels = unname(plot_labels[scenario_order]), las = 2, cex.axis = 0.78
)
mtext("Primary design panel", side = 3, at = 4.5, line = 0.3, cex = 0.85)
mtext("Deliberate stress panel", side = 3, at = 12, line = 0.3, cex = 0.85)
legend(
  "topright",
  legend = c(unname(method_labels), "Nominal 5%"),
  col = c(unname(method_colors), "#b22222"),
  text.col = "black",
  pch = c(unname(method_pch), NA),
  lty = c(NA, NA, 2),
  lwd = c(NA, NA, 1.5),
  bty = "n", cex = 0.85
)
dev.off()

# Rebuild the abrupt-power panel with enough lower margin for the x-axis label.
plot_design_label <- c(
  gaussian_iid = "Gaussian IID",
  hong_var_diag = "Diagonal VAR(1)",
  component_garch = "Component GARCH",
  component_sv = "Component SV",
  component_t5 = "Component t(5)"
)
plot_power_cell <- function(rows, title, show_y, show_x,
                            x_label = "Standardized shift", title_cex = 0.92) {
  rows <- rows[order(rows$delta_standardized), ]
  plot(
    rows$delta_standardized, rows$range_postchange_detection,
    type = "b", pch = 16, col = method_colors[["adjusted_range"]],
    lwd = 1.7, cex = 0.9, ylim = c(0, 1),
    xlab = if (show_x) x_label else "",
    ylab = if (show_y) "Post-change detection probability" else "",
    main = title,
    cex.axis = 0.92, cex.lab = 0.92, cex.main = title_cex
  )
  abline(h = 0.05, col = "grey78", lty = 3)
  lines(
    rows$delta_standardized, rows$shao_postchange_detection,
    type = "b", pch = 17, col = method_colors[["shao"]],
    lwd = 1.7, cex = 0.9
  )
}
abrupt <- power[power$change_shape == "abrupt", ]
directions <- c("sparse_first", "dense_first_two")
designs <- names(plot_design_label)
png(
  file.path(figure_dir, "paired_power_abrupt_size_adjusted_5000.png"),
  width = 2600, height = 1500, res = 180, pointsize = 18
)
layout(
  matrix(c(seq_len(10), rep(11, 5)), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.18)
)
par(mar = c(4.1, 3.7, 2.4, 0.8), oma = c(0.2, 0.2, 0.2, 0.2))
for (direction_index in seq_along(directions)) {
  direction <- directions[direction_index]
  for (j in seq_along(designs)) {
    design <- designs[j]
    rows <- abrupt[abrupt$direction == direction & abrupt$design == design, ]
    direction_text <- if (direction == "sparse_first") {
      "PC1 shift"
    } else {
      "Equal PC1-PC2 shift"
    }
    plot_power_cell(
      rows, paste0(plot_design_label[[design]], "\n", direction_text),
      show_y = j == 1L, show_x = direction_index == length(directions),
      title_cex = 0.82
    )
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = c("Adjusted-range SN", "Quadratic SN", "5% reference"),
  col = c(unname(method_colors), "grey55"),
  pch = c(unname(method_pch), NA),
  lty = c(1, 1, 3), lwd = c(1.7, 1.7, 1.2),
  bty = "n", horiz = TRUE, cex = 0.92
)
dev.off()

# Rebuild the ramp-power panel on the same scale and with the same labels.
ramp <- power[power$change_shape == "ramp", ]
ramp_designs <- c("hong_var_diag", "component_garch", "component_t5")
png(
  file.path(figure_dir, "paired_power_ramp_size_adjusted_5000.png"),
  width = 2600, height = 1050, res = 180, pointsize = 18
)
layout(matrix(c(seq_len(3), rep(4, 3)), nrow = 2, byrow = TRUE),
       heights = c(1, 0.20))
par(mar = c(4.3, 4.0, 2.5, 0.8), oma = c(0.2, 0.2, 0.2, 0.2))
for (j in seq_along(ramp_designs)) {
  design <- ramp_designs[j]
  rows <- ramp[ramp$design == design, ]
  plot_power_cell(
    rows, plot_design_label[[design]], show_y = j == 1L, show_x = TRUE,
    x_label = "Terminal standardized shift"
  )
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = c("Adjusted-range SN", "Quadratic SN", "5% reference"),
  col = c(unname(method_colors), "grey55"),
  pch = c(unname(method_pch), NA),
  lty = c(1, 1, 3), lwd = c(1.7, 1.7, 1.2),
  bty = "n", horiz = TRUE, cex = 0.92
)
dev.off()
