#!/usr/bin/env Rscript
# Reconstruct exactly the paper's daily curves; never substitute return curves.
# Run from the replication-package/repository root. No monitor is re-estimated.
# Required: data.table, base R graphics; ffmpeg for --video=true.
args <- list(data_file = Sys.getenv("SPX_DATA_FILE", ""),
             out_dir = "multimedia", video = "true", frames = "288")
for (arg in commandArgs(trailingOnly = TRUE)) {
  key <- sub("=.*$", "", sub("^--", "", arg))
  val <- sub("^[^=]*=", "", arg)
  if (!key %in% names(args) || !grepl("=", arg, fixed = TRUE)) stop("Unknown argument: ", arg)
  args[[key]] <- val
}
if (!file.exists(args$data_file)) stop("Set SPX_DATA_FILE or --data_file to the licensed SPX.csv file.")
if (!args$video %in% c("true", "false")) stop("--video must be true or false")
n_frames <- as.integer(args$frames)
if (is.na(n_frames) || n_frames < 24L) stop("--frames must be at least 24")
if (args$video == "true" && !nzchar(Sys.which("ffmpeg"))) stop("ffmpeg is required for the MP4.")
source("reproducibility/code/R/intraday_interface.R")
source("reproducibility/code/R/dual_channel_monitor.R")
read_result <- function(s) read.csv(paste0("reproducibility/data/spx_continuous_", s, ".csv"))
paths <- subset(read_result("paths"), design == "continuous")
alarms <- subset(read_result("alarms"), design == "continuous")
cv <- subset(read_result("critical_values"), design == "continuous")
methods <- c("adjusted_range", "quadratic")
labels <- c("(a) Adjusted-range SN", "(b) Quadratic SN (Shao, 2010)")
alarms <- alarms[match(methods, alarms$method), ]
stopifnot(!anyNA(alarms$method), all(cv$m == 748L), all(cv$h == 500L), all(cv$K == 4L))
dates <- as.Date(subset(paths, method == methods[1])$date)
stopifnot(length(dates) == 500L, !anyDuplicated(dates), all(diff(dates) > 0))
for (j in 1:2) {
  p <- subset(paths, method == methods[j])
  first <- which(p$ratio > 1)[1]
  stopifnot(identical(as.Date(p$date), dates), first == alarms$alarm_index[j],
            p$date[first] == alarms$alarm_date[j], !any(p$ratio[format(dates, "%Y") == "2019"] > 1))
}

message("Rebuilding the 500 daily log realized-variance profiles...")
raw <- dc_read_intraday_ohlc(args$data_file)
curves <- dc_build_intraday_volatility_curves(raw, "2019-01-01", "2020-12-31",
          bin_minutes = 5L, min_day_coverage = .95, log_floor = 1e-12, spline_df = 12L)
rm(raw)
stopifnot(identical(curves$dates, dates), ncol(curves$curves) == 78L,
          all(is.finite(curves$curves)))
# A raw-data check against all four saved FPCA coordinates on all 500 dates.
saved <- readRDS("reproducibility/code/empirical/outputs/spx_score_only_range_vs_shao_analysis.rds")
fpca <- saved$fit$adjusted_range$fpca
projected <- dc_project_fpca(fpca, curves$curves)$retained
scores <- read_result("scores")
expected <- as.matrix(scores[match(as.character(dates), scores$date), paste0("PC", 1:4)])
score_error <- max(abs(projected - expected))
stopifnot(is.finite(score_error), score_error < 1e-8)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
x <- 9.5 + 6.5 * curves$grid
z <- t(curves$curves)
zlim <- range(z) + c(-.5, .5)
palette <- hcl.colors(128, "YlGnBu", rev = TRUE)
alarm_col <- "#B2182B"
facet_cols <- function(zz) {
  nr <- nrow(zz); nc <- ncol(zz)
  heights <- (zz[-nr, -nc] + zz[-1, -nc] + zz[-nr, -1] + zz[-1, -1]) / 4
  palette[pmax(1L, pmin(128L, 1L + as.integer(127 * (heights - zlim[1]) / diff(zlim))))]
}
full <- seq_along(dates)
zoom <- which(dates >= as.Date("2020-01-01") & dates <= as.Date("2020-06-30"))
date_label <- function(d) paste(as.integer(format(as.Date(d), "%d")), format(as.Date(d), "%b %Y"))

draw <- function(file, ids = full, theta = 40, phi = 25, width = 2400, height = 1400, paper = FALSE) {
  png(file, width = width, height = height, res = width / if (paper) 7.6 else 12, bg = "white")
  on.exit(dev.off())
  if (paper) layout(matrix(1:3, 3, 1), heights = c(1, 1, .18)) else
    layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), heights = c(1, .19))
  par(oma = c(0, 0, if (paper) .2 else 3.5, 0), family = "sans", cex = if (paper) 1 else .83)
  for (j in 1:2) {
    par(mar = c(2.2, .9, 3.3, .9), xpd = NA)
    zz <- z[, ids, drop = FALSE]
    pmat <- persp(x, ids, zz, theta = theta, phi = phi, expand = .58, xlim = c(9.5, 16),
                  zlim = zlim, col = facet_cols(zz), border = NA, shade = .12,
                  ltheta = -30, axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
    tr <- function(a, b, c) trans3d(a, b, c, pmat)
    floor <- zlim[1]; top <- zlim[2]
    y0 <- min(ids); y1 <- max(ids); x0 <- 9.5; x1 <- 16
    lines(tr(c(x0, x1, x1, x0, x0), c(y0, y0, y1, y1, y0), rep(floor, 5)), col = "#6C737A")
    lines(tr(c(x0, x0), c(y0, y0), c(floor, top)), col = "#6C737A")
    zticks <- pretty(zlim, 4); zticks <- zticks[zticks >= floor & zticks <= top]
    text(tr(rep(x0, length(zticks)), rep(y0, length(zticks)), zticks), labels = zticks,
         pos = 2, cex = .72, offset = .45)
    # Labels remain black. Color is reserved for height and the alarm marker.
    text(tr(c(x0, 12.75, x1), rep(y0, 3), rep(floor, 3)),
         labels = c("09:30", "12:45", "16:00"), pos = 1, cex = .72, offset = .65)
    tick_dates <- if (length(ids) > 200) as.Date(c("2019-01-02", "2019-07-01", "2020-01-02", "2020-07-01", "2020-12-31")) else
      as.Date(c("2020-01-02", "2020-03-02", "2020-05-01", "2020-06-30"))
    ticks <- match(tick_dates, dates)
    text(tr(rep(x1, length(ticks)), ticks, rep(floor, length(ticks))),
         labels = format(tick_dates, "%b %Y"), pos = 4, cex = .68, offset = .45)
    k <- alarms$alarm_index[j]
    polygon(tr(c(x0, x1, x1, x0), rep(k, 4), c(floor, floor, top, top)),
            col = adjustcolor(alarm_col, alpha.f = .12), border = NA)
    lines(tr(c(x0, x1), c(k, k), c(floor, floor)), col = alarm_col, lwd = 2)
    lines(tr(x, rep(k, length(x)), z[, k]), col = alarm_col, lwd = 2.2)
    lines(tr(c(x0, x0), c(k, k), c(floor, top)), col = alarm_col, lwd = 1, lty = 2)
    title(main = labels[j], col.main = "black", cex.main = 1.02, line = 1.65)
    mtext(paste0("First alarm: ", date_label(alarms$alarm_date[j]), "  |  k = ", k),
          side = 3, line = .35, cex = .83, col = "black")
    mtext("Intraday clock time", side = 1, line = .5, cex = .75)
  }
  par(mar = c(.5, 1, .4, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  xx <- seq(.36, .72, length.out = 129)
  rect(xx[-129], .65, xx[-1], .84, col = palette, border = NA)
  text(.34, .74, "Log realized variance", adj = 1, cex = .79)
  text(c(.36, .54, .72), .46, labels = sprintf("%.1f", c(zlim[1], mean(zlim), zlim[2])), cex = .72)
  text(.5, .18, if (paper) "Red plane and profile: first threshold crossing" else
    "Red plane and profile: first threshold crossing; not an estimated change-onset date.", cex = .78)
  if (!paper) {
  mtext("S&P 500 intraday log realized-variance profiles", side = 3, outer = TRUE,
        line = 1.8, cex = 1.15, font = 2)
  mtext(if (length(ids) > 200) "Full monitoring period: 2019-2020 | Identical curves, axes and color scale in both panels" else
    "Detail: January-June 2020 | Alarms still use the full 2019-2020 monitoring horizon", side = 3, outer = TRUE, line = .4, cex = .86)
  }
}

draw(file.path(args$out_dir, "spx_3d_full_horizon.png"))
draw(file.path(args$out_dir, "spx_3d_early_2020.png"), ids = zoom)
draw(file.path(args$out_dir, "spx_3d_early_2020_paper.png"), ids = zoom,
     width = 2000, height = 2100, paper = TRUE)
validation <- data.frame(check = c("curves", "intraday_bins", "max_score_difference", "reference_observations",
   "monitoring_observations", "adjusted_range_alarm", "quadratic_alarm", "trading_observations_between_alarms"),
   value = c(nrow(curves$curves), ncol(curves$curves), format(score_error, digits = 16),
             748, 500, alarms$alarm_date, diff(alarms$alarm_index)))
write.csv(validation, file.path(args$out_dir, "spx_3d_validation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(args$out_dir, "R_sessionInfo.txt"))

if (args$video == "true") {
  frame_dir <- tempfile("spx_3d_frames_")
  dir.create(frame_dir)
  for (i in seq_len(n_frames)) {
    # Two equal-length chapters: full horizon, followed by early-2020 detail.
    chapter_n <- ceiling(n_frames / 2)
    phase <- ((i - 1) %% chapter_n) / max(1, chapter_n - 1)
    theta <- 40 + 24 * sin(2 * pi * phase)
    draw(file.path(frame_dir, sprintf("frame_%04d.png", i)),
         ids = if (i <= chapter_n) full else zoom, theta = theta,
         phi = 25 + 5 * sin(pi * phase), width = 1920, height = 1120)
    if (i %% 24 == 0) message("Rendered ", i, "/", n_frames, " video frames")
  }
  destination <- file.path(args$out_dir, "Movie_S1_SPX_intraday_surface.mp4")
  status <- system2(Sys.which("ffmpeg"), c("-y", "-loglevel", "error", "-framerate", "12",
     "-i", shQuote(file.path(frame_dir, "frame_%04d.png")), "-c:v", "libx264", "-crf", "18",
     "-pix_fmt", "yuv420p", "-movflags", "+faststart", shQuote(destination)))
  if (status != 0L || !file.exists(destination)) stop("Video encoding failed; frames remain at ", frame_dir)
  # Only this run's uniquely created frames are removed after successful encoding.
  unlink(list.files(frame_dir, pattern = "^frame_[0-9]+\\.png$", full.names = TRUE))
  unlink(frame_dir, recursive = FALSE)
  message("Wrote ", destination)
}
print(validation)
