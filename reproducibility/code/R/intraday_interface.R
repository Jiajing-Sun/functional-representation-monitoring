# Interface from minute OHLC data to daily intraday-volatility curves.

dc_read_intraday_ohlc <- function(path) {
  if (!file.exists(path)) stop("Intraday data file does not exist: ", path)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required to read the intraday file efficiently.")
  }
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  required <- c("DateTime", "Open", "Close")
  missing <- setdiff(required, header)
  if (length(missing)) stop("Intraday file is missing columns: ", paste(missing, collapse = ", "))
  data <- data.table::fread(
    path,
    select = required,
    colClasses = list(character = "DateTime", numeric = c("Open", "Close")),
    showProgress = FALSE
  )
  if (!nrow(data)) stop("Intraday file contains no observations.")
  if (any(!is.finite(data$Open)) || any(!is.finite(data$Close)) ||
      any(data$Open <= 0) || any(data$Close <= 0)) {
    stop("Open and Close must be positive finite prices.")
  }
  data
}

dc_interpolate_curve <- function(values) {
  observed <- which(is.finite(values))
  if (!length(observed)) return(values)
  if (length(observed) == 1L) return(rep(values[observed], length(values)))
  stats::approx(observed, values[observed], xout = seq_along(values),
                method = "linear", rule = 2)$y
}

dc_smooth_intraday_curves <- function(curves, grid, spline_df = 12L) {
  spline_df <- as.integer(spline_df)
  if (spline_df < 4L || spline_df >= ncol(curves)) {
    stop("spline_df must be at least four and smaller than the number of bins.")
  }
  basis <- splines::bs(grid, df = spline_df, degree = 3L, intercept = TRUE)
  gram_inverse <- solve(crossprod(basis))
  smoother <- basis %*% gram_inverse %*% t(basis)
  smoothed <- curves %*% smoother
  dimnames(smoothed) <- dimnames(curves)
  smoothed
}

dc_build_intraday_volatility_curves <- function(data, start_date = NULL, end_date = NULL,
                                                bin_minutes = 5L,
                                                min_day_coverage = 0.95,
                                                log_floor = 1e-12,
                                                spline_df = 12L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  required <- c("DateTime", "Open", "Close")
  if (!all(required %in% names(data))) stop("data must contain DateTime, Open, and Close.")
  bin_minutes <- as.integer(bin_minutes)
  if (bin_minutes < 1L || 390L %% bin_minutes != 0L) {
    stop("bin_minutes must divide the 390-minute regular US session.")
  }
  if (!is.numeric(min_day_coverage) || min_day_coverage <= 0 || min_day_coverage > 1) {
    stop("min_day_coverage must lie in (0,1].")
  }

  dt <- data.table::as.data.table(data)[, .(
    DateTime = as.character(DateTime),
    Open = as.numeric(Open),
    Close = as.numeric(Close)
  )]
  dt[, date := data.table::as.IDate(substr(DateTime, 1L, 10L))]
  dt[, clock := substr(DateTime, 12L, 19L)]
  dt[, hour := as.integer(substr(clock, 1L, 2L))]
  dt[, minute := as.integer(substr(clock, 4L, 5L))]
  dt[, minute_index := hour * 60L + minute - (9L * 60L + 31L) + 1L]
  dt <- dt[minute_index >= 1L & minute_index <= 390L]
  if (!is.null(start_date)) dt <- dt[date >= data.table::as.IDate(start_date)]
  if (!is.null(end_date)) dt <- dt[date <= data.table::as.IDate(end_date)]
  if (!nrow(dt)) stop("No regular-session observations remain after date filtering.")

  data.table::setorder(dt, date, minute_index, DateTime)
  dt <- dt[, .SD[.N], by = .(date, minute_index)]
  dt[, previous_close := data.table::shift(Close), by = date]
  dt[, previous_minute := data.table::shift(minute_index), by = date]
  dt[, log_return := data.table::fifelse(
    minute_index == 1L,
    log(Close / Open),
    data.table::fifelse(previous_minute == minute_index - 1L, log(Close / previous_close), NA_real_)
  )]
  dt[, bin := (minute_index - 1L) %/% bin_minutes + 1L]

  n_bins <- 390L %/% bin_minutes
  binned <- dt[, .(
    rv = if (sum(is.finite(log_return)) == bin_minutes) sum(log_return^2) else NA_real_,
    minute_returns = sum(is.finite(log_return))
  ), by = .(date, bin)]
  all_bins <- data.table::CJ(date = sort(unique(dt$date)), bin = seq_len(n_bins), unique = TRUE)
  binned <- binned[all_bins, on = .(date, bin)]
  day_qc <- binned[, .(
    observed_bins = sum(is.finite(rv)),
    total_bins = n_bins,
    coverage = mean(is.finite(rv)),
    interpolated_bins = sum(!is.finite(rv)),
    zero_or_near_zero_bins = sum(is.finite(rv) & rv <= log_floor)
  ), by = date]
  kept_dates <- day_qc[coverage >= min_day_coverage, date]
  if (!length(kept_dates)) stop("No day meets the requested intraday coverage threshold.")
  binned <- binned[date %in% kept_dates]

  curves <- matrix(NA_real_, nrow = length(kept_dates), ncol = n_bins)
  for (day_index in seq_along(kept_dates)) {
    values <- binned[date == kept_dates[day_index]][order(bin)]$rv
    values <- dc_interpolate_curve(values)
    curves[day_index, ] <- log(values + log_floor)
  }
  colnames(curves) <- sprintf("rv_%03d_%03d",
                              (seq_len(n_bins) - 1L) * bin_minutes,
                              seq_len(n_bins) * bin_minutes)
  rownames(curves) <- as.character(kept_dates)
  raw_log_curves <- curves
  grid <- (seq_len(n_bins) - 0.5) / n_bins
  curves <- dc_smooth_intraday_curves(curves, grid = grid, spline_df = spline_df)
  kept_qc <- day_qc[match(kept_dates, date)]

  list(
    curves = curves,
    raw_log_curves = raw_log_curves,
    dates = as.Date(kept_dates),
    grid = grid,
    quality = as.data.frame(day_qc),
    kept_quality = as.data.frame(kept_qc),
    settings = list(
      bin_minutes = bin_minutes,
      n_bins = n_bins,
      min_day_coverage = min_day_coverage,
      log_floor = log_floor,
      spline_df = as.integer(spline_df),
      session = "09:31--16:00"
    )
  )
}

dc_split_intraday_curves <- function(curve_data, fit_dates, calibration_dates, monitor_dates) {
  dates <- curve_data$dates
  select_period <- function(period) {
    which(dates >= as.Date(period[1L]) & dates <= as.Date(period[2L]))
  }
  fit_index <- select_period(fit_dates)
  calibration_index <- select_period(calibration_dates)
  monitor_index <- select_period(monitor_dates)
  if (length(fit_index) < 20L || length(calibration_index) < 20L || !length(monitor_index)) {
    stop("The requested fit/normalization/monitor periods do not contain enough complete days.")
  }
  if (length(intersect(fit_index, calibration_index)) ||
      length(intersect(fit_index, monitor_index)) ||
      length(intersect(calibration_index, monitor_index))) {
    stop("Fit, normalization, and monitoring periods must not overlap.")
  }
  list(
    x_fit = curve_data$curves[fit_index, , drop = FALSE],
    x_calibration = curve_data$curves[calibration_index, , drop = FALSE],
    x_monitor = curve_data$curves[monitor_index, , drop = FALSE],
    fit_dates = dates[fit_index],
    calibration_dates = dates[calibration_index],
    monitor_dates = dates[monitor_index],
    grid = curve_data$grid
  )
}
