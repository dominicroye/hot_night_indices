# nightheat.R
# Functions to compute hot-night duration (HNd) and hot-night excess (HNe)
# from hourly temperature series, based on astronomical sunset/sunrise windows.
# Vectorised implementation using data.table::foverlaps.
#
# Dependencies: lubridate, data.table, suncalc
# Optional (parallel percentile): future, future.apply
#
# Reference:
#   Royé D, Sera F, Tobías A, et al. (2025). Short-term association between
#   hot nights and mortality: a multicountry analysis in 178 locations
#   considering hourly ambient temperature. Environment International.
#   https://doi.org/10.1016/j.envint.2025.109526
#
# Author: Dominic Royé (@dominicroye) · https://dominicroye.github.io
# ─────────────────────────────────────────────────────────────────────────────

library(lubridate)
library(data.table)
library(suncalc)


# ── Internal helpers ──────────────────────────────────────────────────────────

# Standardise input: select columns, rename, adjust time zone
.prep_obs <- function(x, cols, tz) {
  dt <- as.data.table(x)[, .SD, .SDcols = cols]
  setnames(dt, c("Datetime", "Ta"))
  dt[, Datetime := with_tz(Datetime, tz)]
  dt
}

# Compute sunset / sunrise and define each nocturnal window [sunset_i, sunrise_i+1]
.get_sun_times <- function(dates, coord, tz) {
  srs <- as.data.table(getSunlightTimes(
    date = sort(unique(as.Date(dates))),
    lat  = coord[2], lon = coord[1], tz = tz,
    keep = c("sunrise", "sunset")
  ))[, .(date, sunset, sunrise)]
  srs[, night_end := shift(sunrise, type = "lead")]
  srs[, night_hrs := abs(as.numeric(difftime(
    round_date(night_end, "hour"),
    round_date(sunset,    "hour"),
    units = "hours"
  )))]
  srs[!is.na(night_end)]
}

# Join observations to their nocturnal window using foverlaps (interval overlap
# join). Returns only rows that fall inside a night window.
# Note: foverlaps uses closed intervals [start, end]. In hourly series, exact
# coincidence with the sunset/sunrise second is negligible.
.night_join <- function(tab, srs) {
  obs <- copy(tab)[, dt_end := Datetime]   # degenerate interval: point
  setkey(obs, Datetime, dt_end)

  nts <- copy(srs)[, .(date, night_hrs, start = sunset, end = night_end)]
  setkey(nts, start, end)

  foverlaps(obs, nts,
            by.x = c("Datetime", "dt_end"),
            by.y = c("start",    "end"),
            nomatch = NA)[!is.na(date)]
}


# ── Hot-night duration (HNd) ──────────────────────────────────────────────────
#
# Returns a data.table with:
#   date  - date of the sunset that opens the night
#   HNd   - percentage of night hours with Ta > thres
#
# Arguments:
#   x      data.frame / data.table with at least two columns
#   cols   character(2): names of the datetime and temperature columns
#   tz     time-zone string, e.g. "Europe/Madrid"
#   coord  numeric(2): c(longitude, latitude) in decimal degrees (W = negative)
#   thres  temperature threshold (°C)

night_duration <- function(x, cols, tz, coord, thres) {

  tab <- .prep_obs(x, cols, tz)
  srs <- .get_sun_times(tab$Datetime, coord, tz)

  joined <- .night_join(tab, srs)

  above <- joined[Ta > thres, .(n_above = .N), by = .(date, night_hrs)]

  result <- merge(srs[, .(date, night_hrs)], above,
                  by = c("date", "night_hrs"), all.x = TRUE)
  result[is.na(n_above), n_above := 0L]
  result[, HNd := n_above * 100 / pmax(n_above, night_hrs)]

  result[, .(date, HNd)]
}


# ── Hot-night excess (HNe) ────────────────────────────────────────────────────
#
# Returns a data.table with:
#   date  - date of the sunset that opens the night
#   HNe   - sum of (Ta - thres) for hours with Ta > thres (degree-hours)
#
# Arguments: same as nighthour()

night_excess <- function(x, cols, tz, coord, thres) {

  tab <- .prep_obs(x, cols, tz)
  tab[, tadif := Ta - thres]
  srs <- .get_sun_times(tab$Datetime, coord, tz)

  joined <- .night_join(tab, srs)

  excess <- joined[Ta > thres, .(HNe = sum(tadif, na.rm = TRUE)), by = date]

  result <- merge(srs[, .(date)], excess, by = "date", all.x = TRUE)
  result[is.na(HNe), HNe := 0]

  result[, .(date, HNe)]
}


# ── Temperature percentile threshold (fixed or moving) ───────────────────────
#
# Computes the p-th percentile of daily temperature, either as a single fixed
# value per group or as a centred rolling window (moving percentile). The
# result is intended for use as `thres` in nighthour() / nightdegrees().
# Following Royé et al. (2025), the recommended threshold is the moving 95th
# percentile of daily minimum temperature with a 91-day centred window.
#
# Returns a data.table with:
#   date      - date
#   [by cols] - grouping columns (if supplied)
#   pXX       - percentile value (auto-named, e.g. "p95")
#
# Arguments:
#   x       data.frame / data.table with date and temperature columns
#   cols    character(2): names of the date and temperature columns
#   p       percentile as a proportion, e.g. 0.95
#   by      character or NULL: grouping column(s) (station, city, ...)
#   type    "fixed"  -> one value per group (percentile of the full series)
#           "moving" -> centred rolling window of `window` days
#   window  window size in days for type = "moving" (must be odd;
#           even values are rounded up to the next odd integer)
#   workers number of parallel workers (requires {future} and {future.apply});
#           workers = 1 disables parallelism (default)

night_percentile <- function(x,
                             cols,
                             p       = 0.95,
                             by      = NULL,
                             type    = c("fixed", "moving"),
                             window  = 91,
                             workers = 1L) {

  type  <- match.arg(type)
  pname <- paste0("p", round(p * 100))    # output column name, e.g. "p95"

  # Enforce odd window (required for symmetric align = "center")
  if (window %% 2 == 0) {
    window <- window + 1L
    message("night_percentile: `window` rounded up to ", window, " (odd required)")
  }

  dt <- as.data.table(x)[, .SD, .SDcols = c(by, cols)]
  setnames(dt, c(by, "date", "Ta"))
  setorderv(dt, c(by, "date"))

  # NA-safe quantile for rolling apply
  .qfun <- function(v) {
    if (anyNA(v)) return(NA_real_)
    quantile(v, probs = p, names = FALSE)
  }

  if (type == "fixed") {

    result <- dt[, .(date = date,
                     pval = quantile(Ta, probs = p, na.rm = TRUE)),
                 by = by]
    setnames(result, "pval", pname)
    return(result[])

  }

  # ── Moving percentile via frollapply ────────────────────────────────────────
  # frollapply + align = "center" replicates the original lag_array + apply +
  # manual 45-position shift, without building the N x 90 matrix in memory.

  if (workers > 1L) {

    if (!requireNamespace("future",       quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE))
      stop("Install {future} and {future.apply} to use workers > 1")

    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    groups      <- if (is.null(by)) list(dt) else split(dt, by = by)
    result_list <- future.apply::future_lapply(groups, function(g) {
      g[, (pname) := frollapply(Ta, n = window, FUN = .qfun,
                                align = "center", fill = NA)]
      g[, Ta := NULL]
      g
    })
    result <- rbindlist(result_list)

  } else {

    result <- copy(dt)
    if (is.null(by)) {
      result[, (pname) := frollapply(Ta, n = window, FUN = .qfun,
                                     align = "center", fill = NA)]
    } else {
      result[, (pname) := frollapply(Ta, n = window, FUN = .qfun,
                                     align = "center", fill = NA),
             by = by]
    }
    result[, Ta := NULL]

  }

  result[]
}


# ── Example (dummy data) ──────────────────────────────────────────────────────

if (FALSE) {

  set.seed(42)

  coord_scq <- c(-8.54, 42.88)   # Santiago de Compostela (lon, lat)
  tz_scq    <- "Europe/Madrid"

  # Hourly series covering two years (gives the moving percentile enough data)
  datetime_seq <- seq(
    as.POSIXct("2022-01-01 00:00:00", tz = tz_scq),
    as.POSIXct("2023-08-31 23:00:00", tz = tz_scq),
    by = "hour"
  )

  n   <- length(datetime_seq)
  hr  <- hour(datetime_seq)
  doy <- yday(datetime_seq)

  # Synthetic temperature: seasonal cycle + diurnal cycle + Gaussian noise
  Ta <- 14 +
    8  * sin((doy - 172) * 2 * pi / 365) +
    7  * sin((hr  -  14) * 2 * pi / 24 ) +
    rnorm(n, sd = 1.2)

  df_dummy <- data.frame(Datetime = datetime_seq, Ta = round(Ta, 1))

  # ── Step 1: moving P95 of daily minimum temperature (recommended threshold) -
  daily_min <- as.data.table(df_dummy)[
    , .(Ta = min(Ta)), by = .(date = as.Date(Datetime))
  ]

  p95_moving <- night_percentile(
    x      = daily_min,
    cols   = c("date", "Ta"),
    p      = 0.95,
    type   = "moving",
    window = 91
  )

  # Use the median P95 as a representative fixed threshold for this example
  thres_val <- median(p95_moving$p95, na.rm = TRUE)

  # ── Step 2: compute HNd and HNe ───────────────────────────────────────────
  res_dur <- night_duration(
    x     = df_dummy,
    cols  = c("Datetime", "Ta"),
    tz    = tz_scq,
    coord = coord_scq,
    thres = thres_val
  )

  res_exc <- night_excess(
    x     = df_dummy,
    cols  = c("Datetime", "Ta"),
    tz    = tz_scq,
    coord = coord_scq,
    thres = thres_val
  )

  result <- merge(res_dur, res_exc, by = "date")
  print(head(result, 10))

  # ── Multiple stations ──────────────────────────────────────────────────────
  daily_min2 <- rbind(
    data.table(city = "SCQ", daily_min),
    data.table(city = "VGO", daily_min[, Ta := Ta + rnorm(.N, 0.5, 0.3)])
  )

  p95_cities <- night_percentile(
    x      = daily_min2,
    cols   = c("date", "Ta"),
    p      = 0.95,
    by     = "city",
    type   = "moving",
    window = 91
  )

  print(head(p95_cities, 10))

}
