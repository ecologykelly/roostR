#' Detect roost onset time near sunset
#'
#' Finds the time each individual first enters sustained low-activity state
#' near sunset — interpreted as roosting. The algorithm:
#' 1. Restricts detections to a window of ± `window_minutes` around sunset.
#' 2. Flags observations where the rolling median (`vol_col`) falls below
#'    `vol_threshold`.
#' 3. Identifies runs of consecutive low-activity observations using
#'    [data.table::rleid()].
#' 4. Returns the start time of the first run whose duration exceeds
#'    `min_duration` minutes.
#'
#' Default thresholds (`vol_threshold = 3`, `min_duration = 10`) are tuned
#' to song sparrow data where night rolling-median ≈ 1 and day ≈ 5. Adjust
#' for other species or transmitter types.
#'
#' Requires [add_roll_median()] to have been called first.
#'
#' @param data A dataframe of Motus detections with a rolling median column.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param sunset_col Character. POSIXct sunset column (UTC). Default `"Sunset"`.
#' @param vol_col Character. Rolling median activity column from
#'   [add_roll_median()]. Default `"roll_vol"`.
#' @param window_minutes Numeric. Half-width of the search window around
#'   sunset, in minutes. Default `90`.
#' @param vol_threshold Numeric. Rolling median must fall below this value to
#'   qualify as a roost detection. Default `3`.
#' @param min_duration Numeric. Minimum duration (minutes) of sustained low
#'   activity required before calling roost onset. Default `10`.
#'
#' @return A dataframe with one row per individual × day-of-year combination
#'   where a roost onset was detected, containing columns `id_col`, `doy_col`,
#'   and `roost_time` (POSIXct UTC).
#'
#' @examples
#' \dontrun{
#' roost_times <- detect_roost_onset(sp)
#' }
#'
#' @importFrom dplyr group_by mutate filter ungroup transmute
#' @importFrom lubridate minutes
#' @importFrom data.table as.data.table rleid
#' @export
detect_roost_onset <- function(data,
                               id_col = "tagDeployID",
                               time_col = "time",
                               doy_col = "doy",
                               sunset_col = "Sunset",
                               vol_col = "roll_vol",
                               window_minutes = 90,
                               vol_threshold = 3,
                               min_duration = 10) {

  window_df <- data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      window_start = .data[[sunset_col]] - lubridate::minutes(window_minutes),
      window_end   = .data[[sunset_col]] + lubridate::minutes(window_minutes)
    ) %>%
    dplyr::filter(
      .data[[time_col]] >= window_start,
      .data[[time_col]] <= window_end
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      low_vol = .data[[vol_col]] < vol_threshold
    )

  dt <- data.table::as.data.table(window_df)

  dt[, low_run := data.table::rleid(low_vol),
     by = c(id_col, doy_col)]

  runs <- dt[
    low_vol == TRUE,
    .(
      start_time = min(get(time_col)),
      duration_min = as.numeric(
        difftime(max(get(time_col)),
                 min(get(time_col)),
                 units = "mins")
      )
    ),
    by = c(id_col, doy_col, "low_run")
  ]

  roost_onset <- runs[
    duration_min >= min_duration,
    .SD[which.min(start_time)],
    by = c(id_col, doy_col)
  ]

  roost_onset <- roost_onset %>%
    dplyr::transmute(
      !!id_col  := .data[[id_col]],
      !!doy_col := .data[[doy_col]],
      roost_time = start_time
    )

  return(roost_onset)
}


#' Detect roost departure time near sunrise
#'
#' Finds the time each individual first shows sustained activity after roosting,
#' interpreted as departure from the roost. The algorithm:
#' 1. Restricts detections to a window from `window_before` minutes before
#'    sunrise to `window_after` minutes after sunrise.
#' 2. Computes a left-aligned (forward) rolling median over `roll_width`
#'    observations.
#' 3. Flags detections where signal difference exceeds `spike_threshold` AND
#'    the forward rolling median exceeds `median_threshold`.
#' 4. Returns the earliest such detection per individual × day.
#'
#' Default thresholds are tuned to song sparrow data. Adjust `spike_threshold`
#' and `median_threshold` for other species or transmitter types.
#'
#' @param data A dataframe of Motus detections with a signal difference column.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param sunrise_col Character. POSIXct sunrise column (UTC). Default `"Sunrise"`.
#' @param sigdiff_col Character. Signal difference column from
#'   [add_signal_diffs()]. Default `"sig.diff"`.
#' @param window_before Numeric. Minutes before sunrise to start the search
#'   window. Default `90`.
#' @param window_after Numeric. Minutes after sunrise to end the search window.
#'   Default `120`.
#' @param roll_width Integer. Number of observations in the left-aligned
#'   (forward) rolling median. Default `10`.
#' @param spike_threshold Numeric. `sigdiff_col` must exceed this to qualify as
#'   an activity spike. Default `4`.
#' @param median_threshold Numeric. Forward rolling median must exceed this
#'   value simultaneously with the spike. Default `3`.
#'
#' @return A dataframe with one row per individual × day-of-year combination
#'   where a departure was detected, containing columns `id_col`, `doy_col`,
#'   and `leave_roost_time` (POSIXct UTC).
#'
#' @examples
#' \dontrun{
#' leave_times <- detect_roost_departure(sp)
#' }
#'
#' @importFrom dplyr group_by filter arrange mutate ungroup slice_min transmute
#' @importFrom lubridate minutes
#' @importFrom zoo rollapply
#' @export
detect_roost_departure <- function(data,
                                   id_col = "tagDeployID",
                                   time_col = "time",
                                   doy_col = "doy",
                                   sunrise_col = "Sunrise",
                                   sigdiff_col = "sig.diff",
                                   window_before = 90,
                                   window_after = 120,
                                   roll_width = 10,
                                   spike_threshold = 4,
                                   median_threshold = 3) {

  sunrise_window <- data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::filter(
      .data[[time_col]] >= .data[[sunrise_col]] - lubridate::minutes(window_before),
      .data[[time_col]] <= .data[[sunrise_col]] + lubridate::minutes(window_after)
    ) %>%
    dplyr::arrange(.data[[id_col]], .data[[doy_col]], .data[[time_col]]) %>%
    dplyr::ungroup()

  sunrise_window <- sunrise_window %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::mutate(
      forward_med = zoo::rollapply(
        .data[[sigdiff_col]],
        width = roll_width,
        FUN = median,
        align = "left",
        fill = NA
      )
    ) %>%
    dplyr::ungroup()

  sunrise_window <- sunrise_window %>%
    dplyr::mutate(
      departure_signal =
        .data[[sigdiff_col]] > spike_threshold &
        forward_med > median_threshold
    )

  leave_roost <- sunrise_window %>%
    dplyr::filter(departure_signal) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::slice_min(.data[[time_col]], n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      !!id_col  := .data[[id_col]],
      !!doy_col := .data[[doy_col]],
      leave_roost_time = .data[[time_col]]
    )

  return(leave_roost)
}


#' Join roost onset and departure times back to the detection dataframe
#'
#' Performs left joins to attach the per-individual × day-of-year roost onset
#' and departure times (from [detect_roost_onset()] and
#' [detect_roost_departure()]) as columns on the full detection dataframe.
#' Rows with no detected roost onset or departure receive `NA`.
#'
#' @param data A dataframe of Motus detections.
#' @param roost_df Dataframe. Output of [detect_roost_onset()], containing
#'   `id_col`, `doy_col`, and `roost_time`.
#' @param leave_df Dataframe. Output of [detect_roost_departure()], containing
#'   `id_col`, `doy_col`, and `leave_roost_time`.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#'
#' @return The input `data` with two new columns:
#'   - `roost_time`: POSIXct roost onset time (NA if not detected).
#'   - `leave_roost_time`: POSIXct departure time (NA if not detected).
#'
#' @examples
#' \dontrun{
#' roost_times <- detect_roost_onset(sp)
#' leave_times <- detect_roost_departure(sp)
#' sp <- add_roost_times(sp, roost_times, leave_times)
#' }
#'
#' @importFrom dplyr left_join
#' @export
add_roost_times <- function(data,
                            roost_df,
                            leave_df,
                            id_col  = "tagDeployID",
                            doy_col = "doy") {
  data %>%
    dplyr::left_join(roost_df, by = c(id_col, doy_col)) %>%
    dplyr::left_join(leave_df, by = c(id_col, doy_col))
}


#' Convert roost onset and departure times to decimal hour-of-day
#'
#' Adds numeric hour-of-day columns (range 0–24) for roost onset and
#' departure times. Decimal hours are computed as
#' `hour + minute/60 + second/3600` and are useful for plotting and
#' comparing timing across nights and individuals.
#'
#' @param data A dataframe containing roost onset and departure POSIXct columns,
#'   typically after calling [add_roost_times()].
#' @param roost_col Character. Roost onset time column (POSIXct).
#'   Default `"roost_time"`.
#' @param leave_col Character. Departure time column (POSIXct).
#'   Default `"leave_roost_time"`.
#'
#' @return The input dataframe with two new numeric columns:
#'   - `roost_hour`: decimal hour of the roost onset time, or `NA`.
#'   - `leave_roost_hour`: decimal hour of the departure time, or `NA`.
#'
#' @examples
#' \dontrun{
#' sp <- add_roost_hours(sp)
#' }
#'
#' @importFrom dplyr mutate if_else
#' @importFrom lubridate hour minute second
#' @export
add_roost_hours <- function(data,
                            roost_col = "roost_time",
                            leave_col = "leave_roost_time") {
  data %>%
    dplyr::mutate(
      roost_hour = dplyr::if_else(
        !is.na(.data[[roost_col]]),
        lubridate::hour(.data[[roost_col]]) +
          lubridate::minute(.data[[roost_col]]) / 60 +
          lubridate::second(.data[[roost_col]]) / 3600,
        NA_real_
      ),
      leave_roost_hour = dplyr::if_else(
        !is.na(.data[[leave_col]]),
        lubridate::hour(.data[[leave_col]]) +
          lubridate::minute(.data[[leave_col]]) / 60 +
          lubridate::second(.data[[leave_col]]) / 3600,
        NA_real_
      )
    )
}
