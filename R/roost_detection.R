#' Detect roost onset time near sunset
#'
#' Finds the time each individual first enters a sustained low-activity state
#' near sunset — interpreted as roosting. The algorithm:
#' 1. Restricts detections to a window of ± `window_minutes` around sunset.
#' 2. Flags observations where the rolling median (`vol_col`) falls below
#'    `vol_threshold`.
#' 3. Identifies consecutive low-activity runs using [data.table::rleid()].
#' 4. Filters runs meeting `min_duration` minutes of sustained quiet.
#' 5. For each candidate, checks that the following `confirm_min` minutes
#'    contain at least `confirm_frac` low-volatility detections (confirmation
#'    window). Candidates failing this check are discarded.
#' 6. For passing candidates, additionally checks that the period after the
#'    confirmation window (until `window_end`) does not show resumed high
#'    activity. Candidates with post-window activity below `confirm_frac`
#'    low-vol are discarded.
#' 7. Returns the latest surviving confirmed onset per bird-night.
#'
#' Default thresholds are tuned to dark-eyed junco data. Adjust for other
#' species or transmitter types.
#'
#' Requires [add_roll_median()] to have been called first.
#'
#' @param data A dataframe of Motus detections with a rolling median column.
#' @param id_col Character. Individual identifier column. Default `"recvDeployName"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param sunset_col Character. POSIXct sunset column (UTC). Default `"Sunset"`.
#' @param vol_col Character. Rolling median activity column from
#'   [add_roll_median()]. Default `"roll_vol"`.
#' @param window_minutes Numeric. Half-width of the search window around
#'   sunset, in minutes. Default `90`.
#' @param vol_threshold Numeric. Rolling median must fall below this value to
#'   qualify as a low-activity detection. Default `3`.
#' @param min_duration Numeric. Minimum duration (minutes) of sustained low
#'   activity required for a candidate run. Default `10`.
#' @param confirm_min Numeric. Duration (minutes) of the confirmation window
#'   following each candidate onset. Default `30`.
#' @param confirm_frac Numeric. Minimum fraction of detections within the
#'   confirmation window that must be low-volatility for the candidate to be
#'   accepted. Default `0.8`.
#'
#' @return A dataframe with one row per (`id_col`, `doy_col`) combination where
#'   a confirmed roost onset was detected, containing columns `id_col`,
#'   `doy_col`, and `roost_time` (POSIXct UTC).
#'
#' @examples
#' \dontrun{
#' roost_times <- detect_roost_onset(sp)
#' }
#'
#' @importFrom dplyr group_by mutate filter ungroup transmute summarise left_join
#'   select all_of slice_max tibble
#' @importFrom lubridate minutes
#' @importFrom data.table as.data.table rleid
#' @export
detect_roost_onset <- function(data,
                               id_col         = "recvDeployName",
                               time_col       = "time",
                               doy_col        = "doy",
                               sunset_col     = "Sunset",
                               vol_col        = "roll_vol",
                               window_minutes = 90,
                               vol_threshold  = 3,
                               min_duration   = 10,
                               confirm_min    = 30,
                               confirm_frac   = 0.8) {

  # ---- 1. Restrict to search window around sunset and flag low-volatility ----
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
    dplyr::mutate(low_vol = .data[[vol_col]] < vol_threshold)

  # ---- 2. Identify consecutive low-vol runs per bird-night ----
  dt <- data.table::as.data.table(window_df)

  dt[, low_run := data.table::rleid(low_vol),
     by = c(id_col, doy_col)]

  runs <- dt[
    low_vol == TRUE,
    .(
      start_time   = min(get(time_col)),
      duration_min = as.numeric(
        difftime(max(get(time_col)),
                 min(get(time_col)),
                 units = "mins")
      )
    ),
    by = c(id_col, doy_col, "low_run")
  ]

  # ---- 3. Filter to runs meeting minimum duration ----
  candidates <- runs[duration_min >= min_duration] %>%
    dplyr::as_tibble()

  if (nrow(candidates) == 0) return(
    dplyr::tibble(!!id_col := character(), !!doy_col := integer(),
                  roost_time = as.POSIXct(NA))
  )

  # ---- 4. Confirmation window check ----
  # After each candidate onset, the next confirm_min minutes must have at least
  # confirm_frac low-vol detections; otherwise the bird had not truly settled.
  confirm_check <- candidates %>%
    dplyr::left_join(
      window_df %>%
        dplyr::select(dplyr::all_of(c(id_col, doy_col, time_col, "low_vol"))),
      by       = c(id_col, doy_col),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      .data[[time_col]] >= start_time,
      .data[[time_col]] <= start_time + lubridate::minutes(confirm_min)
    ) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]], low_run, start_time) %>%
    dplyr::summarise(
      frac_low_vol = mean(low_vol, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(frac_low_vol),
                  frac_low_vol >= confirm_frac)

  # ---- 4b. Post-confirmation activity check ----
  # If the bird resumes high activity between the confirmation window end and
  # window_end, the candidate is rejected. Candidates with no post-window
  # detections are kept (bird may have drifted out of range after roosting).
  post_check <- confirm_check %>%
    dplyr::left_join(
      window_df %>%
        dplyr::select(dplyr::all_of(c(id_col, doy_col, time_col,
                                      "low_vol", "window_end"))),
      by       = c(id_col, doy_col),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      .data[[time_col]] > start_time + lubridate::minutes(confirm_min),
      .data[[time_col]] <= window_end
    ) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]], low_run, start_time) %>%
    dplyr::summarise(
      frac_low_vol_post = mean(low_vol, na.rm = TRUE),
      .groups = "drop"
    )

  confirm_check <- confirm_check %>%
    dplyr::left_join(post_check,
                     by = c(id_col, doy_col, "low_run", "start_time")) %>%
    dplyr::filter(is.na(frac_low_vol_post) | frac_low_vol_post >= confirm_frac)

  # ---- 5. Take the latest confirmed onset per bird-night ----
  roost_onset <- confirm_check %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::slice_max(start_time, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
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
#' Default thresholds are tuned to dark-eyed junco data. Adjust `spike_threshold`
#' and `median_threshold` for other species or transmitter types.
#'
#' @param data A dataframe of Motus detections with a signal difference column.
#' @param id_col Character. Individual identifier column. Default `"recvDeployName"`.
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
#' @return A dataframe with one row per (`id_col`, `doy_col`) combination where
#'   a departure was detected, containing columns `id_col`, `doy_col`, and
#'   `leave_roost_time` (POSIXct UTC).
#'
#' @examples
#' \dontrun{
#' leave_times <- detect_roost_departure(sp)
#' }
#'
#' @importFrom dplyr group_by filter arrange mutate ungroup slice_min transmute
#' @importFrom lubridate minutes
#' @importFrom stats median
#' @importFrom zoo rollapply
#' @export
detect_roost_departure <- function(data,
                                   id_col           = "recvDeployName",
                                   time_col         = "time",
                                   doy_col          = "doy",
                                   sunrise_col      = "Sunrise",
                                   sigdiff_col      = "sig.diff",
                                   window_before    = 90,
                                   window_after     = 120,
                                   roll_width       = 10,
                                   spike_threshold  = 4,
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
        FUN   = median,
        align = "left",
        fill  = NA
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
#' Performs left joins to attach per-night roost onset and departure times
#' (from [detect_roost_onset()] and [detect_roost_departure()]) as columns on
#' the full detection dataframe. Because roost and departure times are detected
#' on combined-receiver data and represent a single value per night, the join
#' is by `doy_col` only — all receivers on the same night receive the same
#' roost and departure times. Rows with no detected onset or departure receive
#' `NA`.
#'
#' @param data A dataframe of Motus detections.
#' @param roost_df Dataframe. Output of [detect_roost_onset()], containing
#'   `doy_col` and `roost_time`.
#' @param leave_df Dataframe. Output of [detect_roost_departure()], containing
#'   `doy_col` and `leave_roost_time`.
#' @param id_col Character. Individual identifier column (retained for API
#'   consistency). Default `"recvDeployName"`.
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
#' @importFrom dplyr left_join select all_of
#' @export
add_roost_times <- function(data,
                            roost_df,
                            leave_df,
                            id_col  = "recvDeployName",
                            doy_col = "doy") {
  data %>%
    dplyr::left_join(
      roost_df %>% dplyr::select(dplyr::all_of(c(doy_col, "roost_time"))),
      by = doy_col
    ) %>%
    dplyr::left_join(
      leave_df %>% dplyr::select(dplyr::all_of(c(doy_col, "leave_roost_time"))),
      by = doy_col
    )
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
#' @param sunset_col Character. Sunset time column (POSIXct). Default `"Sunset"`.
#' @param sunrise_col Character. Sunrise time column (POSIXct). Default `"Sunrise"`.
#'
#' @return The input dataframe with four new numeric columns:
#'   - `roost_hour`: decimal hour of the roost onset time, or `NA`.
#'   - `leave_roost_hour`: decimal hour of the departure time, or `NA`.
#'   - `sunset_hour`: decimal hour of sunset.
#'   - `sunrise_hour`: decimal hour of sunrise.
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
                            roost_col   = "roost_time",
                            leave_col   = "leave_roost_time",
                            sunset_col  = "Sunset",
                            sunrise_col = "Sunrise") {
  to_hour <- function(x) lubridate::hour(x) + lubridate::minute(x) / 60 +
                          lubridate::second(x) / 3600
  data %>%
    dplyr::mutate(
      roost_hour = dplyr::if_else(
        !is.na(.data[[roost_col]]), to_hour(.data[[roost_col]]), NA_real_
      ),
      leave_roost_hour = dplyr::if_else(
        !is.na(.data[[leave_col]]), to_hour(.data[[leave_col]]), NA_real_
      ),
      sunset_hour  = to_hour(.data[[sunset_col]]),
      sunrise_hour = to_hour(.data[[sunrise_col]])
    )
}
