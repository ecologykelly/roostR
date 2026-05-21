#' Summarise restlessness bouts during the roost interval
#'
#' For each day-of-year, counts restlessness events (signal spikes) that occur
#' within the biological night interval (roost onset on day *d* to departure on
#' day *d+1*).
#'
#' A **spike** is any detection where `spike_col > spike_threshold`. A new
#' **bout** begins when a spike is separated from the previous spike by more
#' than `gap_min` minutes, or when it follows a non-spike observation. Bout
#' duration is measured from the first to the last spike in the bout, with a
#' minimum of `min_bout_min` minutes (giving single-ping spikes a non-zero
#' duration).
#'
#' This function works on data collapsed across receivers: it takes the max
#' signal per (doy, time) before identifying bouts, so receiver interleaving
#' does not inflate bout counts.
#'
#' Default thresholds (`spike_threshold = 4`, `gap_min = 5`,
#' `min_bout_min = 22/60`) are tuned to dark-eyed junco data with a ~22-second
#' ping interval. Requires [add_roost_times()] to have been called.
#'
#' @param data A dataframe of Motus detections with roost timing columns.
#' @param id_col Character. Individual identifier column. Default `"recvDeployName"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param time_col Character. Timestamp column (POSIXct UTC). Default `"time"`.
#' @param roost_col Character. Roost onset time column. Default `"roost_time"`.
#' @param leave_col Character. Departure time column. Default `"leave_roost_time"`.
#' @param spike_col Character. Activity column to threshold. Default `"sig.diff"`.
#' @param spike_threshold Numeric. Values above this are spikes. Default `4`.
#' @param gap_min Numeric. Gaps greater than this many minutes between spikes
#'   begin a new bout. Default `5`.
#' @param min_bout_min Numeric. Minimum bout duration in minutes, applied via
#'   [pmax()] so single-ping bouts have non-zero duration.
#'   Default `22/60` (~22 seconds, one ping interval).
#'
#' @return A dataframe with one row per `doy_col` (nights with a complete roost
#'   interval), containing:
#'   - `doy_col`
#'   - `n_bouts`: number of restlessness bouts.
#'   - `total_restless_min`: total minutes spent in restlessness bouts.
#'   - `max_bout_min`: duration of the longest single bout.
#'
#' @seealso [add_spike_bouts()] to annotate spikes in the full detection
#'   dataframe, [calc_restless_rates()] to normalize by observation time.
#'
#' @examples
#' \dontrun{
#' restless_summary <- calc_restless_all(sp)
#' }
#'
#' @importFrom dplyr filter select distinct mutate group_by summarise ungroup
#'   inner_join bind_rows arrange if_else first last tibble all_of left_join n
#' @export
calc_restless_all <- function(data,
                              id_col          = "tagDeployID",
                              doy_col         = "doy",
                              time_col        = "time",
                              roost_col       = "roost_time",
                              leave_col       = "leave_roost_time",
                              spike_col       = "sig.diff",
                              spike_threshold = 4,
                              gap_min         = 5,
                              min_bout_min    = 22 / 60) {

  required_cols <- c(id_col, doy_col, time_col, roost_col, leave_col, spike_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  id_doy <- c(id_col, doy_col)

  # ---- 1. Build one interval per individual-night ----
  roost_per_night <- data %>%
    dplyr::filter(!is.na(.data[[roost_col]])) %>%
    dplyr::select(dplyr::all_of(id_doy),
                  interval_start = dplyr::all_of(roost_col)) %>%
    dplyr::distinct()

  leave_per_night <- data %>%
    dplyr::filter(!is.na(.data[[leave_col]])) %>%
    dplyr::select(dplyr::all_of(id_doy),
                  interval_end = dplyr::all_of(leave_col)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!doy_col := .data[[doy_col]] - 1L)

  night_intervals <- roost_per_night %>%
    dplyr::left_join(leave_per_night, by = id_doy) %>%
    dplyr::filter(!is.na(interval_start), !is.na(interval_end))

  # ---- 2. Collect detections: evening (doy N) + post-midnight (doy N+1, relabeled N) ----
  night_int_lookup <- night_intervals %>%
    dplyr::select(dplyr::all_of(c(id_doy, "interval_start", "interval_end")))

  night_int_prev <- night_int_lookup %>%
    dplyr::mutate(doy_prev = .data[[doy_col]]) %>%
    dplyr::select(-dplyr::all_of(doy_col))

  dets_same <- data %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col, spike_col))) %>%
    dplyr::inner_join(night_int_lookup, by = id_doy) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col, spike_col)))

  dets_prev <- data %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col, spike_col))) %>%
    dplyr::mutate(doy_prev = .data[[doy_col]] - 1L) %>%
    dplyr::inner_join(night_int_prev, by = c(id_col, "doy_prev")) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::transmute(!!id_col  := .data[[id_col]],
                     !!doy_col := doy_prev,
                     !!time_col := .data[[time_col]],
                     !!spike_col := .data[[spike_col]])

  detections <- dplyr::bind_rows(dets_same, dets_prev) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(id_doy, time_col)))) %>%
    dplyr::summarise(!!spike_col := max(.data[[spike_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(c(id_doy, time_col))))

  if (nrow(detections) == 0) {
    return(dplyr::tibble(
      !!id_col            := data[[id_col]][0],
      !!doy_col           := integer(),
      n_bouts             = integer(),
      total_restless_min  = numeric(),
      max_bout_min        = numeric()
    ))
  }

  # ---- 3. Identify spikes and bouts ----
  detections <- detections %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_doy))) %>%
    dplyr::mutate(
      spike    = .data[[spike_col]] > spike_threshold,
      dt       = as.numeric(difftime(.data[[time_col]],
                                     dplyr::lag(.data[[time_col]]),
                                     units = "mins")),
      new_bout = dplyr::if_else(
        spike & (is.na(dt) | dt > gap_min | !dplyr::lag(spike, default = FALSE)),
        1L, 0L
      ),
      bout_id  = dplyr::if_else(spike, cumsum(new_bout), NA_integer_)
    ) %>%
    dplyr::ungroup()

  # ---- 4. Summarise bouts ----
  restless_summary <- detections %>%
    dplyr::filter(spike) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_doy)), bout_id) %>%
    dplyr::summarise(
      bout_start   = dplyr::first(.data[[time_col]]),
      bout_end     = dplyr::last(.data[[time_col]]),
      duration_min = pmax(
        as.numeric(difftime(bout_end, bout_start, units = "mins")),
        min_bout_min
      ),
      .groups = "drop"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_doy))) %>%
    dplyr::summarise(
      n_bouts            = dplyr::n(),
      total_restless_min = sum(duration_min),
      max_bout_min       = max(duration_min),
      .groups = "drop"
    )

  if (is.numeric(data[[id_col]])) {
    restless_summary <- restless_summary %>%
      dplyr::mutate(!!id_col := as.numeric(.data[[id_col]]))
  }

  return(restless_summary)
}


#' Annotate spike bouts in the full detection dataframe
#'
#' Like [calc_restless_all()], but instead of returning a summary, this
#' function adds spike and bout columns directly to the detection-level
#' dataframe. Useful for plotting individual restlessness events or for
#' further row-level filtering.
#'
#' Spike and bout detection logic is similar to [calc_restless_all()]. Only
#' detections within the biological night interval are flagged; all other rows
#' receive `spike = FALSE` and `bout_id = NA`.
#'
#' Post-midnight rows (UTC hour < 20) are assigned to the previous calendar
#' night via a `night_doy` helper column; this ensures correct interval
#' matching without modifying the original `doy` column.
#'
#' @inheritParams calc_restless_all
#'
#' @return The input dataframe with additional columns:
#'   - `interval_start`, `interval_end`: roost interval bounds (POSIXct).
#'   - `in_roost_interval`: logical, `TRUE` if detection falls within the
#'     roost interval.
#'   - `spike`: logical, `TRUE` if `spike_col > spike_threshold` and within
#'     the roost interval.
#'   - `dt`: minutes since the previous detection in the same night group.
#'   - `new_bout`: `1L` if this spike starts a new bout, `0L` otherwise.
#'   - `bout_id`: integer bout identifier; `NA` for non-spike rows.
#'
#' @seealso [calc_restless_all()] for the summary version.
#'
#' @examples
#' \dontrun{
#' sp <- add_spike_bouts(sp, spike_threshold = 4, gap_min = 5)
#'
#' library(ggplot2)
#' ggplot(sp, aes(x = time, y = sig.diff)) +
#'   geom_point(color = "grey70", size = 0.7) +
#'   geom_point(data = subset(sp, spike), color = "red", size = 1.5) +
#'   geom_hline(yintercept = 4, linetype = 2, color = "blue")
#' }
#'
#' @importFrom dplyr filter select distinct mutate group_by summarise ungroup
#'   left_join if_else lag coalesce rename arrange first
#' @importFrom lubridate hour minute
#' @export
add_spike_bouts <- function(data,
                            id_col          = "recvDeployName",
                            doy_col         = "doy",
                            time_col        = "time",
                            roost_col       = "roost_time",
                            leave_col       = "leave_roost_time",
                            spike_col       = "sig.diff",
                            spike_threshold = 4,
                            gap_min         = 5) {

  # ---- 1. Build one interval per night ----
  roost_per_night <- data %>%
    dplyr::filter(!is.na(.data[[roost_col]])) %>%
    dplyr::select(dplyr::all_of(doy_col),
                  interval_start = dplyr::all_of(roost_col)) %>%
    dplyr::distinct()

  leave_per_night <- data %>%
    dplyr::filter(!is.na(.data[[leave_col]])) %>%
    dplyr::select(dplyr::all_of(doy_col),
                  interval_end = dplyr::all_of(leave_col)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!doy_col := .data[[doy_col]] - 1L)

  night_intervals <- roost_per_night %>%
    dplyr::left_join(leave_per_night, by = doy_col) %>%
    dplyr::filter(!is.na(interval_start), !is.na(interval_end))

  # ---- 2. Assign each row to a night via night_doy, then join intervals ----
  # Post-midnight rows (UTC hour < 20) belong to the previous night.
  # night_doy is used only for the join; doy is never modified.
  data2 <- data %>%
    dplyr::mutate(
      .hour     = lubridate::hour(.data[[time_col]]) +
                  lubridate::minute(.data[[time_col]]) / 60,
      night_doy = dplyr::if_else(.hour < 20, .data[[doy_col]] - 1L, .data[[doy_col]])
    ) %>%
    dplyr::select(-.hour) %>%
    dplyr::left_join(
      night_intervals %>% dplyr::rename(night_doy = !!doy_col),
      by = "night_doy"
    )

  # ---- 3. Compute spikes and bouts inside the roost interval ----
  # Collapse to one row per (night_doy, time) with max signal across receivers
  # so receiver interleaving does not inflate bout counts, then join labels back.
  bout_signals <- data2 %>%
    dplyr::arrange(night_doy, .data[[time_col]]) %>%
    dplyr::group_by(night_doy, .data[[time_col]]) %>%
    dplyr::summarise(
      interval_start = dplyr::first(interval_start),
      interval_end   = dplyr::first(interval_end),
      !!spike_col    := max(.data[[spike_col]], na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    dplyr::group_by(night_doy) %>%
    dplyr::mutate(
      in_roost_interval =
        !is.na(interval_start) &
        !is.na(interval_end)   &
        .data[[time_col]] >= interval_start &
        .data[[time_col]] <= interval_end,

      spike = dplyr::if_else(
        in_roost_interval & .data[[spike_col]] > spike_threshold,
        TRUE, FALSE
      ),

      dt = as.numeric(difftime(.data[[time_col]],
                               dplyr::lag(.data[[time_col]]),
                               units = "mins")),

      new_bout = dplyr::if_else(
        spike & (is.na(dt) | dt > gap_min | !dplyr::lag(spike, default = FALSE)),
        1L, 0L
      ),

      bout_id = dplyr::if_else(spike, cumsum(new_bout), NA_integer_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(night_doy, dplyr::all_of(time_col),
                  in_roost_interval, spike, dt, new_bout, bout_id)

  data2 <- data2 %>%
    dplyr::left_join(bout_signals, by = c("night_doy", time_col)) %>%
    dplyr::mutate(
      in_roost_interval = dplyr::coalesce(in_roost_interval, FALSE),
      spike             = dplyr::coalesce(spike, FALSE)
    ) %>%
    dplyr::select(-night_doy)

  return(data2)
}


#' Summarise antenna bearing during the roost interval
#'
#' For each day-of-year, computes the circular mean and modal bearing of
#' detections falling within the biological night interval. Useful for
#' identifying which direction the roosting animal is oriented relative to
#' the receiver array.
#'
#' The circular mean uses vector averaging over angles converted to radians.
#' The mode returns `NA` when there is no unique most-frequent bearing.
#'
#' @param data A dataframe of Motus detections with roost timing and bearing
#'   columns, typically after calling [add_roost_times()].
#' @param id_col Character. Individual identifier column. Default `"recvDeployName"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param roost_col Character. Roost onset time column. Default `"roost_time"`.
#' @param leave_col Character. Departure time column. Default `"leave_roost_time"`.
#' @param bearing_col Character. Antenna bearing column (degrees, 0–360).
#'   Default `"antBearing"`.
#'
#' @return A dataframe with one row per `doy_col` (nights with a complete roost
#'   interval and bearing data), containing:
#'   - `doy_col`
#'   - `mean_bearing`: circular mean bearing in degrees (0–360).
#'   - `mode_bearing`: modal bearing, or `NA` if no unique mode.
#'
#' @examples
#' \dontrun{
#' bearing_summary <- calc_bearing_summary(sp)
#' }
#'
#' @importFrom dplyr filter select distinct mutate left_join inner_join bind_rows
#'   group_by summarise all_of tibble any_of
#' @export
calc_bearing_summary <- function(data,
                                 id_col      = "recvDeployName",
                                 doy_col     = "doy",
                                 time_col    = "time",
                                 roost_col   = "roost_time",
                                 leave_col   = "leave_roost_time",
                                 bearing_col = "antBearing") {

  circular_mean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    rad <- x * pi / 180
    deg <- atan2(mean(sin(rad)), mean(cos(rad))) * 180 / pi
    (deg + 360) %% 360
  }

  mode_bearing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    tbl   <- table(x)
    modes <- as.numeric(names(tbl[tbl == max(tbl)]))
    if (length(modes) > 1) NA_real_ else modes
  }

  roost_per_night <- data %>%
    dplyr::filter(!is.na(.data[[roost_col]])) %>%
    dplyr::select(dplyr::all_of(doy_col),
                  interval_start = dplyr::all_of(roost_col)) %>%
    dplyr::distinct()

  leave_per_night <- data %>%
    dplyr::filter(!is.na(.data[[leave_col]])) %>%
    dplyr::select(dplyr::all_of(doy_col),
                  interval_end = dplyr::all_of(leave_col)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!doy_col := .data[[doy_col]] - 1L)

  night_intervals <- roost_per_night %>%
    dplyr::left_join(leave_per_night, by = doy_col) %>%
    dplyr::filter(!is.na(interval_start), !is.na(interval_end))

  night_int_lookup <- night_intervals %>%
    dplyr::select(dplyr::all_of(c(doy_col, "interval_start", "interval_end")))

  night_int_prev <- night_int_lookup %>%
    dplyr::mutate(doy_prev = .data[[doy_col]]) %>%
    dplyr::select(-dplyr::all_of(doy_col))

  dets_same <- data %>%
    dplyr::select(-dplyr::any_of(c("interval_start", "interval_end"))) %>%
    dplyr::inner_join(night_int_lookup, by = doy_col) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end)

  dets_prev <- data %>%
    dplyr::select(-dplyr::any_of(c("interval_start", "interval_end"))) %>%
    dplyr::mutate(doy_prev = .data[[doy_col]] - 1L) %>%
    dplyr::inner_join(night_int_prev, by = "doy_prev") %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::mutate(!!doy_col := doy_prev) %>%
    dplyr::select(-doy_prev)

  detections <- dplyr::bind_rows(dets_same, dets_prev) %>%
    dplyr::distinct()

  if (nrow(detections) == 0) {
    return(dplyr::tibble(
      !!doy_col     := integer(),
      mean_bearing  = numeric(),
      mode_bearing  = numeric()
    ))
  }

  detections %>%
    dplyr::group_by(.data[[doy_col]]) %>%
    dplyr::summarise(
      mean_bearing = circular_mean(.data[[bearing_col]]),
      mode_bearing = mode_bearing(.data[[bearing_col]]),
      .groups = "drop"
    )
}


#' Compute restlessness rates normalized by observation time
#'
#' Normalizes restlessness metrics from [calc_restless_all()] by the amount of
#' time the animal was actually observed during the roost interval (from
#' [compute_night_observation()]). This corrects for nights where the receiver
#' only detected the animal for part of the roost period.
#'
#' @param data A dataframe with restlessness and observation columns, typically
#'   after joining outputs of [calc_restless_all()] and
#'   [compute_night_observation()] back to the detection data.
#' @param id_col Character. Individual identifier column. Default `"recvDeployName"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param bouts_col Character. Bout count column from [calc_restless_all()].
#'   Default `"n_bouts"`.
#' @param restless_min_col Character. Total restless time column (minutes).
#'   Default `"total_restless_min"`.
#' @param observed_hr_col Character. Observed time column (hours) from
#'   [compute_night_observation()]. Default `"observed_time_hr"`.
#'
#' @return A dataframe with one row per `doy_col`, containing the input columns
#'   plus:
#'   - `total_restless_hr`: `total_restless_min / 60`.
#'   - `restless_per_obs_hr`: bouts per observed hour (`NA` if
#'     `observed_time_hr` is 0 or `NA`).
#'   - `prop_time_restless`: proportion of observed roost time spent restless
#'     (`NA` if `observed_time_hr` is 0 or `NA`).
#'
#' @examples
#' \dontrun{
#' restless_rates <- calc_restless_rates(sp)
#' }
#'
#' @importFrom dplyr distinct across all_of mutate if_else
#' @export
calc_restless_rates <- function(data,
                                id_col           = "tagDeployID",
                                doy_col          = "doy",
                                bouts_col        = "n_bouts",
                                restless_min_col = "total_restless_min",
                                observed_hr_col  = "observed_time_hr") {
  data %>%
    dplyr::filter(!is.na(.data[[bouts_col]])) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(
      c(id_col, doy_col, bouts_col, restless_min_col, observed_hr_col)
    ))) %>%
    dplyr::mutate(
      total_restless_hr = .data[[restless_min_col]] / 60,

      restless_per_obs_hr = dplyr::if_else(
        !is.na(.data[[observed_hr_col]]) & .data[[observed_hr_col]] > 0,
        .data[[bouts_col]] / .data[[observed_hr_col]],
        NA_real_
      ),

      prop_time_restless = dplyr::if_else(
        !is.na(.data[[observed_hr_col]]) & .data[[observed_hr_col]] > 0,
        total_restless_hr / .data[[observed_hr_col]],
        NA_real_
      )
    )
}
