#' Compute roosting duration and observation coverage per night
#'
#' For each individual × day-of-year, calculates how long the animal was
#' roosting (total interval from onset to next day's departure) and how much
#' of that interval was actually observed by the receiver network.
#'
#' **Biological night interval**: roost onset on day *d* to departure on day
#' *d+1* (paired via [dplyr::lead()]). Nights with no detectable onset or
#' departure are excluded.
#'
#' **Observed time**: continuous detection segments within the roost interval
#' are summed. Gaps larger than `gap_threshold` minutes are treated as periods
#' when the animal moved out of receiver range and are excluded.
#'
#' @param data A dataframe of Motus detections with roost timing columns,
#'   typically after calling [add_roost_times()].
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param roost_col Character. Roost onset time column. Default `"roost_time"`.
#' @param leave_col Character. Departure time column. Default `"leave_roost_time"`.
#' @param sunset_col Character. Sunset column (currently retained for future
#'   fallback logic). Default `"Sunset"`.
#' @param sunrise_col Character. Sunrise column. Default `"Sunrise"`.
#' @param gap_threshold Numeric. Gaps larger than this many minutes between
#'   consecutive detections are treated as out-of-range periods and break an
#'   observation segment. Default `10`.
#'
#' @return A dataframe with one row per individual × day-of-year (nights where
#'   a complete roost interval could be constructed), containing:
#'   - `id_col`, `doy_col`
#'   - `time_roosting_hr`: total roost interval duration in hours.
#'   - `observed_time_hr`: total time within roost interval with detections,
#'     summed across continuous segments.
#'   - `prop_time_observed`: `observed_time_hr / time_roosting_hr`.
#'
#' @examples
#' \dontrun{
#' night_metrics <- compute_night_observation(sp, gap_threshold = 10)
#' summary(night_metrics$prop_time_observed)
#' }
#'
#' @importFrom dplyr distinct arrange group_by mutate lead ungroup filter
#'   select inner_join summarise left_join if_else all_of
#' @export
compute_night_observation <- function(data,
                                      id_col = "tagDeployID",
                                      doy_col = "doy",
                                      time_col = "time",
                                      roost_col = "roost_time",
                                      leave_col = "leave_roost_time",
                                      sunset_col = "Sunset",
                                      sunrise_col = "Sunrise",
                                      gap_threshold = 10) {

  # ---- 1. Construct biological night intervals ----
  night_intervals <- data %>%
    dplyr::distinct(tagDeployID, doy,
                    roost_time, leave_roost_time) %>%
    dplyr::arrange(tagDeployID, doy) %>%
    dplyr::group_by(tagDeployID) %>%
    dplyr::mutate(
      leave_next     = dplyr::lead(leave_roost_time),
      interval_start = roost_time,
      interval_end   = leave_next,
      time_roosting_hr = as.numeric(difftime(interval_end,
                                             interval_start,
                                             units = "hours"))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(time_roosting_hr))

  # ---- 2. Attach interval to detections ----
  detections <- data %>%
    dplyr::select(dplyr::all_of(c(id_col, doy_col, time_col))) %>%
    dplyr::inner_join(
      night_intervals %>%
        dplyr::select(dplyr::all_of(c(id_col, doy_col,
                                      "interval_start",
                                      "interval_end"))),
      by = c(id_col, doy_col)
    ) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::arrange(.data[[id_col]],
                   .data[[doy_col]],
                   .data[[time_col]])

  # ---- 3. Segment continuous detection runs ----
  detections <- detections %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::mutate(
      dt = as.numeric(difftime(.data[[time_col]],
                               dplyr::lag(.data[[time_col]]),
                               units = "mins")),
      new_segment = dplyr::if_else(is.na(dt) | dt > gap_threshold,
                                   1L, 0L),
      segment_id  = cumsum(new_segment)
    ) %>%
    dplyr::ungroup()

  # ---- 4. Compute observed duration per segment ----
  segment_summary <- detections %>%
    dplyr::group_by(.data[[id_col]],
                    .data[[doy_col]],
                    segment_id) %>%
    dplyr::summarise(
      seg_start  = dplyr::first(.data[[time_col]]),
      seg_end    = dplyr::last(.data[[time_col]]),
      seg_dur_hr = as.numeric(difftime(seg_end, seg_start,
                                       units = "hours")),
      .groups = "drop"
    )

  # ---- 5. Sum segments per night ----
  obs_summary <- segment_summary %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::summarise(
      observed_time_hr = sum(seg_dur_hr, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- 6. Combine with total roost duration ----
  night_metrics <- night_intervals %>%
    dplyr::mutate(
      time_roosting_hr = as.numeric(difftime(interval_end,
                                             interval_start,
                                             units = "hours"))
    ) %>%
    dplyr::select(dplyr::all_of(c(id_col, doy_col, "time_roosting_hr"))) %>%
    dplyr::left_join(obs_summary, by = c(id_col, doy_col)) %>%
    dplyr::mutate(
      observed_time_hr = dplyr::if_else(is.na(observed_time_hr),
                                        0,
                                        observed_time_hr),
      prop_time_observed = observed_time_hr / time_roosting_hr
    )

  return(night_metrics)
}
