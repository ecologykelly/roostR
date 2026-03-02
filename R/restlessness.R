#' Summarise restlessness bouts during the roost interval
#'
#' For each individual × day-of-year, counts restlessness events (signal
#' spikes) that occur within the biological night interval (roost onset on
#' day *d* to departure on day *d+1*).
#'
#' A **spike** is any detection where `spike_col > spike_threshold`.
#' A new **bout** begins when a spike is separated from the previous spike by
#' more than `gap_min` minutes, or when it follows a non-spike observation.
#' Bout duration is measured from the first to the last spike in the bout.
#'
#' Default thresholds (`spike_threshold = 4`, `gap_min = 2`) are tuned to
#' dark-eyed junco data. Requires [add_roost_times()] to have been called.
#'
#' @param data A dataframe of Motus detections with roost timing columns.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param time_col Character. Timestamp column (POSIXct UTC). Default `"time"`.
#' @param roost_col Character. Roost onset time column. Default `"roost_time"`.
#' @param leave_col Character. Departure time column. Default `"leave_roost_time"`.
#' @param spike_col Character. Activity column to threshold. Default `"sig.diff"`.
#' @param spike_threshold Numeric. Values above this are spikes. Default `4`.
#' @param gap_min Numeric. Gaps greater than this many minutes between spikes
#'   begin a new bout. Default `2`.
#'
#' @return A dataframe with one row per individual × day-of-year (nights with
#'   a complete roost interval), containing:
#'   - `id_col`, `doy_col`
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
#' @importFrom dplyr distinct across all_of arrange group_by mutate lag ungroup
#'   inner_join filter if_else summarise n
#' @export
calc_restless_all <- function(data,
                              id_col = "tagDeployID",
                              doy_col = "doy",
                              time_col = "time",
                              roost_col = "roost_time",
                              leave_col = "leave_roost_time",
                              spike_col = "sig.diff",
                              spike_threshold = 4,
                              gap_min = 2) {

  required_cols <- c(id_col, doy_col, time_col, roost_col, leave_col, spike_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ---- 1. Construct overnight intervals ----
  night_intervals <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(
      c(id_col, doy_col, roost_col, leave_col)
    ))) %>%
    dplyr::arrange(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      leave_next     = dplyr::lead(.data[[leave_col]]),
      interval_start = .data[[roost_col]],
      interval_end   = leave_next
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(interval_start), !is.na(interval_end))

  # ---- 2. Attach intervals to detections ----
  detections <- data %>%
    dplyr::inner_join(
      night_intervals %>%
        dplyr::select(dplyr::all_of(c(id_col, doy_col,
                                      "interval_start", "interval_end"))),
      by = c(id_col, doy_col)
    ) %>%
    dplyr::filter(
      !!as.name(time_col) >= interval_start,
      !!as.name(time_col) <= interval_end
    ) %>%
    dplyr::arrange(.data[[id_col]], .data[[doy_col]], .data[[time_col]])

  if (nrow(detections) == 0) {
    return(data.frame())
  }

  # ---- 3. Identify spikes and bouts ----
  detections <- detections %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::mutate(
      spike = .data[[spike_col]] > spike_threshold,
      dt    = as.numeric(difftime(.data[[time_col]],
                                  dplyr::lag(.data[[time_col]]),
                                  units = "mins")),
      new_bout = dplyr::if_else(
        spike & (is.na(dt) | dt > gap_min | !dplyr::lag(spike, default = FALSE)),
        1L, 0L
      ),
      bout_id = dplyr::if_else(spike, cumsum(new_bout), NA_integer_)
    ) %>%
    dplyr::ungroup()

  # ---- 4. Summarise bouts ----
  restless_summary <- detections %>%
    dplyr::filter(spike) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]], bout_id) %>%
    dplyr::summarise(
      bout_start   = dplyr::first(.data[[time_col]]),
      bout_end     = dplyr::last(.data[[time_col]]),
      duration_min = as.numeric(difftime(bout_end, bout_start,
                                         units = "mins")),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::summarise(
      n_bouts           = dplyr::n(),
      total_restless_min = sum(duration_min),
      max_bout_min      = max(duration_min),
      .groups = "drop"
    )

  return(restless_summary)
}


#' Annotate spike bouts in the full detection dataframe
#'
#' Like [calc_restless_all()], but instead of returning a summary, this
#' function adds spike and bout columns directly to the detection-level
#' dataframe. Useful for plotting individual restlessness events or for
#' further row-level filtering.
#'
#' Spike and bout detection logic is identical to [calc_restless_all()].
#' Only detections within the biological night interval are flagged; all
#' other rows receive `spike = FALSE` and `bout_id = NA`.
#'
#' @inheritParams calc_restless_all
#'
#' @return The input dataframe with additional columns:
#'   - `interval_start`, `interval_end`: roost interval bounds (POSIXct).
#'   - `in_roost_interval`: logical, `TRUE` if detection falls within the
#'     roost interval.
#'   - `spike`: logical, `TRUE` if `spike_col > spike_threshold` and within
#'     the roost interval.
#'   - `dt`: minutes since the previous detection within the same
#'     individual × day group.
#'   - `new_bout`: `1L` if this spike starts a new bout, `0L` otherwise.
#'   - `bout_id`: integer bout identifier; `NA` for non-spike rows.
#'
#' @seealso [calc_restless_all()] for the summary version.
#'
#' @examples
#' \dontrun{
#' sp <- add_spike_bouts(sp, spike_threshold = 4, gap_min = 2)
#'
#' library(ggplot2)
#' ggplot(sp, aes(x = time, y = sig.diff)) +
#'   geom_point(color = "grey70", size = 0.7) +
#'   geom_point(data = subset(sp, spike), color = "red", size = 1.5) +
#'   geom_hline(yintercept = 4, linetype = 2, color = "blue")
#' }
#'
#' @importFrom dplyr distinct across all_of arrange group_by mutate lead
#'   ungroup filter left_join if_else lag
#' @export
add_spike_bouts <- function(data,
                            id_col = "tagDeployID",
                            doy_col = "doy",
                            time_col = "time",
                            roost_col = "roost_time",
                            leave_col = "leave_roost_time",
                            spike_col = "sig.diff",
                            spike_threshold = 4,
                            gap_min = 2) {

  # ---- 1. Construct overnight intervals ----
  night_intervals <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(
      c(id_col, doy_col, roost_col, leave_col)
    ))) %>%
    dplyr::arrange(.data[[id_col]], .data[[doy_col]]) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      leave_next     = dplyr::lead(.data[[leave_col]]),
      interval_start = .data[[roost_col]],
      interval_end   = leave_next
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(interval_start), !is.na(interval_end))

  # ---- 2. Join intervals back to full dataset ----
  data2 <- data %>%
    dplyr::left_join(
      night_intervals %>%
        dplyr::select(dplyr::all_of(c(id_col, doy_col,
                                      "interval_start", "interval_end"))),
      by = c(id_col, doy_col)
    )

  # ---- 3. Compute spikes and bouts only inside roost interval ----
  data2 <- data2 %>%
    dplyr::arrange(.data[[id_col]], .data[[doy_col]], .data[[time_col]]) %>%
    dplyr::group_by(.data[[id_col]], .data[[doy_col]]) %>%
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
    dplyr::ungroup()

  return(data2)
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
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param bouts_col Character. Bout count column from [calc_restless_all()].
#'   Default `"n_bouts"`.
#' @param restless_min_col Character. Total restless time column (minutes).
#'   Default `"total_restless_min"`.
#' @param observed_hr_col Character. Observed time column (hours) from
#'   [compute_night_observation()]. Default `"observed_time_hr"`.
#'
#' @return A dataframe with one row per individual × day-of-year, containing
#'   the input columns plus:
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
                                id_col = "tagDeployID",
                                doy_col = "doy",
                                bouts_col = "n_bouts",
                                restless_min_col = "total_restless_min",
                                observed_hr_col = "observed_time_hr") {
  data %>%
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
