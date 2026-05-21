#' Compute roosting duration and observation coverage per night
#'
#' For each day-of-year, calculates how long the animal was roosting (from
#' onset on night *d* to departure on day *d+1*) and how much of that interval
#' was actually observed by the receiver network.
#'
#' **Biological night interval**: roost onset on day *d* is paired with
#' departure on day *d+1* via a `doy - 1L` lookup. Nights with no detectable
#' onset or departure are excluded.
#'
#' **Observed time**: detections falling within the roost interval are
#' aggregated into `bin_size_min`-minute bins. The number of occupied bins
#' (bins containing at least one detection) × `bin_size_min / 60` gives the
#' observed time in hours.
#'
#' @param data A dataframe of Motus detections with roost timing columns,
#'   typically after calling [add_roost_times()].
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param doy_col Character. Day-of-year column. Default `"doy"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param roost_col Character. Roost onset time column. Default `"roost_time"`.
#' @param leave_col Character. Departure time column. Default `"leave_roost_time"`.
#' @param bin_size_min Numeric. Width of each time bin in minutes. A bin is
#'   counted as occupied if it contains at least one detection.
#'   Default `2` (appropriate for Motus tags pinging every ~20 sec).
#'
#' @return A dataframe with one row per (`id_col`, `doy_col`) combination where
#'   a complete roost interval could be constructed, containing:
#'   - `id_col`, `doy_col`
#'   - `time_roosting_hr`: total roost interval duration in hours.
#'   - `n_bins_detected`: number of `bin_size_min`-minute bins with detections.
#'   - `observed_time_hr`: `n_bins_detected × bin_size_min / 60`.
#'   - `prop_time_observed`: `observed_time_hr / time_roosting_hr`.
#'
#' @examples
#' \dontrun{
#' night_metrics <- compute_night_observation(sp)
#' summary(night_metrics$prop_time_observed)
#' }
#'
#' @importFrom dplyr distinct filter select mutate left_join bind_rows inner_join
#'   group_by summarise if_else all_of n_distinct transmute
#' @importFrom lubridate floor_date
#' @export
compute_night_observation <- function(data,
                                      id_col       = "tagDeployID",
                                      doy_col      = "doy",
                                      time_col     = "time",
                                      roost_col    = "roost_time",
                                      leave_col    = "leave_roost_time",
                                      bin_size_min = 2) {

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
    dplyr::filter(!is.na(interval_start), !is.na(interval_end)) %>%
    dplyr::mutate(time_roosting_hr = as.numeric(
      difftime(interval_end, interval_start, units = "hours")
    ))

  # ---- 2. Collect detections spanning the overnight period ----
  night_int_lookup <- night_intervals %>%
    dplyr::select(dplyr::all_of(c(id_doy, "interval_start", "interval_end")))

  night_int_prev <- night_int_lookup %>%
    dplyr::mutate(doy_prev = .data[[doy_col]]) %>%
    dplyr::select(-dplyr::all_of(doy_col))

  dets_same <- data %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col))) %>%
    dplyr::inner_join(night_int_lookup, by = id_doy) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col)))

  dets_prev <- data %>%
    dplyr::select(dplyr::all_of(c(id_doy, time_col))) %>%
    dplyr::mutate(doy_prev = .data[[doy_col]] - 1L) %>%
    dplyr::inner_join(night_int_prev, by = c(id_col, "doy_prev")) %>%
    dplyr::filter(.data[[time_col]] >= interval_start,
                  .data[[time_col]] <= interval_end) %>%
    dplyr::transmute(!!id_col  := .data[[id_col]],
                     !!doy_col := doy_prev,
                     !!time_col := .data[[time_col]])

  detections <- dplyr::bind_rows(dets_same, dets_prev) %>%
    dplyr::distinct()

  # ---- 3. Bin detections and count occupied bins ----
  obs_summary <- detections %>%
    dplyr::mutate(
      time_bin = lubridate::floor_date(.data[[time_col]],
                                       unit = paste(bin_size_min, "minutes"))
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_doy))) %>%
    dplyr::summarise(
      n_bins_detected  = dplyr::n_distinct(time_bin),
      observed_time_hr = n_bins_detected * bin_size_min / 60,
      .groups = "drop"
    )

  # ---- 4. Combine with total roost duration and compute proportion observed ----
  night_metrics <- night_intervals %>%
    dplyr::select(dplyr::all_of(c(id_doy, "time_roosting_hr"))) %>%
    dplyr::left_join(obs_summary, by = id_doy) %>%
    dplyr::mutate(
      observed_time_hr   = dplyr::if_else(is.na(observed_time_hr),
                                          0, observed_time_hr),
      prop_time_observed = observed_time_hr / time_roosting_hr
    )

  return(night_metrics)
}
