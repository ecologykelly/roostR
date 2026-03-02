#' Collapse multi-antenna detections to one row per timestamp
#'
#' A Motus tag can be detected simultaneously on up to four antennas at the
#' same receiver. This function reduces those duplicate rows to a single row
#' per individual per timestamp, retaining the maximum signal strength and
#' computing the mean signal across antennas.
#'
#' @param data A dataframe of Motus detections. Must contain columns named by
#'   `id_col`, `time_col`, and `sig_col`.
#' @param id_col Character. Column identifying individual animals (deployment-
#'   specific). Default `"tagDeployID"`. Using `tagDeployID` rather than
#'   `motusTagID` correctly handles cases where a physical tag is redeployed
#'   on a different animal.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param sig_col Character. Signal strength column (dBm). Default `"sig"`.
#'
#' @return The input dataframe collapsed to one row per (`id_col`, `time_col`)
#'   combination. All non-signal columns retain their first value within each
#'   group. Two columns are added or modified:
#'   - `sig`: maximum signal across antennas at that timestamp.
#'   - `sig.mean`: mean signal across antennas at that timestamp.
#'
#' @examples
#' \dontrun{
#' data(sparrow52550)
#' sp <- collapse_motus_time(sparrow52550)
#' }
#'
#' @importFrom dplyr group_by summarise across everything first
#' @export
collapse_motus_time <- function(data,
                                id_col = "tagDeployID",
                                time_col = "time",
                                sig_col = "sig") {
  data %>%
    dplyr::group_by(.data[[id_col]], .data[[time_col]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), dplyr::first),
      !!sig_col := max(.data[[sig_col]], na.rm = TRUE),
      sig.mean = mean(.data[[sig_col]], na.rm = TRUE),
      .groups = "drop"
    )
}


#' Compute signal differences and inter-detection intervals
#'
#' For each individual, calculates the absolute change in signal strength
#' between consecutive detections (`sig.diff`, `sig.diff.mean`) and the
#' elapsed time between them (`time.diff`). These are the primary behavioral
#' activity proxies used in downstream analyses — larger signal changes
#' indicate more movement near the receiver tower.
#'
#' @param data A dataframe of Motus detections, typically the output of
#'   [collapse_motus_time()].
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param sig_col Character. Maximum signal column. Default `"sig"`.
#' @param sig_mean_col Character. Mean signal column. Default `"sig.mean"`.
#'
#' @return The input dataframe with three new columns:
#'   - `sig.diff`: absolute change in maximum signal from the previous
#'     detection. Set to `0` for the first detection per individual.
#'   - `sig.diff.mean`: absolute change in mean signal. Set to `0` for the
#'     first detection.
#'   - `time.diff`: minutes elapsed since the previous detection. `NA` for
#'     the first detection (no prior observation to compare).
#'
#' @examples
#' \dontrun{
#' sp <- collapse_motus_time(sparrow52550)
#' sp <- add_signal_diffs(sp)
#' }
#'
#' @importFrom dplyr group_by arrange mutate ungroup lag
#' @export
add_signal_diffs <- function(data,
                             id_col = "tagDeployID",
                             time_col = "time",
                             sig_col = "sig",
                             sig_mean_col = "sig.mean") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(.data[[time_col]], .by_group = TRUE) %>%
    dplyr::mutate(
      sig.diff = abs(.data[[sig_col]] - dplyr::lag(.data[[sig_col]])),
      sig.diff = ifelse(is.na(sig.diff), 0, sig.diff),

      sig.diff.mean = abs(.data[[sig_mean_col]] - dplyr::lag(.data[[sig_mean_col]])),
      sig.diff.mean = ifelse(is.na(sig.diff.mean), 0, sig.diff.mean),

      time.diff = as.numeric(
        difftime(.data[[time_col]],
                 dplyr::lag(.data[[time_col]]),
                 units = "mins")
      )
    ) %>%
    dplyr::ungroup()
}


#' Flag observation gaps and assign run IDs
#'
#' Marks gaps between consecutive detections and segments the detection series
#' into continuous runs. A "continuous" observation is one where the gap to the
#' previous detection is within `threshold` minutes. The default 2-minute
#' threshold is tuned to Motus transmitters pinging every 19–26 seconds; adjust
#' for other hardware.
#'
#' Requires [add_signal_diffs()] to have been called first (to produce
#' `time.diff`).
#'
#' @param data A dataframe of Motus detections with a `time.diff` column.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param time_col Character. Timestamp column. Default `"time"`.
#' @param time_diff_col Character. Elapsed-time column in minutes, from
#'   [add_signal_diffs()]. Default `"time.diff"`.
#' @param threshold Numeric. Gap threshold in minutes. Detections within this
#'   window are marked as continuous. Default `2`.
#'
#' @return The input dataframe with three new columns:
#'   - `continuous`: `1` if `time.diff <= threshold`, `0` otherwise (first
#'     detection per individual is set to `0`).
#'   - `gap`: logical, `TRUE` if `time.diff > threshold`.
#'   - `run.id`: integer run identifier that increments each time `continuous`
#'     changes value (i.e., each contiguous detection run has a unique ID).
#'
#' @examples
#' \dontrun{
#' sp <- collapse_motus_time(sparrow52550)
#' sp <- add_signal_diffs(sp)
#' sp <- add_continuity_flags(sp, threshold = 2)
#' }
#'
#' @importFrom dplyr group_by arrange mutate ungroup
#' @importFrom data.table rleid
#' @export
add_continuity_flags <- function(data,
                                 id_col = "tagDeployID",
                                 time_col = "time",
                                 time_diff_col = "time.diff",
                                 threshold = 2) {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(.data[[time_col]], .by_group = TRUE) %>%
    dplyr::mutate(
      continuous = ifelse(.data[[time_diff_col]] <= threshold, 1, 0),
      continuous = ifelse(is.na(continuous), 0, continuous),
      gap = .data[[time_diff_col]] > threshold
    ) %>%
    dplyr::mutate(
      run.id = data.table::rleid(continuous)
    ) %>%
    dplyr::ungroup()
}


#' Compute a centered rolling median of signal differences
#'
#' Smooths the signal difference series using a centered rolling median window.
#' This suppresses micro-spikes caused by antenna switching or brief movements
#' while preserving genuine behavioral regime shifts (e.g., settling into roost
#' vs. active flight).
#'
#' At a ~20–25 second ping interval, the default window of `k = 25`
#' observations covers approximately 8–10 minutes. Window size must be odd for
#' symmetric centering. The first and last `floor(k/2)` rows per individual
#' will be `NA` due to the centered alignment.
#'
#' @param data A dataframe of Motus detections with a signal difference column.
#' @param id_col Character. Individual identifier column. Default `"tagDeployID"`.
#' @param time_col Character. Timestamp column. Default `"time"`.
#' @param value_col Character. Column to smooth. Default `"sig.diff.mean"`.
#' @param k Integer. Rolling window width in number of observations. Must be
#'   odd. Default `25` (≈ 8–10 min at a 20 sec ping interval).
#' @param new_col Character. Name of the output smoothed column.
#'   Default `"roll_vol"`.
#'
#' @return The input dataframe with one new column (`new_col`) containing the
#'   centered rolling median. Edge values within `floor(k/2)` observations of
#'   the start or end of each individual's series are `NA`.
#'
#' @examples
#' \dontrun{
#' sp <- collapse_motus_time(sparrow52550)
#' sp <- add_signal_diffs(sp)
#' sp <- add_roll_median(sp, k = 25)
#' }
#'
#' @importFrom dplyr group_by arrange mutate ungroup
#' @importFrom zoo rollmedian
#' @export
add_roll_median <- function(data,
                            id_col = "tagDeployID",
                            time_col = "time",
                            value_col = "sig.diff.mean",
                            k = 25,
                            new_col = "roll_vol") {
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::arrange(.data[[time_col]], .by_group = TRUE) %>%
    dplyr::mutate(
      !!new_col := zoo::rollmedian(
        .data[[value_col]],
        k = k,
        fill = NA,
        align = "center"
      )
    ) %>%
    dplyr::ungroup()
}
