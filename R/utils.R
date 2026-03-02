#' Convert day-of-year and time string to POSIXct
#'
#' Combines a day-of-year integer with a 12-hour time string (e.g. `"6:30 AM"`)
#' to produce a POSIXct datetime. Used to convert sunrise/sunset tables that
#' store times as character strings alongside a day-of-year column.
#'
#' @param time_str Character. Time in `"H:MM AM/PM"` format. Whitespace is
#'   trimmed automatically. Returns `NA` if `time_str` is `NA` or empty.
#' @param doy Integer. Day of year (1–366).
#' @param year Integer. Calendar year used to resolve the date. Default `2024`.
#'
#' @return A character representation of the resulting POSIXct datetime
#'   (suitable for storage in a dataframe column before final conversion), or
#'   `NA` if `time_str` is missing.
#'
#' @examples
#' convert_to_posix("6:30 AM", doy = 120, year = 2024)
#'
#' @export
convert_to_posix <- function(time_str, doy, year = 2024) {
  time_str <- trimws(time_str)
  if (is.na(time_str) || time_str == "") {
    return(NA)
  }
  date_str <- as.Date(doy - 1, origin = paste0(year, "-01-01"))
  datetime_str <- paste(date_str, time_str)
  posix_time <- as.POSIXct(datetime_str,
                            format = "%Y-%m-%d %I:%M %p",
                            tz = "UTC")
  return(as.character(posix_time))
}


#' Wrap UTC hour-of-day values for overnight plotting
#'
#' Shifts hour-of-day values so that a chosen UTC pivot hour becomes the left
#' edge of a plot axis. This is useful for visualizing overnight detections
#' (e.g. dusk to dawn) without a break at midnight. Hours before the pivot
#' are increased by 24 so they appear to the right of midnight on the axis.
#'
#' The default `pivot_hour = 17` corresponds to 17:00 UTC ≈ noon EST, placing
#' noon at the left edge and producing a noon-to-noon view.
#'
#' @param data A dataframe containing one or more numeric hour-of-day columns.
#' @param pivot_hour Numeric. UTC hour that becomes the left edge of the plot
#'   axis. Hours strictly less than `pivot_hour` are shifted +24.
#'   Default `17` (≈ noon Eastern Standard Time).
#' @param hour_cols Character vector. Names of columns to wrap. Columns not
#'   present in `data` are silently skipped. Default wraps `"hour"`,
#'   `"sunset_hour"`, `"sunrise_hour"`, `"roost_hour"`, and
#'   `"leave_roost_hour"`.
#'
#' @return The input dataframe with additional `*_wrap` columns for each
#'   column in `hour_cols` that exists in `data`. For example, `"hour"`
#'   produces a new `"hour_wrap"` column.
#'
#' @examples
#' \dontrun{
#' overnight <- wrap_hours_overnight(sparrow_processed, pivot_hour = 17)
#' }
#'
#' @export
wrap_hours_overnight <- function(data,
                                 pivot_hour = 17,
                                 hour_cols = c("hour",
                                               "sunset_hour",
                                               "sunrise_hour",
                                               "roost_hour",
                                               "leave_roost_hour")) {
  for (col in hour_cols) {
    if (col %in% names(data)) {
      new_col <- paste0(col, "_wrap")
      data[[new_col]] <- ifelse(
        !is.na(data[[col]]) & data[[col]] < pivot_hour,
        data[[col]] + 24,
        data[[col]]
      )
    }
  }
  return(data)
}
