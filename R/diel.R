#' Classify detections as day or night
#'
#' Labels each detection as `"day"` or `"night"` based on per-row sunrise and
#' sunset times. A detection is "day" if its timestamp falls on or after
#' sunrise and before sunset on the same UTC date; all other detections are
#' "night". The resulting column is a factor with levels `c("day", "night")`.
#'
#' All times (`time`, `Sunrise`, `Sunset`) must be POSIXct in UTC.
#'
#' @param data A dataframe of Motus detections containing timestamp,
#'   sunrise, and sunset columns.
#' @param time_col Character. POSIXct timestamp column (UTC). Default `"time"`.
#' @param sunrise_col Character. POSIXct sunrise column (UTC). Default
#'   `"Sunrise"`.
#' @param sunset_col Character. POSIXct sunset column (UTC). Default `"Sunset"`.
#'
#' @return The input dataframe with one new column:
#'   - `diel`: factor with levels `"day"` and `"night"`.
#'
#' @examples
#' \dontrun{
#' # After joining sunrise/sunset times to your detections:
#' sp <- add_day_night(sp)
#' table(sp$diel)
#' }
#'
#' @importFrom dplyr mutate if_else
#' @export
add_day_night <- function(data,
                          time_col    = "time",
                          sunrise_col = "Sunrise",
                          sunset_col  = "Sunset") {
  data %>%
    dplyr::mutate(
      diel = dplyr::if_else(
        .data[[time_col]] >= .data[[sunrise_col]] &
          .data[[time_col]] <  .data[[sunset_col]],
        "day",
        "night"
      ),
      diel = factor(diel, levels = c("day", "night"))
    )
}
