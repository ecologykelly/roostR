#' @keywords internal
"_PACKAGE"

## Suppress R CMD check NOTEs for the pipe operator and data.table walrus operator
#' @importFrom magrittr %>%
#' @importFrom data.table :=
NULL

## Suppress R CMD check NOTEs for NSE column names used as bare symbols inside
## dplyr verbs (mutate, filter, summarise, group_by, arrange, distinct).
utils::globalVariables(c(
  # rlang data pronoun
  ".data",

  # add_signal_diffs
  "sig.diff", "sig.diff.mean",

  # add_continuity_flags
  "continuous",

  # add_day_night
  "diel",

  # detect_roost_onset
  "window_start", "window_end",

  # detect_roost_departure
  "forward_med", "departure_signal",

  # compute_night_observation — hard-coded bare column names
  "tagDeployID", "doy", "roost_time", "leave_roost_time",

  # shared intermediates: compute_night_observation + restlessness
  "leave_next", "interval_start", "interval_end", "time_roosting_hr",
  "dt", "new_segment", "segment_id",
  "seg_start", "seg_end", "seg_dur_hr",

  # add_spike_bouts / calc_restless_all
  "in_roost_interval", "spike", "new_bout", "bout_id",
  "bout_start", "bout_end", "duration_min",

  # calc_restless_rates
  "total_restless_hr"
))
