build_sp_with_roost <- function() {
  data(sparrow52550)
  sp <- collapse_motus_time(sparrow52550)
  sp <- add_signal_diffs(sp)
  sp <- add_continuity_flags(sp, threshold = 2)
  sp <- add_roll_median(sp, window_min = 15)
  sp <- add_day_night(sp)
  roost_times <- detect_roost_onset(sp)
  leave_times <- detect_roost_departure(sp)
  sp <- add_roost_times(sp, roost_times, leave_times)
  sp <- add_roost_hours(sp)
  sp
}

test_that("compute_night_observation returns expected columns", {
  sp  <- build_sp_with_roost()
  out <- compute_night_observation(sp, bin_size_min = 2)
  expect_true(all(c("time_roosting_hr", "n_bins_detected",
                    "observed_time_hr", "prop_time_observed") %in% names(out)))
})

test_that("compute_night_observation prop_time_observed is between 0 and 1", {
  sp  <- build_sp_with_roost()
  out <- compute_night_observation(sp, bin_size_min = 2)
  expect_true(all(out$prop_time_observed >= 0 & out$prop_time_observed <= 1,
                  na.rm = TRUE))
})

test_that("compute_night_observation time_roosting_hr is positive", {
  sp  <- build_sp_with_roost()
  out <- compute_night_observation(sp, bin_size_min = 2)
  expect_true(all(out$time_roosting_hr > 0, na.rm = TRUE))
})
