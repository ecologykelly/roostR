build_sp_full <- function() {
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
  night_metrics <- compute_night_observation(sp, bin_size_min = 2)
  sp <- sp |>
    dplyr::select(-dplyr::any_of(setdiff(names(night_metrics), "doy"))) |>
    dplyr::left_join(night_metrics, by = "doy")
  sp
}

test_that("calc_restless_all returns n_bouts >= 0", {
  sp  <- build_sp_full()
  out <- calc_restless_all(sp, spike_threshold = 4, gap_min = 5,
                            min_bout_min = 22 / 60)
  expect_true("n_bouts" %in% names(out))
  expect_true(all(out$n_bouts >= 0, na.rm = TRUE))
})

test_that("calc_restless_rates returns non-negative rate columns", {
  sp  <- build_sp_full()
  restless <- calc_restless_all(sp, spike_threshold = 4, gap_min = 5,
                                 min_bout_min = 22 / 60)
  sp <- sp |>
    dplyr::select(-dplyr::any_of(setdiff(names(restless), "doy"))) |>
    dplyr::left_join(restless, by = "doy")
  out <- calc_restless_rates(sp)
  expect_true(all(c("prop_time_restless", "restless_per_obs_hr") %in% names(out)))
  expect_true(all(out$prop_time_restless >= 0, na.rm = TRUE))
  expect_true(all(out$restless_per_obs_hr >= 0, na.rm = TRUE))
})

test_that("add_spike_bouts adds logical spike column", {
  sp  <- build_sp_full()
  out <- add_spike_bouts(sp, spike_threshold = 4, gap_min = 5)
  expect_true("spike" %in% names(out))
  expect_type(out$spike, "logical")
})
