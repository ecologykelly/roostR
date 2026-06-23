build_sp <- function() {
  data(sparrow52550)
  sp <- collapse_motus_time(sparrow52550)
  sp <- add_signal_diffs(sp)
  sp <- add_continuity_flags(sp, threshold = 2)
  sp <- add_roll_median(sp, window_min = 15)
  sp <- add_day_night(sp)
  sp
}

test_that("detect_roost_onset returns roost_time column as POSIXct", {
  sp  <- build_sp()
  out <- detect_roost_onset(sp, window_minutes = 90, vol_threshold = 3,
                             min_duration = 10, confirm_min = 30, confirm_frac = 0.8)
  expect_true("roost_time" %in% names(out))
  expect_s3_class(out$roost_time, "POSIXct")
  expect_true(nrow(out) > 0)
})

test_that("detect_roost_departure returns leave_roost_time column as POSIXct", {
  sp  <- build_sp()
  out <- detect_roost_departure(sp, window_before = 90, window_after = 120,
                                 spike_threshold = 4, median_threshold = 3)
  expect_true("leave_roost_time" %in% names(out))
  expect_s3_class(out$leave_roost_time, "POSIXct")
  expect_true(nrow(out) > 0)
})

test_that("add_roost_times adds roost_time and leave_roost_time to detections", {
  sp          <- build_sp()
  roost_times <- detect_roost_onset(sp)
  leave_times <- detect_roost_departure(sp)
  out         <- add_roost_times(sp, roost_times, leave_times)
  expect_true(all(c("roost_time", "leave_roost_time") %in% names(out)))
  expect_equal(nrow(out), nrow(sp))
})

test_that("add_roost_hours adds numeric hour columns in range 0-24", {
  sp          <- build_sp()
  roost_times <- detect_roost_onset(sp)
  leave_times <- detect_roost_departure(sp)
  sp          <- add_roost_times(sp, roost_times, leave_times)
  out         <- add_roost_hours(sp)
  expect_true(all(c("roost_hour", "leave_roost_hour",
                    "sunset_hour", "sunrise_hour") %in% names(out)))
  expect_true(all(out$sunset_hour  >= 0 & out$sunset_hour  <= 24, na.rm = TRUE))
  expect_true(all(out$sunrise_hour >= 0 & out$sunrise_hour <= 24, na.rm = TRUE))
})
