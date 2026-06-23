test_that("collapse_motus_time reduces rows and adds sig.mean", {
  data(sparrow52550)
  out <- collapse_motus_time(sparrow52550)
  expect_true(nrow(out) <= nrow(sparrow52550))
  expect_true("sig.mean" %in% names(out))
})

test_that("add_signal_diffs adds expected columns with valid values", {
  data(sparrow52550)
  sp  <- collapse_motus_time(sparrow52550)
  out <- add_signal_diffs(sp)
  expect_true(all(c("sig.diff", "sig.diff.mean", "time.diff") %in% names(out)))
  expect_true(all(out$sig.diff >= 0, na.rm = TRUE))
  expect_true(all(out$sig.diff.mean >= 0, na.rm = TRUE))
})

test_that("add_continuity_flags adds continuous, gap, and run.id columns", {
  data(sparrow52550)
  sp  <- collapse_motus_time(sparrow52550)
  sp  <- add_signal_diffs(sp)
  out <- add_continuity_flags(sp, threshold = 2)
  expect_true(all(out$continuous %in% c(0, 1)))
  expect_equal(out$continuous[1], 0)
  expect_type(out$gap, "logical")
  expect_type(out$run.id, "integer")
  expect_true(all(diff(out$run.id) >= 0))
})

test_that("add_roll_median adds roll_vol column with non-negative values", {
  data(sparrow52550)
  sp  <- collapse_motus_time(sparrow52550)
  sp  <- add_signal_diffs(sp)
  out <- add_roll_median(sp, window_min = 15)
  expect_true("roll_vol" %in% names(out))
  expect_true(all(out$roll_vol >= 0, na.rm = TRUE))
})
