test_that("add_day_night returns factor with levels day and night", {
  data(sparrow52550)
  sp  <- collapse_motus_time(sparrow52550)
  out <- add_day_night(sp)
  expect_true("diel" %in% names(out))
  expect_s3_class(out$diel, "factor")
  expect_equal(levels(out$diel), c("day", "night"))
  expect_true(all(out$diel %in% c("day", "night")))
})

test_that("add_day_night classifies night detections correctly", {
  data(sparrow52550)
  sp  <- collapse_motus_time(sparrow52550)
  out <- add_day_night(sp)
  night_rows <- out[out$diel == "night", ]
  expect_true(all(night_rows$time < night_rows$Sunrise |
                  night_rows$time >= night_rows$Sunset))
})
