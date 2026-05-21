###############################################################
# roost_workflow.R
# Companion workflow script for the roostR package.
# run Data.Prep.R first 
# Runs the full roost pipeline on per-bird CSV files from junco.clean/ created Data.Prep.R
# and saves intermediate + final results.
#
# Pipeline:
#   Step 1 – Preprocess:  collapse, signal diffs, continuity, rolling median
#   Step 2 – Roost times: detect onset/departure, add to data, add roost hours
#   Step 3 – Night obs:   compute observed time per night
#   Step 4 – Restless:    calc bouts, spike flags, restless rates
#
# Outputs:
#   data/step1_preproc/<bird_id>.rds
#   data/step2_roost/<bird_id>.rds
#   data/step3_metrics/<bird_id>.rds
#   data/step4_restless/<bird_id>.rds
#   data/results/<bird_id>.rds
#   data/results/junco_all.csv  (compiled summary, run Section 6)
#
# K. Williams
###############################################################

#### libraries ####
library(roostR)

library(dplyr)
library(lubridate)
library(readr)
library(zoo)
library(data.table)
library(slider)
library(purrr)
library(ggplot2)
library(viridis)


#### create output folders ####
  dir.create("data/step1_preproc",  recursive = TRUE)
  dir.create("data/step2_roost",    recursive = TRUE)
  dir.create("data/step3_metrics",  recursive = TRUE)
  dir.create("data/step4_restless", recursive = TRUE)
  dir.create("data/results",        recursive = TRUE)
  dir.create("figures/roll_vol",    recursive = TRUE)
  dir.create("figures/overnight",   recursive = TRUE)
  dir.create("figures/roost_timing",recursive = TRUE)
  dir.create("figures/leave_timing",recursive = TRUE)

#### set parameters ####
window_min      <- 15   # rolling median time window in minutes (add_roll_median)
vol_threshold   <- 3    # rolling vol threshold for roost onset
min_duration    <- 10   # min low-vol run (mins) to call roost onset
window_minutes  <- 90   # search window around sunset/sunrise (mins)
spike_threshold <- 4    # sig.diff threshold for restless spikes
gap_min         <- 5    # gap (mins) between spikes to split bouts
min_bout_min    <- 22 / 60  # minimum bout duration (mins); single spikes get one ping interval
bin_size_min    <- 5  # bin width (mins) for observed_time_hr: a bin is counted if >= 1 detection falls in it; observed_time_hr = n_bins * bin_size_min / 60. NB 2 min seems to underestimate still 10 minutes may overestimate and lead to prop time obs > 1
confirm_min     <- 60   # confirmation window duration (mins): after a candidate roost onset, the bird must remain quiet for this long or the candidate is rejected. 
confirm_frac    <- 0.8  # confirmation window fraction: minimum proportion of detections in the confirm_min window that must be low-volatility (roll_vol < vol_threshold) for the onset to be accepted


#### Step 1: Preprocess ####

files <- list.files("data/junco.clean", pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nStep 1:", bird_id)

  junco <- read_csv(f, show_col_types = FALSE) |>
    mutate(time = as_datetime(time))

  junco <- junco |>
    collapse_motus_time() |>
    add_signal_diffs() |>          # grouped by recvDeployName + port
    add_continuity_flags() |>
    add_roll_median(window_min = window_min)

  saveRDS(junco, file.path("data/step1_preproc", paste0(bird_id, ".rds")))
}
cat("\nStep 1 complete.\n")


#### Figures: roll_vol ####
# One PNG per bird: rolling median signal volatility by hour, faceted by doy.
# Color = diel period, shape = receiver (useful when bird detected at >1 tower).
# Figure height scales with number of unique days (~3 in per row of 5 panels).

for (f in list.files("data/step1_preproc", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nFigure roll_vol:", bird_id)

  junco <- readRDS(f) |>
    add_day_night() |>
    mutate(hour = lubridate::hour(time) + lubridate::minute(time) / 60)

  n_doy      <- length(unique(junco$doy))
  fig_height <- max(4, ceiling(n_doy / 5) * 3)

  p <- ggplot(junco, aes(hour, roll_vol, color = diel, shape = recvDeployName)) +
    geom_point(size = 0.3, alpha = 0.4) +
    facet_wrap(~doy) +
    scale_color_manual(values = c(day = "#FDE725FF", night = "#440154FF")) +
    labs(x     = "Hour (UTC)",
         y     = "Rolling median signal difference",
         color = "Period",
         shape = "Receiver") +
    theme_bw()

  ggsave(
    file.path("figures/roll_vol", paste0(bird_id, "_roll_vol.png")),
    plot   = p,
    width  = 14,
    height = fig_height,
    units  = "in",
    dpi    = 150,
    limitsize =  FALSE
  )
}
cat("\nroll_vol figures complete.\n")


#### Step 2: Detect roost onset and departure ####

for (f in list.files("data/step1_preproc", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nStep 2:", bird_id)

  junco <- readRDS(f)

  # Combine signal across receivers before detection.
  # At each unique (doy, time), take max roll_vol and max sig.diff so that
  # activity on any receiver counts. Detection then produces one roost_time
  # and one leave_roost_time per night with no per-receiver collapsing needed.
  junco_combined <- junco |>
    group_by(doy, time) |>
    summarise(
      roll_vol = max(roll_vol, na.rm = TRUE),
      sig.diff = max(sig.diff, na.rm = TRUE),
      Sunset   = first(Sunset),
      Sunrise  = first(Sunrise),
      .groups  = "drop"
    ) |>
    mutate(recvDeployName = "combined")

  roost_times <- detect_roost_onset(
    junco_combined,
    window_minutes = window_minutes,
    vol_threshold  = vol_threshold,
    min_duration   = min_duration,
    confirm_min    = confirm_min,
    confirm_frac   = confirm_frac
  )

  leave_times <- detect_roost_departure(
    junco_combined,
    window_before   = window_minutes,
    window_after    = window_minutes,
    spike_threshold = spike_threshold
  )

  junco <- junco |>
    add_roost_times(roost_times, leave_times) |>
    add_roost_hours()   # adds roost_hour, leave_roost_hour, sunset_hour, sunrise_hour

  saveRDS(junco, file.path("data/step2_roost", paste0(bird_id, ".rds")))
}
cat("\nStep 2 complete.\n")

summary(junco$roost_hour)
summary(junco$leave_roost_hour)

#### Figures: overnight ####
# One PNG per bird: sig.diff by wrapped hour (3:30pm–8:30am EST), faceted by night_doy.
# night_doy = doy of the evening side (doy N). Each panel shows the complete night:
#   evening detections (doy N) + post-midnight detections (doy N+1).
# night_doy is a plot-only label; underlying doy values are never modified.
# Reference lines sourced from the correct side of midnight:
#   sunset + roost onset from evening rows (doy N);
#   sunrise + departure from post-midnight rows (doy N+1, labeled as night N).
# NA wrapped hours (onset/departure not detected) produce no line — existing behavior.
# All receivers on the same panel; shape distinguishes receiver.

for (f in list.files("data/step2_roost", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nFigure overnight:", bird_id)

  junco <- readRDS(f)

  # wrap_hours_overnight adds roost_hour_wrap and leave_roost_hour_wrap (NA where not detected)
  overnight <- wrap_hours_overnight(junco, pivot_hour = 20) |>
    mutate(
      hour              = lubridate::hour(time)    + lubridate::minute(time)    / 60,
      sunset_hour       = lubridate::hour(Sunset)  + lubridate::minute(Sunset)  / 60,
      sunrise_hour      = lubridate::hour(Sunrise) + lubridate::minute(Sunrise) / 60,
      hour_wrap         = ifelse(hour        >= 20, hour,        hour        + 24),
      sunset_hour_wrap  = ifelse(sunset_hour >= 20, sunset_hour, sunset_hour + 24),
      sunrise_hour_wrap = ifelse(sunrise_hour >= 20, sunrise_hour, sunrise_hour + 24),
      night_doy         = if_else(hour < 20, doy - 1L, doy)
    ) |>
    filter(hour_wrap >= 20.5, hour_wrap <= 37.5)

  if (nrow(overnight) == 0) {
    cat(" [no overnight data, skipping]\n")
    next
  }

  # split on midnight so reference lines come from the correct calendar day
  evening <- overnight |> filter(hour >= 20)   # doy N: sunset, roost onset
  morning <- overnight |> filter(hour <  20)   # doy N+1: sunrise, departure

  sunset_lines  <- distinct(evening, night_doy, sunset_hour_wrap)
  roost_lines   <- distinct(evening, night_doy, roost_hour_wrap)
  sunrise_lines <- distinct(morning, night_doy, sunrise_hour_wrap)
  leave_lines   <- distinct(morning, night_doy, leave_roost_hour_wrap)

  p <- ggplot(overnight, aes(hour_wrap, sig.diff, shape = recvDeployName)) +
    geom_point(alpha = 0.4, size = 0.5) +
    geom_vline(data = sunset_lines,
               aes(xintercept = sunset_hour_wrap),
               linetype = 2, color = "#440154FF") +
    geom_vline(data = sunrise_lines,
               aes(xintercept = sunrise_hour_wrap),
               linetype = 2, color = "#95D840FF") +
    geom_vline(data = roost_lines,
               aes(xintercept = roost_hour_wrap),
               color = "#39568CFF") +
    geom_vline(data = leave_lines,
               aes(xintercept = leave_roost_hour_wrap),
               color = "#55C667FF") +
    facet_wrap(~night_doy, ncol = 10) +
    scale_x_continuous(
      breaks = seq(20, 38, by = 2),
      labels = function(x) x %% 24
    ) +
    coord_cartesian(ylim = c(0, 50), xlim = c(20.5, 37.5)) +
    geom_hline(yintercept = 5, color = "#FDE725FF", linetype = 2) +
    labs(x     = "Hour (3:30pm–8:30am EST)",
         y     = "Signal difference",
         shape = "Receiver") +
    ggtitle(bird_id) +
    theme_bw()

  ggsave(
    file.path("figures/overnight", paste0(bird_id, "_overnight.png")),
    plot      = p,
    width     = 14,
    height    = 10,
    units     = "in",
    dpi       = 150,
    limitsize = FALSE
  )
}
cat("\nOvernight figures complete.\n")


#### Figures: roost_timing and leave_timing ####
# One PNG per bird saved to figures/roost_timing/ and figures/leave_timing/.
# Times converted to EST (America/New_York). One point per detected night.
# Dashed line shows the solar event (sunset / sunrise) across the season.

for (f in list.files("data/step2_roost", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nFigure timing:", bird_id)

  junco <- readRDS(f)

  night_est <- junco |>
    mutate(
      roost_time_est   = lubridate::with_tz(roost_time,       "America/New_York"),
      leave_time_est   = lubridate::with_tz(leave_roost_time, "America/New_York"),
      sunset_est       = lubridate::with_tz(Sunset,           "America/New_York"),
      sunrise_est      = lubridate::with_tz(Sunrise,          "America/New_York"),
      roost_hour_est   = lubridate::hour(roost_time_est) +
                         lubridate::minute(roost_time_est) / 60,
      leave_hour_est   = lubridate::hour(leave_time_est) +
                         lubridate::minute(leave_time_est) / 60,
      sunset_hour_est  = lubridate::hour(sunset_est) +
                         lubridate::minute(sunset_est)  / 60,
      sunrise_hour_est = lubridate::hour(sunrise_est) +
                         lubridate::minute(sunrise_est) / 60
    ) |>
    distinct(recvDeployName, doy, roost_hour_est, leave_hour_est,
             sunset_hour_est, sunrise_hour_est)

  # roost time vs sunset
  p_roost <- ggplot(night_est, aes(doy)) +
    geom_line(aes(y = sunset_hour_est), linetype = "dashed") +
    geom_point(aes(y = roost_hour_est)) +
    labs(x = "doy",
         y = "Hour (America/New_York)",
         title = paste(bird_id, "- Roost time vs. sunset")) +
    theme_bw()

  ggsave(
    file.path("figures/roost_timing", paste0(bird_id, "_roost_timing.png")),
    plot   = p_roost,
    width  = 7,
    height = 5,
    units  = "in",
    dpi    = 150
  )

  # leave roost time vs sunrise
  p_leave <- ggplot(night_est, aes(doy)) +
    geom_line(aes(y = sunrise_hour_est), linetype = "dashed") +
    geom_point(aes(y = leave_hour_est)) +
    labs(x = "doy",
         y = "Hour (America/New_York)",
         title = paste(bird_id, "- Leave roost time vs. sunrise")) +
    theme_bw()

  ggsave(
    file.path("figures/leave_timing", paste0(bird_id, "_leave_timing.png")),
    plot   = p_leave,
    width  = 7,
    height = 5,
    units  = "in",
    dpi    = 150
  )
}
cat("\nTiming figures complete.\n")


#### Step 3: Compute night observation ####

for (f in list.files("data/step2_roost", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nStep 3:", bird_id)

  junco <- readRDS(f)

  night_metrics <- compute_night_observation(junco, bin_size_min = bin_size_min)

  junco <- junco |>
    select(-any_of(setdiff(names(night_metrics), "doy"))) |>
    left_join(night_metrics, by = "doy")

  saveRDS(junco, file.path("data/step3_metrics", paste0(bird_id, ".rds")))
}
cat("\nStep 3 complete.\n")

# check
list.files("data/step1_preproc", pattern = "\\.rds$", full.names = TRUE) |>
  map(readRDS) |>
  bind_rows() |>
  filter(!is.na(time.diff)) |>
  group_by(recvDeployName) |>
  summarise(
    pct_gap_over5  = mean(time.diff > 5,  na.rm = TRUE),
    pct_gap_over10 = mean(time.diff > 10, na.rm = TRUE),
    pct_gap_over20 = mean(time.diff > 20, na.rm = TRUE)
  ) |>
  print(n = Inf)


#### Step 4: Restless bouts and rates ####

for (f in list.files("data/step3_metrics", pattern = "\\.rds$", full.names = TRUE)) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  cat("\nStep 4:", bird_id)

  junco <- readRDS(f)

  restless_summary <- calc_restless_all(
    junco,
    spike_threshold = spike_threshold,
    gap_min         = gap_min,
    min_bout_min    = min_bout_min
  )

  junco <- junco |>
    select(-any_of(setdiff(names(restless_summary), "doy"))) |>
    left_join(restless_summary, by = "doy")

  junco <- junco |>
    select(-any_of(c("interval_start", "interval_end", "in_roost_interval",
                     "spike", "dt", "new_bout", "bout_id"))) |>
    add_spike_bouts(spike_threshold = spike_threshold, gap_min = gap_min)

  restless_rates <- calc_restless_rates(junco)

  junco <- junco |>
    select(-any_of(setdiff(names(restless_rates), "doy"))) |>
    left_join(restless_rates, by = "doy")

  saveRDS(junco, file.path("data/step4_restless", paste0(bird_id, ".rds")))
  saveRDS(junco, file.path("data/results",        paste0(bird_id, ".rds")))
}
cat("\nStep 4 complete.\n")

names(junco)

#### Step 5 (optional): Compile per-night summary across all birds ####
# Produces one row per bird per night with summary columns only.
# Run after Step 4 is finished for all birds.

summary_cols <- c(
  "tagDeployID", "recvDeployName", "doy",
  "roost_time", "leave_roost_time",
  "roost_hour", "leave_roost_hour",
  "sunset_hour", "sunrise_hour",
  "time_roosting_hr", "observed_time_hr", "prop_time_observed",
  "n_bouts", "total_restless_min", "max_bout_min",
  "total_restless_hr", "restless_per_obs_hr", "prop_time_restless"
)

result_files <- list.files("data/results", pattern = "\\.rds$", full.names = TRUE)

junco_all <- lapply(result_files, function(f) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  dat <- readRDS(f)

  bearing_summary <- calc_bearing_summary(dat)

  dat |>
    distinct(across(any_of(summary_cols))) |>
    left_join(bearing_summary, by = "doy") |>
    mutate(bird_id = bird_id)
}) |>
  bind_rows()

write_csv(junco_all, "data/results/junco_all.csv")
cat("\nCompiled", nrow(junco_all), "bird-night rows to data/results/junco_all.csv\n")

#### check data and variable summaries before analyses #### 
# Load final results
junco_final <- readRDS("data/results/Sparrow52550.dat.rds")

# Roost near sunset, departure near sunrise (EST ~18-22h roost, ~5-8h departure)
summary(junco_final$roost_hour)
summary(junco_final$leave_roost_hour)
summary(junco_final$prop_time_restless)

# Check key columns are present
names(junco_final)

#2. Verify figures were created
list.files("figures", recursive = TRUE)

#3. Check the compiled summary
junco_all <- read.csv("data/junco_all.csv")  # if Step 5 was run
nrow(junco_all)
head(junco_all)
