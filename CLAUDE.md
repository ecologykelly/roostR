# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

**roostR** converts raw Motus radio-telemetry detections into behavioral metrics
for nocturnally roosting birds. Signal strength changes between consecutive
detections are used as an activity proxy — a stationary roosting bird produces
near-zero changes; a moving bird produces large changes. The package detects
roost onset, roost departure, and restlessness bouts from these patterns.

The source analysis project this package was built from is at
`C:\Users\willi\Desktop\ActivityPatternsinR`.

## Developer commands

All commands are run in RStudio with this project open, or via `devtools` in
any R session pointed at this directory.

```r
devtools::document()       # regenerate NAMESPACE and man/ from roxygen2 comments
devtools::check()          # full R CMD check — must pass before pushing
devtools::install()        # install locally for testing
devtools::build_vignettes() # build vignette HTML manually
```

**First-time setup only:**
```r
source("data-raw/prepare_data.R")  # generates data/sparrow52550.rda
```

## Publishing to GitHub

```r
# In a terminal inside C:\Users\willi\Desktop\roostR:
git init
git add .
git commit -m "Initial roostR package"
# Create a new repo named roostR on github.com, then:
git remote add origin https://github.com/ecologykelly/roostR.git
git push -u origin main
```

After pushing, install from GitHub with:
```r
remotes::install_github("YOUR_USERNAME/roostR")
```

## Package structure

```
R/
├── utils.R           # convert_to_posix, wrap_hours_overnight
├── preprocess.R      # collapse_motus_time, add_signal_diffs,
│                     #   add_continuity_flags, add_roll_median
├── diel.R            # add_day_night
├── roost_detection.R # detect_roost_onset, detect_roost_departure,
│                     #   add_roost_times, add_roost_hours
├── night_metrics.R   # compute_night_observation
├── restlessness.R    # calc_restless_all, add_spike_bouts, calc_restless_rates
└── data.R            # roxygen docs for sparrow52550 dataset
```

## Analysis pipeline and function order

Functions must be called in this sequence — each step's output is required
by the next:

1. `collapse_motus_time()` — one row per individual × timestamp (removes multi-antenna duplicates)
2. `add_signal_diffs()` — adds `sig.diff`, `sig.diff.mean`, `time.diff`
3. `add_continuity_flags()` — adds `continuous`, `gap`, `run.id`
4. `add_roll_median()` — adds `roll_vol` (centered rolling median, default k=25)
5. `add_day_night()` — adds `diel` factor (`"day"` / `"night"`)
6. `detect_roost_onset()` → `detect_roost_departure()` → `add_roost_times()` → `add_roost_hours()`
7. `wrap_hours_overnight()` — for noon-to-noon plotting only
8. `compute_night_observation()` — returns summary df; left-join back to main data
9. `calc_restless_all()` — returns summary df; left-join back to main data
10. `add_spike_bouts()` — annotates spike bouts on the detection-level df
11. `calc_restless_rates()` — returns summary df; left-join back to main data

## Key design decisions

- **`tagDeployID` is the individual key**, not `motusTagID` — handles tags
  redeployed on different animals.
- **All timestamps are UTC.** Sunrise/sunset times must be converted from local
  time before use (see `Activity.Data.Prep.R` in the source project).
- **Overnight intervals span two DOYs** — roost onset is on DOY *d*, departure
  is on DOY *d+1*. Functions pair them with `dplyr::lead(leave_roost_time)`.
- **No bare `library()` calls inside functions** — all dependencies are declared
  in `DESCRIPTION` Imports and called with `package::function()` or via
  `@importFrom` roxygen tags.

## Dependencies

**Imports** (required, in DESCRIPTION):
- `dplyr` — data manipulation throughout
- `lubridate` — POSIXct arithmetic (`minutes()`, `hour()`, `with_tz()`)
- `zoo` — `rollmedian()` in `add_roll_median()`, `rollapply()` in `detect_roost_departure()`
- `data.table` — `rleid()` in `add_continuity_flags()` and `detect_roost_onset()`

**Suggests** (optional, for vignette and user scripts):
`ggplot2`, `viridis`, `tidyr`, `psych`, `knitr`, `rmarkdown`

## Default thresholds (tuned to song sparrow / Motus hardware)

All thresholds are exposed as function parameters and can be overridden.
These defaults were tuned to song sparrow data with ~20 sec ping intervals:

| Parameter | Default | Function | Rationale |
|---|---|---|---|
| Continuity gap | 2 min | `add_continuity_flags` | ~5× ping interval |
| Rolling window k | 25 obs | `add_roll_median` | ≈ 8–10 min smoothing |
| Roost vol threshold | 3 | `detect_roost_onset` | Night median ≈ 1, day ≈ 5 |
| Roost window | ±90 min | `detect_roost_onset` | Around sunset |
| Departure spike threshold | 4 | `detect_roost_departure` | Matches restlessness threshold |
| Departure window | −90/+120 min | `detect_roost_departure` | Around sunrise |
| Night gap threshold | 10 min | `compute_night_observation` | Out-of-range criterion |
| Restlessness spike | 4 | `calc_restless_all`, `add_spike_bouts` | Matched to departure threshold |
| Restlessness gap | 2 min | `calc_restless_all`, `add_spike_bouts` | Bout separation |

## Outstanding tasks

- **`add_day_night()`**: implemented from inline code in `Activity.R`. The user
  has a version in a separate file — replace the body in `R/diel.R` if needed,
  then re-run `devtools::document()` and `devtools::check()`.
- **DESCRIPTION**: update the `email` field in `Authors@R` before publishing.
- **Example data**: run `source("data-raw/prepare_data.R")` once to generate
  `data/sparrow52550.rda` before running `devtools::check()`.
