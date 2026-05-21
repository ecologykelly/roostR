# roostR

**roostR** is an R package for detecting and quantifying nocturnal roosting
behavior in birds tagged with [Motus](https://motus.org/) radio transmitters.
It converts raw signal-strength detections into behavioral metrics: roost onset
time, roost departure time, nightly observation coverage, and restlessness bout
statistics.

## Installation

```r
# Install from GitHub (includes vignette)
remotes::install_github("ecologykelly/roostR", build_vignettes = TRUE)
```

## Quick Start

```r
library(roostR)
library(dplyr)

# Example data: one dark-eyed junco, multiple nights
data(sparrow52550)

# 1. Collapse multi-antenna detections to one row per timestamp
sp <- collapse_motus_time(sparrow52550)

# 2. Compute activity proxies (signal differences, time gaps)
sp <- add_signal_diffs(sp)
sp <- add_continuity_flags(sp)

# 3. Smooth signal and classify diel period
sp <- add_roll_median(sp)
sp <- add_day_night(sp)

# 4. Detect roost timing
roost_times <- detect_roost_onset(sp)
leave_times <- detect_roost_departure(sp)
sp <- add_roost_times(sp, roost_times, leave_times)
sp <- add_roost_hours(sp)

# 5. Night observation metrics — use select(-any_of(...)) before each join
#    to prevent duplicate columns if the script is re-run
night_metrics <- compute_night_observation(sp, gap_threshold = 10)
sp <- sp |>
  select(-any_of(setdiff(names(night_metrics), c("tagDeployID", "doy")))) |>
  left_join(night_metrics, by = c("tagDeployID", "doy"))

# 6. Restlessness bouts
restless_summary <- calc_restless_all(sp, spike_threshold = 4, gap_min = 2)
sp <- sp |>
  select(-any_of(setdiff(names(restless_summary), c("tagDeployID", "doy")))) |>
  left_join(restless_summary, by = c("tagDeployID", "doy"))

# 7. Annotate spikes and compute rates
sp <- sp |>
  select(-any_of(c("interval_start", "interval_end", "in_roost_interval",
                   "spike", "dt", "new_bout", "bout_id"))) |>
  add_spike_bouts(spike_threshold = 4, gap_min = 2)

restless_rates <- calc_restless_rates(sp)
sp <- sp |>
  select(-any_of(setdiff(names(restless_rates), c("tagDeployID", "doy")))) |>
  left_join(restless_rates, by = c("tagDeployID", "doy"))
```

See `vignette("roostR-workflow")` for a full walkthrough with plots.

## Companion Scripts

Two ready-to-run scripts are bundled with the package. Copy them to your working directory to get started:

```r
file.copy(system.file("scripts/Data.Prep.R",     package = "roostR"), ".")
file.copy(system.file("scripts/roost_workflow.R", package = "roostR"), ".")
```

- **`Data.Prep.R`** — prepares raw Motus CSV files: joins sunrise/sunset times, filters low-detection receivers, and writes one clean CSV per bird to `data/junco.clean/`. Run this before `roost_workflow.R`.
- **`roost_workflow.R`** — runs the full roostR pipeline on cleaned CSVs and saves intermediate and final results to `data/` subfolders and `figures/`.

> **Note:** Processing time scales with file size and number of birds. Large detection files (many nights or high ping rates) can take several minutes per bird, particularly in Steps 2 and 4. For datasets of ~50 birds, full pipeline runs can take several hours — consider running each step section independently and allowing it to complete before proceeding, or running R in the background while working on other tasks.

## Input Data Format

roostR expects a dataframe of Motus detections with at minimum:

| Column | Type | Description |
|--------|------|-------------|
| `tagDeployID` | integer | Deployment-specific tag identifier |
| `time` | POSIXct (UTC) | Detection timestamp |
| `doy` | integer | Day of year |
| `sig` | numeric | Signal strength (dBm) |
| `Sunrise` | POSIXct (UTC) | Sunrise time for that day |
| `Sunset` | POSIXct (UTC) | Sunset time for that day |

Astronomical twilight columns (`AT.Start`, `AT.End`, `NT.Start`, `NT.End`,
`CT.Start`, `CT.End`) are optional but used by some downstream analyses.

## Functions

| Function | Purpose |
|----------|---------|
| `convert_to_posix()` | Convert DOY + time string to POSIXct |
| `collapse_motus_time()` | Collapse multi-antenna detections to one row per timestamp |
| `add_signal_diffs()` | Compute signal differences and inter-detection intervals |
| `add_continuity_flags()` | Flag gaps and assign run IDs |
| `add_roll_median()` | Centered rolling median smoothing |
| `add_day_night()` | Classify detections as day or night |
| `detect_roost_onset()` | Detect roost onset time near sunset |
| `detect_roost_departure()` | Detect roost departure time near sunrise |
| `add_roost_times()` | Join roost/departure times back to detection data |
| `add_roost_hours()` | Convert roost times to decimal hour-of-day |
| `wrap_hours_overnight()` | Shift UTC hours for noon-to-noon plotting |
| `compute_night_observation()` | Compute roosting duration and observation coverage |
| `calc_restless_all()` | Count restlessness bouts during roost interval |
| `add_spike_bouts()` | Annotate spike bouts in the detection dataframe |
| `calc_restless_rates()` | Normalize restlessness to per-observed-hour rates |
