##################################
# Data prep.
# data are in data/junco.data = 1 csv per bird
# Workflow:
#   1. Check towers: review detection counts per receiver for each bird
#   2. Set min_detections threshold for receiver filtering
#   3. Process all birds: add doy/year, join sun times (UTC),
#      remove NA receivers and low-detection receivers, reduce columns,
#      save to data/junco.clean/
#   4. Verify cleaned files: identify birds with >1 tower remaining
#      (important - these birds need special handling in roostR)


#### load libraries ####
library(roostR)

library(dplyr)
library(lubridate)
library(readr)


# prep_sun_times() and add_sun_times() are provided by the roostR package


#### columns to keep ####
keep_cols <- c("sig", "noise", "motusTagID", "ambigID", "port", "runLen",
               "motusFilter", "mfgID", "tagDeployID", "fullID",
               "recvDeployName", "antBearing", "speciesEN",
               "time", "year", "doy",
               "Sunrise", "Sunset",
               "AT.Start", "AT.End", "NT.Start", "NT.End", "CT.Start", "CT.End")


#### Step 1: check towers ####
# run this section first to identify birds needing receiver removal

files <- list.files("data-raw", pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("sunrisetimes", files)]

# print sorted tower detection counts for every bird; collect multi-tower birds
multi_tower <- character(0)
for (f in files) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  dat <- read_csv(f, col_select = c("recvDeployName", "ambigID"), show_col_types = FALSE)
  cat("\n---", bird_id, "---\n")
  summary_tbl <- dat |>
    group_by(recvDeployName) |>
    summarise(n = n(), n_ambig = sum(!is.na(ambigID)), .groups = "drop") |>
    arrange(n)
  print(summary_tbl)
  if (length(unique(dat$recvDeployName)) > 1) multi_tower <- c(multi_tower, bird_id)
}
cat("\nBirds with >1 tower:\n")
print(multi_tower)


#### Step 2: receiver filtering criteria ####
# Review Step 1 output before running Step 3.
# Rows with NA receiver are always removed.
# Any receiver with fewer than min_detections rows for a given bird is removed.
# Adjust min_detections as needed.

min_detections <- 100   # change this threshold if needed


#### Step 3: process and save to junco.clean ####
dir.create("data/junco.clean", recursive = TRUE, showWarnings = FALSE)

# only read columns needed from raw files - avoids parsing warnings on unused columns
raw_cols <- c("sig", "noise", "motusTagID", "ambigID", "port", "runLen",
              "motusFilter", "mfgID", "tagDeployID", "fullID",
              "recvDeployName", "antBearing", "speciesEN", "time")

for (f in files) {
  bird_id <- tools::file_path_sans_ext(basename(f))

  junco <- read_csv(f, col_select = all_of(raw_cols), show_col_types = FALSE) |>
    mutate(
      time = as_datetime(time, tz = "UTC"),
      year = year(time),
      doy  = yday(time)
    )

  ref_year <- as.integer(format(min(junco$time, na.rm = TRUE), "%Y"))
  sun <- prep_sun_times("data-raw/sunrisetimes.csv", ref_year = ref_year)

  junco <- junco |> add_sun_times(sun)

  junco <- junco |>
    filter(!is.na(recvDeployName), !is.na(port), !is.na(time)) |>
    group_by(recvDeployName) |>
    filter(n() >= min_detections) |>
    ungroup() |>
    select(any_of(keep_cols))

  write_csv(junco, file.path("data/junco.clean", paste0(bird_id, ".csv")))
}


#### Step 4: verify cleaned files - check for remaining multiple towers ####
# Birds with >1 tower after cleaning will need special handling in roostR
# (e.g. group_by tower before calculating signal differences)
# If additional problems are identified here, adjust min_detections in Step 2
# and re-run Step 3. junco.clean files will be overwritten.

clean_files <- list.files("data/junco.clean", pattern = "\\.csv$", full.names = TRUE)

multi_tower_clean <- character(0)
for (f in clean_files) {
  bird_id <- tools::file_path_sans_ext(basename(f))
  dat <- read_csv(f, col_select = c("recvDeployName", "ambigID"), show_col_types = FALSE)
  if (length(unique(dat$recvDeployName)) > 1) {
    multi_tower_clean <- c(multi_tower_clean, bird_id)
    cat("\n---", bird_id, "---\n")
    print(sort(table(dat$recvDeployName)))
  }
}
cat("\nBirds with >1 tower after cleaning:\n")
print(multi_tower_clean)

