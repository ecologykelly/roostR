# Run this script once to regenerate data/sparrow52550.rda from the source CSVs.
# Requires: usethis, lubridate

# ── Detection data ────────────────────────────────────────────────────────────
sparrow52550 <- read.csv("data-raw/Sparrow52550.dat.csv",
                         stringsAsFactors = FALSE)

# Convert timestamp to POSIXct UTC
sparrow52550$time <- lubridate::as_datetime(sparrow52550$time, tz = "UTC")

# ── Sunrise / sunset times ────────────────────────────────────────────────────
# sunrisetimes.csv has times in local Eastern Time (America/New_York).
# We join on J.day = doy and convert to UTC.

sunrise_raw <- read.csv("data-raw/sunrisetimes.csv",
                        stringsAsFactors = FALSE,
                        fileEncoding     = "UTF-8-BOM")

# Drop empty trailing columns (artifact of extra commas in source file)
sunrise_raw <- sunrise_raw[, vapply(sunrise_raw,
                                    function(x) !all(is.na(x) | x == ""),
                                    logical(1))]

# Derive reference year from detection data to reconstruct full dates
ref_year <- as.integer(format(min(sparrow52550$time, na.rm = TRUE), "%Y"))
dates <- as.Date(sunrise_raw$J.day - 1, origin = paste0(ref_year, "-01-01"))

# Helper: parse "7:48 AM" style local time on a given date → POSIXct UTC
to_utc <- function(time_str, dates, tz_local = "America/New_York") {
  dt_local <- as.POSIXct(
    strptime(paste(dates, time_str), format = "%Y-%m-%d %I:%M %p", tz = tz_local)
  )
  lubridate::with_tz(dt_local, "UTC")
}

time_cols <- c("Sunrise", "Sunset",
               "AT.Start", "AT.End",
               "NT.Start", "NT.End",
               "CT.Start", "CT.End")

for (col in time_cols) {
  if (col %in% names(sunrise_raw)) {
    sunrise_raw[[col]] <- to_utc(sunrise_raw[[col]], dates)
  }
}

# Rename J.day → doy for the join; keep only relevant columns
join_cols  <- c("J.day", intersect(time_cols, names(sunrise_raw)))
sunrise_join <- sunrise_raw[, join_cols]
names(sunrise_join)[names(sunrise_join) == "J.day"] <- "doy"

# Join to detections
sparrow52550 <- merge(sparrow52550, sunrise_join, by = "doy", all.x = TRUE)

# ── Select the 24 documented columns ─────────────────────────────────────────
keep_cols <- c("sig", "noise", "motusTagID", "ambigID", "port", "runLen",
               "motusFilter", "mfgID", "tagDeployID", "fullID",
               "recvDeployName", "antBearing", "speciesEN",
               "time", "year", "doy",
               "Sunrise", "Sunset",
               "AT.Start", "AT.End", "NT.Start", "NT.End", "CT.Start", "CT.End")

sparrow52550 <- sparrow52550[, keep_cols]

usethis::use_data(sparrow52550, overwrite = TRUE)
