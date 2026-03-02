#' Dark-eyed junco Motus detections — example dataset
#'
#' Radio-telemetry detections for a single dark-eyed junco (*Junco hyemalis*)
#' tagged with a Motus transmitter (tag deploy ID 52550). Detections were
#' recorded across multiple nights and span the full range of conditions
#' encountered during nocturnal roosting analysis.
#'
#' This dataset is used throughout the roostR vignette and function examples.
#' It represents the output of the Motus data download pipeline *before* any
#' roostR preprocessing — use it as the starting point for the analysis
#' workflow.
#'
#' @format A data frame with 153,989 rows and 24 columns:
#' \describe{
#'   \item{sig}{Numeric. Signal strength in dBm.}
#'   \item{noise}{Numeric. Background noise level in dBm.}
#'   \item{motusTagID}{Integer. Physical tag model identifier.}
#'   \item{ambigID}{Numeric. Motus ambiguity code; \code{NA} if unambiguous.}
#'   \item{port}{Integer. Receiver antenna port number.}
#'   \item{runLen}{Integer. Number of consecutive detections in this run.}
#'   \item{motusFilter}{Integer. Motus quality-control flag; \code{1} = passed.}
#'   \item{mfgID}{Character. Manufacturer tag identifier.}
#'   \item{tagDeployID}{Integer. Deployment-specific tag identifier (primary
#'     key for individuals across the package).}
#'   \item{fullID}{Character. Full tag identifier string.}
#'   \item{recvDeployName}{Character. Receiver deployment name.}
#'   \item{antBearing}{Numeric. Antenna bearing in degrees.}
#'   \item{speciesEN}{Character. English species name.}
#'   \item{time}{POSIXct (UTC). Detection timestamp.}
#'   \item{year}{Integer. Calendar year of detection.}
#'   \item{doy}{Integer. Day of year (1–366).}
#'   \item{Sunrise}{POSIXct (UTC). Sunrise time for the detection date.}
#'   \item{Sunset}{POSIXct (UTC). Sunset time for the detection date.}
#'   \item{AT.Start}{POSIXct (UTC). Astronomical twilight start (sun 18–20°
#'     below horizon).}
#'   \item{AT.End}{POSIXct (UTC). Astronomical twilight end.}
#'   \item{NT.Start}{POSIXct (UTC). Nautical twilight start (sun 6–12° below
#'     horizon).}
#'   \item{NT.End}{POSIXct (UTC). Nautical twilight end.}
#'   \item{CT.Start}{POSIXct (UTC). Civil twilight start (sun 0–6° below
#'     horizon).}
#'   \item{CT.End}{POSIXct (UTC). Civil twilight end.}
#' }
#'
#' @source Motus Wildlife Tracking System (\url{https://motus.org/}).
#'   Sunrise/sunset and twilight times computed for the receiver location.
"sparrow52550"
