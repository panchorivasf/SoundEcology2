#' Spectral Dynamics Analysis
#'
#' Computes summary statistics based on spectral similarity (consecutive pairwise comparison of PSD vectors) over user-defined time steps for specified frequency bands.
#'
#' @param wave A Wave object representing the audio file to be analyzed.
#' @param timestep Numeric value indicating the time step (in seconds) for segmenting the audio file to calculate the Power Spectral Densitiess.
#' @param rmoffset Logical value indicating whether to remove DC offset from the audio (default is TRUE).
#' @param highpass Numeric value indicating the cutoff frequency for the high-pass filter (in Hz, default is 0 Hz).
#' @param lfc A numeric vector of length 2 indicating the lower frequency range (e.g., c(0, 2) for 0-2 Hz).
#' @param mfc A numeric vector of length 2 indicating the mid-frequency range (e.g., c(2, 10) for 2-10 Hz).
#' @param hfc A numeric vector of length 2 indicating the high-frequency range (e.g., c(10, 16) for 10-16 Hz).
#' @param ufc A numeric vector of length 2 indicating the ultra-frequency range (e.g., c(16, 24) for 16-24 Hz).
#'
#' @return A tibble with summary statistics (standard deviation, variance, mean, range) derived from spectral similarity for each frequency band.
#' The frequency bands are ordered from highest (ufc) to lowest (lfc).
#'
#' @examples
#' # Assuming `wave` is a Wave object loaded in R
#' result <- spectral_dynamics(wave, timestep = 5, lfc = c(0, 2), mfc = c(2, 10), hfc = c(10, 16), ufc = c(16, 24))
#' print(result)
#'
#' @import dplyr
#' @import tidyr
#' @import seewave
#'
#' @export
spectral_dynamics <- function(wave,
                              timestep = 5,
                              rmoffset = TRUE,
                              highpass = 0,
                              lfc = c(0, 2),
                              mfc = c(2, 10),
                              hfc = c(10, 16),
                              ufc = c(16, 24)) {

  # Remove DC offset if specified
  if (rmoffset) {
    wave <- rmoffset(wave, output = "Wave")
  }

  # Apply a high-pass filter if specified
  if (highpass > 0) {
    wave <- ffilter(wave, from = highpass, output = "Wave")
  }

  # Function to extract values according to frequency limits
  extr_range <- function(meanspec_matrix, x_min, x_max) {
    indices <- which(meanspec_matrix[, "x"] >= x_min & meanspec_matrix[, "x"] <= x_max)
    y_values <- meanspec_matrix[indices, "y"]
    return(y_values)
  }

  # Initialize a list to store similarity stats
  similarity_stats <- list()

  # Get the total duration of the audio file
  total_duration <- duration(wave)

  # Ensure total_duration is valid
  if (is.null(total_duration) || total_duration <= 0) {
    stop("The audio file has an invalid duration.")
  }

  # Calculate meanspectra over time windows and compute spectral similarity
  for (i in seq(0, total_duration - timestep, by = timestep)) {
    end_time <- min(i + timestep, total_duration)

    if (end_time <= total_duration) {
      meanspec_current <- meanspec(wave,
                                   from = i, to = end_time,
                                   norm = FALSE,
                                   plot = FALSE,
                                   PSD = TRUE,
                                   wl = 512)

      if (i == 0) {
        prev_meanspec <- meanspec_current
      } else {
        # Calculate similarity for each frequency range
        freq_ranges <- list(lfc = lfc, mfc = mfc, hfc = hfc, ufc = ufc)
        for (range_name in names(freq_ranges)) {
          range <- freq_ranges[[range_name]]
          prev_values <- extr_range(prev_meanspec, range[1], range[2])
          curr_values <- extr_range(meanspec_current, range[1], range[2])
          sim_value <- simspec(prev_values, curr_values)
          similarity_stats[[range_name]] <- c(similarity_stats[[range_name]], sim_value)
        }
        prev_meanspec <- meanspec_current
      }
    }
  }

  # Convert the summary statistics to a tibble
  summary_stats <- lapply(similarity_stats, function(sim_values) {
    tibble(
      sd = round(sd(sim_values),2),
      var = round(var(sim_values),2),
      mean = round(mean(sim_values),2),
      min = round(min(sim_values),2),
      max = round(max(sim_values),2)
    )
  })

  summary_stats_tibble <- bind_rows(summary_stats, .id = "frequency_band")

  # Set the custom order for frequency bands
  summary_stats_tibble <- summary_stats_tibble %>%
    mutate(frequency_band = factor(frequency_band, levels = c("ufc", "hfc", "mfc", "lfc"))) %>%
    arrange(frequency_band)

  return(summary_stats_tibble)
}

