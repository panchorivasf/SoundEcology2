#' Calculate Soundscape Dynamicity
#'
#' Quantifies temporal variability in acoustic energy across frequency bins by comparing 
#' mean power spectra (MPS) between consecutive segments of an audio recording. Useful for 
#' identifying periods of rapid acoustic change in soundscapes.
#'
#' @param wave A `Wave` object (from `tuneR` package) containing the audio data.
#' @param fun Function to apply in spectral calculation (`"mean"`, `"var"`, etc.). 
#'   Passed to `seewave::meanspec(FUN)`.
#' @param segment_l Length (in seconds) of each analysis segment. Longer segments provide 
#'   better frequency resolution but fewer temporal comparisons (default: 5s).
#' @param freq_res Frequency resolution (in Hz) for spectral analysis. Higher values 
#'   increase frequency binning (default: 50Hz).
#' @param plot Logical. If `TRUE` (default), generates a heatmap of spectral changes.
#'
#' @return A list with components:
#' \itemize{
#'   \item `dynamicity`: Scalar metric of overall acoustic variability (mean ΔMPS across all bins).
#'   \item `delta_mps`: Matrix of absolute MPS differences between segments (frequency bins × time).
#'   \item `dynamicity_vector`: Frequency-specific variability (sum of ΔMPS per bin).
#'   \item `frequency`: Vector of frequency bins (kHz).
#'   \item `time`: Vector of segment mid-point times (s).
#' }
#'
#' @details 
#' The function works by:
#' 1. Dividing the audio into equal-length segments.
#' 2. Calculating mean power spectra (MPS) for each segment.
#' 3. Computing absolute differences in MPS between consecutive segments.
#' 4. Summarizing these differences as a "dynamicity" metric.
#'
#' @examples
#' \dontrun{
#' library(tuneR)
#' data(tico)
#' dyn <- dynamicity(tico, segment_l = 1, freq_res = 100)
#' }
#'
#' @export
#' @importFrom tuneR extractWave
#' @importFrom seewave meanspec
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn 
#'   scale_x_continuous scale_y_continuous labs theme_minimal
dynamicity <- function(wave, 
                       fun = "mean",
                       segment_l = 5,
                       freq_res = 50,
                       plot = TRUE) {
  # 1) Calculate mps for each length-second segment
  duration <- length(wave@left) / wave@samp.rate
  segment_length <- segment_l  
  n_segments <- floor(duration / segment_length)
  
  if (n_segments < 2) {
    stop("Select shorter segment_l.")
  }
  
  samples_per_segment <- round(segment_length * wave@samp.rate)
  mps_list <- list()
  freq <- NULL
  
  for (i in 1:n_segments) {
    segment <- extractWave(wave, 
                           from = (i-1) * samples_per_segment + 1, 
                           to = i * samples_per_segment, 
                           xunit = "samples")
    wl <- wave@samp.rate / freq_res  
    mps <- meanspec(segment, 
                    FUN = fun,
                    wl = wl, 
                    norm = FALSE,
                    correction = "energy",
                    plot = FALSE)
    
    mps_list[[i]] <- mps[, 2]  # Store amplitude
    if (is.null(freq)) freq <- mps[, 1]  # Store frequency bins
  }
  
  # 2) Compute mps differences between consecutive segments
  delta_mps <- matrix(nrow = length(freq), ncol = n_segments - 1)
  for (i in 1:(n_segments - 1)) {
    delta_mps[, i] <- abs(mps_list[[i + 1]] - mps_list[[i]])
  }
  
  # 3) Metrics
  dynamicity <- round(sum(delta_mps)/nrow(delta_mps), 4)  # Total variability
  dynamicity_vector <- rowSums(delta_mps)  # Variability per frequency bin
  
  # 4) Plot (if requested)
  if (plot) {
    plot_data <- expand.grid(
      Frequency = freq, 
      Time = seq(2.5, by = 2.5, length.out = n_segments - 1)
    )
    plot_data$Delta <- as.vector(delta_mps)
    
    p <- ggplot(plot_data, aes(x = Time, y = Frequency, fill = Delta)) +
      geom_tile() +
      scale_fill_gradientn(
        colors = c("black", "blue", "yellow", "red"),
        name = "??MPS"
      ) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(
        x = "Time (s)", 
        y = "Frequency (kHz)", 
        title = "Soundscape Dynamicity"
      ) +
      theme_minimal()
    print(p)
  }
  
  # Return results (invisibly)
  return(invisible(list(
    dynamicity = dynamicity,
    delta_mps = delta_mps,  
    dynamicity_vector = dynamicity_vector, 
    frequency = freq, 
    time = seq(segment_l, by = segment_l, length.out = n_segments - 1) 
  )))
}
