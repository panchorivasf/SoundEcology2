#' Calculate Energy Distribution Across Frequency Bands
#'
#' Calculates the proportion of acoustic energy in each frequency band of an
#' audio recording. This function generates a spectrogram and computes the
#' distribution of energy above a specified decibel threshold across
#' user-defined frequency bands.
#'
#' @param wave A Wave object from the tuneR package containing the audio data
#' @param freq.res Numeric. Frequency resolution in Hz (frequency bin width).
#'   Default is 50 Hz
#' @param win.fun Character. Window function for spectrogram generation.
#'   Options include "hanning", "hamming", "blackman", etc. Default is "hanning"
#' @param min.freq Numeric. Minimum frequency in Hz to include in analysis.
#'   Default is 0 Hz
#' @param max.freq Numeric. Maximum frequency in Hz to include in analysis.
#'   Will be adjusted to Nyquist frequency if it exceeds it. Default is 10000 Hz
#' @param n.bands Integer. Number of frequency bands to divide the spectrum into.
#'   Default is 10
#' @param cutoff Numeric. Decibel threshold for including values in calculations.
#'   Values below this threshold are excluded. Default is -60 dB
#' @param norm.spec Logical. Should the spectrogram be normalized? If TRUE,
#'   uses normalized amplitude values. If FALSE, converts to decibels.
#'   Default is FALSE
#' @param noise.red Integer. Noise reduction method: 0 = none, 1 = by rows,
#'   2 = by columns. Default is 0
#' @param rm.offset Logical. Should DC offset be removed from the signal?
#'   Default is TRUE
#' @param db.fs Logical. If TRUE, converts amplitude to dBFS (decibels relative
#'   to full scale) based on bit depth. If FALSE, uses standard power-based
#'   decibel conversion. Only applies when norm.spec = FALSE. Default is TRUE
#'
#' @return A data frame with one row containing:
#'   \itemize{
#'     \item Metadata columns: w_len, w_fun, cutoff, min_f, max_f, n_bands,
#'       norm, noise_red, rm_off, samp_rate, freq_res, nyq, dur, channels
#'     \item Energy proportion columns: Named by channel and frequency range
#'       (e.g., "left_0.000-1.000", "right_0.000-1.000" for stereo files,
#'       or "mono_0.000-1.000" for mono files)
#'   }
#'   Energy proportions are normalized to sum to 1.0 across all bands.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Generating a spectrogram using the specified parameters
#'   \item Dividing the frequency range into n.bands equal-width bands
#'   \item Counting spectrogram cells above the cutoff threshold in each band
#'   \item Calculating the proportion of energy in each band relative to total
#'   \item Normalizing proportions to sum to 1.0
#' }
#'
#' For stereo files, both channels are processed independently and results are
#' returned for each channel.
#'
#' The function uses dBFS (decibels relative to full scale) by default when
#' db.fs = TRUE, which accounts for the bit depth of the recording (16, 24, or
#' 32 bit). This provides a standardized measure across recordings with different
#' bit depths.
#'
#' @examples
#' \dontrun{
#' library(tuneR)
#'
#' # Load an audio file
#' wave <- readWave("recording.wav")
#'
#' # Calculate band energy with default parameters
#' result <- band_energy(wave)
#'
#' # Calculate with custom frequency range and more bands
#' result <- band_energy(wave,
#'                      min.freq = 1000,
#'                      max.freq = 8000,
#'                      n.bands = 20,
#'                      cutoff = -50)
#'
#' # Use normalized spectrogram
#' result <- band_energy(wave, norm.spec = TRUE)
#' }
#'
#' @seealso \code{\link{band_energy_folder}} for batch processing multiple files
#'
#' @export
band_energy <- function(wave,
                        freq.res = 50,
                        win.fun = "hanning",
                        min.freq = 0,
                        max.freq = 10000,
                        n.bands = 10,
                        cutoff = -60,
                        norm.spec = FALSE,
                        noise.red = 0,
                        rm.offset = TRUE,
                        db.fs = TRUE) {
  
  # Validate inputs
  stopifnot(
    "max.freq must be numeric" = is.numeric(max.freq),
    "cutoff must be numeric" = is.numeric(cutoff),
    "n.bands must be positive" = n.bands > 0
  )
  
  # Calculate frequency parameters
  samplingrate <- wave@samp.rate
  nyquist_freq <- samplingrate / 2
  duration <- length(wave@left) / samplingrate
  freq_step <- (max.freq - min.freq) / n.bands
  
  # Adjust max.freq if it exceeds Nyquist
  if (max.freq > nyquist_freq) {
    message(sprintf("max.freq adjusted from %d Hz to Nyquist frequency %d Hz", 
                    max.freq, nyquist_freq))
    max.freq <- nyquist_freq
  }
  
  # Calculate window length
  wlen <- samplingrate / freq.res
  if (wlen %% 2 == 1) wlen <- wlen + 1
  
  # Noise reduction label
  noise_label <- switch(as.character(noise.red),
                        "1" = "rows",
                        "2" = "columns",
                        "none")
  
  # Function to process a single channel
  process_channel <- function(channel_wave, channel_name) {
    
    # Remove DC offset if requested
    if (rm.offset) {
      channel_wave <- seewave::rmoffset(channel_wave, output = "Wave")
    }
    
    # Generate spectrogram
    spec_params <- list(
      wave = channel_wave,
      wl = wlen,
      wn = win.fun,
      plot = FALSE
    )
    
    if (norm.spec) {
      if (noise.red %in% c(1, 2)) {
        spec_params$noisereduction <- noise.red
      }
      specA <- do.call(seewave::spectro, spec_params)$amp
    } else {
      spec_params$norm <- FALSE
      spec_params$dB <- NULL
      spec_params$correction <- 'amplitude'
      if (noise.red %in% c(1, 2)) {
        spec_params$noisereduction <- noise.red
      }
      
      specA <- do.call(seewave::spectro, spec_params)$amp
      
      # Convert to dB
      if (db.fs) {
        amp_max <- switch(as.character(channel_wave@bit),
                          "16" = 32768,
                          "24" = 8388607,
                          "32" = 2147483647,
                          stop("Unsupported bit depth: ", channel_wave@bit))
        specA <- 20 * log10(abs(specA) / amp_max)
      } else {
        specA <- 10 * log10(specA^2)
      }
    }
    
    # Calculate band indices
    freq_per_row <- freq.res
    minspec <- round(min.freq / freq_per_row)
    maxspec <- round(max.freq / freq_per_row)
    
    # Ensure indices are within bounds
    minspec <- max(1, minspec)
    maxspec <- min(nrow(specA), maxspec)
    
    # Extract relevant frequency range
    speclims <- specA[minspec:maxspec, , drop = FALSE]
    total_values <- sum(speclims > cutoff)
    
    # Preallocate results
    band_props <- numeric(n.bands)
    band_ranges <- character(n.bands)
    
    # Calculate proportions for each band (vectorized where possible)
    band_height <- (maxspec - minspec) / n.bands
    
    for (j in seq_len(n.bands)) {
      miny <- minspec + (j - 1) * band_height
      maxy <- minspec + j * band_height - 1
      
      # Ensure indices are integers and within bounds
      miny <- max(1, round(miny))
      maxy <- min(maxspec, round(maxy))
      
      # Extract band subset
      subA <- specA[miny:maxy, , drop = FALSE]
      
      # Count values above threshold
      band_count <- sum(subA > cutoff)
      
      # Calculate proportion
      band_props[j] <- if (total_values == 0) 0 else band_count / total_values
      
      # Create frequency range label
      freq_low <- min.freq + (j - 1) * freq_step
      freq_high <- min.freq + j * freq_step
      band_ranges[j] <- sprintf("%.3f-%.3f", freq_low / 1000, freq_high / 1000)
    }
    
    # Normalize proportions to sum to 1
    if (sum(band_props) > 0) {
      band_props <- band_props / sum(band_props)
    }
    
    # Return as named vector
    names(band_props) <- paste0(channel_name, "_", band_ranges)
    band_props
  }
  
  # Process channels
  if (wave@stereo) {
    left_channel <- tuneR::channel(wave, which = "left")
    right_channel <- tuneR::channel(wave, which = "right")
    
    left_props <- process_channel(left_channel, "left")
    right_props <- process_channel(right_channel, "right")
    
    result_data <- c(left_props, right_props)
    channels_label <- "stereo"
  } else {
    left_channel <- tuneR::channel(wave, which = "left")
    result_data <- process_channel(left_channel, "mono")
    channels_label <- "mono"
  }
  
  # Create metadata
  metadata <- data.frame(
    w_len = wlen,
    w_fun = win.fun,
    cutoff = cutoff,
    min_f = min.freq,
    max_f = max.freq,
    n_bands = n.bands,
    norm = norm.spec,
    noise_red = noise_label,
    rm_off = rm.offset,
    samp_rate = samplingrate,
    freq_res = freq.res,
    nyq = nyquist_freq,
    dur = duration,
    channels = channels_label,
    stringsAsFactors = FALSE
  )
  
  # Combine metadata with band energy values
  result <- cbind(metadata, as.data.frame(t(result_data)))
  
  return(result)
}