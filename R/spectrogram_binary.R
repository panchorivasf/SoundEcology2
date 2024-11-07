#' Generate a Binary Spectrogram
#'
#' This function generates a binary spectrogram from an audio `wave` object. The spectrogram is calculated using the Fast Fourier Transform (FFT) and then converted into a binary matrix based on a specified dB cutoff. Optionally, the function can plot the binary spectrogram.
#'
#' @param wave A `Wave` object from the `tuneR` package representing the audio signal.
#' @param channel A character string indicating which channel to use for stereo audio. Options are `'left'`, `'right'`, or `'mix'` (for combining both channels). If `NULL`, the left channel is used by default.
#' @param freq.res A numeric value specifying the frequency resolution (in Hz) of the spectrogram. Default is 100 Hz.
#' @param cutoff A numeric value specifying the dB cutoff threshold for the binary transformation. Amplitudes below this value will be set to 0. Default is -50 dB.
#' @param plot A logical value indicating whether to plot the binary spectrogram. Default is `TRUE`.
#' @param plot.title An optional character string to use as the title of the plot. Default is `NULL`.
#' @param ggplot A logical value indicating whether to use ggplot2 to plot the spectrogram. If FALSE (default), R's base plot function is used instead, rendering much faster.
#' @param verbose A logical value indicating whether to print additional information (such as the amplitude range and cutoff) to the console. Default is `FALSE`.
#'
#' @return Returns a binary matrix representing the spectrogram, where `1` represents amplitude above the cutoff and `0` represents amplitude below the cutoff.
#' @export
#'
#' @import ggplot2
#' @import reshape2
#' @import tuneR
#' @import seewave
#'
#' @examples
#' \dontrun{
#'   # Load a wave file
#'   wave <- tuneR::readWave("path_to_file.wav")
#'
#'   # Generate and plot a binary spectrogram
#'   binary_spec <- spectrogram_binary(wave, cutoff = -40, verbose = TRUE)
#' }
spectrogram_binary <- function(wave,
                               channel = 'left',
                               freq.res = 100,
                               cutoff = -50,
                               plot = FALSE,
                               plot.title = NULL,
                               verbose = FALSE,
                               ggplot = FALSE) {  # New argument

  require(reshape2)

  # Check if the wave is stereo
  if (wave@stereo) {
    if (channel == 'left') {
      wave <- tuneR::channel(wave, 'left')
    } else if (channel == 'right') {
      wave <- tuneR::channel(wave, 'right')
    } else if (channel == 'mix'){
      wave <- tuneR::mono(wave, "both")
    }
  }

  # Calculate window length according to the sampling rate for a fixed frequency resolution 
  wl = wave@samp.rate / freq.res

  # Force it to be even
  if (wl %% 2 == 1) { wl <- wl + 1 }

  # Remove DC offset
  wave <- rmoffset(wave, output = "Wave")

  # Generate the spectrogram data using soundecology2::spectro2
  spectro_data <- seewave::spectro(wave,
                                   wl = wl,
                                   norm = FALSE,
                                   plot = FALSE,
                                   dB = NULL,
                                   scale = FALSE,
                                   correction = "amplitude")$amp

  raw_range <- range(spectro_data)

  if (verbose) cat("Raw FFT amplitude range:", raw_range, "\n")

  # Calculate maximum possible amplitude based on bit depth
  amp_max <- switch(as.character(wave@bit),
                    `16` = 32767,
                    `24` = 8388607,
                    `32` = 2147483647,
                    stop("Unsupported bit depth"))

  # Convert amplitude to dBFS
  spectro_data <-  20 * log10(abs(spectro_data) / amp_max)

  dBFS_range <- range(spectro_data)

  if (verbose) {
    cat("dBFS dynamic range:", range(spectro_data), "\n")
    cat("Cutoff:", cutoff, "\n")
  }

  # Transform the spectrogram matrix to binary
  spectrogram <- spectro_data > cutoff

  if (plot) {
    # Define frequency ranges
    samp_rate <- wave@samp.rate
    freq_bins <- seq(0, samp_rate / 2, length.out = nrow(spectro_data))
    time <- seq(0, seewave::duration(wave), length.out = ncol(spectrogram))
    freq <- freq_bins / 1000

    binary_spectro_df <- melt(spectrogram)
    binary_spectro_df$Frequency <- freq[binary_spectro_df$Var1]
    binary_spectro_df$Time <- time[binary_spectro_df$Var2]

    if (ggplot) {
      # Use ggplot for plotting
      p <- ggplot(binary_spectro_df, aes(Time, Frequency, fill = value)) +
        geom_tile() +
        scale_fill_manual(values = c("transparent", "black")) +
        labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw() +
        labs(title = plot.title) +
        theme(legend.position = "none",
              panel.grid.major = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.3),
              panel.grid.minor = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.05))

      print(p)
    } else {
      # Use base R plotting
      image(time, freq, t(spectrogram),
            col = c("white", "black"),
            xlab = "Time (s)", ylab = "Frequency (kHz)",
            main = plot.title,
            useRaster = TRUE, axes = FALSE)

      # Add axes
      axis(1, at = pretty(time), labels = pretty(time))
      axis(2, at = pretty(freq), labels = pretty(freq))

      # Add grid lines
      abline(h = pretty(freq), col = "lightgray", lty = "dashed")
      abline(v = pretty(time), col = "lightgray", lty = "dashed")

      # Add box around the plot
      box()

    }
  }

  invisible(list(spectrogram = spectrogram, raw_range = raw_range, dbfs_range = dBFS_range))
}



