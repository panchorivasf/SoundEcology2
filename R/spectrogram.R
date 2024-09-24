#' Plot a Spectrogram with Customizable Parameters
#'
#' This function generates and plots a spectrogram from a Wave object, allowing customization of the frequency resolution, amplitude threshold, and plot aesthetics.  The matrix is obtained with seewave's spectro() function without normalization and the amplitude values are transformed to dBFS. The spectrogram is plotted using `ggplot2` with optional dark mode.
#'
#' @param wave A Wave object (from the tuneR package) representing the audio data to be analyzed.
#' @param freq.res Numeric. The frequency resolution, specified in Hz. This determines the window length for the spectrogram. Default is 100 Hz.
#' @param cutoff Numeric. The minimum amplitude in dB to be displayed in the spectrogram. Values below this threshold are set to NA. Default is -60 dB.
#' @param dark.plot Logical. Whether to apply a dark theme to the spectrogram plot. If `TRUE`, the background and text are dark. Default is `FALSE`.
#' @param plot.title Character. Optional title for the plot. If not provided, no title is displayed. Default is `NULL`.
#'
#' @return A `ggplot2` plot of the spectrogram. The function prints the plot to the graphics device and returns `NULL`.
#' @export
#'
#' @import seewave
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' # Load a Wave file
#' wave <- tuneR::readWave(system.file("extdata", "test.wav", package = "tuneR"))
#'
#' # Generate a spectrogram with default parameters
#' spectrogram(wave)
#'
#' # Generate a spectrogram with a custom frequency resolution and dark plot mode
#' spectrogram(wave, freq.res = 50, cutoff = -50, dark.plot = TRUE, plot.title = "Spectrogram")
#' }
spectrogram <- function(wave, freq.res = 100, cutoff = -60, dark.plot = FALSE, plot.title = NULL){

  total_duration <- seewave::duration(wave)

  # Calculate window length according to the sampling rate for a fixed frequency resolution
  wl = wave@samp.rate / freq.res

  # Force it to be even
  if (wl %% 2 == 1) { wl <- wl + 1 }

  cat("Automatic window length:", wl, "\n")

    # Generate the spectrogram data using soundecology2::spectro2
    spectrogram_matrix <- seewave::spectro(wave,
                                           wl = wl,
                                           norm = FALSE,
                                           plot = FALSE,
                                           dB = NULL,
                                           scale = FALSE,
                                           correction = "amplitude")$amp



    # Calculate maximum possible amplitude based on bit depth
    amp_max <- switch(as.character(wave@bit),
                      `16` = 32767,
                      `24` = 8388607,
                      `32` = 2147483647,
                      stop("Unsupported bit depth"))

    # Convert amplitude to dBFS
    spectrogram_matrix <-  20 * log10(abs(spectrogram_matrix) / amp_max)



  range <- range(spectrogram_matrix)

  range.min <- range[1]
  range.max <- range[2]

  cat("Energy range: ", range.min, " : ", range.max, "\n")
  cat("Nr of frequency bins: ", nrow(spectrogram_matrix), "\n")
  cat("Nr of time frames: ", ncol(spectrogram_matrix), "\n")


  # Create the time axis values
  time_values <- seq(0, total_duration, length.out = ncol(spectrogram_matrix))
  # Create the frequency axis values
  freq_values <- seq(0, by = wave@samp.rate / 2 / nrow(spectrogram_matrix), length.out = nrow(spectrogram_matrix)) / 1000  # Convert to kHz

  # clip the values below the threshold to a new unique value below the threshold
  spectrogram_matrix[spectrogram_matrix < cutoff] <- NA



  # Define frequency ranges
  samp_rate <- wave@samp.rate
  freq_bins <- seq(0, samp_rate / 2, length.out = nrow(spectrogram_matrix))

  time <- seq(0, seewave::duration(wave), length.out = ncol(spectrogram_matrix))
  freq <- freq_bins/1000

  spectro_df <- melt(spectrogram_matrix)
  spectro_df$Frequency <- freq[spectro_df$Var1]
  spectro_df$Time <- time[spectro_df$Var2]


  p <- ggplot(spectro_df, aes(Time, Frequency, fill = value)) +
    geom_tile() +
    scale_fill_continuous(type = "viridis", na.value = "transparent") +
    labs(x = "Time (s)", y = "Frequency (kHz)") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(title=plot.title) +
    theme(legend.position = "none",
          panel.grid.major = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.3),
          panel.grid.minor = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.05))


  # Apply the dark theme logic conditionally
  if (dark.plot) {
    p <- p + theme(
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      axis.line = element_line(color = "white"),
      axis.ticks = element_line(color = "white"),
      panel.grid.major = element_line(color = scales::alpha("white", 0.2), linetype = "solid", linewidth = 0.3),
      panel.grid.minor = element_line(color = scales::alpha("white", 0.2), linetype = "solid", linewidth = 0.05),
      plot.title = element_text(color = "white", face = "bold", size = 14)
    )
  }

  # Display the plot
  print(p)

  return(p)
}
