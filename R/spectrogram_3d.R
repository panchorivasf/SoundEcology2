#' Generate a 3D Spectrogram
#'
#' This function generates a 3D surface plot of a spectrogram from an audio wave object.
#' It uses a fixed frequency resolution and allows setting a minimum amplitude threshold
#' for the plot. Optionally, a cutoff plane can be added to the plot for better visualization
#' of the amplitude floor.
#'
#' @param wave A wave object (imported with the `tuneR` package) containing the audio data to be analyzed.
#' @param freq.res Numeric value indicating the frequency resolution in Hz. Default is 100 Hz.
#' @param min.amp Numeric value representing the minimum amplitude (in dB) to display. Amplitudes lower than this value will be clipped. Default is -60 dB.
#' @param add.cutoff Optional numeric value specifying an amplitude cutoff plane to be added to the plot. If NULL, no cutoff plane is added.
#' @param plot.title Optional string to set a custom title for the plot. If NULL, no title is added.
#'
#' @return A `plotly` object representing the 3D surface plot of the spectrogram.
#' @export
#' 
#' @import plotly
#' @import reshape2
#' @import seewave
#' @import magrittr
#'
#' @examples
#' # Assuming 'wave' is a wave object:
#' spectrogram_3d(wave, freq.res = 100, min.amp = -60, add.cutoff = -40, plot.title = "3D Spectrogram")

spectrogram_3d <- function(wave, freq.res = 100, min.amp = -60, add.cutoff = NULL, plot.title = NULL){
  
  total_duration <- seewave::duration(wave)
  
  # Calculate window length according to the sampling rate for a fixed frequency resolution
  wl = wave@samp.rate / freq.res
  
  # Force it to be even
  if (wl %% 2 == 1) { wl <- wl + 1 }
  
  # Generate the spectrogram data using seewave::spectro
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
  
  # Apply min.amp cutoff to the spectrogram
  spectrogram_matrix[spectrogram_matrix < min.amp] <- min.amp
  
  # Create the time axis values
  time_values <- seq(0, total_duration, length.out = ncol(spectrogram_matrix))
  
  # Create the frequency axis values
  freq_values <- seq(0, by = wave@samp.rate / 2 / nrow(spectrogram_matrix), length.out = nrow(spectrogram_matrix)) / 1000  # Convert to kHz
  
  # Create the color gradient function based on dB values
  color_scale <- list(
    list(0, "white"),
    list(0.1, "#2c7bb6"),
    list(0.2, "#00a6ca"),
    list(0.3, "#00ccbc"),
    list(0.4, "#90eb9d"),
    list(0.5, "#ffff8c"),
    list(0.6, "#f9d057"),
    list(0.7, "#f29e2e"),
    list(0.8, "#e76818"),
    list(1, "#d7191c")
  )
  
  
  # Create a 3D surface plot using plotly
  p <- plot_ly(
    x = ~time_values,
    y = ~freq_values,
    z = ~spectrogram_matrix,
    type = "surface",
    colorscale = color_scale,#list(c(0,0.2, 0.5, 1), c("black", "blue", "yellow", "red")),
    cmin = min.amp,
    cmax = 0,
    colorbar = list(title = "Amplitude (dB)")
  ) %>%
    layout(
      title = plot.title,
      scene = list(
        xaxis = list(title = "Time (s)"),
        yaxis = list(title = "Frequency (kHz)"),
        zaxis = list(title = "Amplitude (dB)")
      )
    )
  
  if(!is.null(add.cutoff)){
   # Add a cutoff plane at min.amp
  p <- p %>% add_trace(
    x = ~c(time_values[1], time_values[length(time_values)], time_values[length(time_values)], time_values[1]),
    y = ~c(freq_values[1], freq_values[1], freq_values[length(freq_values)], freq_values[length(freq_values)]),
    z = ~rep(add.cutoff, 4),
    type = "mesh3d",
    color = I("red"),
    opacity = 0.8,
    name = "Cutoff Plane"
  )
  }
  
  return(p)
}

