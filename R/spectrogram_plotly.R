#' Plot Spectrogram Using Plotly
#'
#' This function generates a spectrogram plot using Plotly, allowing for customization of frequency resolution, amplitude normalization, and plot appearance. 
#' It can display the spectrogram in decibels (dBFS) or linear scale, apply a minimum amplitude threshold, and customize the plot with dark themes and grid settings.
#'
#' @param wave A Wave object representing the audio signal to be analyzed.
#' @param norm Logical. If `TRUE`, the amplitude values are normalized. Defaults to `FALSE`.
#' @param freq.res Numeric. The frequency resolution in Hz, which determines the window length for the spectrogram. Defaults to 100.
#' @param overlap Numeric. The percentage overlap between successive windows in the spectrogram. Defaults to 50.
#' @param min.amp Numeric. The minimum amplitude threshold in decibels (dB) for the spectrogram. Values below this will be clipped. Defaults to -60 dB.
#' @param dark.plot Logical. If `TRUE`, applies a dark theme to the plot. Defaults to `FALSE`.
#' @param grid Logical. If `TRUE`, grid lines are shown in the plot. Defaults to `TRUE`.
#' @param plot.title Character. Optional. A title for the plot. If `NULL`, no title is displayed. Defaults to `NULL`.
#'
#' @return A Plotly object containing the spectrogram heatmap.
#' @export
#'
#' @examples
#' # Load a Wave object
#' data(orni)
#' # Generate a Plotly spectrogram
#' p <- spectrogram_plotly(orni, norm = TRUE, freq.res = 200, min.amp = -50, dark.plot = TRUE, plot.title = "Spectrogram")
#' p

spectrogram_plotly <- function(wave, norm = FALSE, freq.res = 100, overlap = 50, zp = 0, min.amp = -60, dark.plot = FALSE, grid = TRUE, plot.title = NULL){
  
  total_duration <- seewave::duration(wave)
  
  # Calculate window length according to the sampling rate for a fixed frequency resolution
  wl = wave@samp.rate / freq.res
  
  # Force it to be even
  if (wl %% 2 == 1) { wl <- wl + 1 }
  
  if(norm){
    spectrogram_matrix <- seewave::spectro(wave,
                                           wl = wl,
                                           ovlp = overlap,
                                           zp = zp,
                                           norm = TRUE,
                                           plot = FALSE,
                                           dB = "max0",
                                           scale = FALSE)$amp
    
    
  } else {
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
    
    
  }

  # Apply min.amp cutoff
  spectrogram_matrix[spectrogram_matrix < min.amp] <- min.amp
  
  # Create the time axis values
  time_values <- seq(0, total_duration, length.out = ncol(spectrogram_matrix))
  
  # Create the frequency axis values (converted to kHz)
  freq_values <- seq(0, by = wave@samp.rate / 2 / nrow(spectrogram_matrix), length.out = nrow(spectrogram_matrix)) / 1000  # Convert to kHz
  

  # Create a 2D heatmap plot using plotly, with hoverinfo excluding 'Trace 0'
  p <- plot_ly(
    x = ~time_values,
    y = ~freq_values,
    z = ~spectrogram_matrix,
    type = "heatmap",
    colorscale = list(
      c(0, "transparent"),
      c(0.1, "#2c7bb6"),
      c(0.2, "#00a6ca"),
      c(0.4, "#00ccbc"),
      c(0.6, "#90eb9d"),
      c(0.8, "#f9d057"),
      c(1, "#d7191c")
    ),
    zmin = min.amp,
    zmax = 0,
    hoverinfo = "x+y+z", 
    hovertemplate = paste(
      "<b>Energy:</b> %{z:.1f} dB<br><b>Time:</b> %{x:.1f} s<br><b>Frequency:</b> %{y:.1f} kHz<br><extra></extra>"
      
    ),
    showlegend = FALSE,
    colorbar = list(
      title = list(text = "Amplitude (dB)", font = list(color = ifelse(dark.plot, "white", "black"))),
      tickfont = list(color = ifelse(dark.plot, "white", "black"))
    )
  ) %>%
    layout(
      title = plot.title,
      xaxis = list(title = "Time (s)", zeroline = TRUE),
      yaxis = list(title = "Frequency (kHz)", zeroline = TRUE)
    )
  
  # Apply dark plot theme if required
  if (dark.plot) {
    p <- p %>%
      layout(
        paper_bgcolor = "black",
        plot_bgcolor = "black",
        xaxis = list(title = "Time (s)", zeroline = TRUE, color = "white", zerolinecolor = "white"),
        yaxis = list(title = "Frequency (kHz)", zeroline = TRUE, color = "white", zerolinecolor = "white"),
        title = list(font = list(color = "white"))
      )
  }
  
  # Apply grid settings based on grid argument and dark.plot flag
  if (grid) {
    grid_color <- ifelse(dark.plot, "white", "black")
    
    p <- p %>%
      layout(
        xaxis = list(showgrid = TRUE, gridcolor = grid_color),
        yaxis = list(showgrid = TRUE, gridcolor = grid_color)
      )
  } else {
    p <- p %>%
      layout(
        xaxis = list(showgrid = FALSE),
        yaxis = list(showgrid = FALSE)
      )
  }
  
  return(p)
}
