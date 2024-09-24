#' Convert Matrix to Spectrogram Plot
#'
#' This function takes a spectrogram matrix (typically in decibels) and generates a spectrogram plot, with options for customizing the appearance, such as a dark background or minimum amplitude threshold. The spectrogram matrix should represent time-frequency data, where rows correspond to frequency bins and columns to time frames.
#'
#' @param spectrogram_matrix A numeric matrix where rows represent frequency bins and columns represent time frames, typically containing decibel values.
#' @param duration Numeric. The total duration of the audio signal (in seconds).
#' @param samp.rate Numeric. The sampling rate of the audio signal (in Hz).
#' @param min.amp Numeric. The minimum amplitude (in dB) to display. Any values below this threshold will be set to the minimum. Default is -60 dB.
#' @param dark.plot Logical. If TRUE, generates a spectrogram with a dark background. Default is FALSE.
#' @param plot.title Character. The title for the spectrogram plot. Default is NULL.
#'
#' @return A ggplot object representing the spectrogram plot.
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' matrix <- matrix(runif(10000, -80, 0), nrow = 100, ncol = 100)
#' matrix_to_spectrogram(spectrogram_matrix = matrix, duration = 10, samp.rate = 44100)
#' }
matrix_to_spectrogram <- function(spectrogram_matrix,
                                  duration,
                                  samp.rate,
                                  min.amp = -60,
                                  dark.plot = FALSE,
                                  plot.title = NULL){

  total_duration <- duration

  samp_rate <- samp.rate

  range <- range(spectrogram_matrix)

  range.min <- range[1]
  range.max <- range[2]

  cat("Energy range: ", range.min, " : ", range.max, "\n")
  cat("Nr of frequency bins: ", nrow(spectrogram_matrix), "\n")
  cat("Nr of time frames: ", ncol(spectrogram_matrix), "\n")


  # Create the time axis values
  time_values <- seq(0, total_duration, length.out = ncol(spectrogram_matrix))
  # Create the frequency axis values
  freq_values <- seq(0, by = samp_rate / 2 / nrow(spectrogram_matrix), length.out = nrow(spectrogram_matrix)) / 1000  # Convert to kHz


  spectrogram_matrix[spectrogram_matrix < min.amp] <- min.amp


  # Conditionally set na.value based on dark.plot flag
  na_color <- if (dark.plot) "black" else "white"

  # Create the color gradient function based on dB values
  color_func <- scales::col_numeric(palette = c(na_color, "#2c7bb6","#00a6ca","#00ccbc","#90eb9d",
                                                "#ffff8c", "#f9d057", "#f29e2e","#e76818","#d7191c"),
                                    domain = c(min.amp, 0))

  # Create a dataframe for plotting
  plot_data <- as.data.frame(spectrogram_matrix)
  colnames(plot_data) <- time_values
  plot_data <- cbind(Frequency = freq_values, plot_data)
  plot_data <- reshape2::melt(plot_data, id.vars = "Frequency", variable.name = "Time", value.name = "dB")
  plot_data$Time <- as.numeric(as.character(plot_data$Time))

  # Apply the color function to the dB values
  plot_data$Color <- color_func(plot_data$dB)

  # Conditionally set na.value based on dark.plot flag
  na_color <- if (dark.plot) "black" else "white"

  p <- ggplot(plot_data, aes(x = Time, y = Frequency)) +
    geom_tile(aes(fill = Color), color = NA) +
    scale_fill_identity(na.value = na_color) +
    # scale_fill_identity(na.value = "transparent") +
    labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0))

  # Apply the dark theme logic conditionally
  if (dark.plot) {
    p <- p + theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
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
