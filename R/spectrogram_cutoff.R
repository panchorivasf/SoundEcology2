#' Generate a Spectrogram with Energy Cutoff
#'
#' This function generates a spectrogram from a wave object and applies a cutoff to filter out
#' amplitude values below a specified threshold. The resulting spectrogram can be plotted using either
#' ggplot2 or base R plotting functions.
#'
#' @param wave A `Wave` object containing the audio data.
#' @param cutoff A numeric value specifying the dB cutoff threshold for the spectrogram. Values below this threshold will be clipped to `NA`.
#' @param plot A logical value indicating whether to plot the resulting spectrogram. Defaults to `FALSE`.
#' @param plot.title An optional string specifying the title of the plot. Defaults to `NULL`.
#' @param ggplot A logical value indicating whether to use ggplot2 for plotting. If `FALSE`, base R plotting will be used to generate the plot. Defaults to `TRUE`.
#'
#' @return Returns the spectrogram matrix with the cutoff applied. If `plot = TRUE`, a plot of the spectrogram is displayed.
#' @export
#'
#' @importFrom seewave duration
#' @importFrom seewave spectro
#' @import viridis
#'
#' @examples
#' \dontrun{
#' wave <- seewave::synth(d=1, f=8000, cf=1000)
#' spectrogram_cutoff(wave, cutoff = -40, plot = TRUE, plot.title = "Filtered Spectrogram")
#' }
spectrogram_cutoff <- function(wave,
                               cutoff = -60,
                               freq.res = 100,
                               plot = FALSE,
                               plot.title = NULL,
                               ggplot = TRUE,
                               noise.red = NULL){

  total_duration <- seewave::duration(wave)
  samp_rate <- wave@samp.rate
  wl <- samp_rate / freq.res
  if (wl %% 2 == 1) { wl <- wl + 1 }
  cat("Automatic window length:", wl, "\n")

  # Remove DC offset
  wave <- rmoffset(wave, output = "Wave")
  
  if (noise.red == "rows"){
    spectro_res <- seewave::spectro(wave,
                                    wl = wl,
                                    norm = FALSE,
                                    dB = NULL,
                                    plot = FALSE,
                                    scale = FALSE,
                                    correction = "amplitude",
                                    noisereduction = 1)
  } else if (noise.red == "cols"){
    spectro_res <- seewave::spectro(wave,
                                    wl = wl,
                                    norm = FALSE,
                                    dB = NULL,
                                    plot = FALSE,
                                    scale = FALSE,
                                    correction = "amplitude",
                                    noisereduction = 2)
    
  } else if (is.null(noise.red)) {
    # Get the spectrogram matrix
    spectro_res <- seewave::spectro(wave,
                                    wl = wl,
                                    norm = FALSE,
                                    dB = NULL,
                                    plot = FALSE,
                                    scale = FALSE,
                                    correction = "amplitude")
    
  }




  matrix <- spectro_res$amp

  # Calculate amp_max based on bit depth
  amp_max <- if (wave@bit == 16) {
    32768
  } else if (wave@bit == 24) {
    8388607
  } else if (wave@bit == 32) {
    2147483647
  } else {
    stop("Unsupported bit depth")
  }

  # Convert amplitude to dBFS
  spectrogram <- 20 * log10(abs(matrix) / amp_max)
  rm(matrix)

  # Apply the cutoff to filter out values below the threshold
  spectrogram[spectrogram < cutoff] <- NA

  if (plot) {
    # Define frequency ranges
    samp_rate <- wave@samp.rate
    freq_bins <- seq(0, samp_rate / 2, length.out = nrow(spectrogram))
    time <- seq(0, seewave::duration(wave), length.out = ncol(spectrogram))
    freq <- freq_bins / 1000

    if (ggplot) {
      # Convert the spectrogram matrix into a data frame for ggplot
      spectro_df <- melt(spectrogram)
      spectro_df$Frequency <- freq[spectro_df$Var1]
      spectro_df$Time <- time[spectro_df$Var2]

      # Use ggplot for plotting
      p <- ggplot(spectro_df, aes(x = Time, y = Frequency, fill = value)) +
        geom_tile(na.rm = TRUE) +
        scale_fill_viridis_c(na.value = "transparent", option = "D") +  # Use viridis with transparent NA values
        labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw() +
        theme(legend.position = "none",
              panel.grid.major = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.3),
              panel.grid.minor = element_line(color = scales::alpha("black", 0.2), linetype = "solid", linewidth = 0.05))

      print(p)

    } else {
      # Use base R for plotting
      # Use viridis colors and handle NA as transparent
      viridis_colors <- viridis(100, option = "D")

      # Plot the spectrogram using base R
      image(time, freq, t(spectrogram),
            col = viridis_colors, xlab = "Time (s)", ylab = "Frequency (kHz)",
            main = plot.title, useRaster = TRUE, axes = FALSE)

      # Add NA transparency by overwriting the plot with white rectangles for NA values
      na_matrix <- is.na(spectrogram)
      if (any(na_matrix)) {
        rect(time[col(na_matrix)][na_matrix], freq[row(na_matrix)][na_matrix],
             time[col(na_matrix) + 1][na_matrix], freq[row(na_matrix) + 1][na_matrix],
             col = "white", border = NA)
      }

      # Add axes and grid lines
      axis(1, at = pretty(time), labels = pretty(time))
      axis(2, at = pretty(freq), labels = pretty(freq))
      abline(h = pretty(freq), col = "lightgray", lty = "dashed")
      abline(v = pretty(time), col = "lightgray", lty = "dashed")
      box()
    }
  }
  
  if (plot) {
    invisible(list(matrix = spectrogram, plot = p))
  } else {
    invisible(spectrogram)
  }

}
