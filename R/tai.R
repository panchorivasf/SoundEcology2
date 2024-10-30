#' Trill Activity Index
#'
#' This function calculates the trill index of an audio wave object by analyzing the frequency modulation pattern over time. It can operate in either binary or continuous modes and provides options for generating visual plots of the trill activity. The function also identifies potential noise issues in the low- and mid-frequency ranges.
#'
#' @param wave A wave object to be analyzed.
#' @param channel Channel or channels to be analyzed. Options are "left", "right", "each", and "mix".
#' @param cutoff Numeric. The cutoff in decibels for the spectrogram generation.
#' @param n.windows Numeric. Number of time windows to divide the signal into for analysis. Default is 60.
#' @param freq.res Numeric. Frequency resolution (in Hz) for the spectrogram analysis. Default is 100.
#' @param plot Logical. If TRUE, generates a plot of the trill index over time. Default is TRUE.
#' @param plot.title Character. The title to be used for the plot. Default is NULL.
#' @param verbose Logical. If TRUE, provides detailed output during the function's execution. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{summary}: A tibble summarizing TAI statistics, including values for low and mid-frequency noise.
#'   \item \code{spectral}: A 1-column matrix of mean trill index values for each frequency bin
#' }
#'
#' @export
#' @importFrom seewave rmoffset duration
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' wave <- read_audio("example.wav")
#' trill_results <- tai(wave, channel = "left", binary = TRUE, plot = TRUE)
#' }
tai <- function(wave,
                channel = 'each',
                hpf = 0,
                rm.offset = TRUE,
                cutoff = -60,
                n.windows = 120,
                freq.res = 100,
                plot = FALSE,
                plot.title = NULL,
                verbose = FALSE) {


  # Check if the wave is stereo
  if (wave@stereo) {
    if (channel == 'each') {
      wave.left <- tuneR::channel(wave, 'left')
      wave.right <- tuneR::channel(wave, 'right')
    } else if (channel == 'left') {
      wave <- tuneR::channel(wave, 'left')
    } else if (channel == 'right') {
      wave <- tuneR::channel(wave, 'right')
    } else if (channel == 'mix'){
      wave <- tuneR::mono(wave, "both")
    } else {
      stop("Invalid channel selected.")
    }
  }

  total_duration <- seewave::duration(wave)
  samp_rate <- wave@samp.rate

  calculate_index <- function(wave,
                              channel,
                              hpf,
                              rm.offset,
                              cutoff,
                              n.windows,
                              freq.res,
                              plot,
                              plot.title,
                              verbose){

    if(rm.offset){
      wave <- seewave::rmoffset(wave, output = "Wave")
    }

    # Apply high-pass filter
    if (hpf > 0) {
      wave <- seewave::fir(wave, wl = 1024, from = hpf, to = NULL, bandpass = TRUE, output = "Wave")
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }


    spectrogram <- spectrogram_cutoff(wave,
                                      freq.res = freq.res,
                                      plot = FALSE)



    trill_spectral <- matrix(nrow = nrow(spectrogram), ncol = ncol(spectrogram))

    wave.dur = seewave::duration(wave)
    j = round(ncol(spectrogram)/n.windows)

    # Apply scoring mechanism for trill, considering the cutoff
    for (freq_index in 1:nrow(spectrogram)) {
      for (start_col in seq(1, ncol(spectrogram), by = j)) {
        end_col <- min(ncol(spectrogram), start_col + j - 1)
        window <- spectrogram[freq_index, start_col:end_col]

        # Skip window if it contains NA (below cutoff)
        if (any(is.na(window))) {
          trill_spectral[freq_index, start_col:end_col] <- 0  # Set to zero if the window is below cutoff
        } else {
          score <- sum(abs(diff(window)))  # Calculate trill score (you can adjust the metric here)
          trill_spectral[freq_index, start_col:end_col] <- ifelse(length(window) >= floor(j / 2), score, 0)
        }
      }
    }

    # Calculate trill mean for each frequency bin
    trill_mean <- round(rowMeans(trill_spectral, na.rm = TRUE),1)

    if(verbose){
      cat("Trill scores per bin:\n", trill_mean, "\n")
    }

    # Check for possible low-frequency noise in the first 10 frequency bins
    low_freq_bins <- trill_spectral[1:10, ]  # Extract the first 10 frequency bins ~0-1 kHz
    low_freq_means <- colMeans(low_freq_bins, na.rm = TRUE)  # Calculate the mean for each time frame across the low frequency bins

    mid_freq_bins <- trill_spectral[11:20,]
    mid_freq_means <- colMeans(mid_freq_bins, na.rm = TRUE) # Means of time frames across mid-frequency bins (1-2kHz, e.g., traffic noise)

    high_freq_bins <- trill_spectral[21:nrow(trill_spectral),]
    high_freq_mean <- round(mean(colMeans(high_freq_bins, na.rm = TRUE),1)) # Average values above ~2 kHz

    # Count how many time frames have a mean >= 10
    time_frames_with_noise <- sum(low_freq_means >= 10)
    time_frames_with_midf_noise <- sum(mid_freq_means >= 10)

    percent_lowf_noisy_frames <- round((time_frames_with_noise *100)/ ncol(trill_spectral))
    percent_midf_noisy_frames <- round((time_frames_with_midf_noise *100)/ ncol(trill_spectral))

    trill_summary <- tibble(
      index = "trill activity",
      value = round(mean(trill_mean, na.rm = TRUE), 1),
      value.a2k = high_freq_mean,
      lowf.noise = percent_lowf_noisy_frames,
      midf.noise = percent_midf_noisy_frames,
      freq.res = freq.res,
      f.bins = nrow(trill_spectral),
      t.frames = ncol(trill_spectral)
    )

    if (plot) {
      # Prepare data for plotting
      samp_rate <- wave@samp.rate
      freq_bins <- seq(0, samp_rate / 2, length.out = nrow(spectrogram))
      time <- seq(0, seewave::duration(wave), length.out = ncol(spectrogram))
      freq <- freq_bins / 1000

      trill_df <- melt(trill_spectral)

      trill_df$value[trill_df$value == 0] <- NA # To allow transparency of zero values

      trill_df$frequency <- freq[trill_df$Var1]  # Frequency
      trill_df$time <- time[trill_df$Var2]  # Time

      p <- ggplot(trill_df, aes(time, frequency, fill = value)) +
        geom_tile() +
        scale_fill_viridis_c(na.value = "transparent", option = "viridis", name = "Trills") +
        labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
          plot.background = element_rect(fill = "black", color = NA),
          panel.background = element_rect(fill = "black", color = NA),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          axis.ticks = element_line(color = "white"),
          axis.line = element_line(color = "white"),
          legend.position = "right",
          legend.background = element_rect(fill = "black", color = NA),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          plot.title = element_text(color = "white"),
          panel.grid.major = element_line(color = alpha("white", 0.5), linewidth = 0.3),
          panel.grid.minor = element_line(color = alpha("white", 0.2), linewidth = 0.2)
        ) +
        annotate("label", x = max(trill_df$time) * 0.8, y = max(trill_df$frequency) * 0.8,
                 label = paste("TAI: ", trill_summary$value,
                               "\nTAI a2k: ", trill_summary$value.a2k,
                               "\nLowf noise: ", trill_summary$lowf.noise, "%",
                               "\nMidf noise: ", trill_summary$midf.noise, "%",
                               "\nThreshold: ", cutoff, "dB",
                               "\nTime Step: ", round(wave.dur/n.windows, 2), "s"
                 ),
                 hjust = 0, vjust = 1, size = 3, color = "black",
                 label.size = 0.5, label.padding = unit(0.5, "lines"),
                 fontface = "italic", fill = "white")

      print(p)
    }

    return(list(summary = trill_summary, spectral = trill_mean))
  }

  # Calculate the index based on the stereo condition
  if (channel == 'each') {
    if(!wave@stereo){
      stop("Can't select 'each' channel for a mono file.\n")
    }

    if (plot) {
      stop("Plotting is not allowed when calculating TAI over both channels.\n")
    }


    if (verbose) cat("Calculating Trill Activity Index on 2 channels... \n")

    tai_left <- calculate_index(wave.left,
                                hpf = hpf,
                                rm.offset = rm.offset,
                                cutoff = cutoff,
                                n.windows = n.windows,
                                freq.res = freq.res,
                                plot = FALSE,
                                plot.title = NULL,
                                verbose = FALSE)
    tai_right <- calculate_index(wave.right,
                                 hpf = hpf,
                                 rm.offset = rm.offset,
                                 cutoff = cutoff,
                                 n.windows = n.windows,
                                 freq.res = freq.res,
                                 plot = FALSE,
                                 plot.title = NULL,
                                 verbose = FALSE)

    tai_left_summary <- tai_left$summary
    tai_right_summary <- tai_right$summary

    tai_left_spectral <- tai_left$spectral
    tai_right_spectral <- tai_right$spectral

    tai_global <- tibble::tibble(
      index = "tai",
      value_l = tai_left_summary$value,
      value_r = tai_right_summary$value,
      value_avg = round((tai_left_summary$value + tai_right_summary$value) / 2, 3),
      low_freq_noise_l = tai_left_summary$lowf.noise,
      low_freq_noise_r = tai_right_summary$lowf.noise,
      low_freq_noise_avg = round((tai_left_summary$lowf.noise + tai_right_summary$lowf.noise) / 2, 3),
      mid_freq_noise_l = tai_left_summary$midf.noise,
      mid_freq_noise_r = tai_right_summary$midf.noise,
      mid_freq_noise_avg = round((tai_left_summary$midf.noise + tai_right_summary$midf.noise) / 2, 3)
    )

    if(verbose){
      print(tai_global)
    }


    invisible(list(summary = tai_global, spectral_left = tai_left_spectral, spectral_right = tai_right_spectral))




  } else {

    if (verbose) {
      if(channel == 'left'){
        cat("Calculating Trill Activity Index on the left channel... \n")

      }else if(channel == 'right'){
        cat("Calculating Trill Activity Index on the right channel... \n")

      }else if(channel == 'mix'){
        cat("Calculating Trill Activity Index on a mix of the two channels... \n")
      }
    }

    # Calulate NBAI
    tai_global <- calculate_index(wave,
                                  hpf = hpf,
                                  rm.offset = rm.offset,
                                  cutoff = cutoff,
                                  n.windows = n.windows,
                                  freq.res = freq.res,
                                  plot = FALSE,
                                  plot.title = NULL,
                                  verbose = FALSE)

    summary <- tai_global$summary
    spectral <- tai_global$spectral

    if(verbose){
      print(summary)

    }

    invisible(list(summary, spectral))


  }
}
