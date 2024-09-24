#' Broadband Activity Index
#'
#' This function processes an audio signal to detect broadband activity by identifying 'clicks' based on time-frame-wise (i.e., column-wise) amplitude changes in the spectrogram. It computes statistics related to click height, variance, and centroid frequency, and can plot a spectrogram with detected clicks highlighted. The function also classifies whether the signal contains noise or insect based on the variance and centroid frequencies of the clicks.
#'
#' @param wave A `Wave` object containing the audio signal to be analyzed.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param rm.offset Logical. Should the DC offset be removed from the audio signal? Defaults to `TRUE`.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. The amplitude threshold (in dBFS) for removing low-amplitude values in the spectrogram. Default is `-50`.
#' @param click.height Numeric. The minimum height (in frequency bins) for a detected click to be kept. Default is `10`.
#' @param difference Numeric. The maximum difference in amplitude between adjacent frequency bins to be considered part of a single 'click'. Default is `20`.
#' @param gap.allowance Numeric. The size of gaps (in frequency bins) allowed between contiguous parts of a click. Default is `2`. Gaps larger than this value will split clicks.
#' @param plot Logical. Should a spectrogram with highlighted clicks be plotted? Default is `TRUE`.
#' @param dark.plot Logical. Should the plot use a dark theme (black background)? Default is `FALSE`.
#' @param plot.title Character. The title for the plot, if `plot` is `TRUE`. Default is `NULL`.
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.

#' @return A tibble containing the following columns:
#'   - `index`: The name of the index. Useful later when merging data with other indices.
#'   - `value`: The number of clicks detected in the recording.
#'   - `mean`: The mean click height (in frequency bins).
#'   - `variance`: The variance of the click height.
#'   - `sd`: The standard deviation of the click height.

#' @export
#'
#' @importFrom seewave rmoffset
#' @importFrom seewave spectro
#' @importFrom seewave duration
#' @importFrom reshape2 melt
#'
#' @examples bbai(wave)
bbai <- function(wave,
                 channel = 'each',
                 hpf = 0,
                 rm.offset = TRUE,
                 freq.res = 100,
                 cutoff = -60,
                 click.length = 10,
                 difference = 10,
                 gap.allowance = 2,
                 plot = FALSE,
                 dark.plot = FALSE,
                 plot.title = NULL,
                 verbose = TRUE) {


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
                              rm.offset,
                              hpf,
                              freq.res,
                              cutoff){


    # Remove DC offset
    if (rm.offset) {
      wave <- seewave::rmoffset(wave, output = "Wave")
    }


    # Apply high-pass filter
    if (hpf > 0) {
      wave <- fir(wave, from=hpf, to = NULL, bandpass = TRUE, output = "Wave", wl = 1024)
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }

    matrix <- spectrogram_cutoff(wave,
                                 freq.res = freq.res,
                                 cutoff = cutoff)


    # Initialize the number of time frames with clicks and list for click heights
    click_time_frames <- 0
    click_heights <- c()  # Store lengths of clicks in frequency bins
    centroids <- c()      # Store centroids of clicks

    n_freq_bins <- nrow(matrix)
    n_time_frames <- ncol(matrix)

    # Create a logical matrix to store click detection
    click_matrix <- matrix(FALSE, nrow = n_freq_bins, ncol = n_time_frames)

    freq_values <- seq(0, wave@samp.rate / 2, length.out = n_freq_bins) / 1000  # Frequency in kHz

    for (i in 1:n_time_frames) {
      diff_vec <- diff(matrix[, i], na.rm = TRUE)

      if (all(is.na(diff_vec))) next

      is_small_diff <- abs(diff_vec) < difference
      is_small_diff[is.na(is_small_diff)] <- FALSE

      gap_counter <- 0
      contiguous_blocks <- logical(length(is_small_diff) + 1)

      for (j in seq_along(is_small_diff)) {
        if (is_small_diff[j]) {
          contiguous_blocks[j] <- TRUE
          gap_counter <- 0
        } else {
          if (gap_counter < gap.allowance) {
            gap_counter <- gap_counter + 1
            contiguous_blocks[j] <- TRUE
          } else {
            gap_counter <- 0
          }
        }
      }

      rle_result <- rle(contiguous_blocks)
      lengths <- rle_result$lengths
      values <- rle_result$values

      pos <- 1

      for (k in seq_along(lengths)) {
        len <- lengths[k]
        val <- values[k]

        if (val && len > click.length) {
          click_matrix[pos:(pos + len - 1), i] <- TRUE
          click_time_frames <- click_time_frames + 1
          click_heights <- c(click_heights, len)

          click_frequencies <- freq_values[pos:(pos + len - 1)]
          click_amplitudes <- matrix[pos:(pos + len - 1), i]
          centroid <- sum(click_frequencies * click_amplitudes, na.rm = TRUE) / sum(click_amplitudes, na.rm = TRUE)
          centroids <- c(centroids, centroid)
        }
        pos <- pos + len
      }
    }

    # Compute the statistics for the click heights
    if (length(click_heights) > 1) {
      click_sum <- sum(click_heights)
      click_mean <- round(mean(click_heights), 1)
      click_variance <- round(var(click_heights), 1)
      click_sd <- round(sd(click_heights), 1)
    } else {
      click_sum <- NA
      click_mean <- NA
      click_variance <- NA
      click_sd <- NA
    }

    # Analyze centroid frequencies
    if (length(centroids) > 1) {
      mean_centroid <- round(mean(centroids), 1)
      sd_centroid <- round(sd(centroids), 1)
      var_centroid <- round(var(centroids), 1)
    } else {
      mean_centroid <- NA
      sd_centroid <- NA
      var_centroid <- NA
    }


    # Identify click clusters by checking gaps between click frames
    click_times <- which(apply(click_matrix, 2, any))

    if (length(click_times) > 1) {
      click_diffs <- diff(click_times)

      # Calculate the average distance between all clicks (mean.all.click.dist)
      mean_all_click_dist <- round(mean(click_diffs))

    } else {
      # mean_click_dist_in_clust <- NA
      mean_all_click_dist <- NA
      # n_clusters <- NA
    }

    # Calculate the total number of cells in the spectrogram
    total_cells <- n_freq_bins * n_time_frames

    # Step 3: Calculate the proportion of clicks
    broadband_activity <- round((click_sum / total_cells)*100,1)

    n_clicks <- length(click_heights)


    # if (is.na(broadband_activity)) {
    #   broadband_activity <- 0
    #   noise <- FALSE
    #   biophony <- FALSE
    #   # } else if (mean_centroid > 5 || (!is.na(mean_all_click_dist) & mean_all_click_dist < 7 || n_clicks > 1000)) {
    # } else if (mean_centroid > 5 || n_clicks > 1000) {
    #
    #   biophony <- TRUE
    # } else {
    #   biophony <- FALSE
    # }
    #

    # # Check for noise
    # if (click_variance > 100 || var_centroid > 10 || mean_centroid < 3  || broadband_activity < 1) {
    #   noise <- TRUE
    # } else {
    #   noise <- FALSE
    # }
    #

    cat("Broadband Activity: ", broadband_activity, "\n")

    # Update the data frame with the new metric
    summary <- tibble(
      index = "bbai",
      value = broadband_activity,
      nclicks = length(click_heights),
      mean.length = click_mean,
      var.length = click_variance,
      sd.length = click_sd,
      mean.centroid = mean_centroid,
      sd.centroid = sd_centroid,
      var.centroid = var_centroid,
      # n.click.clusters = n_clusters,
      mean.click.dist = mean_all_click_dist
      # clus.mean.click.dist = mean_click_dist_in_clust,
      # biophony = biophony,
      # noise = noise
    )

    # If plot is TRUE, generate the spectrogram with clicks highlighted
    if (plot) {
      time_values <- seq(0, total_duration, length.out = n_time_frames)
      freq_values <- seq(0, by = samp_rate / 2 / n_freq_bins, length.out = n_freq_bins) / 1000

      if (!is.null(click_matrix) && any(click_matrix)) {
        combined_matrix <- matrix
        combined_matrix[click_matrix] <- 0

        plot_data <- as.data.frame(combined_matrix)
        colnames(plot_data) <- time_values
        plot_data <- cbind(Frequency = freq_values, plot_data)
        plot_data <- reshape2::melt(plot_data, id.vars = "Frequency", variable.name = "Time", value.name = "dB")
        plot_data$Time <- as.numeric(as.character(plot_data$Time))
        plot_data$Click <- as.vector(click_matrix)

        color_func <- scales::col_numeric(palette = c(if (dark.plot) "black" else "white", "#2c7bb6","#00a6ca","#00ccbc","#90eb9d",
                                                      "#ffff8c", "#f9d057", "#f29e2e","#e76818","#d7191c"),
                                          domain = c(-60, 0), na.color = "transparent")

        plot_data$Color <- color_func(plot_data$dB)

        p <- ggplot(plot_data, aes(x = Time, y = Frequency)) +
          geom_tile(aes(fill = Color), color = NA) +
          scale_fill_identity(na.value = "transparent") +
          labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
          theme_bw() +
          theme(legend.position = "none")

        p <- p + geom_tile(data = subset(plot_data, Click == TRUE), fill = "red", color = NA) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +

          annotate("label", x = max(plot_data$Time) * 0.8, y = max(plot_data$Frequency) * 0.95,
                   label = paste("BBAI:", summary$value,
                                 "\nN. clicks:", summary$nclicks,
                                 # "\nNoise:", summary$noise,
                                 # "\nInsect:", summary$biophony,
                                 "\nMean.click.dist.:", summary$mean.click.dist, "fr",
                                 # "\nClick length stats:",
                                 "\nMean.length:", summary$mean.length,
                                 "\nVar.length:", summary$var.length,
                                 "\nSD.length:", summary$sd.length,
                                 # "\nClick centroid stats:",
                                 "\nMean.cent:", summary$mean.centroid, "kHz",
                                 "\nVar.cent:", summary$var.centroid, "kHz",
                                 "\nSD.cent:", summary$sd.centroid, "kHz",
                                 # "\nClick dist. stats:",
                                 # "\nN.clusters:", summary$n.click.clusters,
                                 # "\nClus.mean.click.dist.:", summary$clus.mean.click.dist, "fr",
                                 "\nParameters:",
                                 "\nCutoff:", cutoff, "dB",
                                 "\nMin. click h.:", click.length,
                                 "\nMax. amp. diff.:", difference
                   ),
                   hjust = 0, vjust = 1, size = 3.5, color = "black",
                   label.size = 0.5, label.padding = unit(0.5, "lines"),
                   fontface = "italic", fill = "white")

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

        print(p)
      } else {
        warning("No clicks detected, skipping plot generation.")
      }
    }

    return(summary)

  }


  # Calculate the index based on the stereo condition
  if (channel == 'each') {
    if(!wave@stereo){
      stop("Can't select 'each' channel for a mono file.\n")
    }

    if (plot) {
      stop("Plotting is not allowed when calculating BBAI over both channels.\n")
    }


    if (verbose) cat("Calculating Broadband Activity Index on 2 channels... \n")

    bbai_left <- calculate_index(wave.left, rm.offset = rm.offset,
                                 hpf = hpf, freq.res = freq.res,
                                 cutoff = cutoff)
    bbai_right <- calculate_index(wave.right, rm.offset = rm.offset,
                                  hpf = hpf, freq.res = freq.res,
                                  cutoff = cutoff)

    bbai_global <- tibble::tibble(
      index = "bbai",
      value_l = bbai_left$value,
      value_r = bbai_right$value,
      value_avg = round((bbai_left$value + bbai_right$value) / 2, 1),
      n_clicks_l = bbai_left$nclicks,
      n_clicks_r = bbai_right$nclicks,
      n_clicks_avg = round((bbai_left$nclicks + bbai_right$nclicks) / 2, 1),
      mean_length_l = bbai_left$mean.length,
      mean_length_r = bbai_right$mean.length,
      mean_length_avg = round((bbai_left$mean.length + bbai_right$mean.length) / 2, 1),
      var_length_l = bbai_left$var.length,
      var_length_r = bbai_right$var.length,
      var_length_avg = round((bbai_left$var.length + bbai_right$var.length) / 2, 1),
      sd_length_l = bbai_left$sd.length,
      sd_length_r = bbai_right$sd.length,
      sd_length_avg = round((bbai_left$sd.length + bbai_right$sd.length) / 2, 1),
      mean_centroid_l = bbai_left$mean.centroid,
      mean_centroid_r = bbai_right$mean.centroid,
      mean_centroid_avg = round((bbai_left$mean.centroid + bbai_right$mean.centroid) / 2, 1),
      sd_centroid_l = bbai_left$sd.centroid,
      sd_centroid_r = bbai_right$sd.centroid,
      sd_centroid_avg = round((bbai_left$sd.centroid + bbai_right$sd.centroid) / 2, 1),
      var_centroid_l = bbai_left$var.centroid,
      var_centroid_r = bbai_right$var.centroid,
      var_centroid_avg = round((bbai_left$var.centroid + bbai_right$var.centroid) / 2, 1),
      mean_click_dist_l = bbai_left$mean.click.dist,
      mean_click_dist_r = bbai_right$mean.click.dist,
      mean_click_dist_avg = round((bbai_left$mean.click.dist + bbai_right$mean.click.dist) / 2, 1)
    )

    if(verbose){
      print(bbai_global)
    }


    invisible(bbai_global)




  } else {

    if (verbose) {
      if(channel == 'left'){
        cat("Calculating Broadband Activity Index on the left channel... \n")

      }else if(channel == 'right'){
        cat("Calculating Broadband Activity Index on the right channel... \n")

      }else if(channel == 'mix'){
        cat("Calculating Broadband Activity Index on a mix of the two channels... \n")
      }
    }

    # Calulate NBAI
    bbai_global <- calculate_index(wave,
                                   rm.offset = rm.offset,
                                   channel = channel,
                                   hpf = hpf,
                                   freq.res = freq.res,
                                   cutoff = cutoff)

    if(verbose){
      print(bbai_global)

    }

    invisible(bbai_global)


  }

}
