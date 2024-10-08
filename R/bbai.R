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
#' @importFrom seewave duration rmoffset fir
#' @importFrom reshape2 melt
#'
#' @examples bbai(wave)
bbai <- function(wave,
                 channel = 'each',
                 hpf = 0,
                 rm.offset = TRUE,
                 freq.res = 50,
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
      wave <- seewave::fir(wave, wl = 1024, from = hpf, to = NULL, bandpass = TRUE, output = "Wave")
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }

    matrix <- spectrogram_cutoff(wave,
                                 freq.res = freq.res,
                                 cutoff = cutoff)

    # Initialize the number of time frames with clicks and list for click heights
    click_time_frames <- 0
    click_heights <- c()
    centroids <- c()

    n_freq_bins <- nrow(matrix)
    n_time_frames <- ncol(matrix)

    # Create a logical matrix to store click detection
    click_matrix <- matrix(FALSE, nrow = n_freq_bins, ncol = n_time_frames)

    freq_values <- seq(0, wave@samp.rate / 2, length.out = n_freq_bins) / 1000

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
    click_sum <- ifelse(length(click_heights) > 1, sum(click_heights), 0)
    click_mean <- ifelse(length(click_heights) > 1, round(mean(click_heights), 1), 0)
    click_variance <- ifelse(length(click_heights) > 1, round(var(click_heights), 1), 0)
    click_sd <- ifelse(length(click_heights) > 1, round(sd(click_heights), 1), 0)

    # Analyze centroid frequencies
    mean_centroid <- ifelse(length(centroids) > 1, round(mean(centroids), 1), 0)
    sd_centroid <- ifelse(length(centroids) > 1, round(sd(centroids), 1), 0)
    var_centroid <- ifelse(length(centroids) > 1, round(var(centroids), 1), 0)

    # Identify click clusters by checking gaps between click frames
    click_times <- which(apply(click_matrix, 2, any))
    mean_all_click_dist <- ifelse(length(click_times) > 1, round(mean(diff(click_times))), 0)

    # Calculate the total number of cells in the spectrogram
    total_cells <- n_freq_bins * n_time_frames

    # Step 3: Calculate the proportion of clicks
    broadband_activity <- round((click_sum / total_cells) * 100, 1)

    n_clicks <- length(click_heights)

    cat("Broadband Activity: ", broadband_activity, "\n")

    summary <- tibble(
      index = "bbai",
      value = broadband_activity,
      nclicks = n_clicks,
      mean.length = click_mean,
      var.length = click_variance,
      sd.length = click_sd,
      mean.centroid = mean_centroid,
      sd.centroid = sd_centroid,
      var.centroid = var_centroid,
      mean.click.dist = mean_all_click_dist
    )

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
      value_avg = round((bbai_left$value + bbai_right$value) / 2, 2),
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
