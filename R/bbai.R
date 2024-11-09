#' Broadband Activity Index
#'
#' This function processes an audio signal to detect broadband activity by identifying 'clicks' based on time-frame-wise (i.e., column-wise) amplitude changes in the spectrogram. It computes statistics related to click height, variance, and centroid frequency, and can plot a spectrogram with detected clicks highlighted. The function also classifies whether the signal contains noise or insect based on the variance and centroid frequencies of the clicks.
#'
#' @param wave A `Wave` object containing the audio signal to be analyzed.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. The amplitude threshold (in dBFS) for removing low-amplitude values in the spectrogram. Default is `-50`.
#' @param click.height Numeric. The minimum height (in frequency bins) for a detected click to be kept. Default is `10`.
#' @param difference Numeric. The maximum difference in amplitude between adjacent frequency bins to be considered part of a single 'click'. Default is `20`.
#' @param gap.allowance Numeric. The size of gaps (in frequency bins) allowed between contiguous parts of a click. Default is `2`. Gaps larger than this value will split clicks.
#' @param spectrogram Logical. Should a spectrogram with highlighted clicks be plotted? Default is `TRUE`.
#' @param dark.plot Logical. Should the plot use a dark theme (black background)? Default is `FALSE`.
#' @param plot.title Character. The title for the plot, if `plot` is `TRUE`. Default is `NULL`.
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#' @return A tibble and optional spectrogram. 
#' @export
#' @import ggplot2
#' @importFrom tuneR channel mono
#' @importFrom seewave duration fir
#' @importFrom reshape2 melt
#' @importFrom scales col_numeric
#' @importFrom tibble add_column
#'
#' @examples bbai(wave)
bbai <- function(wave,
                 channel = "left",
                 hpf = 0,
                 freq.res = 50,
                 cutoff = -60,
                 click.length = 10,
                 difference = 10,
                 gap.allowance = 2,
                 spectrogram = FALSE,
                 dark.plot = FALSE,
                 plot.title = NULL,
                 verbose = TRUE) {
  
  if (!channel %in% c("left", "right", "mix", "each")) {
    stop("Invalid channel selected. Choose from 'left', 'right', 'mix', or 'each' (for stereo).")
  }
  
  process_channel <- function(wave, channel) {
    if (channel == "left") return(tuneR::channel(wave, "left"))
    if (channel == "right") return(tuneR::channel(wave, "right"))
    if (channel == "mix") return(tuneR::mono(wave, "both"))
  }
  
  if (wave@stereo) {
    # Stereo case: process the selected channel
    if (channel == "each") {
      # Process left and right channels separately for stereo
      wave.left <- process_channel(wave, "left")
      wave.right <- process_channel(wave, "right")
    } else {
      # Process selected channel (left, right, or mix)
      wave <- process_channel(wave, channel)
    }
  } else {
    # Mono case: force `channel` to "left" for unsupported options
    if (channel %in% c("each", "mix", "right")) {
      if (verbose) {
        cat("This is a mono recording. Calculating BBAI over the 'left' channel.\n")
      }
      channel <- "left"
    }
  }
  
  total_duration <- seewave::duration(wave)
  samp_rate <- wave@samp.rate
  
  bbai_mono <- function(wave,
                        hpf,
                        freq.res,
                        cutoff,
                        click.length,
                        difference,
                        gap.allowance,
                        spectrogram,
                        dark.plot,
                        plot.title){
    
    total_duration <- seewave::duration(wave)
    samp_rate <- wave@samp.rate
    
    # Apply high-pass filter
    if (hpf > 0) {
      wave <- seewave::fir(wave, wl = 1024, from = hpf, to = NULL, bandpass = TRUE, output = "Wave")
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }
    
    matrix <- spectrogram_cutoff(wave,
                                 freq.res = freq.res,
                                 cutoff = cutoff,
                                 noise.red = "rows")
    
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
          centroid <- sum(click_frequencies * click_amplitudes, na.rm = TRUE) / sum(click_amplitudes, 
                                                                                    na.rm = TRUE)
          centroids <- c(centroids, centroid)
        }
        pos <- pos + len
      }
    }
    
    # Compute the statistics for the click heights
    click_sum <- ifelse(length(click_heights) > 0, sum(click_heights), 0)
    click_mean <- ifelse(length(click_heights) > 0, round(mean(click_heights), 1), 0)
    click_variance <- ifelse(length(click_heights) > 0, round(var(click_heights), 1), 0)
    click_sd <- ifelse(length(click_heights) > 0, round(sd(click_heights), 1), 0)
    
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
    
    click_frames_prop <- round(click_time_frames / n_time_frames, 1)
    
    click_rate <- round(click_time_frames / total_duration, 1)
    
    n_clicks <- length(click_heights)
    
    cat("Broadband Activity: ", broadband_activity, "\n")
    
    summary <- tibble(
      index = "bbai",
      value = broadband_activity,
      nclicks = n_clicks,
      prop.clicks = click_frames_prop,
      click.rate = click_rate,
      mean.length = click_mean,
      var.length = click_variance,
      sd.length = click_sd,
      mean.centroid = mean_centroid,
      sd.centroid = sd_centroid,
      var.centroid = var_centroid,
      mean.click.dist = mean_all_click_dist
    )
    
    if(spectrogram){
      
      time_values <- seq(0, total_duration, length.out = n_time_frames)
      freq_values <- seq(0, samp_rate / 2, length.out = n_freq_bins) / 1000  # Frequency in kHz
      
      # Create a combined matrix to overlay clicks in red
      combined_matrix <- matrix
      combined_matrix[click_matrix] <- 0  # Set click positions to 0 dB for clarity
      
      # Convert matrix to a dataframe for plotting
      plot_data <- as.data.frame(combined_matrix)
      colnames(plot_data) <- time_values
      plot_data <- cbind(Frequency = freq_values, plot_data)
      
      # Reshape for ggplot
      plot_data <- reshape2::melt(plot_data, id.vars = "Frequency", variable.name = "Time", value.name = "dB")
      plot_data$Time <- as.numeric(as.character(plot_data$Time))
      
      # Mark clicks as a separate layer in the data
      plot_data$Click <- as.vector(click_matrix)
      
      # Set color gradient function for dB values, with transparent `na.value`
      color_func <- scales::col_numeric(palette = c(if (dark.plot) "black" else "white", "#2c7bb6","#00a6ca","#00ccbc","#90eb9d",
                                                    "#ffff8c", "#f9d057", "#f29e2e","#e76818","#d7191c"),
                                        domain = c(cutoff, 0), na.color = "transparent")
      
      plot_data$Color <- color_func(plot_data$dB)
      
      
      # Create the base plot
      p <- ggplot(plot_data, aes(x = Time, y = Frequency)) +
        geom_tile(aes(fill = Color), color = NA) +
        scale_fill_identity(na.value = "transparent") +  # Make NA values transparent
        labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
        theme_bw() +
        theme(legend.position = "none")
      
      # Highlight clicks in red
      p <- p + geom_tile(data = subset(plot_data, Click == TRUE), fill = "red", color = NA) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0))
      
      
      if (dark.plot) {
        p <- p + theme(
          plot.background = element_rect(fill = "black", color = NA),
          panel.background = element_rect(fill = "black", color = NA),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white"),
          axis.line = element_line(color = "white"),
          axis.ticks = element_line(color = "white"),
          panel.grid.major = element_line(color = scales::alpha("white", 0.2), 
                                          linetype = "solid", linewidth = 0.3),
          panel.grid.minor = element_line(color = scales::alpha("white", 0.2), 
                                          linetype = "solid", linewidth = 0.05),
          plot.title = element_text(color = "white", face = "bold", size = 14)
        )
      }
      print(p)
      
      return(list(summary = summary, spectrogram = p))
      
    } else {
      return(list(summary = summary, spectrogram = NULL))
    }
  }
  
  # Calculate the index based on the stereo condition
  if (wave@stereo) {
    if(channel == "each") {
      if (verbose) cat("Calculating Broadband Activity Index on 2 channels... \n")
      bbai_left <- bbai_mono(wave.left,
                             hpf = hpf,
                             freq.res = freq.res,
                             cutoff = cutoff,
                             click.length = click.length,
                             difference = difference,
                             gap.allowance = gap.allowance,
                             spectrogram = spectrogram,
                             dark.plot = dark.plot,
                             plot.title = plot.title)
      bbai_right <- bbai_mono(wave.right,
                              hpf = hpf,
                              freq.res = freq.res,
                              cutoff = cutoff,
                              click.length = click.length,
                              difference = difference,
                              gap.allowance = gap.allowance,
                              spectrogram = spectrogram,
                              dark.plot = dark.plot,
                              plot.title = plot.title)
      
      summary <- tibble::tibble(
        index = "bbai",
        channel = "each",
        value_l = bbai_left$summary$value,
        value_r = bbai_right$summary$value,
        value_avg = round((bbai_left$summary$value + bbai_right$summary$value) / 2, 2),
        n_clicks_l = bbai_left$summary$nclicks,
        n_clicks_r = bbai_right$summary$nclicks,
        n_clicks_avg = round((bbai_left$summary$nclicks + bbai_right$summary$nclicks) / 2, 1),
        mean_length_l = bbai_left$summary$mean.length,
        mean_length_r = bbai_right$summary$mean.length,
        mean_length_avg = round((bbai_left$summary$mean.length + bbai_right$summary$mean.length) / 2, 1),
        var_length_l = bbai_left$summary$var.length,
        var_length_r = bbai_right$summary$var.length,
        var_length_avg = round((bbai_left$summary$var.length + bbai_right$summary$var.length) / 2, 1),
        sd_length_l = bbai_left$summary$sd.length,
        sd_length_r = bbai_right$summary$sd.length,
        sd_length_avg = round((bbai_left$summary$sd.length + bbai_right$summary$sd.length) / 2, 1),
        mean_centroid_l = bbai_left$summary$mean.centroid,
        mean_centroid_r = bbai_right$summary$mean.centroid,
        mean_centroid_avg = round((bbai_left$summary$mean.centroid + bbai_right$summary$mean.centroid) / 2, 1),
        sd_centroid_l = bbai_left$summary$sd.centroid,
        sd_centroid_r = bbai_right$summary$sd.centroid,
        sd_centroid_avg = round((bbai_left$summary$sd.centroid + bbai_right$summary$sd.centroid) / 2, 1),
        var_centroid_l = bbai_left$summary$var.centroid,
        var_centroid_r = bbai_right$summary$var.centroid,
        var_centroid_avg = round((bbai_left$summary$var.centroid + bbai_right$summary$var.centroid) / 2, 1),
        mean_click_dist_l = bbai_left$summary$mean.click.dist,
        mean_click_dist_r = bbai_right$summary$mean.click.dist,
        mean_click_dist_avg = round((bbai_left$summary$mean.click.dist + bbai_right$summary$mean.click.dist) / 2, 1)
      )
      
      print(summary)
      
      if(spectrogram){
        
        invisible(list(summary = summary,
                       spectrogram_l = bbai_left$spectrogram,
                       spectrogram_r = bbai_right$spectrogram))
        
      } else {
        
        invisible(summary)
        
      }
      
    } else if (channel %in% c("mix", "left", "right")) {

      bbai_global <- bbai_mono(wave,
                               hpf = hpf,
                               freq.res = freq.res,
                               cutoff = cutoff,
                               click.length = click.length,
                               difference = difference,
                               gap.allowance = gap.allowance,
                               spectrogram = spectrogram,
                               dark.plot = dark.plot,
                               plot.title = plot.title)
      
      bbai_global$summary <- bbai_global$summary |>
        # add_column(channel = channel, .after = index)
        mutate(channel = channel) |>
        relocate(channel, .after = index)
      
      print(bbai_global$summary)
     
      if(spectrogram){
        invisible(list(summary = bbai_global$summary, 
                       spectrogram = bbai_global$spectrogram))
        
      } else {
        invisible(bbai_global$summary)
      }
      
      
    } 
    
  } else { # Mono file
    bbai_global <- bbai_mono(wave,
                             hpf = hpf,
                             freq.res = freq.res,
                             cutoff = cutoff,
                             click.length = click.length,
                             difference = difference,
                             gap.allowance = gap.allowance,
                             spectrogram = spectrogram,
                             dark.plot = dark.plot,
                             plot.title = plot.title)
    
    bbai_global$summary <- bbai_global$summary |>
      # add_column(channel = channel, .after = index)
      mutate(channel = channel) |>
      relocate(channel, .after = index)
    
    if(verbose){
      print(bbai_global$summary)
    }
    
    if(spectrogram){
      invisible(list(summary = bbai_global$summary, 
                     spectrogram = bbai_global$spectrogram))
      
    } else {
      invisible(bbai_global$summary)
    }
   
    
  }
  
}

