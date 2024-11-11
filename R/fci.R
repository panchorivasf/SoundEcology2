#' Frequency Cover Indices (FCI)
#'
#' This function analyzes the spectral occupancy of 'low', 'mid', 'high', and 'ultra' frequency bands, which limits are defined by the user. It calculates the proportion of cells in each frequency band exceeding a specified dB cutoff and optionally plots the binary spectrogram with frequency bands delimited.
#'
#' @param wave A Wave object containing the audio data to be analyzed.
#' @param channel Character. If the Wave is stereo, select a channel to analyze. Options are "left", "right", and "mix", which combines L and R channels into a single mono Wave. Default is "left".
#' @param cutoff Numeric. The amplitude cutoff in dB above which the frequency bins are considered active. Default is 70 dB.
#' @param freq.res Numeric. The frequency resolution to be used when selecting a window length for the FFT.
#' @param spectrogram Logical. If TRUE, the function will generate and display a plot of the binary spectrogram with the frequency bands marked. Default is TRUE.
#' @param ggplot Logical. If TRUE, the plot is output is a ggplot object. If FALSE (default), the base plot function is used (faster rendering).
#' @param plot.title Character. An optional title for the plot.
#' @param lf.min Numeric. The minimum frequency (in Hz) for the low-frequency band. Default is 1 Hz.
#' @param lf.max Numeric. The maximum frequency (in Hz) for the low-frequency band. Default is 2000 Hz.
#' @param mf.min Numeric. The minimum frequency (in Hz) for the mid-frequency band. Default is 2000 Hz.
#' @param mf.max Numeric. The maximum frequency (in Hz) for the mid-frequency band. Default is 12000 Hz.
#' @param hf.min Numeric. The minimum frequency (in Hz) for the high-frequency band. Default is 12000 Hz.
#' @param hf.max Numeric. The maximum frequency (in Hz) for the high-frequency band. Default is 22000 Hz.
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#' 
#' @return A tibble and an optional spectrogram. 
#' @export
#' @import ggplot2
#' @importFrom tuneR channel mono
#' @importFrom seewave duration fir
#' @importFrom tibble tibble
#' @importFrom dplyr mutate relocate
#'
#' @examples
#' \dontrun{
#' # Assuming 'wave' is a Wave object
#' results <- fci(wave, cutoff = 65, freq.res = 100,  spectrogram = TRUE)
#' print(results)
#' }
fci <- function(wave,
               channel = 'left',
               hpf = 0,
               cutoff = -60,
               freq.res = 50,
               spectrogram = FALSE,
               ggplot = FALSE,
               plot.title = "Frequency Cover Analysis",
               sound.color = "#045E10",
               lf.min = 0,
               lf.max = 1500,
               mf.min = 1500,
               mf.max = 8000,
               hf.min = 8000,
               hf.max = 18000,
               uf.min = 18000,
               uf.max = 24000,
               verbose = FALSE
) {

  nyquist <- wave@samp.rate/2

  if(uf.max > nyquist){
    uf.max <- nyquist
    cat("The selected UF maximum is higher than the Nyquist. UF max has been changed to Nyquist.\n")
  }

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
                              hpf = 0,
                              cutoff = cutoff,
                              freq.res = freq.res,
                              spectrogram = spectrogram,
                              ggplot = ggplot,
                              plot.title = plot.title,
                              sound.color = sound.color,
                              lf.min = lf.min,
                              lf.max = lf.max,
                              mf.min = mf.min,
                              mf.max = mf.max,
                              hf.min = hf.min,
                              hf.max = hf.max,
                              uf.min = uf.min,
                              uf.max = uf.max,
                              verbose = verbose){

    # Apply high-pass filter
    if (hpf > 0) {
      wave <- fir(wave, from=hpf, to = NULL, bandpass = TRUE, output = "Wave", wl = 1024)
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }

    # Generate the amplitude spectrogram
    spectro_data <- spectrogram_binary(wave,
                                       freq.res = freq.res,
                                       cutoff = cutoff,
                                       plot = FALSE)

    if(verbose){
      cat("Amplitude range: ", round(spectro_data$raw_range), "\n")
      cat("dBFS range: ", spectro_data$dbfs_range, "\n")
      cat("Cutoff: ", cutoff, "\n")
    }


    binary_spectro <- spectro_data$spectrogram

    # Define frequency ranges
    samp_rate <- wave@samp.rate
    freq_bins <- seq(0, samp_rate / 2, length.out = nrow(binary_spectro))

    # Define frequency bands
    low_freq_indices <- which(freq_bins >= lf.min & freq_bins <= lf.max)
    mid_freq_indices <- which(freq_bins > mf.min & freq_bins <= mf.max)
    high_freq_indices <- which(freq_bins > hf.min & freq_bins <= hf.max)
    ultra_freq_indices <- which(freq_bins > uf.min & freq_bins <= uf.max)

    # Convert binary spectrogram to a data frame for plotting
    time <- seq(0, seewave::duration(wave), length.out = ncol(binary_spectro))
    freq <- freq_bins/1000

    binary_spectro_df <- data.frame(
      Var1 = rep(1:nrow(binary_spectro), ncol(binary_spectro)),
      Var2 = rep(1:ncol(binary_spectro), each = nrow(binary_spectro)),
      value = as.vector(binary_spectro)
    )
    binary_spectro_df$Frequency <- freq[binary_spectro_df$Var1]
    binary_spectro_df$Time <- time[binary_spectro_df$Var2]

    # Frequency bands in kHz
    lf.min <- lf.min / 1000
    lf.max <- lf.max / 1000
    mf.min <- mf.min / 1000
    mf.max <- mf.max / 1000
    hf.min <- hf.min / 1000
    hf.max <- hf.max / 1000
    uf.min <- uf.min / 1000
    uf.max <- uf.max / 1000
    
    if (verbose) {
    cat(paste("Frequency bands: \n",
              "Low Frequency: ", lf.min, "kHz", "- ", lf.max, " kHz \n",
              "Mid Frequency:  ", mf.min, "kHz", "- ", mf.max, " kHz \n",
              "High Frequency: ", hf.min, "kHz", "- ", hf.max, " kHz \n",
              "Ultra Frequency:", uf.min, "kHz", "-", uf.max, "kHz \n"))
    }

    # Calculate the proportion of cells above the cutoff in each frequency band
    LFC <- round(sum(binary_spectro[low_freq_indices, ]) / (length(low_freq_indices) * ncol(binary_spectro)),6)
    MFC <- round(sum(binary_spectro[mid_freq_indices, ]) / (length(mid_freq_indices) * ncol(binary_spectro)),6)
    HFC <- round(sum(binary_spectro[high_freq_indices, ]) / (length(high_freq_indices) * ncol(binary_spectro)),6)
    UFC <- round(sum(binary_spectro[ultra_freq_indices, ]) / (length(ultra_freq_indices) * ncol(binary_spectro)),6)

    # Results tibble
    summary <- tibble(
      index = c("lfc", "mfc", "hfc", "ufc"),
      value = c(LFC, MFC, HFC, UFC),
      par_cutoff = cutoff,
      par_freq_res = freq.res,
      par_lf_min = lf.min,
      par_lf_max = lf.max,
      par_mf_min = mf.min,
      par_mf_max = mf.max,
      par_hf_min = hf.min,
      par_hf_max = hf.max,
      par_uf_min = uf.min,
      par_uf_max = uf.max
    )

    if (spectrogram) {
      # Position of the labels on the x-axis
      xposition <- seewave::duration(wave) / 100

      if (ggplot) {
        # ggplot version
        p <- ggplot(binary_spectro_df, aes(Time, Frequency, fill = value)) +
          geom_tile() +
          scale_fill_manual(values = c("transparent", sound.color)) +
          labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          geom_hline(yintercept = lf.max, linetype = "dashed", colour = "red") +
          geom_hline(yintercept = mf.max, linetype = "dashed", colour = "red") +
          geom_hline(yintercept = hf.max, linetype = "dashed", colour = "red") +
          geom_hline(yintercept = uf.max, linetype = "dashed", colour = "red") +
          annotate("label", x = xposition, hjust = 0, y = (lf.min + lf.max) / 2,
                   label = paste("LFC:", LFC), fill = "white", alpha = 0.9) +
          annotate("label", x = xposition, hjust = 0, y = (mf.min + mf.max) / 2,
                   label = paste("MFC:", MFC), fill = "white", alpha = 0.9) +
          annotate("label", x = xposition, hjust = 0, y = (hf.min + hf.max) / 2,
                   label = paste("HFC:", HFC), fill = "white", alpha = 0.9) +
          annotate("label", x = xposition, hjust = 0, y = (uf.min + uf.max) / 2,
                   label = paste("UFC:", UFC), fill = "white", alpha = 0.9) +
          theme_light() +
          theme(legend.position = "none")
        
        print(summary)
        print(p)
        
        invisible(list(summary = summary, spectrogram = p))


      } else {
        # Base R version
        image(time, freq, t(binary_spectro),
              col = c("white", sound.color),
              xlab = "Time (s)", ylab = "Frequency (kHz)",
              main = plot.title, axes = FALSE, useRaster = TRUE)

        # Add axes with shorter ticks
        axis(1, at = pretty(time), labels = pretty(time), tck = -0.015, mgp = c(2, 0.5, 0))  # mgp adjusts title and labels closer
        axis(2, at = pretty(freq), labels = pretty(freq), tck = -0.015, mgp = c(2, 0.5, 0))  # mgp adjusts title and labels closer

        # if(grid){
        #   # Add grid lines
        #   abline(h = pretty(freq), col = "lightgray", lty = "dashed")
        #   abline(v = pretty(time), col = "lightgray", lty = "dashed")
        #
        # }

        # Add frequency band lines
        abline(h = c(lf.max, mf.max, hf.max, uf.max), col = "red", lty = "dashed")

        # Add annotations with white background box
        # LFC
        lf.y <- (lf.min + lf.max) / 2
        rect(xleft = xposition - 0.05, xright = xposition + 7, ybottom = lf.y - 0.5, ytop = lf.y + 0.4, col = "white", border = NA,density = NA)
        text(xposition, lf.y, labels = paste("LFC:", LFC), adj = c(0, 0.5), col = "black")

        # MFC
        mf.y <- (mf.min + mf.max) / 2
        rect(xleft = xposition - 0.05, xright = xposition + 7, ybottom = mf.y - 0.5, ytop = mf.y + 0.4, col = "white", border = NA)
        text(xposition, mf.y, labels = paste("MFC:", MFC), adj = c(0, 0.5), col = "black")

        # HFC
        hf.y <- (hf.min + hf.max) / 2
        rect(xleft = xposition - 0.05, xright = xposition + 7, ybottom = hf.y - 0.5, ytop = hf.y + 0.4, col = "white", border = NA)
        text(xposition, hf.y, labels = paste("HFC:", HFC), adj = c(0, 0.5), col = "black")

        # UFC
        uf.y <- (uf.min + uf.max) / 2
        rect(xleft = xposition - 0.05, xright = xposition + 7, ybottom = uf.y - 0.5, ytop = uf.y + 0.4, col = "white", border = NA)
        text(xposition, uf.y, labels = paste("UFC:", UFC), adj = c(0, 0.5), col = "black")

        # Add box around the plot
        box()
        
        return(list(summary = summary, spectrogram = NULL))

      }
    } else {
      return(list(summary = summary, spectrogram = NULL))
    }

  }
  
  if (wave@stereo) {
    
    # Calculate the index based on the stereo condition
    if (channel == 'each') {
      
      if (verbose) cat("Calculating Frequency Cover Index on 2 channels... \n")
      
      fc_left <- calculate_index(wave.left,
                                 hpf = 0,
                                 cutoff = cutoff,
                                 freq.res = freq.res,
                                 spectrogram = spectrogram,
                                 ggplot = ggplot,
                                 plot.title = paste0(plot.title, "_L"),
                                 sound.color = sound.color,
                                 lf.min = lf.min,
                                 lf.max = lf.max,
                                 mf.min = mf.min,
                                 mf.max = mf.max,
                                 hf.min = hf.min,
                                 hf.max = hf.max,
                                 uf.min = uf.min,
                                 uf.max = uf.max,
                                 verbose = verbose)
      fc_right <- calculate_index(wave.right,
                                  hpf = 0,
                                  cutoff = cutoff,
                                  freq.res = freq.res,
                                  spectrogram = spectrogram,
                                  ggplot = ggplot,
                                  plot.title = paste0(plot.title, "_R"),
                                  sound.color = sound.color,
                                  lf.min = lf.min,
                                  lf.max = lf.max,
                                  mf.min = mf.min,
                                  mf.max = mf.max,
                                  hf.min = hf.min,
                                  hf.max = hf.max,
                                  uf.min = uf.min,
                                  uf.max = uf.max,
                                  verbose = verbose)
      
      fc_global <- tibble(
        index = c("lfc", "mfc", "hfc", "ufc"),
        channel = "each")
      
      fc_global <- fc_global |>
        mutate(
          value_l = fc_left$summary$value[match(index, fc_left$summary$index)],
          value_r = fc_right$summary$value[match(index, fc_right$summary$index)],
          value_avg = round((
            fc_left$summary$value[match(index, fc_left$summary$index)] + 
              fc_right$summary$value[match(index, fc_right$summary$index)]
          )/2, 3),
          par_cutoff = cutoff,
          par_freq_res = freq.res,
          par_lf_min = lf.min,
          par_lf_max = lf.max,
          par_mf_min = mf.min,
          par_mf_max = mf.max,
          par_hf_min = hf.min,
          par_hf_max = hf.max
        )

      print(fc_global)
    
    if (spectrogram) {
      invisible(list(summary = fc_global, 
                     spectrogram_l = fc_left$spectrogram,
                     spectrogram_r = fc_right$spectrogam))
      
    } else {
      invisible(fc_global)
    }
    

  }  else if (channel %in% c("left", "right")) {
    # Use one channel from the stereo file or a mix

    if (verbose) {
      if(channel == 'left'){
        cat("Calculating Frequency Cover on the left channel... \n")
      }else if(channel == 'right'){
        cat("Calculating Frequency Cover on the right channel... \n")
      }
    }
 
    results <- calculate_index(wave,
                               hpf = 0,
                               cutoff = cutoff,
                               freq.res = freq.res,
                               spectrogram = spectrogram,
                               ggplot = ggplot,
                               plot.title = plot.title,
                               sound.color = sound.color,
                               lf.min = lf.min,
                               lf.max = lf.max,
                               mf.min = mf.min,
                               mf.max = mf.max,
                               hf.min = hf.min,
                               hf.max = hf.max,
                               uf.min = uf.min,
                               uf.max = uf.max)
    
    results$summary <- results$summary |>
      mutate(channel = channel) |>
      relocate(channel, .after = index)

    print(results$summary)
    
    if (spectrogram) {
      invisible(list(summary = results$summary, 
                     spectrogram = results$spectrogram))
    } else {
    invisible(results$summary)
    }
  } else if(channel == 'mix'){
    cat("Calculating Frequency Cover on a mix of the two channels... \n")
    results <- calculate_index(wave,
                               hpf = 0,
                               cutoff = cutoff,
                               freq.res = freq.res,
                               spectrogram = spectrogram,
                               ggplot = ggplot,
                               plot.title = plot.title,
                               sound.color = sound.color,
                               lf.min = lf.min,
                               lf.max = lf.max,
                               mf.min = mf.min,
                               mf.max = mf.max,
                               hf.min = hf.min,
                               hf.max = hf.max,
                               uf.min = uf.min,
                               uf.max = uf.max)
    
    results$summary <- results$summary |>
      mutate(channel = "mix") |>
      relocate(channel, .after = index)
    
    print(results$summary)
    
    if (spectrogram) {
      invisible(list(summary = results$summary, 
                     spectrogram = results$spectrogram))
    } else {
      invisible(results$summary)
    }
    
    
  }
    
  } else {
    if (verbose) {
        cat("Calculating Frequency Cover on a mono file... \n")
    }
      # Calulate NBAI
      results <- calculate_index(wave,
                                   hpf = 0,
                                   cutoff = cutoff,
                                   freq.res = freq.res,
                                   spectrogram = spectrogram,
                                   ggplot = ggplot,
                                   plot.title = plot.title,
                                   sound.color = sound.color,
                                   lf.min = lf.min,
                                   lf.max = lf.max,
                                   mf.min = mf.min,
                                   mf.max = mf.max,
                                   hf.min = hf.min,
                                   hf.max = hf.max,
                                   uf.min = uf.min,
                                   uf.max = uf.max,
                                   verbose = verbose)
      
      results$summary <- results$summary |>
        mutate(channel = "mono") |>
        relocate(channel, .after = index)
        # add_column(channel = "mono", .after = index)
      
      print(results$summary)

      if (spectrogram) {
        invisible(list(summary = results$summary, 
                       spectrogram = results$spectrogram))
      } else {
        invisible(results$summary)
      }
  }
}


