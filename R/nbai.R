#' Narrow-Band Activity Index
#' @description
#' This index describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary acoustic indices).
#'
#' @param wave A Wave object
#' @param channel Character. Options are "left", "right", "mix", or "each". The "mix" option creates a temporary mono file from a mix of both channels before calculating the index. Default is "left".
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.
#' @param activity.cutoff Numeric. Cutoff percent activity. Only the frequency bands active equal or above this percentage will be considered as "active" in the active band statistics.
#' @param spectrogram Logical. Whether to plot the binary spectrogram used for the analysis. Allowed only when Wave is mono or when one channel is selected from a stereo file.
#' @param dark.plot Logical. If true (default) a the binary spectrogram will have a black background.
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#'
#' @return A list containing: 1) A binary spectrogram (if mono), 2) tibble with the Narrow-Band Activity Index (NBI) summary statistics, and 3) a tibble with NBI spectral, which number of rows equals the number of frequency bins in the analysis.
#' @export
#'
#' @import tuneR
#' @import seewave
#' @import patchwork
#' @import ggplot2
#'
#' @examples nbai(wave, channel = 'left', plot = TRUE, verbose = TRUE)

nbai <- function(wave,
                 channel = 'left',
                 hpf = 250,
                 freq.res = 50,
                 cutoff = -60,
                 activity.cutoff = 10,
                 spectrogram = FALSE,
                 dark.plot = TRUE,
                 plot.title = deparse(substitute(wave)),
                 verbose = TRUE) {
  
  
  if (wave@stereo) {
    if (channel == 'each') {
     cat("Calculating NBAI for each channel... \n")
    } else if (channel == 'left') {
      cat("Calculating NBAI for left channel... \n")
    } else if (channel == 'right') {
      cat("Calculating NBAI for right channel... \n")
    } else if (channel == 'mix'){
      cat("Calculating NBAI for a mix of 2 channels... \n")
    } else {
      stop("Invalid channel selected.")
    }
  } else {
    cat("Calculating NBAI for a mono file... \n")
  }
  
  
  nbai_mono <- function(wave, 
                        hpf = 0, 
                        freq.res = 50, 
                        cutoff = -60, 
                        activity.cutoff = 10, 
                        spectrogram = TRUE, 
                        plot.activity = FALSE, 
                        plot.title = plot.title, 
                        dark.plot = TRUE){
    
    
    # Calculate window length according to the sampling rate for a fixed frequency resolution
    wl = wave@samp.rate / freq.res
    
    # Force it to be even
    if (wl %% 2 == 1) { wl <- wl + 1 }
    
    # Remove DC offset
    wave <- rmoffset(wave, output = "Wave")
    
    # Apply high-pass filter
    if (hpf > 0) {
      wave <- fir(wave, from=hpf, to = NULL, bandpass = TRUE, output = "Wave", wl = 1024)
    } else if (hpf < 0) {
      stop("HPF should be either 0 or a positive number (in Hertz) \n")
    }
    
    binary_spectrogram <- spectrogram_binary(wave,
                                             freq.res = freq.res,
                                             cutoff = cutoff)$spectrogram
    
    # Calculate the percent of active cells (=1) for each frequency bin
    activity_percent <- apply(binary_spectrogram, 1, function(x) {
      sum(x) / length(x) * 100
    })
    
    
    # Count the number of bins where the percent of active cells is >= threshold
    bins_above_threshold <- sum(activity_percent >= activity.cutoff)
    
    # Calculate the total number of frequency bins
    total_bins <- nrow(binary_spectrogram)
    
    # Calculate the percent of bins with activity above user-defined threshold
    percent_bins_above_threshold <- (bins_above_threshold / total_bins) * 100
    
    
    #### Calculate Entropy of activity classes #####
    # Define the 20 class boundaries (5% intervals)
    class_breaks <- seq(0, 100, by = 5)
    
    
    # Classify each value in activity_percent into the N classes
    # Use the right-closed interval to include values like 5.5 into the correct class
    activity_classes <- cut(activity_percent,
                            breaks = class_breaks,
                            include.lowest = TRUE,
                            right = FALSE)
    
    # Count the frequency of values in each class
    class_counts <- table(activity_classes)
    
    # Convert the counts to proportions
    class_proportions <- class_counts / sum(class_counts)
    
    # Calculate Shannon entropy
    shannon_entropy <- -sum(class_proportions * log(class_proportions+0.000001))
    
    #####################
    
    
    # Report results in a tibble
    result <- tibble(
      index = 'nbai',
      value = round(mean(activity_percent), 2),
      w.value = round(value * (shannon_entropy + 1), 2),
      nab = bins_above_threshold,
      pab = round(percent_bins_above_threshold, 2),
      ent = round(shannon_entropy, 2),
      par.hpf = hpf, 
      par.freq.res = freq.res, 
      par.cutoff = cutoff, 
      par.act.cutoff = activity.cutoff, 
    )
    

    if (spectrogram) {
      # Retrieve binary spectrogram
      # binary_spectrogram <- nbai_global$spectrogram
      total_duration <- seewave::duration(wave)
      # Create the time axis values
      time_values <- seq(0, total_duration, length.out = ncol(binary_spectrogram))
      # Create the frequency axis values
      freq_values <- seq(0, by = wave@samp.rate / 2 / nrow(binary_spectrogram), 
                         length.out = nrow(binary_spectrogram)) / 1000  # Convert to kHz
      
      # Create a color gradient based on percent activity
      colors <- scales::col_numeric(palette = c("#2c7bb6","#00a6ca","#00ccbc","#90eb9d",
                                                "#ffff8c", "#f9d057", "#f29e2e","#e76818","#d7191c"),
                                    domain = c(0, 100))(activity_percent)
      
      
      # Create a dataframe for plotting
      plot_data <- as.data.frame(binary_spectrogram)
      colnames(plot_data) <- time_values
      plot_data <- cbind(Frequency = freq_values, plot_data)
      plot_data <- reshape2::melt(plot_data, id.vars = "Frequency", 
                                  variable.name = "Time", value.name = "Binary")
      plot_data$Time <- as.numeric(as.character(plot_data$Time))
      
      # Add activity percent color coding
      
      plot_data$Color <- rep(colors, times = ncol(binary_spectrogram))
      plot_data$Color[plot_data$Binary == 0] <- NA
      
      p <- ggplot(plot_data, aes(x = Time, y = Frequency)) +
        geom_tile(aes(fill = Color), color = NA) +
        scale_fill_identity(na.value = "transparent") +
        labs(x = "Time (s)", y = "Frequency (kHz)", title = plot.title) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_y_continuous(expand = c(0,0)) +
        
        # Add annotation box
        annotate("label", x = max(plot_data$Time) * 0.8, y = max(plot_data$Frequency) * 0.9,
                 label = paste("NBAI:", result$value,
                               "\nNAB: ", result$nab,
                               "\nPAB: ", result$pab,
                               "\nEntropy: ", result$ent,
                               "\nwNBAI: ", result$w.value,
                               "\nEnergy cutoff: ", cutoff, "dB",
                               "\nActivity cutoff: ", activity.cutoff, "%"
                 ),
                 hjust = 0, vjust = 1, size = 3.5, color = "black",
                 label.size = 0.5, label.padding = unit(0.5, "lines"),
                 fontface = "italic", fill = "white")
      
      # # Add a manual gradient color scale legend
      # p <- p + guides(fill = guide_colourbar(barwidth = 10, barheight = 1, title = "Activity (%)",
      #                                        title.position = "top", label.position = "bottom"))
      
      
      # Apply the dark theme logic conditionally outside of the annotate() function
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
      plot(p)
      
      return(list(summary = result, spectral = activity_percent, spectrogram = p))
      
      
    } else {
      
      return(list(summary = result, spectral = activity_percent))
      
    }
    
    # if (plot.activity) {
    #   samp_rate <- wave@samp.rate
    #   nyquist_freq <- samp_rate / 2
    #   freq_values <- seq(0, nyquist_freq, length.out = length(activity_percent))/1000
    #   
    #   # Create a data frame with frequency values and spectral data
    #   df <- data.frame(
    #     frequency = freq_values,
    #     spectral = activity_percent
    #   )
    #   
    #   # Plot the data using ggplot2 and viridis
    #   p <-  ggplot(df, aes(x = freq_values, y = spectral)) +
    #     geom_area() +
    #     labs(title = "Narrow-band Activity", x = "Frequency (kHz)", y = "Percent Activity") +
    #     theme_minimal() +
    #     theme(plot.title = element_text(hjust = 0.5)) +
    #     coord_flip() +
    #     scale_x_continuous(expand = c(0,0)) +
    #     scale_y_continuous(expand = c(0,0), limits = c(0,100))
    #   
    #   plot(p)
    #   
    # }
    # 
    # return(list(summary = result, spectral = activity_percent, spectrogram = p))
  }
  
  
  if (wave@stereo) {
  # Calculate the index based on the stereo condition
  if (channel == 'each') {
    wave.left <- tuneR::channel(wave, 'left')
    wave.right <- tuneR::channel(wave, 'right')
    
    if (verbose) cat("Calculating Narrow Band Activity Index on 2 channels... \n")
    
    nbai_left <- nbai_mono(wave.left, hpf = hpf, freq.res = freq.res, cutoff = cutoff, 
                           activity.cutoff = activity.cutoff, spectrogram = spectrogram, 
                           plot.title = paste(plot.title, "NBAI - Left Channel"))
    nbai_right <- nbai_mono(wave.right, hpf = hpf, freq.res = freq.res, cutoff = cutoff, 
                            activity.cutoff = activity.cutoff, spectrogram = spectrogram,
                            plot.title = paste(plot.title, "NBAI - Right Channel"))
    
    nbai_global <- tibble::tibble(
      index = "nbai",
      channel = 'each',
      value_l = nbai_left$summary$value,
      value_r = nbai_right$summary$value,
      value_avg = round((nbai_left$summary$value + nbai_right$summary$value) / 2, 1),
      w.value_l = nbai_left$summary$w.value,
      w.value_r = nbai_right$summary$w.value,
      w.value_avg = round((nbai_left$summary$w.value + nbai_right$summary$w.value)/2, 1),
      nab_l = nbai_left$summary$nab,
      nab_r = nbai_right$summary$nab,
      nab_avg = round((nbai_left$summary$nab + nbai_right$summary$nab)/2, 1),
      pab_l = nbai_left$summary$pab,
      pab_r = nbai_right$summary$pab,
      pab_avg = round((nbai_left$summary$pab + nbai_right$summary$pab)/2, 1),
      ent_l = nbai_left$summary$ent,
      ent_r = nbai_right$summary$ent,
      ent_avg = round((nbai_left$summary$ent + nbai_right$summary$ent)/2, 1)
    )
    
    if(verbose){
      
      cat(" value: Narrow-band Activity Index (NBAI)\n",
          "w.value: Weighted NBAI \n",
          "nab: Number of active (>", activity.cutoff, "%) frequency bins \n",
          "pab: Percent of active (>", activity.cutoff, "%) bins across the spectrogram \n",
          "ent: Shannon's Entropy\n")
      print(nbai_global)
      
    }
    
    if(spectrogram){
      # if(spectrogram)
      invisible(
        list(summary = nbai_global,
             spectral_left = nbai_left$spectral,
             spectral_right = nbai_right$spectral,
             spectrogram_left = nbai_left$spectrogram,
             spectrogram_right = nbai_right$spectrogram)
      )
    }else{
      invisible(
        list(summary = nbai_global,
             spectral_left = nbai_left$spectral,
             spectral_right = nbai_right$spectral)
      )
    }
    } else if (channel == 'left') {
    
    wave <- tuneR::channel(wave, 'left')
    
    
    nbai_left <- nbai_mono(wave, 
                             hpf = hpf, 
                             freq.res = freq.res, 
                             cutoff = cutoff, 
                             activity.cutoff = activity.cutoff, 
                             spectrogram = spectrogram, 
                             plot.title = plot.title, 
                             dark.plot = dark.plot)
    
    nbai_left$summary <- nbai_left$summary |> 
      mutate(channel = "left") |> 
      relocate(channel, .after = index)
    
    if(verbose){
      
      cat(" value: Narrow-band Activity Index (NBAI)\n",
          "w.value: Weighted NBAI \n",
          "nab: Number of active (>", activity.cutoff, "%) frequency bins \n",
          "pab: Percent of active (>", activity.cutoff, "%) bins across the spectrogram \n",
          "ent: Shannon's Entropy\n")
      print(nbai_left$summary)
      
    }
    
    
    if(spectrogram){
      invisible(list(summary = nbai_left$summary,
                     spectral = nbai_left$spectral,
                     spectrogram = nbai_left$spectrogram))
    } else {
      invisible(list(summary = nbai_left$summary,
                     spectral = nbai_left$spectral))
      
    } 
    
    } else if (channel == 'right') {
      
      wave <- tuneR::channel(wave, 'right')
      
      
      nbai_right <- nbai_mono(wave, 
                             hpf = hpf, 
                             freq.res = freq.res, 
                             cutoff = cutoff, 
                             activity.cutoff = activity.cutoff, 
                             spectrogram = spectrogram, 
                             plot.title = plot.title, 
                             dark.plot = dark.plot)
      
      nbai_right$summary <- nbai_right$summary |> 
        mutate(channel = "right") |> 
        relocate(channel, .after = index)
      
      
      if(verbose){
        cat(" value: Narrow-band Activity Index (NBAI)\n",
            "w.value: Weighted NBAI \n",
            "nab: Number of active (>", activity.cutoff, "%) frequency bins \n",
            "pab: Percent of active (>", activity.cutoff, "%) bins across the spectrogram \n",
            "ent: Shannon's Entropy\n")
        print(nbai_right$summary)
        
      }
      
      if(spectrogram){
        invisible(list(summary = nbai_right$summary,
                       spectral = nbai_right$spectral,
                       spectrogram = nbai_right$spectrogram))
      } else {
        invisible(list(summary = nbai_right$summary,
                       spectral = nbai_right$spectral))
        
      }
      
    
    } else if (channel == 'mix') {
      
      wave <- tuneR::mono(wave, "both")
      
      
      nbai_mix <- nbai_mono(wave, 
                              hpf = hpf, 
                              freq.res = freq.res, 
                              cutoff = cutoff, 
                              activity.cutoff = activity.cutoff, 
                              spectrogram = spectrogram, 
                              plot.title = plot.title, 
                              dark.plot = dark.plot)
      
      nbai_mix$summary <- nbai_mix$summary |> 
        mutate(channel = "mix") |> 
        relocate(channel, .after = index)
      
      if(verbose){
        cat(" value: Narrow-band Activity Index (NBAI)\n",
            "w.value: Weighted NBAI \n",
            "nab: Number of active (>", activity.cutoff, "%) frequency bins \n",
            "pab: Percent of active (>", activity.cutoff, "%) bins across the spectrogram \n",
            "ent: Shannon's Entropy\n")
        print(nbai_mix$summary)
        
      }
      
      if(spectrogram){
        invisible(list(summary = nbai_mix$summary,
                       spectral = nbai_mix$spectral,
                       spectrogram = nbai_mix$spectrogram))
      } else {
        invisible(list(summary = nbai_mix$summary,
                       spectral = nbai_mix$spectral))
        
      }
    }
  } else {
    # MONO RECORDING
    nbai_monofile <- nbai_mono(wave, 
                             hpf = hpf, 
                             freq.res = freq.res, 
                             cutoff = cutoff, 
                             activity.cutoff = activity.cutoff, 
                             spectrogram = spectrogram, 
                             plot.title = plot.title, 
                             dark.plot = dark.plot)
    
    nbai_monofile$summary <- nbai_monofile$summary |> 
      mutate(channel = "mono") |> 
      relocate(channel, .after = index)
    
    
    if(verbose){
      cat(" value: Narrow-band Activity Index (NBAI)\n",
          "w.value: Weighted NBAI \n",
          "nab: Number of active (>", activity.cutoff, "%) frequency bins \n",
          "pab: Percent of active (>", activity.cutoff, "%) bins across the spectrogram \n",
          "ent: Shannon's Entropy\n")
      print(nbai_monofile$summary)
      
    }
    
    if(spectrogram){
      invisible(list(summary = nbai_monofile$summary,
                     spectral = nbai_monofile$spectral,
                     spectrogram = nbai_monofile$spectrogram))
    } else {
      invisible(list(summary = nbai_monofile$summary,
                     spectral = nbai_monofile$spectral))
      
    }
    
    
  }
    
  
}

