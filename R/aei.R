#' Acoustic Evenness Index
#'
#' @description
#' Acoustic Evenness Index from Villanueva-Rivera \emph{et al.} 2011.
#' The AEI is calculated by dividing the spectrogram into frequency bands (default 10),
#' taking the proportion of energy in each band above an energy threshold (default -50 dBFS),
#' and then calculating the Gini Coefficient from those proportions.
#' The new version allows the user to choose between different ways to compute
#' the proportions before calculating the Gini, among other new parameter options (see Details)
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param w.len the window length to compute the spectrogram (i.e., FFT window size).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram.
#' @param max.freq maximum frequency to compute the spectrogram.
#' @param n.bands number of bands to split the spectrogram.
#' @param cutoff dB threshold to calculate energy proportions (if norm.spec = FALSE, set to 5 or above).
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param props logical; if set to TRUE, the function stores the energy proportion values for each frequency band and channel. Default = TRUE.
#' @param prop.den numeric; indicates how the energy proportion is calculated.
#'
#' @return A tibble (data frame) with the AEI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @importFrom tuneR readWave
#' @import seewave
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @examples aei(tropicalsound)
aei <- function(wave,
                w.len = 512,
                w.fun = "hanning",
                min.freq = 0,
                max.freq = 10000,
                n.bands = 10,
                cutoff = 5,
                norm.spec = FALSE,
                noise.red = 2,
                rm.offset = TRUE,
                props = TRUE,
                prop.den = 1){




  # Store the frequency step (band "height") # NEW 09/25/2023 Francisco Rivas
  freq_step <- (max.freq - min.freq)/n.bands


  cutoff <- as.numeric(cutoff)

  #test arguments
  if (is.numeric(as.numeric(max.freq))){
    max.freq <- as.numeric(max.freq)
  } else{
    stop(" max.freq is not a number.")
  }

  if (is.numeric(as.numeric(cutoff))){
    cutoff <- as.numeric(cutoff)
  } else{
    stop(" cutoff is not a number.")
  }

  if (is.numeric(as.numeric(freq_step))){
    freq_step <- as.numeric(freq_step)
  } else{
    stop(" freq_step is not a number.")
  }


  #function that gets the proportion of values over a db
  # value in a specific band of frequencies.
  # Frequency is in Hz
  # getscore<-function(spectrum, minf, maxf, db, freq_row){
  #   miny<-round((minf)/freq_row)
  #   maxy<-round((maxf)/freq_row)
  #
  #   subA=spectrum[miny:maxy,]
  #
  #   index1<-length(subA[subA>db])/length(subA)
  #
  #   return(index1)
  # }
  # Function that gets the proportion of values higher than the
  # db threshold in a specific frequency band. The frequencies are in Hz
  getscore <- function(spectrum, minf, maxf, db, freq_row){
    miny<-round((minf)/freq_row) # the minimum frequency of the frequency band
    maxy<-round((maxf)/freq_row) # the maximum frequency of the frequency band

    subA = spectrum[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)


    minspec <- round(0/freq_row) # lower end of the spectrogram defined by min.freq
    maxspec <- round(max.freq/freq_row) # upper end of the spectrogram defined by max.freq
    speclims <- spectrum[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq

    # Calculate the proportion of cells in f that are higher than the dB threshold
    if(prop.den == 1){ #original AEI proportion calculation (within frequency band)
      index1 <- length(subA[subA>db]) / length(subA)

    }else if(prop.den == 2){
      # Alternative 2: over the user-defined spectrogram range
      # (cells above the energy threshold across the spectrogram)
      index1 <- length(subA[subA>db]) / length(speclims[speclims>db])

    }else if(prop.den == 3){
      # Alternative 3: over the whole spectrogram
      # (over all the pixels above the threshold up to the Nyquist frequency)
      index1 <- length(subA[subA>db]) / length(spectrum[spectrum>db])
    }
    return(index1)
  }


  # Save the denominator used in the proportion calculation
  if(prop.den == 1){
    prop.denom <- "within band"
  }else if(prop.den == 2){
    prop.denom <- "max.freq"
  }else if(prop.den == 3){
    prop.denom <- "nyquist"
  }

  # Add information about noise reduction procedure
  if(noise.red == 1){
    noise <- "rows"
  } else if(noise.red == 2){
    noise <- "columns"
  } else {
    noise <- "none"
  }

  #Some general values
  #Get sampling rate
  samplingrate <- wave@samp.rate
  duration <- length(wave@left)/samplingrate

  #Get Nyquist frequency in Hz
  nyquist_freq <- samplingrate/2

  #window length for the spectro and spec functions
  #to keep each row every 10Hz
  #Frequencies and seconds covered by each
  # freq_per_row = 10
  # w.len = samplingrate/freq_per_row
  # Calculate frequency resolution (i.e., frequency bin width)
  freq_per_row = samplingrate/w.len

  # Adding 1 if w.len is an odd number (new behavior in seewave)
  # fix by JSueur
  if(w.len%%2 == 1) {w.len <- w.len+1}


  #Stereo file
  if (wave@stereo == TRUE) {

    cat("Calculating AEI on a stereo file... \n")

    left<-channel(wave, which = c("left"))
    right<-channel(wave, which = c("right"))
    rm(wave)

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- rm.offsetset(left)
      right <- rm.offsetset(right)
    }

    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    }


    # #matrix of values
    # specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE)$amp
    # specA_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE)$amp
    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len,
                              wn = w.fun, noise.reduction = noise.red,
                              plot = FALSE)$amp
        specA_right <- spectro(right, f = samplingrate, wl = w.len,
                               wn = w.fun, noise.reduction = noise.red,
                               plot = FALSE)$amp
      }else if (noise.red == 0) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len,
                              wn = w.fun, noise.reduction = NULL,
                              plot = FALSE)$amp
        specA_right <- spectro(right, f = samplingrate, wl = w.len,
                               wn = w.fun, noise.reduction = NULL,
                               plot = FALSE)$amp
      }


      rm(left, right)


    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")

      # if(!is.null(noise.red)){
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noise.reduction = noise.red)$amp
        specA_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
                               norm=FALSE,dB=NULL,unit="power",
                               noise.reduction = noise.red)$amp
      }else if(noise.red == 0){
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power")$amp
        specA_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
                               norm=FALSE,dB=NULL,unit="power")$amp
      }


      # Transform to decibels
      specA_left <- 10*log10(specA_left^2)
      specA_right <- 10*log10(specA_right^2)

      rm(left, right)

    }



    # rm(left,right)

    if (max.freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max.freq <- nyquist_freq
    }

    # Set the frequency bands (Fran's comment)
    # Freq<-seq(from = 0, to = max.freq - freq_step, by = freq_step)
    Freq <- seq(from = min.freq, to = max.freq - freq_step, by = freq_step)


    #LEFT CHANNEL
    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), cutoff, freq_per_row)
    }


    left_vals=Score

    #RIGHT CHANNEL
    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_right, Freq[j], (Freq[j] + freq_step), cutoff, freq_per_row)
    }

    right_vals = Score

    # 		cat(" ==============================================\n")
    # 		cat(paste(" Results (with a dB threshold of ", cutoff, ")\n\n", sep=""))

    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))

    # 		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
    # 		cat("Frequency range (Hz), left channel proportion, right channel proportion\n")
    # for (j in seq(length(Freq), 1, by = -1)) {
    #   # 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), ",", round(right_vals[j],6), "\n", sep=""))
    #   left_bandvals_return[j] = round(left_vals[j], 6)
    #   right_bandvals_return[j] = round(right_vals[j], 6)
    #   left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
    #   right_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
    # }

    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      right_bandvals_return[j] = round(right_vals[j], 6)
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
      right_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }

    # 		cat("\n Plot of proportions in each band: \n\n")
    # 		cat("  Left channel\n")
    # 		cat("   Freq. range (Hz) |--------------------|\n")

    #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
    # for (j in seq(length(Freq), 1, by = -1)) {
    #   this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
    #   this_row_size <- nchar(this_row_name)
    #   this_row_space <- 17 - this_row_size
    #
    #   this_row_spaces = ""
    #
    #   for (f in seq(1,this_row_space,by = 1)) {
    #     this_row_spaces = paste(this_row_spaces, " ", sep = "")
    #   }
    #
    #   # 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
    #   # 			temp_val=round(left_vals[j],2)*20
    #   # 			if (temp_val>0){
    #   # 				for (i in 1:temp_val) {
    #   # 					cat("*")
    #   # 				}
    #   # 			}
    #   # 			cat("\n")
    #   # 			rm(temp_val)
    # }

    # 		cat("\n  Right channel\n")
    # 		cat("   Freq. range (Hz) |--------------------|\n")

    #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
    # for (j in seq(length(Freq), 1, by = -1)) {
    #   this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
    #   this_row_size <- nchar(this_row_name)
    #   this_row_space <- 17 - this_row_size
    #
    #   this_row_spaces = ""
    #
    #   for (f in seq(1,this_row_space, by = 1)) {
    #     this_row_spaces = paste(this_row_spaces, " ", sep="")
    #   }
    #
    #   # 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
    #
    #   # 			temp_val=round(right_vals[j],2)*20
    #   # 			if (temp_val>0){
    #   # 				for (i in 1:temp_val) {
    #   # 					cat("*")
    #   # 				}
    #   # 			}
    #   # 			cat("\n")
    #   # 			rm(temp_val)
    # }

    #cat("\n")
    # cat("Acoustic Evenness Index: \n")
    # cat(paste("   Left channel: ", round(Gini(left_vals), 6), "\n", sep=""))
    # cat(paste("   Right channel: ", round(Gini(right_vals), 6), "\n\n", sep=""))
    left_gini_return = round(Gini(left_vals), 6)
    right_gini_return = round(Gini(right_vals), 6)

    aeiOutputStereo <- tibble(value_l = left_gini_return,
                              value_r = right_gini_return)

    aeiOutputStereo <- aeiOutputStereo %>%
      add_column(value_avg = ((aeiOutputStereo$value_l+aeiOutputStereo$value_r)/2), .after = "value_r")

    # Add metadata columns
    aeiOutputStereo <- aeiOutputStereo %>%
      add_column(w.len = w.len,
                 w.fun = w.fun,
                 dbth = cutoff,
                 minf = min.freq,
                 maxf = max.freq,
                 n.bands = n.bands,
                 norm = norm.spec,
                 noise.red = noise,
                 rm.offset = rm.offset,
                 prop.den = prop.denom,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "stereo")

    aeiOutputStereo <- aeiOutputStereo %>%
      add_column(index = "aei", .before = "value_l")

    if(props == TRUE){
      # Add data on proportions for each frequency band per channel:
      proportions = tibble(band_kHz = left_bandrange_return,
                           left = left_bandvals_return,
                           right = right_bandvals_return)
      # Wider format
      proportions <- pivot_wider(proportions, names_from = "band_kHz",
                                 values_from = c("left", "right"))

      # Append to the original data frame
      aeiOutputStereo <- bind_cols(aeiOutputStereo, proportions)

      cat("Reporting AEI for 2 channels, metadata and energy proportions per frequency band: \n")

      # aeiOutputStereo <- aeiOutputStereo %>%
      #   add_column(index = "aei", .before = "value_l")

      # return(aeiOutputStereo)

    } else {

      cat("Reporting AEI for 2 channels and metadata:\n")

    }

    # aeiOutputStereo <- aeiOutputStereo %>%
    #   add_column(index = "adi", .before = "value_l")

    return(aeiOutputStereo)




  } else
  {
    # MONO
    cat("Calculating AEI on a mono file... \n")


    left<-channel(wave, which = c("left"))
    rm(wave)

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- rm.offsetset(left)
    }


    if(noise.red == 1){
      cat("Applying noise reduction filter (subtract median amplitude) to each row...\n")
    }else if(noise.red == 2){
      cat("Applying a noise reduction filter (subtract median amplitude) to each column...\n")
    }

    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                              noise.reduction = noise.red)$amp
      }else if (noise.red == 3) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE)$amp
      }
      rm(left)

    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noise.reduction = noise.red)$amp
      }else if (noise.red == 3) {
        specA_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noise.reduction = noise.red)$amp
      }
      # Transform to decibels
      specA_left <- 10*log10(specA_left^2)

      rm(left)

    }

    if (max.freq>nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max.freq <- nyquist_freq
    }

    Freq <- seq(from = min.freq, to = max.freq - freq_step, by = freq_step)

    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), cutoff, freq_per_row)
    }

    left_vals = Score


    left_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))

    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
    }

    # 		cat("\n Plot of proportions in each band: \n\n")
    # 		cat("   Freq. range (Hz) |--------------------|\n")

    #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
    # for (j in seq(length(Freq), 1, by = -1)) {
    #   this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
    #   this_row_size <- nchar(this_row_name)
    #   this_row_space <- 17 - this_row_size
    #
    #   this_row_spaces = ""
    #
    #   for (f in seq(1, this_row_space, by = 1)) {
    #     this_row_spaces = paste(this_row_spaces, " ", sep = "")
    #   }
    #
    #   # 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
    #   # 			temp_val=round(left_vals[j],2)*20
    #   # 			if (temp_val>0){
    #   # 				for (i in 1:temp_val) {
    #   # 					cat("*")
    #   # 				}
    #   # 			}
    #   # 			cat("\n")
    #   # 			rm(temp_val)
    # }

    #cat("\n")
    # cat("  Acoustic Evenness Index: ")
    # cat(paste(round(Gini(left_vals), 6), "\n", sep = ""))
    left_gini_return = round(Gini(left_vals), 6)
    # right_gini_return = NA

    aeiOutputMono <- tibble(value = left_gini_return)

    # Add metadata columns
    aeiOutputMono <- aeiOutputMono %>%
      add_column(w.len = w.len,
                 w.fun = w.fun,
                 dbth = cutoff,
                 minf = min.freq,
                 maxf = max.freq,
                 n.bands = n.bands,
                 norm = norm.spec,
                 noise.red = noise,
                 rm.offset = rm.offset,
                 prop.den = prop.denom,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "mono")

    aeiOutputMono <- aeiOutputMono %>%
      add_column(index = "aei", .before = "value")

    if(props == TRUE){
      # Add data on proportions for each frequency band per channel:
      proportions = tibble(band_kHz = left_bandrange_return,
                           left = left_bandvals_return)
      # Wider format
      proportions <- pivot_wider(proportions, names_from = "band_kHz",
                                 values_from = c("left"))

      # Append to the original data frame
      aeiOutputMono <- bind_cols(aeiOutputMono, proportions)

      cat("Reporting AEI for 1 channel, metadata and energy proportions per frequency band. \n")

      # aeiOutputMono <- aeiOutputMono %>%
      #   add_column(index = "aei", .before = "value")

      # return(aeiOutputMono)

    } else {

      cat("Reporting AEI for 1 channel and metadata.\n")

    }
    #
    # aeiOutputMono <- aeiOutputMono %>%
    #   add_column(index = "aei", .before = "value")

    # return(aeiOutputMono)



  }
  return(aeiOutputMono)
}
