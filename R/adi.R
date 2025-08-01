#' Acoustic Diversity Index
#' @description
#' Acoustic Diversity Index from Villanueva-Rivera et al. 2011.
#' The ADI is calculated by dividing the spectrogram into frequency bands (default 10),
#' taking the proportion of energy in each band above an energy threshold (default -50 dBFS),
#' and then calculating the Shannon's Diversity from those proportions.
#' The new version allows the user to choose between different ways to compute
#' the proportions before calculating the Shannon, among other new parameter options (see Details).
#'
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, 
#' therefore, the window length to be used (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", 
#' "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram.
#' @param max.freq maximum frequency to compute the spectrogram.
#' @param n.bands number of bands to split the spectrogram.
#' @param cutoff dB threshold to calculate energy proportions.
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended 
#' because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied 
#' to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each 
#' column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param props logical; if set to TRUE, the function stores the energy proportion values for each frequency band 
#' and channel. Default = TRUE.
#' @param prop.den numeric; indicates how the energy proportion is calculated.
#' @param db.fs logical; if TRUE, the amplitude scale is expressed as decibels Full Scale (dBFS). Only used when norm = FALSE. 
#' @param use.vegan logical; if TRUE, the diversity() function from the vegan package is called to compute Shannon's entropy. This is 
#' only intended for when the user wants to obtain a result equivalent to the original soundecology package. Default = FALSE.
#' @return A tibble (data frame) with the ADI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#'
#' @importFrom tuneR readWave channel
#' @importFrom seewave spectro rmoffset
#' @importFrom tibble tibble add_column
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom dplyr bind_cols
#'
#' @details
#' Options for the 'prop.den' parameter: 1 = The original calculation from the 
#' "soundecology" package is applied. 
#' The denominator of the proportion equals to all the cells in the same 
#' frequency band. 2 = A "true Shannon" proportion is calculated, 
#' where the "whole community" equals to the cells above the decibel threshold 
#' across the spectrogram (up to 'max_freq').
#' Another important update is that now the spectrogram is not normalized by 
#' default, which made recordings#' with different signal-to-noise ratio not comparable.
#' @examples
#' data(tropicalsound)
#' adi(tropicalsound)
adi <- function(wave,
                 freq.res = 50,
                 win.fun = "hanning",
                 min.freq = 0,
                 max.freq = 10000,
                 n.bands = 10,
                 cutoff = -60,
                 norm.spec = FALSE,
                 noise.red = 0,
                 rm.offset = TRUE,
                 props = TRUE,
                 prop.den = 2,
                 db.fs = TRUE,
                 use.vegan = FALSE){
  
  # Store the frequency step (band "height") # NEW 09/25/2023 Francisco Rivas
  freq_step <- (max.freq - min.freq)/n.bands
  
  # Store the decibel threshold as numeric
  cutoff <- as.numeric(cutoff)
  
  # Test arguments
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
  
  
  # Function that gets the proportion of values higher than the
  # db threshold in a specific frequency band. The frequencies are in Hz
  getscore <- function(spectrum, minf, maxf, db, freq_row, min.freq, max.freq){
    miny<-round((minf)/freq_row) # the minimum frequency of the frequency band
    maxy<-round((maxf)/freq_row) # the maximum frequency of the frequency band
    
    subA = spectrum[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)
    
    minspec <- round(min.freq/freq_row) # lower end of the spectrogram defined by min.freq
    maxspec <- round(max.freq/freq_row) # upper end of the spectrogram defined by max.freq
    speclims <- spectrum[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq
    
    # Calculate the proportion of cells in f that are higher than the dB threshold
    if(prop.den == 1){ #original ADI proportion calculation (within frequency band)
      
      index1 <- length(subA[subA>db]) / length(subA)
      
    } else if(prop.den == 2){
      # Calculate the total number of values in the ENTIRE spectrogram (within the specified frequency range) that exceed the threshold
      total_values <- length(speclims[speclims > db])
      
      # Calculate the number of values in the CURRENT band that exceed the threshold
      index1 <- length(subA[subA > db])
      
      # Calculate the proportion of values in the current band relative to the total values in the spectrogram
      if (total_values == 0) { # avoid division by zero
        index1 <- 0
      } else {
        index1 <- index1 / total_values
      }
    }

    # else if(prop.den == 3){
    #   # True Shannon's over the whole spectrogram
    #   # (over all the pixels above the threshold up to the Nyquist frequency)
    #   index1 <- length(subA[subA>db]) / length(spectrum[spectrum>db])
    # }
    
    return(index1)
  }
  
  # Function to normalize the proportions so that they add exactly to zero
  norm_props <- function(Score) {
    total <- sum(Score)
    if (total == 0) {
      return(Score) 
    }
    return(Score / total)
  }
  
  
  # Save the denominator used in the proportion calculation
  if(prop.den == 1){
    propdenom <- "within band"
  }else if(prop.den == 2){
    propdenom <- "whole spec"
  }
  # else if(prop.den == 3){
  #   propdenom <- "nyquist"
  # }
  
  
  #Get sampling rate
  samplingrate <- wave@samp.rate
  
  duration <- length(wave@left)/samplingrate
  
  #Get Nyquist frequency in Hz
  nyquist_freq <- samplingrate/2
  
  # Add information about noise reduction procedure
  if(noise.red == 1){
    # cat("Applying noise reduction filter to each row...\n")
    noise <- "rows"
  } else if (noise.red == 2){
    # cat("Applying noise reduction filter to each column...\n")
    noise <- "columns"
  } else {
    noise <- "none"
  }
  
  
  # In the original ADI from "soundecology", the frequency bin width was fixed to 10 Hz,
  # which is equivalent of a window length of 4,800 samples for a recording with
  # a sampling rate of 48 kHz. This creates a blurred spectrogram, distorting the
  # real distribution of energy in time.
  
  # window length for the spectro and spec functions (legacy)
  # to keep each row every 10Hz  (legacy)
  # Frequencies and seconds covered by each (legacy)
  # freq_per_row = 10 (legacy)
  # wlen = samplingrate/freq_per_row (legacy)
  
  wlen = samplingrate/freq.res
  
  # Adding 1 if wlen is an odd number (new behavior in seewave)
  # fix by JSueur
  if(wlen%%2 == 1) {wlen <- wlen+1}
  
  # Calculate frequency resolution (i.e., frequency bin width)
  # freq_per_row = samplingrate/wlen
  freq_per_row = freq.res
  
  
  #Stereo file
  if (wave@stereo == TRUE) {
    
    cat("Calculating ADI on a stereo file... \n")
    
    # Add information about noise reduction procedure
    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    } 
    
    left<-channel(wave, which = c("left"))
    right<-channel(wave, which = c("right"))
    
    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left, output = "Wave")
      right <- rmoffset(right, output = "Wave")
    }
    
    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){
      
      cat("Using normalized spectrograms.\n\n")
      
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- seewave::spectro(left, wl = wlen,
                                       wn = win.fun,
                                       noisereduction = noise.red,
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right, wl = wlen,
                                        wn = win.fun,
                                        noisereduction = noise.red,
                                        plot = FALSE)$amp
      }else if (noise.red == 0) {
        
        # Use the original parameters for window length and function
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       wn = win.fun,
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right,
                                        wl = wlen,
                                        wn = win.fun,
                                        plot = FALSE)$amp
      }
      
      
      rm(left, right)
      cat("Left range: ", range(specA_left), "\n")
      cat("Right range: ", range(specA_right), "\n")
      
      
    }else{
      # Without normalizing the spectrogram. Procedure extracted from Xu et al
      cat("No spectrogram normalization...\n\n")
      
      # if(!is.null(noise.red)){
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       norm=FALSE,
                                       dB=NULL,
                                       noisereduction = noise.red,
                                       correction = 'amplitude',
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right,
                                        wl = wlen,
                                        norm=FALSE,
                                        dB=NULL,
                                        noisereduction = noise.red,
                                        correction = 'amplitude',
                                        plot = FALSE)$amp
      }else if(noise.red == 0){
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       norm=FALSE,
                                       dB=NULL,
                                       correction = 'amplitude',
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right,
                                        wl = wlen,
                                        norm=FALSE,
                                        dB=NULL,
                                        correction = 'amplitude',
                                        plot = FALSE)$amp
      }
      
      if(db.fs==TRUE){
        
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
        specA_left <- 20 * log10(abs(specA_left) / amp_max)
        specA_right <- 20 * log10(abs(specA_right) / amp_max)
        
        
        
      }else{
        
        # Transform to decibels
        specA_left <- 10*log10(specA_left^2)
        specA_right <- 10*log10(specA_right^2)
        
        
      }
      
      rm(left)
      rm(right)
      
      cat("Left dB range: ", range(specA_left), "\n")
      cat("Right dB range: ", range(specA_right), "\n")
      
    }
    
    cat("Spectrograms extracted.\n")
    
    # The maximum frequency must be equal or smaller than the Nyquist frequency
    if (max.freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
      max.freq <- nyquist_freq
    }
    
    # Set the frequency bands (Fran's comment)
    # Freq <- seq(from = min.freq, to = max.freq - freq_step, by = freq_step)
    
    Freq <- seq(from = 0, to = max.freq - freq_step, by = freq_step)
    
    
    
    #LEFT CHANNEL
    minspec <- round(min.freq/freq_per_row) # lower end of the spectrogram defined by min.freq
    maxspec <- round(max.freq/freq_per_row) # upper end of the spectrogram defined by max.freq
    speclims <- specA_left[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq
    
    total_values <- length(speclims[speclims > cutoff]) # Total values above threshold in the whole spectrogram
    
    Score <- rep(0, n.bands) # Initialize Score vector with 0s
    
    # Calculate the proportion-of-energy score, looping over each frequency band
    for (j in 1:n.bands) {
      
      miny <- minspec + (j-1) * round((maxspec - minspec)/n.bands) # Calculate starting index of band
      maxy <- minspec + j * round((maxspec - minspec)/n.bands) -1 # Calculate ending index of band. -1 is crucial
      
      # If maxy exceeds maxspec, then force maxy=maxspec
      if(maxy > maxspec){
        maxy <- maxspec
      }
      
      subA = specA_left[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)
      
      index1 <- length(subA[subA > cutoff])
      
      if (total_values == 0) {
        index1 <- 0
      } else {
        index1 <- index1 / total_values
      }
      
      Score[j] <- index1
    }
    
    # Normalize Score values 
    Score <- norm_props(Score)
    
    left_vals <- Score
    
    
    if(use.vegan == TRUE){
      
      Score_left <- vegan::diversity(Score, index = "shannon")
      
      
    }else{
      
      Score_left <- -sum(Score * log(Score), na.rm = TRUE)
      
      
    }
    
    
    
    #RIGHT CHANNEL
    minspec <- round(min.freq/freq_per_row) # lower end of the spectrogram defined by min.freq
    maxspec <- round(max.freq/freq_per_row) # upper end of the spectrogram defined by max.freq
    speclims <- specA_right[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq
    
    total_values <- length(speclims[speclims > cutoff]) # Total values above threshold in the whole spectrogram
    
    Score <- rep(0, n.bands) # Initialize Score vector with 0s
    
    # Calculate the proportion-of-energy score, looping over each frequency band
    for (j in 1:n.bands) {
      
      miny <- minspec + (j-1) * round((maxspec - minspec)/n.bands) # Calculate starting index of band
      maxy <- minspec + j * round((maxspec - minspec)/n.bands) -1 # Calculate ending index of band. -1 is crucial
      
      #Ensure that the maxy does not exceed the maxspec, if it does, then maxy=maxspec
      if(maxy > maxspec){
        maxy <- maxspec
      }
      
      subA = specA_right[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)
      
      index1 <- length(subA[subA > cutoff])
      
      if (total_values == 0) {
        index1 <- 0
      } else {
        index1 <- index1 / total_values
      }
      
      Score[j] <- index1
    }
    
    # Normalize Score values  
    Score <- norm_props(Score)
    
    right_vals <- Score
    
    
    if(use.vegan == TRUE){
      
      Score_right <- vegan::diversity(Score, index = "shannon")
      
    }else{
      
      Score_right <- -sum(Score * log(Score), na.rm = TRUE)
      
    }
    
    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))
    
    
    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = left_vals[j]
      right_bandvals_return[j] = right_vals[j]
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
      right_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }
    
    left_adi_return = round(Score_left, 3)
    right_adi_return = round(Score_right, 3)
    
    
    adiOutputStereo <- tibble(index = "adi",
                              value_l = left_adi_return,
                              value_r = right_adi_return)
    
    adiOutputStereo <- adiOutputStereo |>
      add_column(value_avg = ((adiOutputStereo$value_l+adiOutputStereo$value_r)/2), 
                 .after = "value_r")
    
    
    
    # Add metadata columns
    adiOutputStereo <- adiOutputStereo |>
      add_column(w_len = wlen,
                 w_fun = win.fun,
                 cutoff = cutoff,
                 min_f = min.freq,
                 max_f = max.freq,
                 n_bands = n.bands,
                 norm = norm.spec,
                 noise_red = noise,
                 rm_off = rm.offset,
                 prop_den = propdenom,
                 samp_rate = samplingrate,
                 freq_res = freq_per_row,
                 nyq = nyquist_freq,
                 dur = duration,
                 channels = "stereo")
    
    
    
    if(props == TRUE){
      # Add data on proportions for each frequency band per channel:
      proportions = tibble(band_kHz = left_bandrange_return,
                           left = left_bandvals_return,
                           right = right_bandvals_return)
      # Wider format
      proportions <- pivot_wider(proportions, names_from = "band_kHz",
                                 values_from = c("left", "right"))
      
      # Append to the original data frame
      adiOutputStereo <- bind_cols(adiOutputStereo, proportions)
      
      cat("Reporting ADI for 2 channels, metadata and energy proportions per frequency band. \n")
      
      return(adiOutputStereo)
      
    } else {
      
      cat("Reporting ADI for 2 channels and metadata.\n")
      
    }
    
    return(adiOutputStereo)
    
    
  } else
    
  {
    # MONO
    cat("Calculating ADI on a mono file... \n")
    
    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    } 
    
    left<-channel(wave, which = c("left"))
    
    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left, output = "Wave")
    }
    
    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){
      
      cat("Using normalized spectrograms.\n\n")
      
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       noisereduction = noise.red,
                                       plot = FALSE)$amp
      }else if (noise.red == 0) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       plot = FALSE)$amp
      }
      rm(left)
      
      cat("dB range: ", range(specA_left), "\n")
      
    }else{
      
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       norm=FALSE,
                                       dB=NULL,
                                       noisereduction = noise.red,
                                       plot = FALSE)$amp
      }else if (noise.red == 0) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       norm = FALSE,
                                       dB = NULL,
                                       plot = FALSE)$amp
      }
      
      # Transform to decibels
      if(db.fs==TRUE){
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
        specA_left <- 20 * log10(abs(specA_left) / amp_max)
        
      }else{
        # Transform to decibels
        specA_left <- 10*log10(specA_left^2)
      }
      
      rm(left)
      cat("dB range: ", range(specA_left), "\n")
      
    }
    
    
    if (max.freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max.freq <- nyquist_freq
    }
    
    minspec <- round(min.freq/freq_per_row) # lower end of the spectrogram defined by min.freq
    maxspec <- round(max.freq/freq_per_row) # upper end of the spectrogram defined by max.freq
    speclims <- specA_left[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq
    
    total_values <- length(speclims[speclims > cutoff]) # Total values above threshold in the whole spectrogram
    
    Score <- rep(0, n.bands) # Initialize Score vector with 0s
    
    # Calculate the proportion-of-energy score, looping over each frequency band
    for (j in 1:n.bands) {
      
      miny <- minspec + (j-1) * round((maxspec - minspec)/n.bands) # Calculate starting index of band
      maxy <- minspec + j * round((maxspec - minspec)/n.bands) -1 # Calculate ending index of band. -1 is crucial
      
      #Ensure that the maxy does not exceed the maxspec, if it does, then maxy=maxspec
      if(maxy > maxspec){
        maxy <- maxspec
      }
      
      subA = specA_left[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)
      
      index1 <- length(subA[subA > cutoff])
      
      if (total_values == 0) {
        index1 <- 0
      } else {
        index1 <- index1 / total_values
      }
      
      Score[j] <- index1
    }
    
    
    Score <- norm_props(Score)
    
    
    left_vals = Score
    
    if(use.vegan == TRUE){
      
      Score_left <- vegan::diversity(Score, index = "shannon")
      
    }else{
      
      Score_left <- -sum(Score * log(Score), na.rm = TRUE)
      
      
    }
    
    left_adi_return = round(Score_left, 3)
    
    left_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    
    
    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = left_vals[j]
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }
    
    left_adi_return = round(Score_left, 3)
    
    adiOutputMono <- tibble(index = "adi",
                            value = left_adi_return)
    
    # Add metadata columns
    adiOutputMono <- adiOutputMono |>
      add_column(w_len = wlen,
                 w_fun = win.fun,
                 cutoff = cutoff,
                 min_f = min.freq,
                 max_f = max.freq,
                 n_bands = n.bands,
                 norm = norm.spec,
                 noise_red = noise,
                 rm_off = rm.offset,
                 prop_den = propdenom,
                 samp_rate = samplingrate,
                 freq_res = freq_per_row,
                 nyq = nyquist_freq,
                 dur = duration,
                 channels = "mono")
    
    if(props == TRUE){
      # Add data on proportions for each frequency band per channel:
      proportions = tibble(band_kHz = left_bandrange_return,
                           left = left_bandvals_return)
      # Wider format
      proportions <- pivot_wider(proportions, names_from = "band_kHz",
                                 values_from = c("left"))
      
      # Append to the original data frame
      adiOutputMono <- bind_cols(adiOutputMono, proportions)
      
      cat("Reporting ADI for 1 channel, metadata and energy proportions per frequency band. \n")
      
    }else {
      
      cat("Reporting ADI for 1 channel and metadata.\n")
      
    }
    
    return(adiOutputMono)
    
  }
  
}
