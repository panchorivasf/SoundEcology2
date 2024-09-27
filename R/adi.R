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
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram.
#' @param max.freq maximum frequency to compute the spectrogram.
#' @param n.bands number of bands to split the spectrogram.
#' @param cutoff dB threshold to calculate energy proportions.
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param props logical; if set to TRUE, the function stores the energy proportion values for each frequency band and channel. Default = TRUE.
#' @param prop.den numeric; indicates how the energy proportion is calculated.
#' @param use.vegan logical; if TRUE, the function uses the \emph{diversity} function from the \emph{vegan} package to calculate Shannon's Diversity.
#' @param db.fs logical; if TRUE, the amplitude scale is expressed as decibels Full Scale (dBFS). Only used when norm = FALSE.
#'
#' @return A tibble (data frame) with the ADI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#'
#' @importFrom tuneR readWave
#' @import seewave
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import vegan
#'
#' @details
#' Options for the 'prop.den' parameter: 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion equals to all the cells in the same frequency band. 2 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max_freq'). 3 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.
#' Another important update is that now the spectrogram is not normalized by default, which made recordings
#' with different signal-to-noise ratio not comparable.
#' @examples
#' data(tropicalsound)
#' adi(tropicalsound)

adi <- function(wave,
                freq.res = 50,
                w.fun = "hanning",
                min.freq = 0,
                max.freq = 10000,
                n.bands = 10,
                cutoff = -60,
                norm.spec = FALSE,
                noise.red = 0,
                rm.offset = TRUE,
                props = FALSE,
                prop.den = 1,
                use.vegan = FALSE,
                db.fs = TRUE){



  # Store the frequency step (band "height") # NEW 09/25/2023 Francisco Rivas
  freq_step <- (max.freq - min.freq)/n.bands

  # Test arguments
  # Check if the maximum frequency, decibel threshold,
  # and frequency step arguments are numbers:
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
  getscore <- function(spectrum, minf, maxf, db, freq_row){
    # miny<-round((minf)/freq_row) # the minimum frequency of the frequency band
    # maxy<-round((maxf)/freq_row) # the maximum frequency of the frequency band
    miny <- floor(minf / freq_row) # more efficient
    maxy <- ceiling(maxf / freq_row) # more efficient

    subA = spectrum[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)



    # Calculate the proportion of cells in f that are higher than the dB threshold
    if(prop.den == 1){ #original ADI proportion calculation (within frequency band)

      # index1 <- length(subA[subA>db]) / length(subA)
      index1 <- mean(subA > db)


    }else if(prop.den == 2){
      minspec <- round(0/freq_row) # lower end of the spectrogram defined by min.freq
      maxspec <- ceiling(max.freq/freq_row) # upper end of the spectrogram defined by max.freq
      speclims <- spectrum[minspec:maxspec,] # the spectrogram with range defined by min.freq and max.freq

      # Alternative 2:'true Shannon', over the user-defined spectrogram range
      index1 <- length(subA[subA>db]) / length(speclims[speclims>db])

    }else if(prop.den == 3){
      # True Shannon's over the whole spectrogram
      # (over all the pixels above the threshold up to the Nyquist frequency)
      index1 <- length(subA[subA>db]) / length(spectrum[spectrum>db])
    }
    return(index1)
  }

  # Save the denominator used in the proportion calculation
  if(prop.den == 1){
    propdenom <- "within band"
  }else if(prop.den == 2){
    propdenom <- "max.freq"
  }else if(prop.den == 3){
    propdenom <- "nyquist"
  }

  #Get sampling rate
  samplingrate <- wave@samp.rate

  duration <- length(wave@left)/samplingrate

  #Get Nyquist frequency in Hz
  nyquist_freq <- samplingrate/2

  # Add information about noise reduction procedure
  if(noise.red == 1){
    noise <- "rows"
  } else if (noise.red == 2){
    noise <- "columns"
  } else {
    noise <- "none"
  }


  # NOTE: In the original ADI from "soundecology", the frequency bin width was hard-coded
  # to 10 Hz, which is equivalent of a window length of 4,800 points for a recording with
  # a sampling rate of 48 kHz. This creates a blurred spectrogram, distorting the real
  # distribution of energy in time.


  wlen = samplingrate/freq.res

  # Adding 1 if wlen is an odd number (new behavior in seewave)
  # fix by JSueur
  if(wlen%%2 == 1) {wlen <- wlen+1}

  # Calculate frequency resolution (i.e., frequency bin width)
  freq_per_row = freq.res

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


  #Stereo file
  if (wave@stereo == TRUE) {

    cat("Calculating ADI on a stereo file... \n")

    # Add information about noise reduction procedure
    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
      # noise <- "rows"
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
      # noise <- "columns"
    } else {
      # noise <- "none"
    }

    left<-channel(wave, which = c("left"))
    right<-channel(wave, which = c("right"))

    rm(wave)

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- seewave::rmoffset(left, output = "Wave")
      right <- seewave::rmoffset(right, output = "Wave")
    }


    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noise.red == 1 || noise.red == 2) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       wn = w.fun,
                                       noisereduction = noise.red,
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right,
                                        wl = wlen,
                                        wn = w.fun,
                                        noisereduction = noise.red,
                                        plot = FALSE)$amp
        rm(left, right)


      }else if (noise.red == 0) {

        # Use the original parameters for window length and function
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       wn = w.fun,
                                       plot = FALSE)$amp
        specA_right <- seewave::spectro(right,
                                        wl = wlen,
                                        wn = w.fun,
                                        plot = FALSE)$amp

        rm(left, right)

      }


      cat("Left range: ", range(specA_left), "\n")
      cat("Right range: ", range(specA_right), "\n")


    }else{
      # Without normalizing the spectrogram. Procedure extracted from Xu et al
      cat("No spectrogram normalization...\n\n")

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


        rm(left, right)


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

        rm(left, right)

      }

      if(db.fs==TRUE){ # Added by Francisco Rivas, September 2024.


        # Transform raw amplitude to dBFS
        specA_left <- 20 * log10(abs(specA_left) / amp_max)
        specA_right <- 20 * log10(abs(specA_right) / amp_max)



      }else{

        # Transform to decibels
        specA_left <- 10*log10(specA_left^2)
        specA_right <- 10*log10(specA_right^2)


      }


      cat("Left dB range: ", range(specA_left), "\n")
      cat("Right dB range: ", range(specA_right), "\n")

    }

    cat("Spectrograms extracted.\n")

    # The maximum frequency must be equal or smaller than the Nyquist frequency
    if (max.freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
      max.freq <- nyquist_freq
    }

    # Set the frequency bands
    # Freq <- seq(from = min.freq, to = max.freq - freq_step, by = freq_step)
    Freq <- seq(from = 0, to = max.freq - freq_step, by = freq_step)



    #LEFT CHANNEL

    # Score for the left channel
    Score <- sapply(Freq, function(f) {
      getscore(specA_left, f, (f + freq_step), cutoff, freq_per_row)
    })

    rm(specA_left)

    # Calculate Score1
    Score1 <- sum(Score * log(Score + 1e-7))


    # Store the score vector as 'left_vals' (as a copy) (Fran's comment)
    left_vals = Score




    if(use.vegan){

      Score_left <- vegan::diversity(Score, index = "shannon")


    }else{

      Score_left <- -sum(Score * log(Score + 0.000001))



    }


    # Score for the left channel
    Score <- sapply(Freq, function(f) {
      getscore(specA_right, f, (f + freq_step), cutoff, freq_per_row)
    })

    rm(specA_right)

    # Calculate Score1
    Score1 <- sum(Score * log(Score + 1e-7))

    right_vals = Score


    if(use.vegan){

      Score_right <- vegan::diversity(Score, index = "shannon")


    }else{

      # Use a vectorized version to improve efficiency
      Score_right <- -sum(Score * log(Score+ 0.000001))

    }



    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))


    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      right_bandvals_return[j] = round(right_vals[j], 6)
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
      right_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }

    left_adi_return = round(Score_left, 3)
    right_adi_return = round(Score_right, 3)



    adiOutputStereo <- tibble(index = "adi",
                              value_l = left_adi_return,
                              value_r = right_adi_return)

    adiOutputStereo <- adiOutputStereo %>%
      add_column(value_avg = ((adiOutputStereo$value_l+adiOutputStereo$value_r)/2), .after = "value_r")



    # Add metadata columns
    adiOutputStereo <- adiOutputStereo %>%
      add_column(w_len = wlen,
                 w_fun = w.fun,
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


  } else {
    # MONO
    cat("Calculating ADI on a mono file... \n")

    # Add information about noise reduction procedure
    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    } else {
    }

    left<-channel(wave, which = c("left"))

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- seewave::rmoffset(left, output = "Wave")
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
        rm(left)

      }else if (noise.red == 0) {
        specA_left <- seewave::spectro(left,
                                       wl = wlen,
                                       norm = FALSE,
                                       dB = NULL,
                                       plot = FALSE)$amp
        rm(left)

      }

      # # Transform to decibels
      # if(db.fs==TRUE){
      #   # Calculate amp_max based on bit depth
      #   amp_max <- if (wave@bit == 16) {
      #     32768
      #   } else if (wave@bit == 24) {
      #     8388607
      #   } else if (wave@bit == 32) {
      #     2147483647
      #   } else {
      #     stop("Unsupported bit depth")
      #   }

        # Convert amplitude to dBFS
        specA_left <- 20 * log10(abs(specA_left) / amp_max)

      # }else{
        # Transform to decibels
        specA_left <- 10*log10(specA_left^2)
      # }

      rm(left)
      cat("dB range: ", range(specA_left), "\n")

    }


    if (max.freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max.freq <- nyquist_freq
    }

    Freq <- seq(from = min.freq, to = max.freq - freq_step, by = freq_step)

    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), cutoff, freq_per_row)
    }

    rm(specA_left)

    left_vals = Score

    if(prop.den == 1){

      Score_left <- vegan::diversity(Score, index = "shannon")


    }else{

      Score_left <- -sum(Score * log(Score + 0.000001))

    }

    left_adi_return = round(Score_left, 3)
    left_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))


    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }

    left_adi_return = round(Score_left, 3)

    adiOutputMono <- tibble(index = "adi",
                            value = left_adi_return)


    # Add metadata columns
    adiOutputMono <- adiOutputMono %>%
      add_column(w_len = wlen,
                 w_fun = w.fun,
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
