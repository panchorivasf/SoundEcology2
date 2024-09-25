#' Bioacoustic Index
#' @description
#' Inspired by the "Bioacoustic Index" from the paper:Boelman NT, Asner GP, Hart PJ, Martin RE. 2007.
#' Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne
#' remote sensing. Ecol Applications 17(8):2137-44. Based on Matlab code provided by NT Boelman.
#' Boelman et al. 2007 used min.freq=2000, max.freq=8000, w.len=512.
#' Several parts where changed, in particular log math, so this won't be
#' directly comparable to the original code in the paper.
#'
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param w.len the window length to compute the spectrogram (i.e., FFT window size).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq miminum frequency to use when calculating the value, in Hertz. Default = NA.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0 (Default), noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.

#' @return A tibble (data frame) with the BI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @import tuneR
#' @import seewave
#' @examples bi(tropicalsound)

bi <- function(wave,
               w.len = 512,
               w.fun = "hanning",
               min.freq = 2000,
               max.freq = 8000,
               norm.spec = FALSE,
               noise.red = 0,
               rm.offset = TRUE){

  #test arguments
  if (is.numeric(as.numeric(min.freq))){
    min.freq <- as.numeric(min.freq)
  } else{
    stop(" min.freq is not a number.")
  }

  if (is.numeric(as.numeric(max.freq))){
    max.freq <- as.numeric(max.freq)
  } else{
    stop(" max.freq is not a number.")
  }

  if (is.numeric(as.numeric(w.len))){
    w.len <- as.numeric(w.len)
  } else{
    stop(" w.len is not a number.")
  }

  # Add information about noise reduction procedure
  if(noise.red == 1){
    cat("Applying noise reduction filter (subtract median amplitude) to each row...\n")
    noise <- "rows"
  } else if(noise.red == 2){
    cat("Applying a noise reduction filter (subtract median amplitude) to each column...\n")
    noise <- "columns"
  } else if(noise.red == 0) {
    noise <- "none"
  }

  #Get sampling rate
  samplingrate <- wave@samp.rate
  duration <- length(wave@left)/samplingrate
  #Get Nyquist frequency in Hz
  nyquist_freq <- samplingrate/2

  # The maximum frequency must be equal or smaller than the Nyquist frequency
  if (max.freq > nyquist_freq) {
    cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
    max.freq <- nyquist_freq
  }


  freq_per_row = samplingrate/w.len

  # freq_per_row = 10
  # w.len = samplingrate/freq_per_row



  if (max.freq > nyquist_freq) {
    cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz.\n\n", sep = ""))
    #break
  }

  #Stereo file
  if (wave@stereo == TRUE) {
    cat("Calculating BI on a stereo file... \n")

    left <- channel(wave, which = c("left"))
    right <- channel(wave, which = c("right"))
    # rm(wave)

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



    # #Get values
    # spec_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE, dB = "max0")$amp
    # spec_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE, dB = "max0")$amp
    # #Clear from memory
    # rm(left, right)



    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              wn = w.fun,
                              noise.reduction = noise.red,
                              plot = FALSE)$amp
        specA_right <- spectro(right,
                               f = samplingrate,
                               wl = w.len,
                               wn = w.fun,
                               noise.reduction = noise.red,
                               plot = FALSE)$amp
      }else if (noise.red == 0) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              wn = w.fun,
                              noise.reduction = NULL,
                              plot = FALSE)$amp
        specA_right <- spectro(right,
                               f = samplingrate,
                               wl = w.len,
                               wn = w.fun,
                               noise.reduction = NULL,
                               plot = FALSE)$amp
      }


      rm(left, right)


    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")

      # if(!is.null(noise.red)){
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              plot = FALSE,
                              norm=FALSE,
                              dB=NULL,
                              noise.reduction = noise.red)$amp
        specA_right <- spectro(right,
                               f = samplingrate,
                               wl = w.len,
                               plot = FALSE,
                               norm=FALSE,
                               dB=NULL,
                               noise.reduction = noise.red)$amp
      }else if(noise.red == 0){
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              plot = FALSE,
                              norm=FALSE,
                              dB=NULL)$amp
        specA_right <- spectro(right,
                               f = samplingrate,
                               wl = w.len,
                               plot = FALSE,
                               norm=FALSE,
                               dB=NULL)$amp
      }

      amp_max <- switch(as.character(wave@bit),
                        `16` = 32767,
                        `24` = 8388607,
                        `32` = 2147483647,
                        stop("Unsupported bit depth"))

      # Convert amplitude to dBFS
      specA_left <- 20*log10(specA_left/amp_max)
      # specA_left <- 10*log10(specA_left^2)
      specA_right <- 20*log10(specA_right/amp_max)
      # specA_right <- 10*log10(specA_right^2)

      rm(left)
      rm(right)

    }

    # Get average amplitude in time
    specA_left <- apply(specA_left, 1, meandB)
    specA_right <- apply(specA_right, 1, meandB)

    # bin width
    rows_width = length(specA_left) / nyquist_freq

    min_row = min.freq * rows_width
    max_row = max.freq * rows_width

    #Select rows
    specA_left_segment <- specA_left[min_row:max_row]
    specA_right_segment <- specA_right[min_row:max_row]

    freq_range <- max.freq - min.freq
    freqs <- seq(from = min.freq, to = max.freq, length.out = length(specA_left_segment))

    specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)
    specA_right_segment_normalized <- specA_right_segment - min(specA_right_segment)

    #left_area <- trapz(freqs, specA_left_segment_normalized)
    left_area <- sum(specA_left_segment_normalized * rows_width)
    right_area <- sum(specA_right_segment_normalized * rows_width)

    #cat("\n")
    # cat("Bioacoustic Index:\n")

    # cat("   Left channel: ")
    # cat(left_area)
    # cat("\n   Right channel: ")
    # cat(right_area)
    # cat("\n\n")

    biOutputStereo <- tibble(
      index = "bi",
      value_l = left_area,
      value_r = right_area)

    # Add average
    biOutputStereo <- biOutputStereo %>%
      add_column(value_avg = ((biOutputStereo$value_l+biOutputStereo$value_r)/2), .after = "value_r")

    biOutputStereo <- biOutputStereo %>%
      add_column(w.len = w.len,
                 w.fun = w.fun,
                 minf = min.freq,
                 maxf = max.freq,
                 norm = norm.spec,
                 noise.red = noise,
                 rm.offset = rm.offset,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "stereo")

    return(biOutputStereo)

  } else
  { # MONO
    cat("Calculating BI on a mono file... \n")

    #Get left channel
    left<-channel(wave, which = c("left"))

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- seewave::rmoffsetset(left, output = "Wave")
    }

    # # Add information about noise reduction procedure
    # if(noise.red == 1){
    #   cat("Applying noise reduction filter (subtract median amplitude) to each row...\n")
    #   noise <- "rows"
    # } else if(noise.red == 2){
    #   cat("Applying a noise reduction filter (subtract median amplitude) to each column...\n")
    #   noise <- "columns"
    # } else if(noise.red == 0) {
    #   noise <- "none"
    # }
    #


#
#     #Get values
#     # cat("\n Calculating index. Please wait... \n\n")
#     spec_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE, dB = "max0")$amp
#     #Clear from memory
#     rm(left)

    # Generate normalized spectrogram if norm.spec = TRUE
    if(norm.spec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              wn = w.fun,
                              noise.reduction = noise.red,
                              plot = FALSE)$amp
        # specA_right <- spectro(right, f = samplingrate, wl = w.len,
        #                        wn = w.fun, noise.reduction = noise.red,
        #                        plot = FALSE)$amp
      }else if (noise.red == 0) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              wn = w.fun,
                              noise.reduction = NULL,
                              plot = FALSE)$amp
        # specA_right <- spectro(right, f = samplingrate, wl = w.len,
        #                        wn = w.fun, noise.reduction = NULL,
        #                        plot = FALSE)$amp
      }


      rm(left)


    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")

      # if(!is.null(noise.red)){
      if(noise.red == 1 || noise.red == 2) {
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              plot = FALSE,
                              norm=FALSE,
                              dB=NULL,
                              noise.reduction = noise.red)$amp
        # specA_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
        #                        norm=FALSE,dB=NULL,unit="power",
        #                        noise.reduction = noise.red)$amp
      }else if(noise.red == 0){
        specA_left <- spectro(left,
                              f = samplingrate,
                              wl = w.len,
                              plot = FALSE,
                              norm=FALSE,
                              dB=NULL,unit="power")$amp
        # specA_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
        #                        norm=FALSE,dB=NULL,unit="power")$amp
      }

      amp_max <- switch(as.character(wave@bit),
                        `16` = 32767,
                        `24` = 8388607,
                        `32` = 2147483647,
                        stop("Unsupported bit depth"))

      # Convert amplitude to dBFS
      specA_left <- 20*log10(specA_left/amp_max)
      # Transform to decibels
      # specA_left <- 10*log10(specA_left^2)
      # specA_right <- 10*log10(specA_right^2)

      rm(left)
      # rm(right)

    }

    #Get average in time
    specA_left <- apply(specA_left, 1, meandB)

    #How much Hz are covered per row
    rows_width = length(specA_left) / nyquist_freq

    min_row = min.freq * rows_width
    max_row = max.freq * rows_width

    #Select rows
    specA_left_segment <- specA_left[min_row:max_row]
    freq_range <- max.freq - min.freq
    freqs <- seq(from = min.freq, to = max.freq, length.out = length(specA_left_segment))

    specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)

    #left_area <- trapz(freqs, specA_left_segment_normalized)
    left_area <- sum(specA_left_segment_normalized * rows_width)

    # cat("Bioacoustic Index: ")
    # cat(left_area)
    # cat("\n\n")
    # right_area <- NA

    biOutputMono <- tibble(index = "bi",
                           value = left_area)


    biOutputMono <- biOutputMono %>%
      add_column(w.len = w.len,
                 w.fun = w.fun,
                 minf = min.freq,
                 maxf = max.freq,
                 norm = norm.spec,
                 noise.red = noise,
                 rm.offset = rm.offset,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "mono")


    return(biOutputMono)


  }
  # invisible(list(left_area = left_area, right_area = right_area))
}
