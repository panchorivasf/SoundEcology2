#' Acoustic Complexity Index
#'
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param wlen the window length to compute the spectrogram (i.e., FFT window size).
#' @param wfun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min_freq minimum frequency to use when calculating the value, in Hertz. Default = 0.
#' @param max_freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param j the cluster size, in seconds. Default = NA (Duration of the audio file).
#' @param noisered numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rmoff logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @description
#' Acoustic Complexity Index (ACI) from Pieretti, et al. 2011. The ACI is based
#' on the "observation that many biotic sounds, such as bird songs,
#' are characterized by an intrinsic variability of intensities,
#' while some types of human generated noise (such as car passing or airplane transit)
#' present very constant intensity values" (Pieretti, et al. 2011).
#' @details
#' The index was tested to the ACItot calculated using SoundscapeMeter v 1.0.14.05.2012,
#' courtesy of A. Farina. The results are accumulative.
#' Very long samples will return comparatively larger values for ACI.
#' The current version (\emph{soundecology2}) normalizes the output (i.e., "j" equals the
#' duration of the audio file) to make it equivalent to the default results in seewave's version.
#'
#' Reference: N. Pieretti, A. Farina, D. Morri. 2011. A new methodology to infer
#' the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators 11: 868-873.
#'
#'
#' @return A tibble (data frame) with the ACI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @importFrom tuneR readWave
#' @import seewave
#' @import tibble
#'
#' @examples aci(tropicalsound)
aci <- function(wave,
                wlen = 512,
                wfun = "hanning",
                min_freq = NA,
                max_freq = NA,
                j = NA,
                noisered = 2,
                rmoff = TRUE
                ){


  #test arguments
  if (is.na(max_freq)){
    max_freq <- wave@samp.rate / 2
    # cat(paste("\n max_freq not set, using value of:", max_freq, "\n\n"))
  }

  if (is.na(min_freq)){
    min_freq <- 0
    # cat(paste("\n min_freq not set, using value of:", min_freq, "\n\n"))
  }

  if (is.numeric(as.numeric(min_freq))){
    min_freq <- as.numeric(min_freq)
  } else{
    stop(" min_freq is not a number.")
  }

  if (is.numeric(as.numeric(max_freq))){
    max_freq <- as.numeric(max_freq)
  } else{
    stop(" max_freq is not a number.")
  }



  if (is.numeric(as.numeric(wlen))){
    wlen <- as.numeric(wlen)
  } else{
    stop(" wlen is not a number.")
  }



  #Some general values
  #Get sampling rate
  samplingrate <- wave@samp.rate
  duration <- length(wave@left)/samplingrate


  if(is.na(j)){
    j <- duration
  }

  if (is.numeric(as.numeric(j))){
    j <- as.numeric(j)
  } else{
    stop(" j is not a number.")
  }


  #function that gets the difference of values
  get_d <- function(spectrum, freq_row, min_col, max_col){
    D = 0
    for (k in min_col:(max_col - 1)) {
      D = D + abs(spectrum[freq_row,k] - spectrum[freq_row,k + 1])
    }

    return(D)
  }


  #Get Nyquist frequency in Hz
  nyquist_freq <- (samplingrate/2)
  if (max_freq>nyquist_freq) {
    cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
    max_freq <- nyquist_freq
    #break
  }

  # #window length for the spectro and spec functions
  # wlen = fft_w

  # Adding 1 if wlen is an odd number (new behavior in seewave)
  # fix by JSueur
  if(wlen%%2 == 1) {wlen <- wlen+1}

  # Calculate frequency resolution (i.e., frequency bin width)
  freq_per_row = samplingrate/wlen



  #Stereo file
  if (wave@stereo == TRUE) {
    cat("Calculating ACI on a stereo file... \n")


    left <- channel(wave, which = c("left"))
    right <- channel(wave, which = c("right"))

    # Remove DC offset
    if(rmoff == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left)
      right <- rmoffset(right)
    }


    if(noisered == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noisered == 2){
      cat("Applying noise reduction filter to each column...\n")
    }

    if(noisered == 1 || noisered == 2) {
    spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                         norm = TRUE, dB = NULL, scale = FALSE, wn = wfun,
                         noisereduction = noisered)
    }else if (noisered == 0) {
      spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                           norm = TRUE, dB = NULL, scale = FALSE, wn = wfun)
    }



    specA_left <- spec_left$amp

    min_freq1k = min_freq/1000
    max_freq1k = max_freq/1000

    which_min_freq <- which(abs(spec_left$freq - min_freq1k)==min(abs(spec_left$freq - min_freq1k)))
    which_max_freq <- which(abs(spec_left$freq - max_freq1k)==min(abs(spec_left$freq - max_freq1k)))

    if (which_min_freq <1){
      which_min_freq = 1
    }

    if (which_max_freq > dim(specA_left)[1]){
      which_max_freq = dim(specA_left)[1]-1
    }

    # 		cat(which_min_freq)
    # 		cat(",")
    # 		cat(which_max_freq)
    # 		cat(",")
    # 		cat(dim(specA_left))
    specA_left <- spec_left$amp[which_min_freq:which_max_freq,]
    rm(spec_left)


    if(noisered == 1 || noisered == 2) {
    spec_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE,
                          norm = TRUE, dB = NULL, scale = FALSE, wn = wfun,
                          noisereduction = noisered)
    }else if (noisered == 3) {
      spec_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE,
                            norm = TRUE, dB = NULL, scale = FALSE, wn = wfun)
    }



    specA_right <- spec_right$amp[which_min_freq:which_max_freq,]

    rm(spec_right)

    rm(left,right)

    # 		specA_rows <- dim(specA_left)[1]
    # 		specA_cols <- dim(specA_left)[2]
    #
    # 		freq_per_row <- specA_rows/nyquist_freq
    #
    # 		max_row <- round(max_freq * freq_per_row)
    #
    # 		specA_left <- specA_left[1:max_row,]
    # 		specA_right <- specA_right[1:max_row,]
    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]

    fl <- rep(NA, specA_rows)
    delta_fl <- ( max_freq - min_freq ) / specA_rows
    delta_tk <- (length(wave@left)/wave@samp.rate) / specA_cols

    #m <- floor(duration / j)
    #q <- specA_rows
    no_j <- floor(duration / j)

    #Number of values, in each row, for each j period (no. of columns)
    I_per_j <- floor(j/delta_tk)

    ACI_left_vals <- rep(NA, no_j)
    ACI_fl_left_vector <- rep(NA, no_j)
    ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))

    ACI_right_vals <- rep(NA, no_j)
    ACI_fl_right_vector <- rep(NA, no_j)
    ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))

    #Left channel
    #For each frequency bin fl
    for (q_index in 1:specA_rows) {

      #For each j period of time
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j

        D <- get_d(specA_left, q_index, min_col, max_col)
        sum_I <- sum(specA_left[q_index,min_col:max_col])
        ACI_left_vals[j_index] <- D / sum_I
        ACI_left_matrix[q_index, j_index] <- D / sum_I
      }

      ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
    }

    ACI_tot_left <- sum(ACI_fl_left_vector)
    # ACI_tot_left <- as.numeric(ACI_tot_left)

    #Right channel
    #For each frequency bin fl
    for (q_index in 1:specA_rows) {

      #For each j period of time
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j

        D <- get_d(specA_right, q_index, min_col, max_col)
        sum_I <- sum(specA_right[q_index, min_col:max_col])
        ACI_right_vals[j_index] <- D / sum_I
        ACI_right_matrix[q_index, j_index] <- D / sum_I
      }

      ACI_fl_right_vector[q_index] <- sum(ACI_right_vals)
    }

    ACI_tot_right <- sum(ACI_fl_right_vector)
    # ACI_tot_right <- as.numeric(ACI_tot_right)

    ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 60)#, 2)
    ACI_tot_right_by_min <- round((ACI_tot_right/duration) * 60)#, 2)

    # cat(paste("  Acoustic Complexity Index (total):\n", "   Left channel: ", sep=""))
    # cat(ACI_tot_left)
    # cat(paste("\n", "   Right channel: ", sep=""))
    # cat(ACI_tot_right)
    # cat("\n\n")
    # if (duration > 60){
    #   cat(paste("  Acoustic Complexity Index (by minute):\n", "   Left channel: ", sep=""))
    #   cat(ACI_tot_left_by_min)
    #   cat(paste("\n", "   Right channel: ", sep=""))
    #   cat(ACI_tot_right_by_min)
    #   cat("\n\n")
    # }

    # if(noisered == 1){
    #   noise <- "rows"
    # } else if(noisered == 2){
    #   noise <- "columns"
    # } else {
      # noise <- "none"
    # }


    aciOutputStereo <- tibble(value_l = ACI_tot_left,
                              value_r = ACI_tot_right)

    aciOutputStereo <- aciOutputStereo %>%
      add_column(value_avg = ((aciOutputStereo$value_l+aciOutputStereo$value_r)/2), .after = "value_r")


    # Add metadata columns
    aciOutputStereo <- aciOutputStereo %>%
      add_column(wlen = wlen,
                 wfun = wfun,
                 j = j,
                 minf = min_freq,
                 maxf = max_freq,
                 noisered = "none",
                 rmoff = rmoff,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "stereo")

    aciOutputStereo <- aciOutputStereo %>%
      add_column(index = "aci", .before = "value_l")

    return(aciOutputStereo)



  } else
  {
    cat("Calculating ACI on a mono file... \n")

    left<-channel(wave, which = c("left"))


    # Remove DC offset
    if(rmoff == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left)
      # right <- rmoffset(right)
    }

    if(noisered == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noisered == 2){
      cat("Applying noise reduction filter to each column...\n")
    }


    #matrix of values
    # cat("\n Calculating index. Please wait... \n\n")
    if(noisered == 1 || noisered == 2) {
      spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                           norm = TRUE, dB = NULL, scale = FALSE, wn = wfun,
                           noisereduction = noisered)
    }else if (noisered == 0) {
      spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                           norm = TRUE, dB = NULL, scale = FALSE, wn = wfun)
    }

    # spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, wn = wfun)

    specA_left <- spec_left$amp

    min_freq1k = min_freq/1000
    max_freq1k = max_freq/1000

    which_min_freq <- which(abs(spec_left$freq - min_freq1k)==min(abs(spec_left$freq - min_freq1k)))
    which_max_freq <- which(abs(spec_left$freq - max_freq1k)==min(abs(spec_left$freq - max_freq1k)))

    specA_left <- specA_left[which_min_freq:which_max_freq,]
    rm(spec_left)

    rm(left)

    #LEFT CHANNEL
    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]
    #
    # 		freq_per_row <- specA_rows/nyquist_freq
    #
    # 		max_row <- round(max_freq * freq_per_row)
    #
    # 		specA_left <- specA_left[1:max_row,]
    # 		specA_rows <- dim(specA_left)[1]

    fl <- rep(NA, specA_rows)
    delta_fl <- ( max_freq - min_freq ) / specA_rows
    delta_tk <- (length(wave@left)/wave@samp.rate) / specA_cols

    no_j <- floor(duration / j)
    #q <- specA_rows
    #m <- floor(duration / j)

    #Number of values, in each row, for each j period (no. of columns)
    I_per_j <- floor(j/delta_tk)

    ACI_left_vals <- rep(NA, no_j)
    ACI_fl_left_vector <- rep(NA, no_j)
    ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))

    ACI_right_vals <- rep(NA, no_j)
    ACI_fl_right_vector <- rep(NA, no_j)
    ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))

    #Left channel
    #For each frequency bin fl
    for (q_index in 1:specA_rows) {

      #For each j period of time
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j

        D <- get_d(specA_left, q_index, min_col, max_col)
        sum_I <- sum(specA_left[q_index, min_col:max_col])
        ACI_left_vals[j_index] <- D / sum_I
        ACI_left_matrix[q_index, j_index] <- D / sum_I
      }

      ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
    }

    ACI_tot_left <- sum(ACI_fl_left_vector)
    ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 60, 2)

    ACI_tot_right <- NA
    ACI_tot_right_by_min <- NA

    # cat("  Acoustic Complexity Index (total): ")
    # cat(ACI_tot_left)
    # cat("\n\n")
    # if (duration > 60){
    #   cat("  Acoustic Complexity Index (by minute): ")
    #   cat(ACI_tot_left_by_min)
    #   cat("\n\n")
    # }

  # if(noisered == 1){
  #   noise <- "rows"
  # } else if(noisered == 2){
  #   noise <- "columns"
  # } else {
  #   noise <- "none"
  # }


    aciOutputMono <- tibble(value = ACI_tot_left)

     # Add metadata columns
    aciOutputMono <- aciOutputMono %>%
      add_column(wlen = wlen,
                 wfun = wfun,
                 j = j,
                 minf = min_freq,
                 maxf = max_freq,
                 noisered = "none",
                 rmoff = rmoff,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "mono")

    aciOutputMono <- aciOutputMono %>%
      add_column(index = "aci", .before = "value")

    return(aciOutputMono)


  }

  # invisible(list(AciTotAll_left = ACI_tot_left, AciTotAll_right = ACI_tot_right,
  #                AciTotAll_left_bymin = ACI_tot_left_by_min, AciTotAll_right_bymin = ACI_tot_right_by_min,
  #                #AciIfTotAll_left=ACIif_tot_left, AciIfTotAll_right=ACIif_tot_right,
  #                aci_fl_left_vals = ACI_fl_left_vector, aci_fl_right_vals = ACI_fl_right_vector,
  #                #aci_if_left_vals=ACI_if_left_vector, aci_if_right_vals=ACI_if_right_vector,
  #                aci_left_matrix = ACI_left_matrix, aci_right_matrix = ACI_right_matrix))
}
