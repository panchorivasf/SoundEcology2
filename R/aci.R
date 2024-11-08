#' Acoustic Complexity Index
#'
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param freq.res numeric. The frequency resolution to use (Hz per bin) which will determine the window length for the FFT (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to use when calculating the value, in Hertz. Default = 0.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param j the cluster size, in seconds. Default = NA (Duration of the audio file).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
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
#' @importFrom tuneR readWave channel
#' @importFrom seewave rmoffset spectro
#' @importFrom tibble tibble add_column
#'
#' @examples 
#' \dontrun{
#' aci(tropicalsound)}
aci <- function(wave,
                freq.res = 50,
                win.fun = "hanning",
                min.freq = NA,
                max.freq = NA,
                j = NA,
                noise.red = 0,
                rm.offset = TRUE
                ){
  #test arguments
  if (is.na(max.freq)){
    max.freq <- wave@samp.rate / 2
  }

  if (is.na(min.freq)){
    min.freq <- 0
  }

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

  if (is.numeric(as.numeric(freq.res))){
    freq.res <- as.numeric(freq.res)
  } else{
    stop(" freq.res is not a number.")
  }

  # Extract general values
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

  if (max.freq>nyquist_freq) {
    cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max.freq, "Hz. The value of max.freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
    max.freq <- nyquist_freq
    #break
  }

  w.len <- samplingrate/freq.res

  # #window length for the spectro and spec functions
  # w.len = fft_w

  # Adding 1 if w.len is an odd number (new behavior in seewave)
  # fix by JSueur
  if(w.len%%2 == 1) {w.len <- w.len+1}

  # Calculate frequency resolution (i.e., frequency bin width)
  freq_per_row = freq.res



  #Stereo file
  if (wave@stereo == TRUE) {
    cat("Calculating ACI on a stereo file... \n")


    left <- channel(wave, which = c("left"))
    right <- channel(wave, which = c("right"))

    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- seewave::rmoffset(left, output = "Wave")
      right <- seewave::rmoffset(right, output = "Wave")
    }


    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    }

    if(noise.red == 1 || noise.red == 2) {
    spec_left <- seewave::spectro(left,
                         f = samplingrate,
                         wl = w.len,
                         plot = FALSE,
                         norm = TRUE,
                         dB = NULL,
                         scale = FALSE,
                         wn = win.fun,
                         noisereduction = noise.red)
    }else if (noise.red == 0) {
      spec_left <- seewave::spectro(left,
                           f = samplingrate,
                           wl = w.len,
                           plot = FALSE,
                           norm = TRUE,
                           dB = NULL,
                           scale = FALSE,
                           wn = win.fun)
    }



    specA_left <- spec_left$amp

    min.freq1k = min.freq/1000
    max.freq1k = max.freq/1000

    which_min.freq <- which(abs(spec_left$freq - min.freq1k)==min(abs(spec_left$freq - min.freq1k)))
    which_max.freq <- which(abs(spec_left$freq - max.freq1k)==min(abs(spec_left$freq - max.freq1k)))

    if (which_min.freq <1){
      which_min.freq = 1
    }

    if (which_max.freq > dim(specA_left)[1]){
      which_max.freq = dim(specA_left)[1]-1
    }

    specA_left <- spec_left$amp[which_min.freq:which_max.freq,]
    rm(spec_left)


    if(noise.red == 1 || noise.red == 2) {
    spec_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
                          norm = TRUE, dB = NULL, scale = FALSE, wn = win.fun,
                          noise.reduction = noise.red)
    }else if (noise.red == 0) {
      spec_right <- spectro(right, f = samplingrate, wl = w.len, plot = FALSE,
                            norm = TRUE, dB = NULL, scale = FALSE, wn = win.fun)
    }



    specA_right <- spec_right$amp[which_min.freq:which_max.freq,]

    rm(spec_right)

    rm(left,right)

    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]

    fl <- rep(NA, specA_rows)
    delta_fl <- ( max.freq - min.freq ) / specA_rows
    delta_tk <- (length(wave@left)/wave@samp.rate) / specA_cols

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

    ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 60)
    ACI_tot_right_by_min <- round((ACI_tot_right/duration) * 60)


    aciOutputStereo <- tibble(value_l = ACI_tot_left,
                              value_r = ACI_tot_right)

    aciOutputStereo <- aciOutputStereo |>
      add_column(value_avg = ((aciOutputStereo$value_l+aciOutputStereo$value_r)/2), .after = "value_r")


    # Add metadata columns
    aciOutputStereo <- aciOutputStereo |>
      add_column(w.len = w.len,
                 win.fun = win.fun,
                 j = j,
                 minf = min.freq,
                 maxf = max.freq,
                 noise.red = "none",
                 rm.offset = rm.offset,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "stereo")

    aciOutputStereo <- aciOutputStereo |>
      add_column(index = "aci", .before = "value_l")

    return(aciOutputStereo)



  } else
  {
    cat("Calculating ACI on a mono file... \n")

    left<-channel(wave, which = c("left"))


    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- seewave::rmoffset(left, output = "Wave")
    }

    if(noise.red == 1){
      cat("Applying noise reduction filter to each row...\n")
    } else if (noise.red == 2){
      cat("Applying noise reduction filter to each column...\n")
    }


    #matrix of values
    if(noise.red == 1 || noise.red == 2) {
      spec_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                           norm = TRUE, dB = NULL, scale = FALSE, wn = win.fun,
                           noise.reduction = noise.red)
    }else if (noise.red == 0) {
      spec_left <- spectro(left, f = samplingrate, wl = w.len, plot = FALSE,
                           norm = TRUE, dB = NULL, scale = FALSE, wn = win.fun)
    }


    specA_left <- spec_left$amp

    min.freq1k = min.freq/1000
    max.freq1k = max.freq/1000

    which_min.freq <- which(abs(spec_left$freq - min.freq1k)==min(abs(spec_left$freq - min.freq1k)))
    which_max.freq <- which(abs(spec_left$freq - max.freq1k)==min(abs(spec_left$freq - max.freq1k)))

    specA_left <- specA_left[which_min.freq:which_max.freq,]
    rm(spec_left)

    rm(left)

    #LEFT CHANNEL
    specA_rows <- dim(specA_left)[1]
    specA_cols <- dim(specA_left)[2]


    fl <- rep(NA, specA_rows)
    delta_fl <- ( max.freq - min.freq ) / specA_rows
    delta_tk <- (length(wave@left)/wave@samp.rate) / specA_cols

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

    aciOutputMono <- tibble(value = ACI_tot_left)

     # Add metadata columns
    aciOutputMono <- aciOutputMono |>
      add_column(freq.res = freq.res,
                 w.len = w.len,
                 win.fun = win.fun,
                 j = j,
                 minf = min.freq,
                 maxf = max.freq,
                 noise.red = "none",
                 rm.offset = rm.offset,
                 samp = samplingrate,
                 freqres = freq_per_row,
                 nyq = nyquist_freq,
                 duration = duration,
                 channels = "mono")

    aciOutputMono <- aciOutputMono |>
      add_column(index = "aci", .before = "value")

    return(aciOutputMono)


  }

}
