#' Normalized Difference Soundscape Index
#' @description
#' Normalized Difference Soundscape Index (NDSI) from REAL and Kasten, et al. 2012. The NDSI seeks to "estimate the level of anthropogenic disturbance on the soundscape by computing the ratio of human-generated (anthrophony) to biological (biophony) acoustic components found in field collected sound samples" (Kasten, et al. 2012).
#'
#' @param wave an object of class Wave imported with the \emph{readWave} function of the \emph{tuneR} package.
#' @param w.len numeric. The window length for the FFT (sampling rate / frequency resolution).
#' @param anthro.min minimum value of the range of frequencies of the anthrophony.
#' @param anthro.max maximum value of the range of frequencies of the anthrophony.
#' @param bio.min minimum value of the range of frequencies of the biophony.
#' @param bio.max maximum value of the range of frequencies of the biophony.
#' @param rm.offset logical. Whether to remove the DC offset.
#'
#' @return a wide format tibble with NDSI values per channel (if stereo), parameters used and audio metadata
#' @export
#'
#' @import seewave
#' @import tuneR
#' @import oce
#' @import pracma
#' @import ineq
#'
#' @examples ndsi(tropicalsound)
ndsi <- function(wave,
                w.len = 512,
                anthro.min = 1000,
                anthro.max = 2000,
                bio.min = 2000,
                bio.max = 11000,
                rm.offset = TRUE){

  #test arguments
  if (is.numeric(as.numeric(w.len))){
    w.len <- as.numeric(w.len)
  } else{
    stop(" freq.res is not a number.")
  }

  if (is.numeric(as.numeric(anthro.min))){
    anthro.min <- as.numeric(anthro.min)
  } else{
    stop(" anthro.min is not a number.")
  }

  if (is.numeric(as.numeric(anthro.max))){
    anthro.max <- as.numeric(anthro.max)
  } else{
    stop(" anthro.max is not a number.")
  }

  if (is.numeric(as.numeric(bio.min))){
    bio.min <- as.numeric(bio.min)
  } else{
    stop(" bio.min is not a number.")
  }

  if (is.numeric(as.numeric(bio.max))){
    bio.max <- as.numeric(bio.max)
  } else{
    stop(" bio.max is not a number.")
  }
  
  #Get sampling rate
  samplingrate <- wave@samp.rate
  duration <- length(wave@left)/wave@samp.rate
  freq.res <- samplingrate/w.len



  # Adding 1 if w.len is an odd number (new behavior in seewave)
  # fix by JSueur
  if(w.len%%2 == 1) {w.len <- w.len+1}


  #Some general values
  hz_interval = anthro.max - anthro.min


  #Get Nyquist frequency in Hz
  nyquist_freq <- (samplingrate/2)

  #Check errors
  if (bio.max > nyquist_freq) {
    stop(paste("The maximum frequency of biophony (", bio.max, " Hz) can not be higher than the maximum frequency of the file (", nyquist_freq, " Hz)\n\n Change the value of bio.max to less than ", nyquist_freq, "\n\n", sep = ""))
  }

  if (anthro.max > bio.min) {
    stop(paste("The maximum frequency of anthrophony (", anthro.max, " Hz) can not be higher than the minimum frequency of biophony (", bio.min, " Hz)\n\n Change the value of anthro.max to equal or less than bio.min\n\n", sep = ""))
  }

  if (anthro.max < anthro.min) {
    stop(paste("The minimum frequency of anthrophony (", anthro.min, " Hz) can not be higher than the maximum frequency of anthrophony (", anthro.max, " Hz)\n\n Change the value of anthro.min to less than anthro.max\n\n", sep = ""))
  }

  if (bio.max < bio.min) {
    stop(paste("The minimum frequency of biophony (", bio.min, " Hz) can not be higher than the maximum frequency of biophony (", bio.max, " Hz)\n\n Change the value of anthro.min to less than anthro.max\n\n", sep = ""))
  }


  #Stereo file
  if (wave@stereo == TRUE) {

    left <- channel(wave, which = c("left"))
    right <- channel(wave, which = c("right"))
    
    # Remove DC offset
    if(rm.offset == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left, output = "Wave")
      right <- rmoffset(right, output = "Wave")
    }
    
    rm(wave)

    cat("\n Calculating NDSI on a stereo file... \n")

    #LEFT CHANNEL
    left1 <- as.vector(cutw(left, from = 0, to = length(left@left) / left@samp.rate))
    left2 <- data.frame(matrix(NA, nrow = samplingrate, ncol = floor(duration)))

    for (i in 0:(floor(duration)-1)){
      j <- i+1
      start1 <- (i * samplingrate)+ 1
      end <- start1 + samplingrate - 1
      left2[, j] <- left1[start1:end]
    }

    left3 <- data.frame(matrix(NA, nrow = w.len/2, ncol = floor(duration)))

    left4 <- apply(left2, 2, pwelch, fs = samplingrate, nfft = w.len, plot = FALSE)

    for (i in 1:floor(duration)){
      left3[, i] <- left4[[i]]$spec
    }

    specA_left <- apply(left3, 1, mean)
    specA_rows <- length(specA_left)

    freq_per_row <- specA_rows/nyquist_freq

    anthro_vals_range <- anthro.max - anthro.min
    bio_vals_range <- bio.max - bio.min
    bio_bins <- round(bio_vals_range/hz_interval)

    anthro_bins <- rep(NA, round(anthro_vals_range/hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range/hz_interval))

    anthro.min_row <- round(anthro.min * freq_per_row)
    anthro.max_row <- round(anthro.max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range/length(bio_bins))
    bio.min_row <- round(bio.min * freq_per_row)
    bio.max_row <- bio.min_row + bio_step_range


    #Get the area for each bin of anthrophony and biophony
    #Anthrophony
    for (i in 1:length(anthro_bins)){
      anthro_bins[i] <- trapz(specA_left[anthro.min_row:anthro.max_row])
    }

    #Biophony
    for (i in 1:length(bio_bins)){

      if (bio.max_row >= specA_rows){
        bio.max_row <- specA_rows
      }

      bio_bins[i] <- trapz(specA_left[bio.min_row:bio.max_row])

      bio.min_row <- bio.min_row + bio_step_range
      bio.max_row <- bio.max_row + bio_step_range
    }

    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    #Normalize
    freqbins = freqbins / norm(as.matrix(freqbins), "F")

    #All bins
    freqbins.SumAll <- sum(freqbins)
    #All biophony bins
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    #Single anthrophony bin
    freqbins.Anthro <- freqbins[1]

    #Result
    NDSI_left <- (freqbins.SumBio - freqbins.Anthro) / (freqbins.SumBio + freqbins.Anthro)
    biophony_left <- freqbins.SumBio
    anthrophony_left <- freqbins.Anthro

    #Right channel
    right1 <- as.vector(cutw(right, from = 0, to = length(right@left) / right@samp.rate))

    right2 <- data.frame(matrix(NA, nrow = samplingrate, ncol = floor(duration)))

    for (i in 0:(floor(duration)-1)){
      j <- i+1
      start1 <- (i * samplingrate)+ 1
      end <- start1 + samplingrate - 1
      right2[, j] <- right1[start1:end]
    }

    right3 <- data.frame(matrix(NA, nrow = w.len/2, ncol = floor(duration)))

    right4 <- apply(right2, 2, pwelch, fs = samplingrate, nfft = w.len, plot = FALSE)

    for (i in 1:floor(duration)){
      right3[, i] <- right4[[i]]$spec
    }

    specA_right <- apply(right3, 1, mean)
    specA_rows <- length(specA_right)

    freq_per_row <- specA_rows / nyquist_freq

    anthro_vals_range <- anthro.max - anthro.min
    bio_vals_range <- bio.max - bio.min
    bio_bins <- round(bio_vals_range/hz_interval)

    anthro_bins <- rep(NA, round(anthro_vals_range / hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range / hz_interval))

    anthro.min_row <- round(anthro.min * freq_per_row)
    anthro.max_row <- round(anthro.max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range/length(bio_bins))
    bio.min_row <- round(bio.min * freq_per_row)
    bio.max_row <- bio.min_row + bio_step_range

    #Anthrophony
    for (i in 1:length(anthro_bins)){
      anthro_bins[i] <- trapz(specA_right[anthro.min_row:anthro.max_row])
    }

    #Biophony
    for (i in 1:length(bio_bins)){

      if (bio.max_row >= specA_rows){
        bio.max_row <- specA_rows
      }

      bio_bins[i] <- trapz(specA_right[bio.min_row:bio.max_row])

      bio.min_row <- bio.min_row + bio_step_range
      bio.max_row <- bio.max_row + bio_step_range
    }

    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    freqbins = freqbins / norm(as.matrix(freqbins), "F")


    freqbins.SumAll <- sum(freqbins)
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    freqbins.Anthro <- freqbins[1]

    NDSI_right <- (freqbins.SumBio - freqbins.Anthro) / (freqbins.SumBio + freqbins.Anthro)
    biophony_right <- freqbins.SumBio
    anthrophony_right <- freqbins.Anthro

    ndsiOutputStereo <- tibble(index = "ndsi",
                               value_l = NDSI_left,
                               value_r = NDSI_right,
                               bio_l = biophony_left,
                               bio_r = biophony_right,
                               anthro_l = anthrophony_left,
                               anthro_r = anthrophony_right,
                               anthro.min = anthro.min,
                               anthro.max = anthro.max,
                               bio.min = bio.min,
                               bio.max = bio.max,
                               freq.res = freq.res,
                               w.len = w.len,
                               samprate = samplingrate,
                               freqres = freq_per_row,
                               nyquist = nyquist_freq,
                               duration = duration,
                               channels = "stereo")

    ndsiOutputStereo <- ndsiOutputStereo %>%
      add_column(value_avg = ((ndsiOutputStereo$value_l+ndsiOutputStereo$value_r)/2), .after = "value_r") %>%
      add_column(bio_avg = ((ndsiOutputStereo$bio_l+ndsiOutputStereo$bio_r)/2), .after = "bio_r") %>%
      add_column(anthro_avg = ((ndsiOutputStereo$anthro_l+ndsiOutputStereo$anthro_r)/2), .after = "anthro_r")


    return(ndsiOutputStereo)

  } else
    # MONO FILE
  {
    # cat("\n This is a mono file.\n")

    #LEFT CHANNEL
    left<-channel(wave, which = c("left"))
    rm(wave)

    cat("\n Calculating NDSI on a mono file... \n")

    left1 <- as.vector(cutw(left, from = 0, to = length(left@left) / left@samp.rate))
    left2 <- data.frame(matrix(NA, nrow = samplingrate, ncol = floor(duration)))

    #divide in 1-second segment to avoid a long for loop in pwelch
    for (i in 0:(floor(duration)-1)){
      j <- i+1
      start1 <- (i * samplingrate)+ 1
      end <- start1 + samplingrate - 1
      left2[, j] <- left1[start1:end]
    }

    left3 <- data.frame(matrix(NA, nrow = w.len/2, ncol = floor(duration)))

    left4 <- apply(left2, 2, pwelch, fs = samplingrate, nfft = w.len, plot = FALSE)

    for (i in 1:floor(duration)){
      left3[, i] <- left4[[i]]$spec
    }

    specA_left <- apply(left3, 1, mean)

    specA_rows <- length(specA_left)

    freq_per_row <- specA_rows / nyquist_freq

    anthro_vals_range <- anthro.max - anthro.min
    bio_vals_range <- bio.max - bio.min
    bio_bins <- round(bio_vals_range / hz_interval)

    anthro_bins <- rep(NA, round(anthro_vals_range / hz_interval))
    bio_bins <- rep(NA, round(bio_vals_range / hz_interval))

    anthro.min_row <- round(anthro.min * freq_per_row)
    anthro.max_row <- round(anthro.max * freq_per_row)
    bio_step_range <- freq_per_row * (bio_vals_range / length(bio_bins))
    bio.min_row <- round(bio.min * freq_per_row)
    bio.max_row <- bio.min_row + bio_step_range

    #Anthrophony
    for (i in 1:length(anthro_bins)){
      anthro_bins[i] <- trapz(specA_left[anthro.min_row:anthro.max_row])
    }

    #Biophony
    for (i in 1:length(bio_bins)){

      if (bio.max_row >= specA_rows){
        bio.max_row <- specA_rows
      }

      bio_bins[i] <- trapz(specA_left[bio.min_row:bio.max_row])

      bio.min_row <- bio.min_row + bio_step_range
      bio.max_row <- bio.max_row + bio_step_range
    }

    freqbins <- rep(NA, sum(length(anthro_bins), length(bio_bins)))
    freqbins <- c(anthro_bins, bio_bins)
    freqbins = freqbins / norm(as.matrix(freqbins), "F")


    freqbins.SumAll <- sum(freqbins)
    freqbins.SumBio <- sum(freqbins[2:length(freqbins)])
    freqbins.Anthro <- freqbins[1]

    NDSI_left <- (freqbins.SumBio - freqbins.Anthro) / (freqbins.SumBio + freqbins.Anthro)
    biophony_left <- freqbins.SumBio
    anthrophony_left <- freqbins.Anthro
    biophony_right <- NA
    anthrophony_right <- NA

    #Right channel
    NDSI_right = NA

    ndsiOutputMono <- tibble(index = "ndsi",
                             value = NDSI_left,
                             bio = biophony_left,
                             anthro = anthrophony_left,
                             anthro.min = anthro.min,
                             anthro.max = anthro.max,
                             bio.min = bio.min,
                             bio.max = bio.max,
                             freq.res = freq.res,
                             w.len = w.len,
                             samprate = samplingrate,
                             freqres = freq_per_row,
                             nyquist = nyquist_freq,
                             duration = duration,
                             channels = "mono")

    return(ndsiOutputMono)

  }

}
