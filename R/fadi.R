#' Frequency-dependent Acoustic Diversity Index
#' @description
#' The Frequency-dependent Acoustic Diversity Index by Xu et al. (2023) obtains 
#' a floating noise profile before calculating the Acoustic Diversity Index and 
#' it doesn't use normalized spectrogram. Alternatively it can take a noise 
#' sample to reduce noise from the analyzed files.
#' @param soundfile A wave object imported with readWave().
#' @param noise_file An R object of class Wave containing noise-only information 
#' if needed. Default = NULL.
#' @param NEM Numeric. Options are 1 or 2.
#' When NEM = 1, floating thresholds are estimated based on noise_file.
#' When NEM = 2, floating thresholds are calculated based on sound file using an
#' automatic noise level estimation method (median of each row in the 
#' spectrogram). Default = 2.
#' @param min_freq Minimum frequency in Hertz when calculating the global 
#' threshold. Default = 200.
#' @param max_freq Maximum frequency in Hertz when calculating the FADI value. 
#' Default = 10000.
#' @param threshold_fixed A negative number in dB for calculating the global 
#' threshold. Default = −50.
#' @param freq_step Bandwidth of each frequency band, in Hertz. Default = 1000.
#' @param gamma A positive number in dB for calculating the floating thresholds. 
#' Default = 13.
#' @param props Logical; if TRUE, the energy proportion values for each 
#' frequency ban and channel are added to the output tibble. Default = TRUE.
#'
#' @return A tibble with the FADI value per channel, energy proportions, 
#' metadata, and parameters used.
#' @export
#'
#' @importFrom tuneR readWave
#' @import seewave
#' @import tidyverse
#'
#' @examples fadi(tropicalsound)
#'
#' @details
#' Modified version of the Frequency-dependent Acoustic Diversity Index by Xu 
#' et al. (2023).FADI was introduced in: 
#' https://www.sciencedirect.com/science/article/pii/S1470160X23010828.
#' This version returns a wide format (one row per audio file) tibble as output 
#' instead of a nested list.
#' To see the original version as in the paper, use the 
#' \emph{frequency_dependent_acoustic_diversity()} function.
fadi <- function(soundfile,
                 noise_file= NULL,
                 NEM = 2,
                 min_freq = 200,
                 max_freq = 10000,
                 threshold_fixed = -50,
                 freq_step = 1000,
                 gamma = 13,
                 props = TRUE){

  # require(tuneR)
  # require(seewave)
  # require(tidyverse)
  #start input parameters validation

  threshold_fixed <- as.numeric(threshold_fixed)
  NEM <- as.numeric(NEM)

  #test arguments
  if (is.numeric(as.numeric(max_freq))){
    max_freq <- as.numeric(max_freq)
  } else{
    stop(" max_freq is not a number.")
  }

  if (is.numeric(as.numeric(threshold_fixed))){
    threshold_fixed<- as.numeric(threshold_fixed)
  } else{
    stop(" threshold_fixed is not a number.")
  }

  if (is.numeric(as.numeric(freq_step))){
    freq_step <- as.numeric(freq_step)
  } else{
    stop(" freq_step is not a number.")
  }

  if (NEM!=1 & NEM!=2){
    stop("Parameter Error. NEM should be either 1 or 2")
  }
  #finish input parameters validation

  #histogram based noise level estimation embedded here
  noise_estimation<-function(spec,max_f,freq_row)
  {
    maxy <- round((max_f)/freq_row)
    noise_db<- rep(NA, maxy)
    for (j in 1:maxy) {
      specA_hist<-hist(spec[j, ],40, plot = FALSE)
      counts<-unlist(specA_hist[2])
      mids<-unlist(specA_hist[4])
      noise_db[j]<-mids[which.max(counts)]
    }
    return(noise_db)
  }

  #floating threshold calculation at each frequency bin
  db_threshold <- function(spec,noise,min_f,max_f,freq_row,threshold_db,alpha)
  {
    miny <- round((min_f)/freq_row)
    maxy <- round((max_f)/freq_row)
    the_global_threshold=max(spec[miny:maxy, ])+threshold_db
    db_threshold_floating<- rep(NA, maxy)
    db_threshold_narrowband<- rep(NA, maxy)
    for (j in 1:maxy) {
      db_threshold_floating[j] = noise[j]+alpha
      db_threshold_narrowband[j]=max(db_threshold_floating[j],the_global_threshold)
    }
    return(db_threshold_narrowband)
  }

  #proportion calculation in each frequency band
  getscore <- function(spectrum, minf,maxf, db, freq_row)
  {
    miny<-round((minf)/freq_row)
    maxy<-round((maxf)/freq_row)
    index1<-0
    for(m in seq(miny+1,maxy,by=1)){
      for (n in 1:ncol(spectrum)){
        if(spectrum[m,n]>db[m]){
          index1<-index1+1
        }
      }
    }
    index1 <- index1 / length(spectrum[seq(miny+1,maxy,by=1),])
    return(index1)
  }


  #start main processing procedure
  freq_per_row=10    #interval in Hertz between adjacent frequency bins
  samplingrate <- soundfile@samp.rate
  nyquist_freq <- samplingrate/2

  wlen = samplingrate/freq_per_row

  if(wlen%%2 == 1) {wlen <- wlen+1}

  # IF STEREO
  if (soundfile@stereo == TRUE) {
    # cat("\n Calculating FADI on a stereo file...\n")
    left<-channel(soundfile, which = c("left"))
    right<-channel(soundfile, which = c("right"))


    #matrix of values
    # cat("\n Calculating index. Please wait... \n\n")
    specA_left <- seewave::spectro(left, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
    specA_right <- seewave::spectro(right, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp

    specA_left<-specA_left^2
    specA_left<-10*log10(specA_left)
    specA_right<-specA_right^2
    specA_right<-10*log10(specA_right)

    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
      max_freq <- nyquist_freq
    }

    if (NEM==1){     #if noise_file should be provided
      if (is.null(noise_file)){     #parameter validation
        cat(paste("\n WARNING: \n No input noise file .\n Then histogram-based noise-estimation algorithms will be used to estimate the noise of the recording. \n\n"))
        NEM=2
      }else{
        if(length(noise_file)/noise_file@samp.rate<5){
          stop("\n Warning: \n The input noise file length is insufficient. \n\n")
        }

        noise_left<-channel(noise_file, which = c("left"))
        noise_right<-channel(noise_file, which = c("right"))
        noise_specA_left <- seewave::spectro(noise_left, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
        noise_specA_right <- seewave::spectro(noise_right, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp

        noise_specA_left<-noise_specA_left^2
        noise_specA_right<-noise_specA_right^2

        noise_specA_left<-apply(noise_specA_left,1,mean)
        noise_db_left<-10*log10(noise_specA_left)

        noise_specA_right<-apply(noise_specA_right,1,mean)
        noise_db_right<-10*log10(noise_specA_right)
      }
    }

    if (NEM==2){  #if noise_file is omitted

      if(length(soundfile)/soundfile@samp.rate < 30){
        stop("\n WARNING: \n The input soundfile length is insufficient. \n\n")
      }

      noise_db_left<-noise_estimation(specA_left,max_freq,freq_per_row)
      noise_db_right<-noise_estimation(specA_right,max_freq,freq_per_row)
    }

    Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)

    #left channel
    threshold_left=db_threshold(specA_left,noise_db_left,min_freq,max_freq,freq_per_row,threshold_fixed,gamma)

    Score <- rep(NA, length(Freq))
    Scorez <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Scorez[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), threshold_left, freq_per_row)
    }

    Score= Scorez/sum(Scorez)

    left_vals = Score

    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 0.0000001))
    }

    Shannon_left <- -Score1   #FADI result for left channel

    #right channel
    threshold_right=db_threshold(specA_right,noise_db_right,min_freq,max_freq,freq_per_row,threshold_fixed,gamma)
    Score <- rep(NA, length(Freq))
    Scorez <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Scorez[j] = getscore(specA_right, Freq[j], (Freq[j] + freq_step), threshold_right, freq_per_row)
    }
    Score= Scorez/sum(Scorez)

    right_vals = Score

    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 0.0000001))
    }

    Shannon_right <- -Score1  #FADI result for right channel

    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    right_bandrange_return <- rep(NA, length(Freq))


    for (j in seq(length(Freq), 1, by = -1)) {
      # 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), ",", round(right_vals[j],6), "\n", sep=""))
      left_bandvals_return[j] = round(left_vals[j], 6)
      right_bandvals_return[j] = round(right_vals[j], 6)
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), " kHz", sep = "")
      right_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), " kHz", sep = "")
    }



    #precision control
    left_adi_return = round(Shannon_left, 6)
    right_adi_return = round(Shannon_right, 6)



    fadiOutputStereo <- tibble(value_l = left_adi_return,
                               value_r = right_adi_return)

    fadiOutputStereo <- fadiOutputStereo %>%
      add_column(value_avg = ((fadiOutputStereo$value_l+fadiOutputStereo$value_r)/2), .after = "value_r")



    # Print the FADI per channel
    cat(" Frequency-dependent Acoustic Diversity Index: \n")
    cat(paste(" Left channel: ", left_adi_return, "\n", sep = ""))
    cat(paste(" Right channel: ", right_adi_return, "\n", sep = ""))



    if(NEM==1){
      noisered <- "noise file"
    }else if(NEM==2){
      noisered <- "dynamic"
    }



    # Add metadata columns
    fadiOutputStereo <- fadiOutputStereo %>%
      add_column(wlen = wlen,
                 wfun = "hanning",
                 dbth = threshold_fixed,
                 minf = min_freq,
                 maxf = max_freq,
                 nbands = (max_freq/freq_step),
                 norm = "FALSE",
                 noisered = noisered,
                 rmoff = "false",
                 propden = "within band",
                 samp = samplingrate,
                 nyq = nyquist_freq,
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
      fadiOutputStereo <- bind_cols(fadiOutputStereo, proportions)

      cat("Reporting FADI for 2 channels plus metadata and energy proportions per frequency band. \n")
    } else {

      cat("Reporting FADI for 2 channels plus metadata.")

    }

    fadiOutputStereo <- fadiOutputStereo %>%
      add_column(index = "fadi", .before = "value_l") %>%
      add_column(gamma = gamma, .after = "noisered")



    return(fadiOutputStereo)


  }
  else {      #if input soundfile is monochannel
    # cat("\n Calculating FADI on a mono file...\n")


    if(NEM==1){
      noisered <- "noise file"
    }else if(NEM==2){
      noisered <- "dynamic"
    }



    #matrix of values
    # cat("\n Calculating index. Please wait... \n\n")
    specA_left <- seewave::spectro(soundfile, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
    specA_left<-specA_left^2
    specA_left<-10*log10(specA_left) # transform to decibels

    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max_freq <- nyquist_freq
    }

    if (NEM==1){
      if (is.null(noise_file)){
        cat(paste("\n WARNING: \n No input noise file .\n Then histogram-based noise-estimation algorithms will be usde to estimate the noise of the data. \n\n"))
        NEM=2
      }else{
        if(length(noise_file)/noise_file@samp.rate<5){
          stop("\n Warning: \n The input noise signal length is insufficient. \n\n")
        }

        noise_left<-noise_file
        noise_specA_left <- seewave::spectro(noise_left, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp

        noise_specA_left<-noise_specA_left^2

        noise_specA_left<-apply(noise_specA_left,1,mean)
        noise_db_left<-10*log10(noise_specA_left)
      }
    }
    if (NEM==2){

      if(length(soundfile)/soundfile@samp.rate < 30){
        stop("\n WARNING: \n The input signal length is insufficient. \n\n")
      }

      noise_db_left<-noise_estimation(specA_left,max_freq,freq_per_row)
    }

    Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)


    threshold_left=db_threshold(specA_left,noise_db_left,min_freq,max_freq,freq_per_row,threshold_fixed,gamma)

    Score <- rep(NA, length(Freq))
    Scorez <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Scorez[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), threshold_left, freq_per_row)
    }
    Score= Scorez/sum(Scorez)

    left_vals = Score

    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (Score[i] * log(Score[i] + 0.0000001))
    }


    Shannon_left <- -Score1
    Shannon_right <- NA


    #
    #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
    left_bandvals_return <- rep(NA, length(Freq))
    right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    # left_bandrange_return <- left_bandrange_return/1000
    right_bandrange_return <- rep(NA, length(Freq))
    # right_bandrange_return <- right_bandrange_return/1000
    for (j in seq(length(Freq), 1, by = -1)) {
      # 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), "\n", sep=""))
      left_bandvals_return[j] = round(left_vals[j], 6)
      # Get the frequency bands in kHz
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), " kHz", sep = "")
    }

    left_adi_return = round(Shannon_left, 6)


    bands = left_bandrange_return
    left_proportions = left_bandvals_return
    # right_proportions = right_bandvals_return


    fadiOutputMono <- tibble(value = left_adi_return)
    # Add metadata columns
    fadiOutputMono <- fadiOutputMono %>%
      add_column(wlen = wlen,
                 wfun = "hanning",
                 dbth = threshold_fixed,
                 minf = min_freq,
                 maxf = max_freq,
                 nbands = (max_freq/freq_step),
                 norm = "FALSE",
                 noisered = noisered,
                 rmoff = "false",
                 propden = "within band",
                 samp = samplingrate,
                 nyq = nyquist_freq,
                 channels = "mono")





    # Print the FADI per channel
    # cat(" Frequency-dependent Acoustic Diversity Index: \n")
    # cat(paste(" FADI: ", left_adi_return, "\n", sep = ""))




    if(props == TRUE){
      # Add data on proportions for each frequency band per channel:
      proportions = tibble(band_kHz = left_bandrange_return,
                           left = left_bandvals_return)
      # Wider format
      proportions <- pivot_wider(proportions, names_from = "band_kHz",
                                 values_from = c("left"))

      # Append to the original data frame
      fadiOutputMono <- bind_cols(fadiOutputMono, proportions)

      # cat("Reporting FADI for 1 channel plus metadata and energy proportions per frequency band. \n")
    } else {

      # cat("Reporting FADI for 1 channel plus metadata.")

    }

    fadiOutputMono <- fadiOutputMono %>%
      add_column(index = "fadi", .before = "value") %>%
      add_column(gamma = gamma, .after = "noisered")



    return(fadiOutputMono)



  }

}
