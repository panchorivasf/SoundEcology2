adi <- function(wave, wlen = 512,
                wfun = "hanning",
                min_freq = 0,
                max_freq = 10000,
                nbands = 10,
                db_threshold = 5,
                normspec = FALSE,
                noisered = 2,
                rmoff = TRUE,
                props = TRUE,
                entropy = 1){


  require(tuneR)
  require(seewave)
  require(tidyr)
  require(dplyr)

  # Store the frequency step (band "height") # NEW 09/25/2023 Francisco Rivas
  freq_step <- (max_freq - min_freq)/nbands


  # Store the decibel threshold as numeric
  db_threshold <- as.numeric(db_threshold)

  # Test arguments
  # Check if the maximum frequency, decibel threshold,
  # and frequency step arguments are numbers:
  if (is.numeric(as.numeric(max_freq))){
    max_freq <- as.numeric(max_freq)
  } else{
    stop(" max_freq is not a number.")
  }

  if (is.numeric(as.numeric(db_threshold))){
    db_threshold <- as.numeric(db_threshold)
  } else{
    stop(" db_threshold is not a number.")
  }

  if (is.numeric(as.numeric(freq_step))){
    freq_step <- as.numeric(freq_step)
  } else{
    stop(" freq_step is not a number.")
  }

  # Function that gets the proportion of values higher than the
  # db threshold in a specific frequency band. The frequencies are in Hz
  getscore <- function(spectrum, minf, maxf, db, freq_row){
    miny<-round((minf)/freq_row) # the minimum frequency of the frequency band
    maxy<-round((maxf)/freq_row) # the maximum frequency of the frequency band

    subA = spectrum[miny:maxy,] # a subset f of the amplitude matrix (i.e. a single frequency band)


    minspec <- round(0/freq_row) # lower end of the spectrogram defined by min_freq
    maxspec <- round(max_freq/freq_row) # upper end of the spectrogram defined by max_freq
    speclims <- spectrum[minspec:maxspec,] # the spectrogram with range defined by min_freq and max_freq

    # Calculate the proportion of cells in f that are higher than the dB threshold
    if(entropy == 1){ #original ADI proportion calculation (within frequency band)
      index1 <- length(subA[subA>db]) / length(subA)

    }else if(entropy == 2){ #alternative 2:'true Shannon', over the user-defined spectrogram range
      # (cells above the energy threshold across the spectrogram)
      index1 <- length(subA[subA>db]) / length(speclims[speclims>db])

    }else if(entropy == 3){
      # True Shannon's over the whole spectrogram
      # (over all the pixels above the threshold up to the Nyquist frequency)
      index1 <- length(subA[subA>db]) / length(spectrum[spectrum>db])
    }
    return(index1)
  }


  # Save the denominator used in the proportion calculation
  if(entropy == 1){
    propden <- "within band"
  }else if(entropy == 2){
    propden <- "max_freq"
  }else if(entropy == 3){
    propden <- "nyquist"
  }


  #Get sampling rate
  samplingrate <- wave@samp.rate

  #Get Nyquist frequency in Hz
  nyquist_freq <- samplingrate/2

  # In the original ADI from "soundecology", the frequency bin width was fixed to 10 Hz,
  # which is equivalent of a window length of 4,800 samples for a recording with
  # a sampling rate of 48 kHz. This creates a blurred spectrogram, distorting the
  # real distribution of energy in time.

  # window length for the spectro and spec functions (legacy)
  # to keep each row every 10Hz  (legacy)
  # Frequencies and seconds covered by each (legacy)
  # freq_per_row = 10 (legacy)
  # wlen = samplingrate/freq_per_row (legacy)

  # Calculate frequency resolution (i.e., frequency bin width)
  freq_per_row = samplingrate/wlen

  # Adding 1 if wlen is an odd number (new behavior in seewave)
  # fix by JSueur
  if(wlen%%2 == 1) {wlen <- wlen+1}

  #Stereo file
  if (wave@stereo == TRUE) {

    cat("2 channels \n")

    left<-channel(wave, which = c("left"))
    right<-channel(wave, which = c("right"))
    rm(wave)

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


    # Generate normalized spectrogram if normspec = TRUE
    if(normspec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noisered == 1 || noisered == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen,
                              wn = wfun, noisereduction = noisered,
                              plot = FALSE)$amp
        specA_right <- spectro(right, f = samplingrate, wl = wlen,
                               wn = wfun, noisereduction = noisered,
                               plot = FALSE)$amp
      }else if (noisered == 3) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen,
                              wn = wfun, noisereduction = NULL,
                              plot = FALSE)$amp
        specA_right <- spectro(right, f = samplingrate, wl = wlen,
                               wn = wfun, noisereduction = NULL,
                               plot = FALSE)$amp
      }


      rm(left, right)


    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")

      # if(!is.null(noisered)){
      if(noisered == 1 || noisered == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noisereduction = noisered)$amp
        specA_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE,
                               norm=FALSE,dB=NULL,unit="power",
                               noisereduction = noisered)$amp
      }else if(noisered == 3){
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power")$amp
        specA_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE,
                               norm=FALSE,dB=NULL,unit="power")$amp
      }


      # Transform to decibels
      specA_left <- 10*log10(specA_left^2)
      specA_right <- 10*log10(specA_right^2)

      rm(left, right)

    }

    # The maximum frequency must be equal or smaller than the Nyquist frequency
    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
      max_freq <- nyquist_freq
    }

    # Set the frequency bands (Fran's comment)
    Freq <- seq(from = min_freq, to = max_freq - freq_step, by = freq_step)


    #LEFT CHANNEL

    # Create an 'empty' (NA-filled) vector to store the Score (Fran's comment)
    Score <- rep(NA, length(Freq))

    # Calculate the proportion-of-energy score (i.e., the proportion of cells which value is higher than the dB threshold),
    # looping over each frequency band (Fran's comment)
    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), db_threshold, freq_per_row)
    }

    # Store the score vector as 'left_vals' (as a copy) (Fran's comment)
    left_vals = Score

    # Create a new empty numeric object (Fran's comment)
    Score1 = 0
    # Loop over the score vector and calculate the Shannon's Entropy

    for (i in 1:length(Freq)) {
      Score1 = Score1 + (-Score[i] * log(Score[i]+ 0.000001))
    }

    Score_left = Score1


    #RIGHT CHANNEL

    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_right, Freq[j], (Freq[j] + freq_step), db_threshold, freq_per_row)
    }

    right_vals = Score

    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (-Score[i] * log(Score[i]+ 0.000001))
    }

    Score_right = Score1

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


    adiOutputStereo <- tibble(value_l = left_adi_return,
                              value_r = right_adi_return)

    adiOutputStereo <- adiOutputStereo %>%
      add_column(value_avg = ((adiOutputStereo$value_l+adiOutputStereo$value_r)/2), .after = "value_r")

    # Add information about noise reduction procedure
    if(noisered == 1){
      noise <- "rows"
    } else if(noisered == 2){
      noise <- "columns"
    } else {
      noise <- "none"
    }

    # Add metadata columns
    adiOutputStereo <- adiOutputStereo %>%
      add_column(wlen = wlen,
                 wfun = wfun,
                 dbth = db_threshold,
                 minf = min_freq,
                 maxf = max_freq,
                 nbands = nbands,
                 norm = normspec,
                 noisered = noise,
                 rmoff = rmoff,
                 propden = propden,
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
      adiOutputStereo <- bind_cols(adiOutputStereo, proportions)

      cat("Reporting ADI for 2 channels plus metadata and energy proportions per frequency band. \n")
    } else {

      cat("Reporting ADI for 2 channels plus metadata.")

    }

    adiOutputStereo <- adiOutputStereo %>%
      add_column(index = "adi", .before = "value_l") %>%

      return(adiOutputStereo)


  } else

  {
    # MONO
    cat("1 channel \n")

    left<-channel(wave, which = c("left"))
    rm(wave)

    # Remove DC offset
    if(rmoff == TRUE){
      cat("Removing DC offset...\n")
      left <- rmoffset(left)
    }


    if(noisered == 1){
      cat("Applying noise reduction filter (subtract median amplitude) to each row...\n")
    }else if(noisered == 2){
      cat("Applying a noise reduction filter (subtract median amplitude) to each column...\n")
    }

    # Generate normalized spectrogram if normspec = TRUE
    if(normspec == TRUE){

      cat("Using normalized spectrograms.\n\n")

      if(noisered == 1 || noisered == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                              noisereduction = noisered)$amp
      }else if (noisered == 3) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE)$amp
      }
      rm(left)

    }else{
      # Without normalizing the spectrogram
      cat("Using raw amplitude values (no spectrogram normalization)...\n\n")
      if(noisered == 1 || noisered == 2) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noisereduction = noisered)$amp
      }else if (noisered == 3) {
        specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,
                              norm=FALSE,dB=NULL,unit="power",
                              noisereduction = noisered)$amp
      }
      # Transform to decibels
      specA_left <- 10*log10(specA_left^2)

      rm(left)

    }


    if (max_freq > nyquist_freq) {
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep=""))
      max_freq <- nyquist_freq
    }

    Freq<-seq(from = min_freq, to = max_freq - freq_step, by = freq_step)

    Score <- rep(NA, length(Freq))

    for (j in 1:length(Freq)) {
      Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), db_threshold, freq_per_row)
    }

    left_vals = Score

    Score1 = 0
    for (i in 1:length(Freq)) {
      Score1 = Score1 + (-Score[i] * log(Score[i]+0.00001))
    }

    Score_left = Score1

    left_adi_return = round(Score_left, 3)



    left_bandvals_return <- rep(NA, length(Freq))
    # right_bandvals_return <- rep(NA, length(Freq))
    left_bandrange_return <- rep(NA, length(Freq))
    # right_bandrange_return <- rep(NA, length(Freq))


    for (j in seq(length(Freq), 1, by = -1)) {
      left_bandvals_return[j] = round(left_vals[j], 6)
      # right_bandvals_return[j] = round(right_vals[j], 6)
      left_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
      # right_bandrange_return[j] = paste((Freq[j]/1000), "-", ((Freq[j]/1000) + (freq_step/1000)), sep = "")
    }

    left_adi_return = round(Score_left, 3)
    # right_adi_return = round(Score_right, 3)

    adiOutputMono <- tibble(value = left_adi_return)


    # Add information about noise reduction procedure
    if(noisered == 1){
      noise <- "rows"
    } else if(noisered == 2){
      noise <- "columns"
    } else if(noisered == 3) {
      noise <- "none"
    }


    # Add metadata columns
    adiOutputMono <- adiOutputMono %>%
      add_column(wlen = wlen,
                 wfun = wfun,
                 dbth = db_threshold,
                 minf = min_freq,
                 maxf = max_freq,
                 nbands = nbands,
                 norm = normspec,
                 noisered = noise,
                 rmoff = rmoff,
                 propden = propden,
                 samp = samplingrate,
                 nyq = nyquist_freq,
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

      cat(" Reporting ADI for 1 channel plus metadata and energy proportions per frequency band. \n")

    }else {

      cat("Reporting ADI for 1 channel plus metadata.")

    }
    adiOutputMono <- adiOutputMono %>%
      add_column(index = "adi", .before = "value")

    return(adiOutputMono)

  }

}

