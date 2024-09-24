frequency_dependent_acoustic_diversity <- function(soundfile,noisefile=NULL,NEM=2,min_freq=200, max_freq = 10000, threshold_fixed = -50, freq_step = 1000, gamma = 13)
{
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
    stop("Parameter Error")
  }
  #finish input parameters validation 
  
  #histogram based noise level estimation embedded here
  noise_estimation<-function(spec,max_f,freq_row)
  {          
    maxy <- round((max_f)/freq_row)
    noise_db<- rep(NA, maxy) 
    for (j in 1:maxy) {
      specA_hist<-hist(spec[j, ],40)        
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
  if (soundfile@stereo == TRUE) {  
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    left<-channel(soundfile, which = c("left"))  
    right<-channel(soundfile, which = c("right"))
      
      
    #matrix of values
    cat("\n Calculating index. Please wait... \n\n")
    specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
    specA_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
      
    specA_left<-specA_left^2
    specA_left<-10*log10(specA_left)
    specA_right<-specA_right^2
    specA_right<-10*log10(specA_right)
      
    if (max_freq > nyquist_freq) {  
      cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
      max_freq <- nyquist_freq    
     }
      
     if (NEM==1){     #if noisefile should be provided
        if (is.null(noisefile)){     #parameter validation
          cat(paste("\n WARNING: \n No input noise file .\n Then histogram-based noise-estimation algorithms will be usde to estimate the noise of the recording. \n\n"))
          NEM=2
        }else{
          if(length(noisefile)/noisefile@samp.rate<5){
            stop("\n Warning: \n The input noise file length is insufficient. \n\n")
          }
          
        noise_left<-channel(noisefile, which = c("left"))
        noise_right<-channel(noisefile, which = c("right"))
        noise_specA_left <- spectro(noise_left, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
        noise_specA_right <- spectro(noise_right, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
          
        noise_specA_left<-noise_specA_left^2
        noise_specA_right<-noise_specA_right^2
          
        noise_specA_left<-apply(noise_specA_left,1,mean)
        noise_db_left<-10*log10(noise_specA_left)
          
        noise_specA_right<-apply(noise_specA_right,1,mean)
        noise_db_right<-10*log10(noise_specA_right)
        } 
      }
      
    if (NEM==2){  #if noisefile is omitted
        
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
      
    # 		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
    # 		cat("Frequency range (Hz), left channel proportion, right channel proportion\n")
    for (j in seq(length(Freq), 1, by = -1)) {
      # 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), ",", round(right_vals[j],6), "\n", sep=""))
      left_bandvals_return[j] = round(left_vals[j], 6)   
      right_bandvals_return[j] = round(right_vals[j], 6)  
      left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")  
      right_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "") 
    }
      
      # 		cat("\n Plot of proportions in each band: \n\n")
      # 		cat("  Left channel\n")
      # 		cat("   Freq. range (Hz) |--------------------|\n")
      
      #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
    for (j in seq(length(Freq), 1, by = -1)) {
      this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
      this_row_size <- nchar(this_row_name)
      this_row_space <- 17 - this_row_size
        
      this_row_spaces = ""
        
      for (f in seq(1, this_row_space, by = 1)) {
        this_row_spaces = paste(this_row_spaces, " ", sep = "")
      }
    }
      
      #precision control
      left_adi_return = round(Shannon_left, 6)   
      right_adi_return = round(Shannon_right, 6)  
 
      cat("  Frequency-dependent Acoustic Diversity Index: \n")
      cat(paste("   Left channel: ", left_adi_return, "\n", sep = ""))
      cat(paste("   Right channel: ", right_adi_return, "\n", sep = ""))
  }
   else {      #if input soundfile is monochannel
      cat("\n This is a mono file.\n")  
      
      #matrix of values
      cat("\n Calculating index. Please wait... \n\n")
      specA_left <- spectro(soundfile, f = samplingrate, wl = wlen, plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
      specA_left<-specA_left^2
      specA_left<-10*log10(specA_left)
      
      if (max_freq > nyquist_freq) {
        cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep=""))
        max_freq <- nyquist_freq
      }
      
      if (NEM==1){  
        if (is.null(noisefile)){
          cat(paste("\n WARNING: \n No input noise file .\n Then histogram-based noise-estimation algorithms will be usde to estimate the noise of the data. \n\n"))
          NEM=2
        }else{
          if(length(noisefile)/noisefile@samp.rate<5){
            stop("\n Warning: \n The input noise signal length is insufficient. \n\n")
          }
          
          noise_left<-noisefile
          noise_specA_left <- spectro(noise_left, f=samplingrate, wl = wlen,plot = FALSE,norm=FALSE,dB=NULL,unit="power")$amp
          
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
      
      # 		cat(" ==============================================\n")
      # 		cat(paste(" Results (with a dB threshold of ", db_threshold, ")\n\n", sep=""))
      # 		
      # 		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
      # 		cat("Frequency range (Hz), proportion\n")
      
      #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
      left_bandvals_return <- rep(NA, length(Freq))
      right_bandvals_return <- rep(NA, length(Freq))
      left_bandrange_return <- rep(NA, length(Freq))
      right_bandrange_return <- rep(NA, length(Freq))
      for (j in seq(length(Freq), 1, by = -1)) {
        # 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), "\n", sep=""))
        left_bandvals_return[j] = round(left_vals[j], 6)
        left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
      }
      
      # 		cat("\n Plot of proportions in each band: \n\n")
      # 		cat("   Freq. range (Hz) |--------------------|\n")
      
      #printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
      for (j in seq(length(Freq), 1, by = -1)) {
        this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
        this_row_size <- nchar(this_row_name)
        this_row_space <- 17 - this_row_size
        
        this_row_spaces = ""
        
        for (f in seq(1, this_row_space, by = 1)) {
          this_row_spaces = paste(this_row_spaces, " ", sep = "")
        }
        
      }
      
      
      cat("\n  Frequency-dependent Acoustic Diversity Index: ")
      right_adi_return = NA
    
        cat(paste(round(Shannon_left, 6), "\n", sep = ""))
        left_adi_return = round(Shannon_left, 6)
    
  }
  invisible(list(fadi_left = left_adi_return, fadi_right = right_adi_return, left_band_values = left_bandvals_return, right_band_values = right_bandvals_return, left_bandrange_values = left_bandrange_return, right_bandrange_values = right_bandrange_return))
  #return output
}

