# Calculate the Acoustic Diversity Index for all the files in a folder, with extended parameter options.
# It uses parallel processing with all but one of the available cores.
# Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024

# ARGUMENTS:
# folder: Path to the folder where the WAV files are located
# wlen: window length to compute the spectrogram
# wfun: window function (filter to deal with spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
# min_freq: minimum frequency to compute the spectrogram
# max_freq: maximum frequency to compute the spectrogram
# nbands: number of bands to split the spectrogram
# db_threshold = dB threshold to calculate energy proportions (if normspec = FALSE, set to 5 or above)
# normspec: logical. Whether to normalize the spectrogram (not recommended) or not (normalized spectrograms with different SNR are not comparable).
# noisered: numeric. Noise reduction (subtract median from the amplitude values); 1=rows, 2=columns, 3=none.
# rmoff: logical. Whether to remove DC offset before computing ADI (recommended) or not.
# props: logical. Whether to store the energy proportion values for each frequency band and channel (default).
# entropy: numeric. Indicates how the energy proportion is calculated.

# Options for entropy:
# 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion represents all the cells in the same frequency band.
# 2 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max_freq')
# 3 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.

adi_folder <- function (folder,
                        save_csv = FALSE,
                        csv_name = "adi_results.csv",
                        wlen = 512,
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

  require(doParallel)
  require(foreach)
  require(parallel)
  require(tuneR)
  require(tidyverse)
  require(seewave)
  require(lubridate)
  # require(soundecology2)

  cat("Evaluating the job...\n\n")


  setwd(folder)
  audiolist <- list.files(path=folder, pattern = ".wav|.WAV")

  fileName <- tibble(file_name = audiolist)
  nFiles <- length(audiolist)

  # Check if the first sound is stereo or mono
  sound <- readWave(audiolist[1])
  type <- ifelse(sound@stereo, "stereo", "mono")
  rm(sound)

  args_list <- list(wlen = wlen, wfun = wfun, min_freq = min_freq,
                    max_freq = max_freq, nbands = nbands, db_threshold = db_threshold,
                    normspec = normspec, noisered = noisered, rmoff = rmoff,
                    props = props, entropy = entropy)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  adi1 <- adi(sound1, args_list$wlen, args_list$wfun, args_list$min_freq,
                   args_list$max_freq, args_list$nbands, args_list$db_threshold,
                   args_list$normspec, args_list$noisered, args_list$rmoff,
                   args_list$props, args_list$entropy)

  tibble(file_name = "filename") %>% bind_cols(adi1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(adi1)

  # Declare the number of cores to be used (all but one of the available cores)
  cores <- detectCores() - 1 # Leave one core free
  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Setup parallel back-end
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)

  cat("Process start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", cores, "cores... \n")


  # Start loop
  results <- foreach(file = audiolist, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate ADI4 index and keep its default output columns
                       adi <- soundecology2::adi(sound, wlen = args_list$wlen, wfun = args_list$wfun,
                                       min_freq = args_list$min_freq, max_freq = args_list$max_freq,
                                       nbands = args_list$nbands, db_threshold = args_list$db_threshold,
                                       normspec = args_list$normspec, noisered = args_list$noisered,
                                       rmoff = args_list$rmoff, props = args_list$props,
                                       entropy = args_list$entropy)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(adi)


                     }


  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)


  stopCluster(cl)

  if(save_csv == TRUE){
    write.csv(resultsWithMetadata, csv_name, row.names = FALSE)
  }

  cat(paste("Analysis complete!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n"))

  return(resultsWithMetadata)

}

# Example
# path <-  "F:\\2024_eclipse_hydrophone\\test"
# adi_test <- adi_folder(folder=path)


