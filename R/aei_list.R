#' Acoustic Evenness Index - list input
#' @description
#' Calculates the Acoustic Evenness Index for all the files in a list, with extended parameter options.
#' It uses parallel processing with all but one of the available cores.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com)  April 2024.
#'
#' @param audiolist a list with the audio files to import.
#' @param save_csv logical. Whether to save a csv in the working directory.
#' @param csv_name character vector. When 'save_csv' is TRUE, optionally provide a file name.
#' @param wlen window length to compute the spectrogram
#' @param wfun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min_freq minimum frequency to compute the spectrogram
#' @param max_freq maximum frequency to compute the spectrogram
#' @param nbands number of bands to split the spectrogram
#' @param db_threshold dB threshold to calculate energy proportions (if normspec = FALSE, set to 5 or above)
#' @param normspec logical. Whether to normalize the spectrogram (not recommended) or not (normalized spectrograms with different SNR are not comparable).
#' @param noisered numeric. noise reduction (subtract median from the amplitude values); 1=rows, 2=columns, 3=none.
#' @param rmoff logical. Whether to remove DC offset before computing aei (recommended) or not.
#' @param props logical. Whether to store the energy proportion values for each frequency band and channel (default) or not.
#' @param propden numeric. Indicates how the energy proportion is calculated by manipulating the denominator.
#'
#' @return a tibble (data frame) with the aei values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @details
#' Options for propden:
#' 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion equals to all the cells in the same frequency band.
#' 2 = The "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max_freq')
#' 3 = The "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a list of audio files before importing them into R.
#' The working directory should be set to the folder containing the files.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import tuneR
#' @import tidyverse
#' @import seewave
#' @import lubridate

#'
#' @examples
#' files <- list.files(pattern=".wav|.WAV")
#' aei_list(files[1:5])
aei_list <- function (audiolist,
                        save_csv = FALSE,
                        csv_name = "aei_results.csv",
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
                        propden = 1){

  
  #  Quiet function from SimDesign package to run functions without printing
  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      tmpf <- tempfile()
      sink(tmpf)
      on.exit({sink(); file.remove(tmpf)})
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }
  
  cat("Evaluating the job...\n\n")


  fileName <- tibble(file_name = audiolist)
  nFiles <- length(audiolist)


  args_list <- list(wlen = wlen, wfun = wfun, min_freq = min_freq,
                    max_freq = max_freq, nbands = nbands, db_threshold = db_threshold,
                    normspec = normspec, noisered = noisered, rmoff = rmoff,
                    props = props, propden = propden)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  aei1 <- quiet(aei(sound1, args_list$wlen, args_list$wfun, args_list$min_freq,
                    args_list$max_freq, args_list$nbands, args_list$db_threshold,
                    args_list$normspec, args_list$noisered, args_list$rmoff,
                    args_list$props, args_list$propden))

  tibble(file_name = "filename") %>% bind_cols(aei1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(aei1)

  # Declare the number of cores to be used (all but one of the available cores)
  cores <- detectCores() - 1 # Leave one core free
  # Limit the number of cores to the number of files, if 'cores' was initially a higher number
  if (cores > nFiles){
    cores <- nFiles
  }

  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Setup parallel back-end
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)

  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", cores, "cores... \n")


  # Start loop
  results <- foreach(file = audiolist, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate AEI and keep its default output columns
                       aei <- aei(sound, wlen = args_list$wlen, wfun = args_list$wfun,
                                  min_freq = args_list$min_freq, max_freq = args_list$max_freq,
                                  nbands = args_list$nbands, db_threshold = args_list$db_threshold,
                                  normspec = args_list$normspec, noisered = args_list$noisered,
                                  rmoff = args_list$rmoff, props = args_list$props,
                                  propden = args_list$propden)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(aei)


                     }


  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)


  stopCluster(cl)

  if(save_csv == TRUE){
    write.csv(resultsWithMetadata, csv_name, row.names = FALSE)
  }

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(resultsWithMetadata)

}
