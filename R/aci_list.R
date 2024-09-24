#' Acoustic Complexity Index - list input
#'
#' @param audiolist a list of audio files to import.
#' @param save_csv logical. Whether to save a csv in the working directory.
#' @param csv_name character vector. When 'save_csv' is TRUE, optionally provide a file name.
#' @param wlen the window length to compute the spectrogram (i.e., FFT window size).
#' @param wfun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min_freq minimum frequency to use when calculating the value, in Hertz. Default = 0.
#' @param max_freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param j the cluster size, in seconds. Default = NA (Duration of the audio file).
#' @param noisered numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rmoff logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.

#' @return A tibble (data frame) with the ACI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import tuneR
#' @import tidyverse
#' @import seewave
#' @import lubridate
#'
#' @details
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a folder of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' files <- list.files(pattern=".wav|.WAV")
#' aci_list(files[1:5])

aci_list <- function (audiolist,
                        save_csv = FALSE,
                        csv_name = "aci_results.csv",
                        wlen = 512,
                        wfun = "hanning",
                        min_freq = NA,
                        max_freq = NA,
                        j = NA,
                        noisered = 2,
                        rmoff = TRUE){

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

  args_list <- list(wlen = wlen,
                    wfun = wfun,
                    min_freq = min_freq,
                    max_freq = max_freq,
                    j = j,
                    noisered = noisered,
                    rmoff = rmoff
  )

  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  aci1 <- quiet(aci(sound1,
                    args_list$wlen,
                    args_list$wfun,
                    args_list$min_freq,
                    args_list$max_freq,
                    args_list$j,
                    args_list$noisered,
                    args_list$rmoff))

  tibble(file_name = "filename") %>% bind_cols(aci1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(aci1)

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

                       # Calculate ACI and keep its default output columns
                       aci <- aci(sound,
                                  args_list$wlen,
                                  args_list$wfun,
                                  args_list$min_freq,
                                  args_list$max_freq,
                                  args_list$j,
                                  args_list$noisered,
                                  args_list$rmoff)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(aci)

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
