#' Frequency-dependent Acoustic Diversity Index - folder input
#' @description
#' The Frequency-dependent Acoustic Diversity Index by Xu et al. (2023) obtains a floating noise profile
#' before calculating the Acoustic Diversity Index and it doesn't use normalized spectrogram.
#' Alternatively it can take a noise sample to reduce noise in the analyzed files.
#' @param folder a path to the folder with audio files to import.
#' @param save_csv logical. Whether to save a csv in the working directory.
#' @param csv_name character vector. When 'save_csv' is TRUE, optionally provide a file name.
#' @param noise_file An R object of class Wave containing noise-only information if needed. Default = NULL.
#' @param NEM Numeric. Options are 1 or 2.
#' When NEM = 1, floating thresholds are estimated based on noise_file.
#' When NEM = 2, floating thresholds are calculated based on sound file using an
#' automatic noise level estimation method (median of each row in the spectrogram). Default = 2.
#' @param min_freq Minimum frequency in Hertz when calculating the global threshold. Default = 200.
#' @param max_freq Maximum frequency in Hertz when calculating the FADI value. Default = 10000.
#' @param threshold_fixed A negative number in dB for calculating the global threshold. Default = −50.
#' @param freq_step Bandwidth of each frequency band, in Hertz. Default = 1000.
#' @param gamma A positive number in dB for calculating the floating thresholds. Default = 13.
#' @param props Logical; if TRUE, the energy proportion values for each frequency band
#' and channel are added to the output tibble. Default = TRUE.
#'
#' @return A tibble with the FADI value per channel, energy proportions, metadata, and parameters used.
#' @export
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import tuneR
#' @import tidyverse
#' @import seewave
#' @import lubridate
#'
#' @details
#' Modified version of the Frequency-dependent Acoustic Diversity Index by Xu et al. (2023).
#' FADI was introduced in: https://www.sciencedirect.com/science/article/pii/S1470160X23010828.
#' This version returns a wide format (one row per audio file) tibble as output instead of a nested list.
#' To see the original version as in the paper, use the \emph{frequency_dependent_acoustic_diversity()} function.

#' @examples
#' fadi_folder(folder=pathB, "fadi_hydro_b.csv")

fadi_folder <- function (folder,
                         save_csv = FALSE,
                         csv_name = "fadi_results.csv",
                         noise_file=NULL,
                         NEM=2,
                         min_freq=200,
                         max_freq = 10000,
                         threshold_fixed = -50,
                         freq_step = 1000,
                         gamma = 13,
                         props=TRUE){

  
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
  
  cat("Evaluating the job...\n")

  setwd(folder)
  audiolist <- list.files(path=folder, pattern = ".wav|.WAV")

  fileName <- tibble(file_name = audiolist)
  nFiles <- length(audiolist)


  args_list <- list(noise_file=noise_file,NEM=NEM,min_freq=min_freq,
                    max_freq = max_freq, threshold_fixed = threshold_fixed,
                    freq_step = freq_step, gamma = gamma, props=props)



  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")
  fadi1 <- quiet(fadi(sound1, args_list$noise_file, args_list$NEM,
                     args_list$min_freq, args_list$max_freq, args_list$threshold_fixed,
                     args_list$freq_step, args_list$gamma, args_list$props))

  tibble(file_name = "filename") %>% bind_cols(fadi1)

  rm(sound1)
  rm(fadi1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime

  # Assess number of cores to be used (all but one of the available cores)
  cores <- detectCores() - 1
  # If the number of cores is larger than the number if files...
  if (cores > nFiles){
    cores <- nFiles
  }

  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(cores)
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + estimatedTotalTime
  # Setup parallel back-end
  cl <- makeCluster(cores[1]) # Leave one core free
  registerDoParallel(cl)


  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", cores, "cores... \n")

  # Start loop
  results <- foreach(file = audiolist, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate FADI and keep its default output columns
                       fadi <- fadi(soundfile=sound,
                                         args_list$noise_file,
                                         args_list$NEM,
                                         args_list$min_freq,
                                         args_list$max_freq,
                                         args_list$threshold_fixed,
                                         args_list$freq_step,
                                         args_list$gamma,
                                         args_list$props)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(fadi)
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

