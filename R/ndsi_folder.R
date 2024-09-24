#' Normalized Difference Soundscape Index - folder input
#' @description
#' Normalized Difference Soundscape Index (NDSI) from REAL and Kasten,
#' et al. 2012. The NDSI seeks to "estimate the level of anthropogenic
#' disturbance on the soundscape by computing the ratio of human-generated
#' (anthrophony) to biological (biophony) acoustic components found in field
#' collected sound samples" (Kasten, et al. 2012).
#' This version is optimized to work with a path to the folder containing the audio files.
#' @param folder a path to the folder containing the audio files.
#' @param save_csv logical. Whether to save a csv in the working directory.
#' @param csv_name character vector. When 'save_csv' is TRUE, optionally provide a file name.
#' @param fft_w FFT window size.
#' @param anthro_min minimum value of the range of frequencies of the anthrophony.
#' @param anthro_max maximum value of the range of frequencies of the anthrophony.
#' @param bio_min minimum value of the range of frequencies of the biophony.
#' @param bio_max maximum value of the range of frequencies of the biophony.
#'
#' @return a wide format tibble with NDSI values per channel (if stereo), parameters used and audio metadata
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
#' @examples
#' path <- readClipboard() #use this to paste the folder path from the clipboard.
#' ndsi_folder(path)
ndsi_folder <- function (folder,
                       save_csv = FALSE,
                       csv_name = "ndsi_results.csv",
                       fft_w = 1024,
                       anthro_min = 1000,
                       anthro_max = 2000,
                       bio_min = 2000,
                       bio_max = 11000){

  
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

  setwd(folder)
  audiolist <- list.files(path=folder, pattern = ".wav|.WAV")

  fileName <- tibble(file_name = audiolist)
  nFiles <- length(audiolist)

  args_list <- list(fft_w = fft_w,
                    anthro_min = anthro_min,
                    anthro_max = anthro_max,
                    bio_min = bio_min,
                    bio_max = bio_max)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  ndsi1 <- quiet(ndsi(sound1,
                args_list$fft_w,
                args_list$anthro_min,
                args_list$anthro_max,
                args_list$bio_min,
                args_list$bio_max))


   tibble(file_name = "filename") %>% bind_cols(ndsi1)


  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(ndsi1)

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

                       library(soundecology2)
                       library(seewave)
                       library(tuneR)

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate NDSI and keep its default output columns
                       ndsi <- ndsi(sound,
                                    args_list$fft_w,
                                    args_list$anthro_min,
                                    args_list$anthro_max,
                                    args_list$bio_min,
                                    args_list$bio_max)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(ndsi)


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
