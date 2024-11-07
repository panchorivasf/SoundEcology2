#' Normalized Difference Soundscape Index - list imput
#' @description
#' Normalized Difference Soundscape Index (NDSI) from REAL and Kasten,
#' et al. 2012. The NDSI seeks to "estimate the level of anthropogenic
#' disturbance on the soundscape by computing the ratio of human-generated
#' (anthrophony) to biological (biophony) acoustic components found in field
#' collected sound samples" (Kasten, et al. 2012). This version is optimized to work with lists of files.
#' @param audio.list a list of audio files from the working directory.
#' @param folder a path to the folder containing the audio files.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param w.len numeric. The window length for the FFT.
#' @param anthro.min minimum value of the range of frequencies of the anthrophony.
#' @param anthro.max maximum value of the range of frequencies of the anthrophony.
#' @param bio.min minimum value of the range of frequencies of the biophony.
#' @param bio.max maximum value of the range of frequencies of the biophony.
#' @param rm.offset logical. Whether to remove the DC offset.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#'
#' @return a wide format tibble with NDSI values per channel (if stereo), parameters used and audio metadata
#' @export
#' @details
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a list of audio files before importing them in
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
#' files <- list_waves(path/to/folder)
#' ndsi_list(files[1:5])
ndsi_list <- function (audio.list,
                      folder = NULL,
                      save.csv = TRUE,
                      csv.name = "ndsi_results.csv",
                      w.len = 50,
                      anthro.min = 1000,
                      anthro.max = 2000,
                      bio.min = 2000,
                      bio.max = 11000,
                      rm.offset = TRUE,
                      n.cores = -1){

  if(is.null(folder)){
    folder <- getwd()
  }
  
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


  fileName <- tibble(file_name = audio.list)
  nFiles <- length(audio.list)

  args_list <- list(w.len = w.len,
                    anthro.min = anthro.min,
                    anthro.max = anthro.max,
                    bio.min = bio.min,
                    bio.max = bio.max,
                    rm.offset = rm.offset)

  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audio.list[1], from = 0, to = 2 , units ='seconds')
  type <- ifelse(sound1@stereo, "stereo", "mono")
  ndsi1 <- quiet(ndsi(sound1,
              args_list$w.len,
              args_list$anthro.min,
              args_list$anthro.max,
              args_list$bio.min,
              args_list$bio.max,
              args_list$rm.offset))

  tibble(file_name = "filename") %>% bind_cols(ndsi1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(ndsi1)

  # Declare the number of cores to be used
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1
  }else{
    num_cores <- n.cores
  }

  if(nFiles < num_cores){
    num_cores <- nFiles
  }

  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Setup parallel back-end
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)

  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")

  # Start loop
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate NDSI and keep its default output columns
                       ndsi2 <- ndsi(sound,
                                    args_list$w.len,
                                    args_list$anthro.min,
                                    args_list$anthro.max,
                                    args_list$bio.min,
                                    args_list$bio.max,
                                    args_list$rm.offset)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(ndsi2)


                     }


  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)


  stopCluster(cl)

  if(save.csv == TRUE){
    write.csv(resultsWithMetadata, csv.name, row.names = FALSE)
  }

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(resultsWithMetadata)


}


