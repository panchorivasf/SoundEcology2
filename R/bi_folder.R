#' Bioacoustic Index - folder input
#' @description
#' Inspired by the "Bioacoustic Index" from the paper:Boelman NT, Asner GP, Hart PJ, Martin RE. 2007.
#' Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne
#' remote sensing. Ecol Applications 17(8):2137-44. Based on Matlab code provided by NT Boelman.
#' Boelman et al. 2007 used min.freq=2000, max.freq=8000, w.len=512.
#' Several parts where changed, in particular log math, so this won't be
#' directly comparable to the original code in the paper.
#'
#' @param folder a path to the folder with audio files to import.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param w.len the window length to compute the spectrogram (i.e., FFT window size).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq miminum frequency to use when calculating the value, in Hertz. Default = NA.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0 (Default), noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).

#' @return A tibble (data frame) with the BI values for each channel (if stereo), metadata, and the parameters used for the calculation.
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
#' bi_folder(path/to/folder)

bi_folder <- function (folder,
                       save.csv = FALSE,
                       csv.name = "bi_results.csv",
                       w.len = 512,
                       w.fun = "hanning",
                       min.freq = 2000,
                       max.freq = 8000,
                       norm.spec = FALSE,
                       noise.red = 0,
                       rm.offset = TRUE,
                       n.cores = -1){


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

  args_list <- list(w.len = w.len,
                    w.fun = w.fun,
                    min.freq = min.freq,
                    max.freq = max.freq,
                    norm.spec = norm.spec,
                    noise.red = noise.red,
                    rm.offset = rm.offset)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audiolist[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  bi1 <- quiet(bi(sound1,
                  args_list$w.len,
                  args_list$w.fun,
                  args_list$min.freq,
                  args_list$max.freq,
                  args_list$norm.spec,
                  args_list$noise.red,
                  args_list$rm.offset))

  tibble(file_name = "filename") %>% bind_cols(bi1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(bi1)

  # Assess the number of cores to be used
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
  results <- foreach(file = audiolist, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate BI and keep its default output columns
                       bi <- bi(sound,
                                  args_list$w.len,
                                  args_list$w.fun,
                                  args_list$min.freq,
                                  args_list$max.freq,
                                  args_list$norm.spec,
                                  args_list$noise.red,
                                  args_list$rm.offset)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(bi)


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
