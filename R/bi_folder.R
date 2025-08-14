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
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param save.csv logical. Whether to save a CSV output.
#' @param save.to character. Path to where the output CSV will be saved. Default
#' is NULL (save in working directory).
#' @param csv.name character. When 'save.csv' is TRUE, optionally provide a file name.
#' @param w.len numeric. The window length to compute the spectrogram (i.e., FFT window size).
#' @param w.fun character. The window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", 
#' "hamming", "hanning", or "rectangle".
#' @param min.freq miminum frequency to use when calculating the value, in Hertz. Default = NA.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param norm.spec logical; if TRUE, the spectrogram is normalized, scaled by its maximum value (not recommended 
#' because normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied 
#' to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each 
#' column similarly. If set to 0 (Default), noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param n.cores The number of cores to use for parallel processing. Default is
#' -1 to use all but one core.
#' 
#' @return A tibble (data frame) with the BI values for each channel (if stereo), metadata, and the parameters used 
#' for the calculation.
#' @export
#'
#' @import doParallel 
#' @import foreach
#' @import seewave
#' @importFrom parallel detectCores makeCluster
#' @importFrom tuneR readWave
#' @importFrom dplyr bind_cols tibble
#'
#' @details
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a folder of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' bi_folder(path/to/folder)
#' 
bi_folder <- function (folder = NULL,
                       recursive = FALSE,
                       start = 0,
                       end = 1,
                       unit = "minutes",
                       list = NULL,
                       save.csv = TRUE,
                       save.to = NULL,
                       csv.name = "bi_results.csv",
                       w.len = 512,
                       w.fun = "hanning",
                       min.freq = 2000,
                       max.freq = 8000,
                       norm.spec = FALSE,
                       noise.red = 0,
                       rm.offset = TRUE,
                       n.cores = -1){
  cat("Working on it...\n")
  
  args_list <- list(w.len = w.len,
                    w.fun = w.fun,
                    min.freq = min.freq,
                    max.freq = max.freq,
                    norm.spec = norm.spec,
                    noise.red = noise.red,
                    rm.offset = rm.offset)
  
  original_wd <- getwd()
  
  if(is.null(folder)) {
    folder <- getwd()
  }
  
  if(is.null(save.to)){
    save.to <- folder
  }
  
  if(!dir.exists(save.to)){
    dir.create(save.to)
  }
  
  setwd(folder)
  
  if(is.null(list)){
    audio.list <- list_waves(recursive = recursive)
    
  } else {
    audio.list <- list
  }
  n.files <- length(audio.list)
  
  # Declare the number of cores to be used 
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1 # Leave one core free 
  }else{
    num_cores <- n.cores
  } 
  if (num_cores > n.files){
    num_cores <- n.files  # Limit the number of cores to the number of files
  }
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)

  # Evaluate the duration of the analysis if n.files > 10
  if(n.files>10){
    cat("Evaluating the job...\n\n")
    
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1],
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    bi1 <- quiet(do.call(bi, c(list(sound1), args_list)))
    
    tibble(file_name = "filename")  |>  bind_cols(bi1 )
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(bi1 )
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * n.files) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  } else {
    sound1 <- readWave(audio.list[1], from = 0, to = 2 , units ='seconds')
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  cat("Analyzing", n.files, type, "files using", num_cores, "cores... \n")


  # Start loop
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "seewave", "dplyr")) %dopar% {
                       
                       sound <- tryCatch({
                         readWave(file,
                                  from = start,
                                  to = end,
                                  units = unit)
                       }, error = function(e) {
                         message(paste("Error reading file:", file, "Skipping to the next file."))
                         return(NULL) 
                       })
                       if (is.null(sound)) {
                         return(NULL)
                       }
                       
                       bi_result <- quiet(do.call(bi, c(list(sound), args_list)))
                       
                       tibble(file_name = file)  |> 
                         bind_cols(bi_result)
                       
                     }

  stopCluster(cl)
  setwd(original_wd)
  
  # Combine results with metadata and return
  results <- addMetadata(results)
  
  if(save.csv == TRUE){
    
    sensor <- unique(results$sensor_id)
    
    results$datetime <- format(results$datetime, 
                               "%Y-%m-%d %H:%M:%S")
    
    # Export results to CSV
    if (length(sensor) == 1){
      write.csv(results, file = paste0(save.to, "/", sensor,"_", 
                                       csv.name, ".csv"), row.names = FALSE)
    } else {
      write.csv(results, file = paste0(save.to, "/", csv.name, ".csv"), 
                row.names = FALSE)
    }
    
  }

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(results)

}
