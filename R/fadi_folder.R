#' Frequency-dependent Acoustic Diversity Index - folder input
#' @description
#' The Frequency-dependent Acoustic Diversity Index by Xu et al. (2023) obtains a floating noise profile
#' before calculating the Acoustic Diversity Index and it doesn't use normalized spectrogram.
#' Alternatively it can take a noise sample to reduce noise in the analyzed files.
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
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide 
#' a file name.
#' @param noise.file An R object of class Wave containing noise-only information 
#' if needed. Default = NULL.
#' @param NEM Numeric. Options are 1 or 2.
#' When NEM = 1, floating thresholds are estimated based on noise.file.
#' When NEM = 2, floating thresholds are calculated based on sound file using an
#' automatic noise level estimation method (median of each row in the 
#' spectrogram). Default = 2.
#' @param min.freq Minimum frequency in Hertz when calculating the global 
#' threshold. Default = 200.
#' @param max.freq Maximum frequency in Hertz when calculating the FADI value. 
#' Default = 10000.
#' @param threshold.fixed A negative number in dB for calculating the global 
#' threshold. Default = −50.
#' @param freq.step Bandwidth of each frequency band, in Hertz. Default = 1000.
#' @param gamma A positive number in dB for calculating the floating thresholds. 
#' Default = 13.
#' @param props Logical; if TRUE, the energy proportion values for each 
#' frequency band and channel are added to the output tibble. Default = TRUE.
#' @param n.cores The number of cores to use for parallel processing. Default is 
#' `n.cores = -1` to use all but one core. 
#'
#' @return A tibble with the FADI value per channel, energy proportions, 
#' metadata, and parameters used.
#' 
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

fadi_folder <- function (folder = NULL,
                         recursive = FALSE,
                         list = NULL,
                         start = 0,
                         end = 1,
                         unit = "minutes",
                         save.csv = TRUE,
                         save.to = NULL,
                         csv.name = "fadi_results",
                         noise.file = NULL,
                         NEM = 2,
                         min.freq = 200,
                         max.freq = 10000,
                         threshold.fixed = -50,
                         freq.step = 1000,
                         gamma = 13,
                         props = FALSE,
                         n.cores = -1){
  cat("Working on it...\n")
  
  args_list <- list(noise.file=noise.file,
                    NEM=NEM,
                    min.freq=min.freq,
                    max.freq = max.freq, 
                    threshold.fixed = threshold.fixed,
                    freq.step = freq.step, 
                    gamma = gamma, 
                    props=props)


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
  
  nFiles <- length(audio.list)
  
  
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
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  if(nFiles>10){
    cat("Evaluating the job...\n")

    startTime <- Sys.time()

    sound1 <- readWave(audio.list[1],
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    fadi1 <- quiet(fadi(sound1,args_list))

    tibble(file_name = "filename") |> bind_cols(fadi1)

    rm(sound1)
    rm(fadi1)

    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    timePerFile <- timePerFile + 3 # Add overhead

    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + estimatedTotalTime

    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  } else {
    sound1 <- readWave(audio.list[1], from = start, to = end , units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  
  # Start loop
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "dplyr", "seewave")) %dopar% {

                       sound <- tryCatch({
                         readWave(file,
                                  from = start,
                                  to = end,
                                  units = unit)
                       }, error = function(e) {
                         message(paste("Error reading file:", file, "Skipping to the next file."))
                         return(NULL) 
                       })
                       
                       if (is.null(sound)) return(NULL)
                       
                       fadi_result <- quiet(do.call(fadi, c(list(sound), args_list)))
                       
                       tibble(file_name = file)  |> 
                         bind_cols(fadi_result)

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

