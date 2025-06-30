#' Normalized Difference Soundscape Index - folder input
#' @description
#' Normalized Difference Soundscape Index (NDSI) from Kasten et al.,
#' et al. 2012. The NDSI seeks to "estimate the level of anthropogenic
#' disturbance on the soundscape by computing the ratio of human-generated
#' (anthrophony) to biological (biophony) acoustic components found in field
#' collected sound samples" (Kasten, et al. 2012).
#' This version is optimized to work with a path to the folder containing the audio files.
#' @param folder a path to the folder containing the audio files.
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param list An optional list (subset) of files in the folder to analyze. 
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name. 
#' Default name is "ndsi_results.csv"
#' @param w.len numeric. Window length for the FFT (sampling rate / frequency resolution).
#' @param anthro.min minimum value of the range of frequencies of the anthrophony.
#' @param anthro.max maximum value of the range of frequencies of the anthrophony.
#' @param bio.min minimum value of the range of frequencies of the biophony.
#' @param bio.max maximum value of the range of frequencies of the biophony.
#' @param rm.offset logical. Whether to remove the DC offset.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all 
#' but one core. Default is NULL (single-core processing).
#'
#' @return a wide format tibble with NDSI values per channel (if stereo), parameters used and audio 
#' metadata
#' @export
#'
#' @import foreach
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom tuneR readWave
#' @importFrom dplyr tibble bind_cols
#' @import seewave
#'
#' @examples
#' \dontrun{
#' path <- readClipboard() #use this to paste the folder path from the clipboard.
#' ndsi_folder(path)
#' }

ndsi_folder <- function (folder = NULL,
                         recursive = FALSE,
                         list = NULL,
                         start = 0,
                         end = 1,
                         unit = "minutes",
                         save.csv = TRUE,
                         csv.name = "ndsi_results.csv",
                         w.len = 512,
                         anthro.min = 1000,
                         anthro.max = 2000,
                         bio.min = 2000,
                         bio.max = 11000,
                         rm.offset = TRUE,
                         n.cores = -1){
  
  cat("Working on it...\n")
  
  args_list <- list(w.len = w.len,
                    anthro.min = anthro.min,
                    anthro.max = anthro.max,
                    bio.min = bio.min,
                    bio.max = bio.max,
                    rm.offset = rm.offset)
  
  if(is.null(folder)){
    folder <- getwd()
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
    num_cores <- parallel::detectCores() - 1 
  }else{
    num_cores <- n.cores
  } 
  if (num_cores > n.files){
    num_cores <- n.files  
  }
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  if(n.files>10){
    cat("Evaluating the job...\n\n")
    
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1],
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    ndsi1 <- quiet(do.call(ndsi, c(list(sound1), args_list)))
    
    tibble(file_name = "filename")  |>  bind_cols(ndsi1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(ndsi1)
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * n.files) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  } else {
    sound1 <- readWave(audio.list[1], 
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  
  cat("Analyzing", n.files, type, "files using", num_cores, "cores... \n")
  
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
                       if (is.null(sound)) {
                         return(NULL)
                       }
                       
                       ndsi_result <- quiet(do.call(ndsi, c(list(sound), args_list)))
                       
                       tibble(file_name = file) |>
                         bind_cols(ndsi_result)
                     }
  stopCluster(cl)
  
  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)
  
  if(save.csv == TRUE){
    sensor <- unique(resultsWithMetadata$sensor_id)
    
    resultsWithMetadata$datetime <- format(resultsWithMetadata$datetime, "%Y-%m-%d %H:%M:%S")
    # Export results to CSV
    if (length(sensor) == 1){
      write.csv(resultsWithMetadata, file = paste0(sensor,"_",csv.name), row.names = FALSE)
    } else {
      write.csv(resultsWithMetadata, file = csv.name, row.names = FALSE)
    }
    
  }
  
  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(resultsWithMetadata)
  
}
