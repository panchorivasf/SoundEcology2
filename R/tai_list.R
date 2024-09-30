#' Calculate Trill Activity Index for all the Files in a List
#'
#' This function calculates the trill index of an audio wave object by analyzing the frequency modulation pattern over time. It can operate in either binary or continuous modes and provides options for generating visual plots of the trill activity. The function also identifies potential noise issues in the low- and mid-frequency ranges.
#'
#' @param audio.list a list of audio files to analyze
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param wave A wave object to be analyzed.
#' @param channel Channel or channels to be analyzed. Options are "left", "right", "each", and "mix".
#' @param cutoff Numeric. The cutoff in decibels for the spectrogram generation.
#' @param n.windows Numeric. Number of time windows to divide the signal into for analysis. Default is 60.
#' @param freq.res Numeric. Frequency resolution (in Hz) for the spectrogram analysis. Default is 100.
#' @param plot Logical. If TRUE, generates a plot of the trill index over time. Default is TRUE.
#' @param plot.title Character. The title to be used for the plot. Default is NULL.
#' @param verbose Logical. If TRUE, provides detailed output during the function's execution. Default is FALSE.
#' @return A tibble summarizing TAI statistics, including values for low and mid-frequency noise.

#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom lubridate seconds
#'
#' @examples 
#' files.list <- list_waves(path/to/folder)
#' tai_results <- tai_list(files.list)
tai_list <- function(audio.list,
                     channel = 'each',
                     hpf = 0,
                     rm.offset = TRUE,
                     cutoff = -60,
                     n.windows = 120,
                     freq.res = 100,
                     plot = FALSE,
                     plot.title = NULL,
                     verbose = TRUE,
                     output.csv = "tai_results.csv",
                     n.cores = -1) {
  
  
  cat("Evaluating the job...\n")
  
  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      tmpf <- tempfile()
      sink(tmpf)
      on.exit({sink(); file.remove(tmpf)})
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }
  
  # setwd(folder)
  # files <- list.files(path=folder, pattern = ".wav|.WAV")
  
  filename <- tibble(filename = audio.list)
  nFiles <- length(audio.list)
  
  
  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()
  
  sound1 <- readWave(audio.list[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")
  
  tai1 <- quiet(tai(sound1, channel = channel))
  
  tibble::tibble(file_name = "filename") %>% 
    dplyr::bind_cols(tai1)
  
  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(lubridate::seconds(2.2))
  
  rm(sound1)
  rm(tai1)
  
  # Determine number of cores to use
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1
  }else{
    num_cores <- n.cores
  }
  
  # Release unused cores
  if(nFiles < num_cores){
    num_cores <- nFiles
  }
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Set up parallel processing
  
  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  
  
  # Define parallel computation
  results <- foreach(file = audio.list, .packages = c("tuneR", "seewave", "tibble")) %dopar% {
    
    filename <- basename(file)  # Get file name without path
    
    # Read audio file
    audio <- readWave(file)
    
    # Initialize an empty tibble for the results
    result_list <- list()
    
    tai <- tai(audio,
               channel = channel,
               hpf = hpf,
               rm.offset = rm.offset,
               cutoff = cutoff,
               n.windows = n.windows,
               freq.res = freq.res,
               plot = FALSE,
               plot.title = NULL,
               verbose = FALSE)$summary
    
    result_list <- list(
      tibble(file_name = filename, channel = channel, tai)
    )
    
    
    return(do.call(rbind, result_list))
  }
  
  # Combine all the results into a single tibble
  combined_results <- do.call(rbind, results)
  
  combined_results <- addMetadata(combined_results)
  
  # Export results to CSV
  write.csv(combined_results, file = output.csv, row.names = FALSE)
  
  # Stop parallel cluster
  stopCluster(cl)
  
  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(combined_results)
}



