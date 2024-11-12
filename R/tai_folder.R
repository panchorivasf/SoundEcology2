#' Calculate Trill Activity Index for all the Files in a Folder
#'
#' This function calculates the trill index of an sound wave object by analyzing the frequency modulation pattern over time. It can operate in either binary or continuous modes and provides options for generating visual plots of the trill activity. The function also identifies potential noise issues in the low- and mid-frequency ranges.
#'
#' @param folder Character. The path to the folder containing the WAV files to analyze.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param wave A wave object to be analyzed.
#' @param channel Channel or channels to be analyzed. Options are "left", "right", "each", and "mix".
#' @param cutoff Numeric. The cutoff in decibels for the spectrogram generation.
#' @param n.windows Numeric. Number of time windows to divide the signal into for analysis. Default is 60.
#' @param freq.res Numeric. Frequency resolution (Hz per bin) for the spectrogram analysis. Default is 100.
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
#' \dontrun{
#' path <- "path/to/folder"
#' tai_results <- tai_folder(path)
#' }

tai_folder <- function(folder = NULL,
                       list = NULL,
                       channel = 'each',
                       hpf = 0,
                       cutoff = -60,
                       n.windows = 120,
                       freq.res = 100,
                       plot = FALSE,
                       plot.title = NULL,
                       verbose = TRUE,
                       output.csv = "tai_results.csv",
                       n.cores = -1) {
  cat("Working on it...\n")
  
  args_list <- list(channel = channel,
                    hpf = hpf,
                    cutoff = cutoff,
                    n.windows = n.windows,
                    freq.res = freq.res,
                    plot = plot,
                    plot.title = plot.title,
                    verbose = verbose)
  
  if(is.null(folder)){
    folder <- getwd()
  }
  setwd(folder)
  
  if(is.null(list)){
    audio.list <- quiet(list_waves())
    
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
    
    sound1 <- readWave(audio.list[1])
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    tai1 <- quiet(do.call(tai, c(list(sound1), args_list))$summary)
    
    tibble(file_name = "filename")  |>  bind_cols(tai1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(tai1)
    
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
  
  # Define parallel computation
  results <- foreach(file = audio.list, .packages = c("tuneR", "seewave", "dplyr")) %dopar% {
    
    filename <- basename(file)  
    
    sound <- tryCatch({
      readWave(file)
    }, error = function(e) {
      message(paste("Error reading file:", file, "Skipping to the next file."))
      return(NULL) 
    })
    if (is.null(sound)) {
      return(NULL)
    }
    
    # Initialize an empty tibble for the results
    result_list <- list()
    
    tai_result <- quiet(do.call(tai, c(list(sound), args_list))$summary)
    
    result_list <- list(
      tibble(file_name = filename, tai_result)
    )
    
    
    return(do.call(rbind, result_list))
  }
  
  # Combine all the results into a single tibble
  combined_results <- do.call(rbind, results)
  
  combined_results <- addMetadata(combined_results)
  combined_results$datetime <- format(combined_results$datetime, "%Y-%m-%d %H:%M:%S")
  
  # Export results to CSV
  write.csv(combined_results, file = output.csv, row.names = FALSE)
  
  # Stop parallel cluster
  stopCluster(cl)
  
  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(combined_results)
}
