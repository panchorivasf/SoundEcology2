#' Calculate Trill Activity Index for all the Files in a Folder
#'
#' This function calculates the trill index of an sound wave object by analyzing the frequency modulation pattern over time. It can operate in either binary or continuous modes and provides options for generating visual plots of the trill activity. The function also identifies potential noise issues in the low- and mid-frequency ranges.
#'
#' @param folder Character. The path to the folder containing the WAV files to analyze.
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param start Numeric. Where to start reading the Wave. 
#' @param end Numeric. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param hpf Numeric. Whether to use a high-pass filter. Default is 0.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param save.csv logical. Whether to save a CSV output.
#' @param save.to character. Path to where the output CSV will be saved. Default
#' is NULL (save in working directory).
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param cutoff Numeric. The cutoff in decibels for the spectrogram generation.
#' @param n.windows Numeric. Number of time windows to divide the signal into for analysis. Default is 60.
#' @param freq.res Numeric. Frequency resolution (Hz per bin) for the spectrogram analysis. Default is 100.
#' @param plot Logical. If TRUE, generates a plot of the trill index over time. Default is TRUE.
#' @param plot.title Character. The title to be used for the plot. Default is NULL.
#' @param n.cores The number of cores to use for parallel processing. Default is
#' 1 to use all but one core. 
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
                       recursive = FALSE,
                       list = NULL,
                       start = 0,
                       end = 1,
                       unit = "minutes",
                       channel = 'each',
                       save.csv = TRUE,
                       save.to = NULL,
                       csv.name = "tai_results.csv",
                       hpf = 0,
                       cutoff = -60,
                       n.windows = 120,
                       freq.res = 100,
                       plot = FALSE,
                       plot.title = NULL,
                       n.cores = -1,
                       verbose = TRUE) {
  cat("Working on it...\n")
  
  args_list <- list(channel = channel,
                    hpf = hpf,
                    cutoff = cutoff,
                    n.windows = n.windows,
                    freq.res = freq.res,
                    plot = plot,
                    plot.title = plot.title,
                    verbose = verbose)
  
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
    audio.list <- quiet(list_waves(recursive = recursive))
    
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
    sound1 <- readWave(audio.list[1], 
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  
  cat("Analyzing", n.files, type, "files using", num_cores, "cores... \n")
  
  # Define parallel computation
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "dplyr", "seewave")) %dopar% {
                       
                       sound <- tryCatch({
                         readWave(file,
                                  from = start,
                                  to = end,
                                  units = unit)
                       }, error = function(e) {
                         message(paste("Error reading file:", file, 
                                       "Skipping to the next file."))
                         return(NULL) 
                       })
                       if (is.null(sound)) {
                         return(NULL)
                       }
                       
                       tai_result <- quiet(do.call(tai, c(list(sound), args_list)))
                       
                       tibble(file_name = file) |>
                         bind_cols(tai_result)
                     }
  
  stopCluster(cl)
  setwd(original_wd)
  
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
  
  # # Combine all the results into a single tibble
  # combined_results <- do.call(rbind, results)
  # 
  # combined_results <- addMetadata(combined_results)
  # combined_results$datetime <- format(combined_results$datetime, "%Y-%m-%d %H:%M:%S")
  # 
  # # Export results to CSV
  # write.csv(combined_results, file = csv.name, row.names = FALSE)
  # 

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(results)
}
