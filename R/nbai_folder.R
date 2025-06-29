#' Narrow-Band Activity Index for Files in a Folder
#' @description
#' NBAI describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. 
#' This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary
#' acoustic indices).
#'
#' @param folder Character. The path to the folder containing the wave files to analyze.
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param list An optional list (subset) of files in the folder to analyze. 
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right"
#'  to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), 
#'  results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals 
#' of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, 
#' therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.
#' @param activity.cutoff Numeric. Cutoff percent activity. Only the frequency bands active equal or above this 
#' percentage will be considered as "active" in the active band statistics.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core.
#' Default is NULL (single-core processing).
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#'
#' @return A list containing: 1) A binary spectrogram (if mono), 2) tibble with the Narrow-Band Activity Index (NBI) 
#' summary statistics, and 3) a tibble with NBI spectral, which number of rows equals the number of frequency bins 
#' in the analysis.
#' @export
#'
#' @import doParallel 
#' @import foreach
#' @import seewave
#' @importFrom parallel detectCores makeCluster
#' @importFrom tuneR readWave
#' @importFrom dplyr bind_cols tibble
#'
#' @examples 
#' \dontrun{
#' nbai(wave, channel = 'left', plot = TRUE, verbose = TRUE)
#' }

nbai_folder <- function(folder = NULL,
                        recursive = FALSE,
                        list = NULL,
                        start = 0,
                        end = 1,
                        unit = "minutes",
                        channel = "each",
                        hpf = 0,
                        freq.res = 50,
                        cutoff = -60,
                        activity.cutoff = 10, 
                        output.csv = "nbai_results.csv",
                        n.cores = -1,
                        verbose = TRUE) {
  cat("Working on it...\n")
  
  args_list <- list(channel = channel,
                    hpf = hpf,
                    freq.res = freq.res,
                    cutoff = cutoff,
                    activity.cutoff = activity.cutoff, 
                    verbose = verbose)
  
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
    
    nbai1 <- quiet(do.call(nbai, c(list(sound1), args_list)))
    
    tibble(file_name = "filename")  |>  bind_cols(nbai1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(nbai1)
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * n.files) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  } else {
    sound1 <- readWave(audio.list[1], from = start, to = end , units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  
  cat("Analyzing", n.files, type, "files using", num_cores, "cores... \n")

  results <- foreach(file = audio.list,
                     .packages = c("tuneR", "seewave", "dplyr")) %dopar% {
                       filename <- basename(file)  
                       sound <- readWave(file,
                                         from = start,
                                         to = end,
                                         units = unit)
                       result_list <- list()
                       
                       nbai_result <- quiet(do.call(nbai, c(list(sound), args_list)))$summary
                       
                       result_list <- list(tibble(file_name = filename, nbai_result))
                       
                       return(do.call(rbind, result_list))
                     }
  

  combined_results <- do.call(rbind, results)
  combined_results <- addMetadata(combined_results)
  
  # Export results to CSV
  write.csv(combined_results, file = output.csv, row.names = FALSE)
  
  # Stop parallel cluster
  stopCluster(cl)
  
  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(combined_results)
}






