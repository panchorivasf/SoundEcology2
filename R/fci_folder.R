#' Frequency Cover Indices on Multiple WAV Files in a Folder
#'
#' This function performs a frequency cover indices (FCI) on all WAV files within a specified folder.
#' It processes each sound file in parallel (if specified) and saves the results in a CSV file.
#'
#' @param folder The folder containing the WAV files to analyze.
#' @param recursive Logical. Whether to search in subfolders. Default is FALSE.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name The name of the CSV file where results will be saved. Default is "frequency_cover_results.csv".
#' @param channel The channel to analyze: 'left', 'right', 'mix' (combine both), or 'each' (process left and right channels separately). Default is 'each'.
#' @param hpf High-pass filter cutoff frequency in Hz. If 0, no high-pass filter is applied. Default is 0.
#' @param cutoff The amplitude threshold (in dB) below which frequencies will be considered inactive. Default is -60.
#' @param freq.res Frequency resolution of the spectrogram in Hz. Default is 100.
#' @param plot Logical. Whether to generate plots for each file. Default is FALSE.
#' @param ggplot Logical. Whether to use ggplot2 for plotting (if plot is TRUE). Default is FALSE.
#' @param plot.title Title for the plot if plotting is enabled. Default is "Frequency Cover Analysis".
#' @param sound.color The color to use for active frequencies in the plot. Default is "#045E10" (a dark green shade).
#' @param lf.min The minimum frequency (in Hz) for the low-frequency band. Default is 0.
#' @param lf.max The maximum frequency (in Hz) for the low-frequency band. Default is 1500.
#' @param mf.min The minimum frequency (in Hz) for the mid-frequency band. Default is 1500.
#' @param mf.max The maximum frequency (in Hz) for the mid-frequency band. Default is 8000.
#' @param hf.min The minimum frequency (in Hz) for the high-frequency band. Default is 8000.
#' @param hf.max The maximum frequency (in Hz) for the high-frequency band. Default is 18000.
#' @param uf.min The minimum frequency (in Hz) for the ultra-high-frequency band. Default is 18000.
#' @param uf.max The maximum frequency (in Hz) for the ultra-high-frequency band. Default is 24000.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#'
#' @return A tibble containing the frequency cover analysis results for each file.
#' @export
#' @import foreach
#' @import seewave
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom tuneR readWave
#' @importFrom dplyr tibble bind_cols
#' @examples
#' \dontrun{
#' # Run frequency cover analysis on all WAV files in the folder "sound_files"
#' fci_folder("sound_files", csv.name = "results.csv", channel = "left", n.cores = 4)
#' }

fci_folder <- function(folder = NULL,
                       recursive = FALSE,
                       list = NULL,
                       start = 0,
                       end = 1,
                       unit = "minutes",
                       save.csv = TRUE,
                       csv.name = "fci_results.csv",
                       channel = 'each',
                       hpf = 0,
                       cutoff = -60,
                       freq.res = 50,
                       spectrogram = FALSE,
                       ggplot = FALSE,
                       plot.title = "Frequency Cover Analysis",
                       sound.color = "#045E10",
                       lf.min = 0,
                       lf.max = 1500,
                       mf.min = 1500,
                       mf.max = 8000,
                       hf.min = 8000,
                       hf.max = 18000,
                       uf.min = 18000,
                       uf.max = 24000,
                       n.cores = -1) {
  cat("Working on it...\n")
  
  args_list <- list(channel = channel,
                    hpf = hpf,
                    cutoff = cutoff,
                    freq.res = freq.res,
                    spectrogram = FALSE,
                    ggplot = FALSE,
                    plot.title = NULL,
                    sound.color = "#045E10",
                    lf.min = lf.min,
                    lf.max = lf.max,
                    mf.min = mf.min,
                    mf.max = mf.max,
                    hf.min = hf.min,
                    hf.max = hf.max,
                    uf.min = uf.min,
                    uf.max = uf.max,
                    verbose = FALSE)
  
  if(is.null(folder)){
    folder <- getwd()
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
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  if (nFiles > 10) {
    cat("Evaluating the job...\n")
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1], 
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    fci1 <- quiet(do.call(fci, c(list(sound1), channel = 'each')))
    # 
    # fci1 <- quiet(fci(sound1, channel = 'each'))
    tibble(file_name = "filename") |> bind_cols(fci1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 3
    
    rm(sound1)
    rm(fci1)
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    # Set up parallel processing
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  } else {
    sound1 <- readWave(audio.list[1], from = start,
                       to = end, units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  
  # Define parallel computation
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "dplyr", "seewave")) %dopar% {
                       
                       # Try to read the sound file, handle errors gracefully
                       sound <- tryCatch({
                         readWave(file, 
                                  from = start,
                                  to = end,
                                  units = unit)
                       }, error = function(e) {
                         message(paste("Error reading file:", file, "Skipping to the next file."))
                         return(NULL) 
                       })
                       
                       # Skip processing if the sound is NULL (i.e., readWave failed)
                       if (is.null(sound)) {
                         return(NULL)
                       }
                       
                       fci_result <- quiet(do.call(fci, c(list(sound), args_list)))
                       
                       tibble(file_name = file)|> 
                         bind_cols(fci_result)
                       
                     }
  stopCluster(cl)
  
  resultsWithMetadata <- addMetadata(results)
  
  if (save.csv){
    
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

