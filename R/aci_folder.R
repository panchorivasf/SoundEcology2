#' Acoustic Complexity Index - Batch process
#'
#' @param folder a path to the folder with audio files to import.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param freq.res numeric. The frequency resolution to use (Hz per bin) which will determine 
#' the window length for the FFT (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", 
#' "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to use when calculating the value, in Hertz. Default = 0.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param j the cluster size, in seconds. Default = NA (Duration of the audio file).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, 
#' noise reduction is applied to each row by subtracting the median from the amplitude values. 
#' If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. 
#' Default = TRUE.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. 
#' Default is NULL (single-core processing).


#' @return A tibble (data frame) with the ACI values for each channel (if stereo), metadata, and 
#' the parameters used for the calculation.
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
#' Optimized to facilitate working with a folder of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' aci_folder(path/to/folder)

aci_folder <- function (folder = NULL,
                        list = NULL,
                        start = 0,
                        end = 1,
                        unit = "minutes",
                        save.csv = TRUE,
                        csv.name = "aci_results.csv",
                        freq.res = 50,
                        win.fun = "hanning",
                        min.freq = NA,
                        max.freq = NA,
                        j = NA,
                        noise.red = 2,
                        rm.offset = TRUE,
                        n.cores = -1){
  cat("Working on it...\n")
  
  args_list <- list(freq.res = freq.res,
                    win.fun = win.fun,
                    min.freq = min.freq,
                    max.freq = max.freq,
                    j = j,
                    noise.red = noise.red,
                    rm.offset = rm.offset)
  
  if(is.null(folder)){
    folder <- getwd()
  }
  setwd(folder)
  
  if(is.null(list)){
    audio.list <- list_waves()
    
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
    
    aci1 <- quiet(do.call(aci, c(list(sound1), args_list)))
    
    tibble(file_name = "filename")  |>  bind_cols(aci1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(aci1)
    
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
                       
                       # Calculate ACI and keep its default output columns
                       aci_result <- quiet(do.call(aci, c(list(sound), args_list)))
                       
                       tibble(file_name = file)  |> 
                         bind_cols(aci_result)
                       
                     }
  stopCluster(cl)
  
  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)
  
  if (save.csv) {
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
