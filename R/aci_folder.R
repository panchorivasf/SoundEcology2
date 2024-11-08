#' Acoustic Complexity Index - folder input
#'
#' @param folder a path to the folder with audio files to import.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param freq.res numeric. The frequency resolution to use (Hz per bin) which will determine the window length for the FFT (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to use when calculating the value, in Hertz. Default = 0.
#' @param max.freq maximum frequency to use when calculating the value, in Hertz. Default = NA (Nyquist).
#' @param j the cluster size, in seconds. Default = NA (Duration of the audio file).
#' @param noise.red numeric; controls the application of noise reduction. If set to 1, noise reduction is applied to each row by subtracting the median from the amplitude values. If set to 2, noise reduction is applied to each column similarly. If set to 0, noise reduction is not applied.
#' @param rm.offset logical; if set to TRUE, the function will remove DC offset before computing ADI. Default = TRUE.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).


#' @return A tibble (data frame) with the ACI values for each channel (if stereo), metadata, and the parameters used for the calculation.
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
#' Optimized to facilitate working with a folder of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' aci_folder(path/to/folder)

aci_folder <- function (folder = NULL,
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
  
  audio.list <- list_waves()
  nFiles <- length(audio.list)
  
  # Declare the number of cores to be used 
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1 # Leave one core free 
  }else{
    num_cores <- n.cores
  } 
  if (num_cores > nFiles){
    num_cores <- nFiles  # Limit the number of cores to the number of files
  }
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  # Evaluate the duration of the analysis if nFiles > 10
  if(nFiles>10){
    cat("Evaluating the job...\n\n")
    
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1])
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    aci1 <- quiet(do.call(aci, c(list(sound1), args_list)))
    
    tibble(file_name = "filename")  |>  bind_cols(aci1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + as.numeric(seconds(2.2))
    
    rm(sound1)
    rm(aci1)
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
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
  
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  
  # Start loop
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {
                       
                       # Try to read the sound file, handle errors gracefully
                       sound <- tryCatch({
                         readWave(file)
                       }, error = function(e) {
                         message(paste("Error reading file:", file, "Skipping to the next file."))
                         return(NULL) # Skip this iteration and continue with the next file
                       })
                       
                       # Skip processing if the sound is NULL (i.e., readWave failed)
                       if (is.null(sound)) {
                         return(NULL)
                       }
                       
                       # Calculate ACI and keep its default output columns
                       aci1 <- quiet(do.call(aci, c(list(sound), args_list)))
                       
                       # Combine the results for each file into a single row
                       tibble(file_name = file)  |> 
                         bind_cols(aci)
                       
                     }
  
  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)
  
  stopCluster(cl)
  
  if(save.csv == TRUE){
    resultsWithMetadata$datetime <- format(resultsWithMetadata$datetime, "%Y-%m-%d %H:%M:%S")
    write.csv(resultsWithMetadata, csv.name, row.names = FALSE)
  }
  
  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(resultsWithMetadata)
  
}
