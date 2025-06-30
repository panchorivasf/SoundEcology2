#' Acoustic Evenness Index - folder input
#' @description
#' Calculates the Acoustic Evenness Index for all the files in a folder, with extended parameter options.
#' It uses parallel processing with all but one of the available cores.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com)  April 2024.
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
#' @param frew.res the frequency resolution  (Hz per bin) to use. From this value the window length for the FFT will be calculated (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram
#' @param max.freq maximum frequency to compute the spectrogram
#' @param n.bands number of bands to split the spectrogram
#' @param cutoff dB threshold to calculate energy proportions (if norm.spec = FALSE, set to 5 or above)
#' @param norm.spec logical. Whether to normalize the spectrogram (not recommended) or not (normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric. noise reduction (subtract median from the amplitude values); 1=rows, 2=columns, 3=none.
#' @param rm.offset logical. Whether to remove DC offset before computing aei (recommended) or not.
#' @param props logical. Whether to store the energy proportion values for each frequency band and channel (default) or not.
#' @param prop.den numeric. Indicates how the energy proportion is calculated by manipulating the denominator.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#'
#' @return a tibble (data frame) with the aei values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @details
#' Options for prop.den:
#' 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion equals to all the cells in the same frequency band.
#' 2 = The "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max.freq')
#' 3 = The "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a list of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @import doParallel 
#' @import foreach
#' @import seewave
#' @importFrom parallel detectCores makeCluster
#' @importFrom tuneR readWave
#' @importFrom dplyr bind_cols tibble
#'
#' @examples
#' aei_folder("path/to/folder")
aei_folder <- function (folder = NULL,
                        list = NULL,
                        recursive = FALSE,
                        start = 0,
                        end = 1,
                        unit = "minutes",
                        save.csv = TRUE,
                        csv.name = "aei_results.csv",
                        freq.res = 50,
                        win.fun = "hanning",
                        min.freq = 0,
                        max.freq = 10000,
                        n.bands = 10,
                        cutoff = -60,
                        norm.spec = FALSE,
                        noise.red = 0,
                        rm.offset = TRUE,
                        props = TRUE,
                        prop.den = 1,
                        db.fs = TRUE,
                        n.cores = -1){
  cat("Working on it...\n")
  
  args_list <- list(freq.res = freq.res, win.fun = win.fun, min.freq = min.freq,
                    max.freq = max.freq, n.bands = n.bands, cutoff = cutoff,
                    norm.spec = norm.spec, noise.red = noise.red, rm.offset = rm.offset,
                    props = props, prop.den = prop.den,  db.fs = db.fs)
  
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
  
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1  
  }else{
    num_cores <- n.cores
  }
  if(n.files < num_cores){
    num_cores <- n.files
  }
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  
  if (n.files > 10){
    cat("Evaluating the job...\n\n")
    
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1],
                       from = start,
                       to = end,
                       units = unit)
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    aei1 <- quiet(do.call(aei, c(list(sound1), args_list)))
    
    tibble(file_name = "filename") |> bind_cols(aei1)
    
    timePerFile <-  Sys.time() - startTime
    timePerFile <- timePerFile + 2.2
    
    rm(sound1)
    rm(aei1)
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * n.files) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  } else {
    sound1 <- readWave(audio.list[1], to = 2, units = "seconds")
    type <- ifelse(sound1@stereo, "stereo", "mono")
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
                         message(paste("Error reading file:", 
                                       file, "Skipping to the next file."))
                         return(NULL)
                       })
                       
                       if (is.null(sound)) return(NULL)
                       
                       
                       # Calculate AEI and keep its default output columns
                       aei_result <- quiet(do.call(aei, c(list(sound), args_list)))
                       
                       # Combine the results for each file into a single row
                       tibble(file_name = file) |>
                         bind_cols(aei_result)
                       
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
