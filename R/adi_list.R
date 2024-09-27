#' Calculate the Acoustic Diversity Index on the Files in a List
#'
#' @description
#' Calculates the Acoustic Diversity Index for all the files in a list, with extended parameter options.
#' It uses parallel processing with all but one of the available cores.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com)  April 2024.
#'
#' @param audio.list a list of audio files to analyze.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram
#' @param max.freq maximum frequency to compute the spectrogram
#' @param n.bands number of bands to split the spectrogram
#' @param cutoff dB threshold to calculate energy proportions (if norm.spec = FALSE, set to 5 or above)
#' @param norm.spec logical. Whether to normalize the spectrogram (not recommended) or not (normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric. noise reduction (subtract median from the amplitude values); 1=rows, 2=columns, 3=none.
#' @param rm.offset logical. Whether to remove DC offset before computing ADI (recommended) or not.
#' @param props logical. Whether to store the energy proportion values for each frequency band and channel (default) or not.
#' @param prop.den numeric. Indicates how the energy proportion is calculated.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).

#'
#' @return a tibble (data frame) with the ADI values for each channel (if stereo), metadata, and the parameters used for the calculation.
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
#' Options for prop.den:
#' 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion equals to all the cells in the same frequency band.
#' 2 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max.freq')
#' 3 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a list of audio files before importing them into R.
#' The working directory should be set to the folder containing the files.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' files <- list.files(pattern=".wav|.WAV")
#' adi_list(files[1:5])

adi_list <- function (audio.list,
                      save.csv = TRUE,
                      csv.name = "adi_results.csv",
                      freq.res = 50,
                      w.fun = "hanning",
                      min.freq = 0,
                      max.freq = 10000,
                      n.bands = 10,
                      cutoff = -60,
                      norm.spec = FALSE,
                      noise.red = 0,
                      rm.offset = TRUE,
                      props = FALSE,
                      prop.den = 1,
                      db.fs = TRUE,
                      n.cores = -1){


  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      tmpf <- tempfile()
      sink(tmpf)
      on.exit({sink(); file.remove(tmpf)})
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }





  fileName <- tibble(file_name = audio.list)
  nFiles <- length(audio.list)

  
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1  # Detect available cores
  }else{
    num_cores <- n.cores
  }
  
  if(nFiles < num_cores){
    num_cores <- nFiles
  }
  
  # Setup parallel back-end
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  args_list <- list(freq.res = freq.res,
                    w.fun = w.fun,
                    min.freq = min.freq,
                    max.freq = max.freq,
                    n.bands = n.bands,
                    cutoff = cutoff,
                    norm.spec = norm.spec,
                    noise.red = noise.red,
                    rm.offset = rm.offset,
                    props = props,
                    prop.den = prop.den,
                    db.fs = db.fs)
  
  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()
  
  if(nFiles>10){
    cat("Evaluating the job...\n\n")
    
    
    sound1 <- readWave(audio.list[1])
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    adi1 <- quiet(adi(sound1, args_list$freq.res, args_list$w.fun, args_list$min.freq,
                      args_list$max.freq, args_list$n.bands, args_list$cutoff,
                      args_list$norm.spec, args_list$noise.red, args_list$rm.offset,
                      args_list$props, args_list$prop.den, args_list$db.fs))
    
    tibble(file_name = "filename") %>% bind_cols(adi1)
    
    # Assess how long it takes to parse 1 file
    timePerFile <-  Sys.time() - startTime
    # Add overhead per file
    timePerFile <- timePerFile + as.numeric(seconds(2.2))
    
    rm(sound1)
    rm(adi1)
    
    
    
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  }
  
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  

  # Start loop
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Import the sounds
                       sound <- readWave(file)

                       # Calculate ADI and keep its default output columns
                       adi_result <- adi(sound, freq.res = args_list$freq.res, w.fun = args_list$w.fun,
                                         min.freq = args_list$min.freq, max.freq = args_list$max.freq,
                                         n.bands = args_list$n.bands, cutoff = args_list$cutoff,
                                         norm.spec = args_list$norm.spec, noise.red = args_list$noise.red,
                                         rm.offset = args_list$rm.offset, props = args_list$props,
                                         prop.den = args_list$prop.den, db.fs = args_list$db.fs)

                       # Log for debugging
                       print(paste("Processing:", file))
                       print(paste("Left channel value:", adi_result$value_l))
                       print(paste("Right channel value:", adi_result$value_r))

                       # Combine the results for each file into a single row
                       result <- tibble(file_name = file) %>%
                         bind_cols(adi_result)

                     }


  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)


  stopCluster(cl)

  if(save.csv == TRUE){
    write.csv(resultsWithMetadata, csv.name, row.names = FALSE)
  }

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(resultsWithMetadata)

}
