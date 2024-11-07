#' Acoustic Evenness Index - folder input
#' @description
#' Calculates the Acoustic Evenness Index for all the files in a folder, with extended parameter options.
#' It uses parallel processing with all but one of the available cores.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com)  April 2024.
#'
#' @param folder a path to the folder with audio files to import.
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param frew.res the frequency resolution  (Hz per bin) to use. From this value the window length for the FFT will be calculated (sampling rate / frequency resolution).
#' @param w.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
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
#' @import parallel
#' @import tuneR
#' @import tidyverse
#' @import seewave
#' @import lubridate


#'
#' @examples
#' aei_folder("path/to/folder")
aei_folder <- function (folder,
                        save.csv = TRUE,
                        csv.name = "aei_results.csv",
                        freq.res = 50,
                        w.fun = "hanning",
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
  
  if(is.null(folder)){
    folder <- getwd()
  }


  cat("Evaluating the job...\n\n")

  #  Quiet function from SimDesign package to run functions without printing
  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      tmpf <- tempfile()
      sink(tmpf)
      on.exit({sink(); file.remove(tmpf)})
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }

  setwd(folder)
  audiolist <- list.files(path=folder, pattern = ".wav|.WAV")

  fileName <- tibble(file_name = audiolist)
  nFiles <- length(audiolist)


  args_list <- list(freq.res = freq.res, w.fun = w.fun, min.freq = min.freq,
                    max.freq = max.freq, n.bands = n.bands, cutoff = cutoff,
                    norm.spec = norm.spec, noise.red = noise.red, rm.offset = rm.offset,
                    props = props, prop.den = prop.den, db.fs = db.fs)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audio.list[1], from = 0, to = 2 , units ='seconds')
  type <- ifelse(sound1@stereo, "stereo", "mono")

  aei1 <- quiet(aei(sound1, args_list$freq.res, args_list$w.fun, args_list$min.freq,
                    args_list$max.freq, args_list$n.bands, args_list$cutoff,
                    args_list$norm.spec, args_list$noise.red, args_list$rm.offset,
                    args_list$props, args_list$prop.den, args_list$db.fs))

  tibble(file_name = "filename") %>% bind_cols(aei1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(aei1)

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


  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Setup parallel back-end
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)

  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")



  # Start loop
  results <- foreach(file = audiolist, .combine = rbind,
                     .packages = c("tuneR", "tidyverse", "seewave")) %dopar% {

                       # Try to read the sound file
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


                       # Calculate AEI and keep its default output columns
                       aei <- aei(sound, freq.res = args_list$freq.res, w.fun = args_list$w.fun,
                                  min.freq = args_list$min.freq, max.freq = args_list$max.freq,
                                  n.bands = args_list$n.bands, cutoff = args_list$cutoff,
                                  norm.spec = args_list$norm.spec, noise.red = args_list$noise.red,
                                  rm.offset = args_list$rm.offset, props = args_list$props,
                                  prop.den = args_list$prop.den, db.fs = args_list$db.fs)

                       # Combine the results for each file into a single row
                       tibble(file_name = file) %>%
                         bind_cols(aei)


                     }


  # Combine results with metadata and return
  resultsWithMetadata <- addMetadata(results)


  stopCluster(cl)

  if(save.csv == TRUE){
    resultsWithMetadata$datetime <- format(resultsWithMetadata$datetime, "%Y-%m-%d %H:%M:%S")
    write.csv(resultsWithMetadata, csv.name, row.names = FALSE)
  }

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  cat("Acoustic Evenness Index:")

  return(resultsWithMetadata)

}

