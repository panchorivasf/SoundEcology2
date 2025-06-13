#' Calculate the Acoustic Diversity Index on the Files in a Folder
#' @description
#' Calculates the Acoustic Diversity Index for all the files in a folder, with extended parameter options.
#' It uses parallel processing with all but one of the available cores.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com)  April 2024.
#'
#' @param folder a path to the folder with audio files to import.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param save.csv logical. Whether to save a csv in the working directory.
#' @param csv.name character vector. When 'save.csv' is TRUE, optionally provide a file name.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param win.fun window function (filter to handle spectral leakage); "bartlett", "blackman", "flattop", "hamming", "hanning", or "rectangle".
#' @param min.freq minimum frequency to compute the spectrogram
#' @param max.freq maximum frequency to compute the spectrogram
#' @param n.bands number of bands to split the spectrogram
#' @param cutoff dB threshold to calculate energy proportions (if normspec = FALSE, set to 5 or above)
#' @param norm.spec logical. Whether to normalize the spectrogram (not recommended) or not (normalized spectrograms with different SNR are not comparable).
#' @param noise.red numeric. noise reduction (subtract median from the amplitude values); 0=none, 1=rows, 2=columns.
#' @param rmoff logical. Whether to remove DC offset before computing ADI (recommended) or not.
#' @param props logical. Whether to store the energy proportion values for each frequency band and channel (default) or not.
#' @param prop.den numeric. Indicates how the energy proportion is calculated.
#' @param db.fs logical; if TRUE, the amplitude scale is expressed as decibels Full Scale (dBFS). Only used when norm = FALSE.
#' @param use.vegan logical; if TRUE, the diversity() function from the vegan package is called to compute Shannon's entropy. Default = FALSE.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#'
#' @return a tibble (data frame) with the ADI values for each channel (if stereo), metadata, and the parameters used for the calculation.
#' @export
#' @import doParallel 
#' @import foreach
#' @import parallel
#' @import seewave
#' @importFrom tuneR readWave
#' @importFrom dplyr bind_cols tibble
#' 
#' @details
#' Options for propden:
#' 1 = The original calculation from the "soundecology" package is applied. The denominator of the proportion equals to all the cells in the same frequency band.
#' 2 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the spectrogram (up to 'max.freq')
#' 3 = A "true Shannon" proportion is calculated, where the "whole population across species" equals the cells above the decibel threshold across the whole spectrogram (up to the Nyquist frequency. This might return a smaller range of values.
#' It uses parallel processing with all but one of the available cores.
#' Optimized to facilitate working with a list of audio files before importing them into R.
#' Modifications by Francisco Rivas (frivasfu@purdue.edu // fcorivasf@gmail.com) April 2024
#'
#' @examples
#' adi_folder("path/to/folder")
adi_folder <- function(folder = NULL,
                       list = NULL,
                       recursive = FALSE,
                       save.csv = TRUE,
                       csv.name = "adi_results.csv",
                       freq.res = 50,
                       win.fun = "hanning",
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
                       use.vegan = FALSE,
                       n.cores = -1) {
  cat("Working on it...\n")
  
  # Store the arguments
  args_list <- list(freq.res = freq.res,
                    win.fun = win.fun,
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
  
  if(is.null(folder)) {
    folder <- getwd()
  }
  setwd(folder)
  
  if(is.null(list)) {
    audio.list <- list.files(pattern = "\\.wav$",
                             recursive = recursive,
                             full.names = TRUE)
  } else {
    audio.list <- list
  }
  
  nFiles <- length(audio.list)
  
  if(nFiles == 0) {
    stop("No WAV files found in the specified directory")
  }
  
  # Setup parallel processing
  if(is.null(n.cores)) {
    num_cores <- 1
  } else if(n.cores == -1) {
    num_cores <- parallel::detectCores() - 1  
  } else {
    num_cores <- n.cores
  }
  if(nFiles < num_cores) {
    num_cores <- nFiles
  }
  
  cl <- makeCluster(num_cores[1])
  registerDoParallel(cl)
  
  # Initialize skipped files log
  skipped_files <- character(0)
  
  # Start processing
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "dplyr", "seewave"),
                     .errorhandling = "pass") %dopar% {
                       
                       tryCatch({
                         
                         # Clean the file path right at the start
                         file <- sub("^\\./", "", file)
                         
                         # Attempt to read file
                         sound <- tryCatch({
                           readWave(file)
                         }, error = function(e) {
                           message(paste("Error reading file:", file, "-", e$message))
                           skipped_files <<- c(skipped_files, file)
                           return(NULL)
                         })
                         
                         if(is.null(sound)) return(NULL)
                         
                         # Validate the wave object
                         if(length(sound@left) == 0) {
                           message(paste("Empty audio file:", file))
                           skipped_files <<- c(skipped_files, file)
                           return(NULL)
                         }
                         
                         # Calculate ADI
                         adi_result <- tryCatch({
                           quiet(do.call(adi, c(list(sound), args_list)))
                         }, error = function(e) {
                           message(paste("Error processing file:", file, "-", e$message))
                           skipped_files <<- c(skipped_files, file)
                           return(NULL)
                         })
                           
                           if(is.null(adi_result)) return(NULL)
                           
                           # Return successful result
                           tibble(file_name = file) |> bind_cols(adi_result)
                           
                       }, error = function(e) {
                         message(paste("Unexpected error with file:", file, "-", e$message))
                         skipped_files <<- c(skipped_files, file)
                         return(NULL)
                       })
                     }
  
  stopCluster(cl)
  
  # Report skipped files
  if(length(skipped_files) > 0) {
    cat("\nThe following files were skipped due to errors:\n")
    cat(paste(skipped_files, collapse = "\n"))
    cat("\n\n")
  }
  
  # Process results
  if(is.null(results) || nrow(results) == 0) {
    stop("No files could be processed successfully")
  }
  
  resultsWithMetadata <- addMetadata(results)
  sensor_id <- unique(resultsWithMetadata$sensor_id)
  
  if(save.csv) {
    resultsWithMetadata$datetime <- format(resultsWithMetadata$datetime, "%Y-%m-%d %H:%M:%S")
    
    if(recursive) {
      parent_dir <- dirname(folder)
      csv_path <- file.path(parent_dir, paste0(sensor_id, "_", csv.name))
    } else {
      csv_path <- paste0(sensor_id, "_", csv.name)
    }
    
    write.csv(resultsWithMetadata, file = csv_path, row.names = FALSE)
  }
  
  cat(paste("Done!\nProcessed", nFiles - length(skipped_files), "of", nFiles, "files\n"))
  cat(paste("Time of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))
  
  return(resultsWithMetadata)
}

