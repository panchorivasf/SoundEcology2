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
  
  start_time <- Sys.time()
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
  original_wd <- getwd()
  setwd(folder)
  
  if(is.null(list)) {
    audio.list <- list.files(pattern = "\\.wav$",
                             recursive = recursive,
                             full.names = TRUE)
    # Normalize paths
    audio.list <- normalizePath(audio.list, winslash = "/")
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
  
  # Initialize progress tracking
  processed_count <- 0
  skipped_files <- character(0)
  cat(sprintf("Processing %d files...\n", nFiles))
  
  # Start processing
  results <- foreach(file = audio.list, .combine = rbind,
                     .packages = c("tuneR", "dplyr", "seewave"),
                     .errorhandling = "pass") %dopar% {
                       
                       clean_file <- basename(file)  # Store just filename without path
                       full_path <- file
                       
                       tryCatch({
                         sound <- tryCatch({
                           readWave(full_path)
                         }, error = function(e) {
                           message(paste("Error reading file:", clean_file, "-", e$message))
                           skipped_files <<- c(skipped_files, clean_file)
                           return(NULL)
                         })
                         
                         if(is.null(sound)) return(NULL)
                         
                         adi_result <- tryCatch({
                           quiet(do.call(adi, c(list(sound), args_list)))
                         }, error = function(e) {
                           message(paste("Error processing file:", clean_file, "-", e$message))
                           skipped_files <<- c(skipped_files, clean_file)
                           return(NULL)
                         })
                           
                           if(is.null(adi_result)) return(NULL)
                           
                           # Update progress (thread-safe)
                           processed_count <<- processed_count + 1
                           if (processed_count %% max(1, floor(nFiles/10)) == 0) {
                             cat(sprintf("..%d/%d (%.0f%%) \n", 
                                         processed_count, nFiles, 
                                         processed_count/nFiles*100))
                           }
                           
                           tibble(file_name = clean_file) |> 
                             bind_cols(adi_result)
                           
                       }, error = function(e) {
                         message(paste("Unexpected error with file:", clean_file, "-", e$message))
                         skipped_files <<- c(skipped_files, clean_file)
                         return(NULL)
                       })
                     }
  
  stopCluster(cl)
  setwd(original_wd)  # Restore original working directory
  
  # Report skipped files
  if(length(skipped_files) > 0) {
    cat("\nSkipped", length(skipped_files), "files due to errors:\n")
    cat(paste(head(skipped_files, 10), collapse = "\n"))
    if(length(skipped_files) > 10) cat("\n... (", length(skipped_files)-10, "more)")
    cat("\n\n")
  }
  
  if(is.null(results) || nrow(results) == 0) {
    stop("No files could be processed successfully")
  }
  
  resultsWithMetadata <- addMetadata(results)
  sensor_id <- unique(resultsWithMetadata$sensor_id)
  
  if(save.csv) {
    resultsWithMetadata$datetime <- format(resultsWithMetadata$datetime, "%Y-%m-%d %H:%M:%S")
    
    csv_path <- if(recursive) {
      file.path(dirname(folder), paste0(sensor_id, "_", csv.name))
    } else {
      file.path(folder, paste0(sensor_id, "_", csv.name))
    }
    
    # Ensure directory exists
    if(!dir.exists(dirname(csv_path))) {
      dir.create(dirname(csv_path), recursive = TRUE)
    }
    
    write.csv(resultsWithMetadata, file = csv_path, row.names = FALSE)
    cat("Results saved to:", normalizePath(csv_path), "\n")
  }
  
  cat(paste(
    sprintf("\nDone! Processed %d/%d files (%.1f%%)", 
            nFiles - length(skipped_files), nFiles,
            100*(nFiles - length(skipped_files))/nFiles),
    paste("\nTime elapsed:", format(round(Sys.time() - start_time, 1))),
    "\n"
  ))
  
  return(resultsWithMetadata)
}

