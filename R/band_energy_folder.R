#' Batch Calculate Band Energy for Multiple Audio Files
#'
#' Processes multiple audio files in parallel to calculate the distribution of
#' acoustic energy across frequency bands for each file. Results can be exported
#' to CSV with optional metadata extraction.
#'
#' @param folder Character. Path to folder containing audio files. If NULL,
#'   uses current working directory. Default is NULL
#' @param list Character vector. Optional list of specific file paths to process.
#'   If NULL, searches for all WAV files in folder. Default is NULL
#' @param recursive Logical. Should subdirectories be searched for audio files?
#'   Default is FALSE
#' @param start Numeric. Starting point for audio extraction. Default is 0
#' @param end Numeric. Ending point for audio extraction. Default is 1
#' @param unit Character. Time unit for start and end parameters. Options are
#'   "samples", "seconds", "minutes", "hours". Default is "minutes"
#' @param save.csv Logical. Should results be saved to a CSV file? Default is TRUE
#' @param save.to Character. Directory path where CSV should be saved. If NULL,
#'   saves to folder directory. Default is NULL
#' @param csv.name Character. Base name for output CSV file (without .csv extension).
#'   Default is "band_energy_results"
#' @param freq.res Numeric. Frequency resolution in Hz. Passed to band_energy().
#'   Default is 50 Hz
#' @param win.fun Character. Window function for spectrogram. Passed to band_energy().
#'   Default is "hanning"
#' @param min.freq Numeric. Minimum frequency in Hz. Passed to band_energy().
#'   Default is 0 Hz
#' @param max.freq Numeric. Maximum frequency in Hz. Passed to band_energy().
#'   Default is 10000 Hz
#' @param n.bands Integer. Number of frequency bands. Passed to band_energy().
#'   Default is 10
#' @param cutoff Numeric. Decibel threshold. Passed to band_energy().
#'   Default is -60 dB
#' @param norm.spec Logical. Normalize spectrogram? Passed to band_energy().
#'   Default is FALSE
#' @param noise.red Integer. Noise reduction method (0, 1, or 2). Passed to
#'   band_energy(). Default is 0
#' @param rm.offset Logical. Remove DC offset? Passed to band_energy().
#'   Default is TRUE
#' @param db.fs Logical. Use dBFS conversion? Passed to band_energy().
#'   Default is TRUE
#' @param n.cores Integer. Number of CPU cores to use for parallel processing.
#'   If -1, uses all available cores minus 1. If NULL or 1, runs sequentially.
#'   Default is -1
#' @param verbose Logical. Should diagnostic information be printed during
#'   processing? Includes a test run on the first file. Default is TRUE
#'
#' @return A data frame with one row per successfully processed audio file,
#'   containing:
#'   \itemize{
#'     \item file_name: Full path to the audio file
#'     \item Metadata columns from band_energy()
#'     \item Energy proportion columns for each frequency band
#'     \item Additional metadata columns if addMetadata() function is available
#'       (e.g., sensor_id, datetime extracted from filename)
#'   }
#'
#' @details
#' This function automates band energy calculation across multiple files using
#' parallel processing for efficiency. Key features include:
#' \itemize{
#'   \item Parallel processing using foreach and doParallel packages
#'   \item Automatic testing on first file before batch processing
#'   \item Robust error handling - files that fail are skipped with warnings
#'   \item Optional metadata extraction via addMetadata() function
#'   \item Progress reporting and timing information
#'   \item Flexible file selection via folder scanning or explicit file list
#' }
#'
#' The function requires the list_waves() helper function to scan for audio files.
#' If an addMetadata() function exists in the environment, it will be applied to
#' extract additional metadata from filenames (e.g., sensor IDs, timestamps).
#'
#' @section Parallel Processing:
#' When n.cores > 1, files are processed in parallel. The band_energy() function
#' must be available in the global environment. For debugging, set n.cores = 1
#' to see detailed error messages for each file.
#'
#' @section Output Files:
#' When save.csv = TRUE, output filename is constructed as:
#' \itemize{
#'   \item If single sensor: "[sensor_id]_[csv.name].csv"
#'   \item If multiple sensors or no sensor_id: "[csv.name].csv"
#' }
#'
#' @examples
#' \dontrun{
#' # Process all WAV files in current directory
#' results <- band_energy_folder()
#'
#' # Process files in specific folder with custom parameters
#' results <- band_energy_folder(
#'   folder = "path/to/audio",
#'   recursive = TRUE,
#'   n.bands = 20,
#'   max.freq = 8000,
#'   n.cores = 4
#' )
#'
#' # Process specific files
#' my_files <- c("file1.wav", "file2.wav", "file3.wav")
#' results <- band_energy_folder(list = my_files, n.cores = -1)
#'
#' # Process first 30 seconds of each file
#' results <- band_energy_folder(
#'   start = 0,
#'   end = 30,
#'   unit = "seconds"
#' )
#'
#' # Debug mode - sequential processing with verbose output
#' results <- band_energy_folder(n.cores = 1, verbose = TRUE)
#' }
#'
#' @seealso \code{\link{band_energy}} for single file processing
#'
#' @export
band_energy_folder <- function(folder = NULL,
                               list = NULL,
                               recursive = FALSE,
                               start = 0,
                               end = 1,
                               unit = "minutes",
                               save.csv = TRUE,
                               save.to = NULL,
                               csv.name = "band_energy_results",
                               freq.res = 50,
                               win.fun = "hanning",
                               min.freq = 0,
                               max.freq = 10000,
                               n.bands = 10,
                               cutoff = -60,
                               norm.spec = FALSE,
                               noise.red = 0,
                               rm.offset = TRUE,
                               db.fs = TRUE,
                               n.cores = -1,
                               verbose = FALSE) {
  
  start_time <- Sys.time()
  cat("Working on it...\n")
  
  # Store the arguments for band_energy
  args_list <- list(freq.res = freq.res,
                    win.fun = win.fun,
                    min.freq = min.freq,
                    max.freq = max.freq,
                    n.bands = n.bands,
                    cutoff = cutoff,
                    norm.spec = norm.spec,
                    noise.red = noise.red,
                    rm.offset = rm.offset,
                    db.fs = db.fs)
  
  original_wd <- getwd()
  
  # Set working directories
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
  
  # Get list of audio files
  if(is.null(list)) {
    audio.list <- list_waves(recursive = recursive)
  } else {
    audio.list <- list
  }
  
  nFiles <- length(audio.list)
  
  if(nFiles == 0) {
    stop("No WAV files found in the specified directory")
  }
  
  cat(sprintf("Found %d audio files\n", nFiles))
  
  # TEST ON FIRST FILE BEFORE PARALLEL PROCESSING
  if(verbose) {
    cat("\nTesting on first file...\n")
    test_file <- audio.list[1]
    cat(sprintf("Test file: %s\n", test_file))
    
    tryCatch({
      test_sound <- readWave(test_file, from = start, to = end, units = unit)
      cat(sprintf("  - Successfully read audio: %d samples, %d Hz, %s\n", 
                  length(test_sound@left), test_sound@samp.rate,
                  ifelse(test_sound@stereo, "stereo", "mono")))
      
      test_result <- do.call(band_energy, c(list(wave = test_sound), args_list))
      cat(sprintf("  - Successfully calculated band energy: %d columns\n", 
                  ncol(test_result)))
      cat("  - Column names:", paste(head(colnames(test_result), 5), collapse = ", "), "...\n")
      cat("✓ Test successful!\n\n")
      
      rm(test_sound, test_result)
      gc(verbose = FALSE)
    }, error = function(e) {
      cat("✗ Test FAILED with error:\n")
      cat(sprintf("  %s\n\n", e$message))
      cat("Full error:\n")
      print(e)
      stop("Cannot proceed - fix the error above first")
    })
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
  
  cat(sprintf("Processing %d files using %d cores...\n", nFiles, num_cores))
  
  cl <- parallel::makeCluster(num_cores[1])
  doParallel::registerDoParallel(cl)
  
  # Export the band_energy function to cluster
  parallel::clusterExport(cl, "band_energy", envir = environment())
  
  # Process files in parallel
  results <- foreach::foreach(file = audio.list, 
                              .combine = rbind,
                              .packages = c("tuneR", "seewave"),
                              .errorhandling = "remove",
                              .verbose = FALSE) %dopar% {
                                
                                tryCatch({
                                  # Read audio file
                                  sound <- tuneR::readWave(file, 
                                                           from = start,
                                                           to = end,
                                                           units = unit)
                                  
                                  # Calculate band energy
                                  band_result <- do.call(band_energy, c(list(wave = sound), args_list))
                                  
                                  # Add filename to results
                                  result <- data.frame(file_name = file, stringsAsFactors = FALSE)
                                  result <- cbind(result, band_result)
                                  
                                  # Clean up
                                  rm(sound, band_result)
                                  
                                  return(result)
                                  
                                }, error = function(e) {
                                  # Return NULL on error (will be removed by .errorhandling = "remove")
                                  return(NULL)
                                })
                              }
  
  parallel::stopCluster(cl)
  setwd(original_wd)  
  
  # Check if any files were processed successfully
  if(is.null(results) || nrow(results) == 0) {
    cat("\n✗ ERROR: No files could be processed successfully\n")
    cat("Possible issues:\n")
    cat("  1. band_energy() function not exported to parallel workers\n")
    cat("  2. Audio files are corrupted or unreadable\n")
    cat("  3. Parameter mismatch between test and parallel processing\n")
    cat("  4. Memory issues on parallel workers\n")
    cat("\nTry setting n.cores = 1 to see detailed error messages\n")
    stop("No files could be processed successfully")
  }
  
  # Add metadata if addMetadata function exists
  if(exists("addMetadata")) {
    results <- addMetadata(results)
  }
  
  # Count successful vs skipped files
  n_processed <- nrow(results)
  n_skipped <- nFiles - n_processed
  
  # Save results to CSV if requested
  if(save.csv == TRUE){
    
    # Format datetime if it exists
    if("datetime" %in% colnames(results)) {
      results$datetime <- format(results$datetime, "%Y-%m-%d %H:%M:%S")
    }
    
    # Generate filename
    if("sensor_id" %in% colnames(results)) {
      sensor <- unique(results$sensor_id)
      if (length(sensor) == 1){
        csv_file <- file.path(save.to, paste0(sensor, "_", csv.name, ".csv"))
      } else {
        csv_file <- file.path(save.to, paste0(csv.name, ".csv"))
      }
    } else {
      csv_file <- file.path(save.to, paste0(csv.name, ".csv"))
    }
    
    write.csv(results, file = csv_file, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", csv_file))
  }
  
  # Print summary
  cat(paste(
    sprintf("\n✓ Done! Processed %d/%d files (%.1f%%)", 
            n_processed, nFiles, 100 * n_processed / nFiles),
    sprintf("Time elapsed: %s", format(round(Sys.time() - start_time, 1))),
    "\n"
  ))
  
  if(n_skipped > 0) {
    cat(sprintf("⚠ Skipped %d files due to errors\n", n_skipped))
  }
  
  return(results)
}