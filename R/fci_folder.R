#' Frequency Cover Indices on Multiple WAV Files in a Folder
#'
#' This function performs a frequency cover indices (FCI) on all WAV files within a specified folder.
#' It processes each sound file in parallel (if specified) and saves the results in a CSV file.
#'
#' @param folder The folder containing the WAV files to analyze.
#' @param output.csv The name of the CSV file where results will be saved. Default is "frequency_cover_results.csv".
#' @param channel The channel to analyze: 'left', 'right', 'mix' (combine both), or 'each' (process left and right channels separately). Default is 'left'.
#' @param rm.offset Logical. Whether to remove the DC offset from the signal. Default is TRUE.
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
#'
#' @examples
#' \dontrun{
#' # Run frequency cover analysis on all WAV files in the folder "sound_files"
#' fci_folder("sound_files", output.csv = "results.csv", channel = "left", n.cores = 4)
#' }

fci_folder <- function(folder,
                      output.csv = "fci_results.csv",
                      channel = 'each',
                      rm.offset = TRUE,
                      hpf = 0,
                      cutoff = -60,
                      freq.res = 100,
                      plot = FALSE,
                      ggplot = FALSE,
                      plot.title = "Frequency Cover Indices",
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


  cat("Evaluating the job...\n")

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
  files <- list.files(path=folder, pattern = ".wav|.WAV")

  filename <- tibble(filename = files)
  nFiles <- length(files)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(files[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  fci1 <- quiet(fci(sound1, channel = 'mix', plot = FALSE))
  tibble::tibble(file_name = "filename") %>% bind_cols(fci1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(fci1)


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

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Set up parallel processing

  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")


  # Define parallel computation
  results <- foreach(file = files, .packages = c("soundecology2", "tuneR", "seewave", "tibble")) %dopar% {

    filename <- basename(file)  # Get file name without path

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


    # Initialize an empty tibble for the results
    result_list <- list()

    results <- fci(wave = sound,
                  channel = channel,
                  rm.offset = rm.offset,
                  hpf = hpf,
                  cutoff = cutoff,
                  freq.res = freq.res,
                  plot = FALSE,
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


    result_list <- list(
      tibble(file_name = filename, channel = channel, results)
    )


    return(do.call(rbind, result_list))
  }

  # Combine all the results into a single tibble
  combined_results <- do.call(rbind, results)

  combined_results <- addMetadata(combined_results)

  # Export results to CSV
  write.csv(combined_results, file = output.csv, row.names = FALSE)

  # Stop parallel cluster
  stopCluster(cl)

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(combined_results)
}


