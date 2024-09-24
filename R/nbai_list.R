#' Narrow-Band Activity Index for Files in a List
#' @description
#' NBAI describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary acoustic indices).
#'
#' @param files Character. A list containing the names of the wave files to analyze.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.
#' @param activity.cutoff Numeric. Cutoff percent activity. Only the frequency bands active equal or above this percentage will be considered as "active" in the active band statistics.
#' @param plot.binary.spec Logical. Whether to plot the binary spectrogram used for the analysis. Allowed only when Wave is mono or when one channel is selected from a stereo file.
#' @param dark.plot Logical. If true (default) a the binary spectrogram will have a black background.
#' @param plot.activity Logical. If true, a barplot depicting the percent activity in each frequency bin will be returned. Default is FALSE.
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#'
#' @return A list containing: 1) A binary spectrogram (if mono), 2) tibble with the Narrow-Band Activity Index (NBI) summary statistics, and 3) a tibble with NBI spectral, which number of rows equals the number of frequency bins in the analysis.
#' @export
#'
#' @import tuneR
#' @import seewave
#' @import patchwork
#' @import ggplot2
#' @import dplyr
#'
#' @examples nbai(wave, channel = 'left', plot = TRUE, verbose = TRUE)

nbai_list <- function(files, channel = "each", output.csv = "nbai_results.csv", ncores = NULL) {
  library(tibble)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(tuneR)
  library(seewave)
  library(tidyverse)

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

  filename <- tibble(filename = files)
  nFiles <- length(files)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(files[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  nbai1 <- quiet(nbai(sound1, channel = 'mix'))
  tibble::tibble(file_name = "filename") %>% bind_cols(nbai1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(nbai1)


  if(is.null(ncores)){
    num_cores <- 1
  }else if(ncores == -1){
    num_cores <- parallel::detectCores() - 1  # Detect available cores
  }else{
    num_cores <- ncores
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

    # Read audio file
    audio <- readWave(file)

    # Initialize an empty tibble for the results
    result_list <- list()

    # Handle different channel selections
    if (channel == "each") {
      # Process each channel separately
      nbai_left <- nbai(audio, cutoff = cutoff, plot.title = filename, channel = "left")$summary
      nbai_right <- nbai(audio, cutoff = cutoff, plot.title = filename, channel = "right")$summary

      nbai <- nbai_left %>%
        dplyr::select(-c(value))

      # Combine the results into one row with 'value_l' and 'value_r' columns
      result_list <- list(tibble(
        file_name = filename,
        value_l = nbai_left$value,
        value_r = nbai_right$value,
        value_avg = round((nbai_left$value + nbai_right$value) / 2, 3),
        nbai
      ))


    } else {
      # Process selected channel ("left", "right", or "mix")
      nbai <- nbai(audio, cutoff = cutoff, plot.title = filename, channel = channel)$summary
      result_list <- list(
        tibble(file_name = filename, channel = channel, nbai)
      )
    }

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



