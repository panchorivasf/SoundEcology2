#' Narrow-Band Activity Index for Files in a Folder
#' @description
#' NBAI describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary acoustic indices).
#'
#' @param folder Character. The path to the folder containing the wave files to analyze.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.
#' @param activity.cutoff Numeric. Cutoff percent activity. Only the frequency bands active equal or above this percentage will be considered as "active" in the active band statistics.
#' @param plot.binary.spec Logical. Whether to plot the binary spectrogram used for the analysis. Allowed only when Wave is mono or when one channel is selected from a stereo file.
#' @param dark.plot Logical. If true (default) a the binary spectrogram will have a black background.
#' @param plot.activity Logical. If true, a barplot depicting the percent activity in each frequency bin will be returned. Default is FALSE.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#'
#' @return A list containing: 1) A binary spectrogram (if mono), 2) tibble with the Narrow-Band Activity Index (NBI) summary statistics, and 3) a tibble with NBI spectral, which number of rows equals the number of frequency bins in the analysis.
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols select
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom tuneR readWave
#' @importFrom lubridate seconds
#' @importFrom utils write.csv
#'
#' @examples nbai(wave, channel = 'left', plot = TRUE, verbose = TRUE)

nbai_folder <- function(folder,
                        channel = "each",
                        hpf = 250,
                        freq.res = 50,
                        cutoff = -60,
                        activity.cutoff = 10,
                        output.csv = "nbai_results.csv",
                        n.cores = -1) {
  
  if(is.null(folder)){
    folder <- getwd()
  }

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

  nbai1 <- quiet(nbai(sound1, channel = channel))
  tibble::tibble(file_name = "filename") %>% bind_cols(nbai1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(nbai1)


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
  results <- foreach(file = files, .packages = c("tuneR", "seewave", "tibble")) %dopar% {

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

    # Process selected channel
    nbai <- nbai(sound,
                 channel = channel,
                 hpf = hpf,
                 freq.res = freq.res,
                 cutoff = cutoff,
                 activity.cutoff = activity.cutoff)


    result_list <- list(
      tibble(file_name = filename, channel = channel, nbai)
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



