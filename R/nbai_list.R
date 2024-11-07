#' Narrow-Band Activity Index for audio.list in a List
#' @description
#' NBAI describes the relative amount of narrow-band persistent sound activity, like that of Cicadas and Orthopterans. This index can be used to evaluate insect activity and their influence on other soundscape metrics (e.g., summary acoustic indices).
#'
#' @param audio.list Character. A list containing the names of the wave audio.list to analyze.
#' @param folder Character. The path to the folder containing the wave files to analyze.
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. Cutoff threshold defining the sounds that will be analyzed, in dBFS.
#' @param activity.cutoff Numeric. Cutoff percent activity. Only the frequency bands active equal or above this percentage will be considered as "active" in the active band statistics.
#' @param output.csv = Character. Name for the csv output.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#' @param verbose Logical. If TRUE, processing info and results will be printed in the console. 
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
#' @examples nbai_list(wave.list, channel = 'left', plot = TRUE)

nbai_list <- function(audio.list,
                      folder = NULL,
                      channel = "each",
                      hpf = 0,
                      freq.res = 50,
                      cutoff = -60,
                      activity.cutoff = 10, 
                      output.csv = "nbai_results.csv",
                      n.cores = -1,
                      verbose = TRUE) {
  
  if(is.null(folder)){
    folder <- getwd()
  }
  
  setwd(folder)
  filename <- tibble(filename = audio.list)
  naudio.list <- length(audio.list)

  # Set up parallel workers
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1
  }else{
    num_cores <- n.cores
  }
  
  if(naudio.list < num_cores){
    num_cores <- naudio.list
  }
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  

  # Evaluate the expected time of completion when files > 10
  if (length(audio.list) > 10) {
    
    quiet <- function(..., messages=FALSE, cat=FALSE){
      if(!cat){
        tmpf <- tempfile()
        sink(tmpf)
        on.exit({sink(); file.remove(tmpf)})
      }
      out <- if(messages) eval(...) else suppressMessages(eval(...))
      out
    }
    
  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audio.list[1], from = 0, to = 2 , units ='seconds')
  type <- ifelse(sound1@stereo, "stereo", "mono")

  nbai1 <- quiet(nbai(sound1, channel = 'mix'))
  tibble::tibble(file_name = "filename") %>% bind_cols(nbai1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(lubridate::seconds(2.2))

  rm(sound1)
  rm(nbai1)
  # Estimate total time accounting for parallel processing
  estimatedTotalTime <- (timePerFile * naudio.list) / as.numeric(num_cores)
  # Add overhead time
  adjustedTotalTime <- estimatedTotalTime
  # Calculate the end time
  expectedCompletionTime <- Sys.time() + adjustedTotalTime
  # Set up parallel processing

  cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
  cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
  cat("Analyzing", naudio.list, "files using", num_cores, "cores... \n")
  } else {
    cat("Analyzing", naudio.list, "files using", num_cores, "cores... \n")
    
  }

  # Define parallel computation
  results <- foreach(file = audio.list, .packages = c("tuneR", "seewave", "tibble")) %dopar% {

    filename <- basename(file)  # Get file name without path
    # Read audio file
    audio <- readWave(file)

    # Initialize an empty tibble for the results
    result_list <- list()
    
    nbai <- nbai(audio, 
                 channel = channel,
                 hpf = hpf,
                 freq.res = freq.res,
                 cutoff = cutoff,
                 activity.cutoff = activity.cutoff,
                 verbose = verbose)$summary
    
    result_list <- list(tibble(file_name = filename, 
             nbai)
    )
    

    return(do.call(rbind, result_list))
  }

  # Combine all the results into a single tibble
  combined_results <- do.call(rbind, results)

  combined_results <- addMetadata(combined_results)
  
  # combined_results <- apply(combined_results,2,as.character)

  # Export results to CSV
  write.csv(combined_results, file = output.csv, row.names = FALSE)

  # Stop parallel cluster
  stopCluster(cl)

  cat(paste("Done!\nTime of completion:", format(Sys.time(), "%H:%M:%S"), "\n\n"))

  return(combined_results)
}



