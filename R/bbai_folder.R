#' Calculate BBAI for all the Files in a Folder
#'
#' This function processes an sound signal to detect broadband activity by identifying 'clicks' based on 
#' time-frame-wise (i.e., column-wise) amplitude changes in the spectrogram. It computes statistics related 
#' to click height, variance, and centroid frequency, and can plot a spectrogram with detected clicks highlighted. 
#' The function also classifies whether the signal contains noise or insect based on the variance and centroid 
#' frequencies of the clicks.
#'
#' @param folder Character. The path to a folder with the wave files to analyze.
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @param start numerical. Where to start reading the Wave. 
#' @param end numerical. Where to end reading the Wave.
#' @param unit character. Unit of measurement for 'start' and 'end'. Options are
#' 'samples', 'seconds', 'minutes', 'hours'. Default is 'minutes'.
#' @param list An optional list (subset) of files in the folder to analyze. If provided, 
#' files outside the list will be excluded. 
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either 
#' "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". 
#' If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless 
#' signals of interest are below that threshold.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency 
#' bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. The amplitude threshold (in dBFS) for removing low-amplitude values in the spectrogram. 
#' Default is `-50`.
#' @param click.height Numeric. The minimum height (in frequency bins) for a detected click to be kept. 
#' Default is `10`.
#' @param difference Numeric. The maximum difference in amplitude between adjacent frequency bins to be 
#' considered part of a single 'click'. Default is `20`.
#' @param gap.allowance Numeric. The size of gaps (in frequency bins) allowed between contiguous parts of a click.
#'  Default is `2`. Gaps larger than this value will split clicks.
#' @param spectrogram Logical. Should a spectrogram with highlighted clicks be plotted? Default is `TRUE`.
#' @param dark.plot Logical. Should the plot use a dark theme (black background)? Default is `FALSE`.
#' @param plot.title Character. The title for the plot, if `plot` is `TRUE`. Default is `NULL`.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. 
#' Default is NULL (single-core processing).
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.
#' @param output.csv Character. Name for the csv file. Default is "bbai_results.csv".
#' @param n.cores Numeric. Number of cores to be used in parallel. Use -1 (Default) to use all but one. 
#' 
#' @return A tibble.

#' @export
#'
#' @import doParallel
#' @import seewave
#' @import foreach 
#' @importFrom tuneR readWave
#' @importFrom dplyr bind_cols tibble
#' @importFrom parallel detectCores makeCluster
#'
#' @examples 
#' \dontrun{
#' bbai_folder(path/to/folder)
#' }  

bbai_folder <- function(folder = NULL,
                        recursive = FALSE,
                        list = NULL,
                        start = 0,
                        end = 1,
                        unit = "minutes",
                        channel = 'each',
                        hpf = 0,
                        freq.res = 50,
                        cutoff = -60,
                        click.length = 10,
                        difference = 10,
                        gap.allowance = 2,
                        spectrogram = FALSE,
                        dark.plot = FALSE,
                        plot.title = "",
                        verbose = FALSE,
                        output.csv = "bbai_results.csv",
                        n.cores = -1) { 
  
  cat("Working on it...\n")
  
  args_list <- list(channel = channel,
                    hpf = hpf,
                    freq.res = freq.res,
                    cutoff = cutoff,
                    click.length = click.length,
                    difference = difference,
                    gap.allowance = gap.allowance,
                    spectrogram = spectrogram,
                    dark.plot = dark.plot,
                    plot.title = plot.title,
                    verbose = verbose)
  
  if(is.null(folder)){
    folder <- getwd()
  }
  setwd(folder)
  
  if(is.null(list)){
    audio.list <- list_waves(recursive = recursive)
  } else {
    audio.list <- list
  }
  
  nFiles <- length(audio.list)
  
  if (nFiles == 0) {
    stop("No valid audio files found after filtering.")
  }
  
  if(is.null(n.cores)){
    num_cores <- 1
  }else if(n.cores == -1){
    num_cores <- parallel::detectCores() - 1
  }else{
    num_cores <- n.cores
  }
  if(nFiles < num_cores){
    num_cores <- nFiles
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  if (length(audio.list) > 10) {
    cat("Evaluating the job...\n")
    
    startTime <- Sys.time()
    
    sound1 <- readWave(audio.list[1])
    type <- ifelse(sound1@stereo, "stereo", "mono")
    
    bbai1 <- quiet(do.call(bbai, c(list(sound1), args_list)))
    tibble::tibble(file_name = "filename") |> bind_cols(bbai1)
    
    timePerFile <-  Sys.time() - startTime
    timePerFile <- timePerFile + 3
    
    rm(sound1)
    rm(bbai1)
    # Estimate total time accounting for parallel processing
    estimatedTotalTime <- (timePerFile * nFiles) / as.numeric(num_cores)
    # Add overhead time
    adjustedTotalTime <- estimatedTotalTime
    # Calculate the end time
    expectedCompletionTime <- Sys.time() + adjustedTotalTime
    # Set up parallel processing
    
    cat("Start time:", format(Sys.time(), "%H:%M"), "\n")
    cat("Expected time of completion:", format(expectedCompletionTime, "%H:%M"),"\n\n")
    
  } else {
    sound1 <- readWave(audio.list[1], from = 0, to = 2 , units ='seconds')
    type <- ifelse(sound1@stereo, "stereo", "mono")
    rm(sound1)
  }
  cat("Analyzing", nFiles, type, "files using", num_cores, "cores... \n")
  
  # Define parallel computation
  results <- foreach(file = audio.list,
                     .packages = c("tuneR", "seewave", "tibble")) %dopar% {
                       filename <- basename(file)
                       sound <- readWave(file,
                                         from = start,
                                         to = end,
                                         units = unit
                                         )
                       result_list <- list()
                       
                       results <- quiet(do.call(bbai, c(list(sound), args_list)))
                       
                       result_list <- list(tibble(file_name = filename, results))
                       
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







