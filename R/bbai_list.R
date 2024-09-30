#' Calculate BBAI for all the Files in a List
#'
#' This function processes an audio signal to detect broadband activity by identifying 'clicks' based on time-frame-wise (i.e., column-wise) amplitude changes in the spectrogram. It computes statistics related to click height, variance, and centroid frequency, and can plot a spectrogram with detected clicks highlighted. The function also classifies whether the signal contains noise or insect based on the variance and centroid frequencies of the clicks.
#'
#' @param audio.list a list of audio files to analyze
#' @param channel Character. If Wave is stereo and you want to use only one channel, pass either "left" or "right" to this argument. If you want to analyze a mix of both channels, select "mix". If NULL (default), results are returned for each channel.
#' @param hpf Numeric. High-pass filter. The default (500 Hz) should be used always for consistency unless signals of interest are below that threshold.
#' @param rm.offset Logical. Should the DC offset be removed from the audio signal? Defaults to `TRUE`.
#' @param freq.res Numeric. Frequency resolution in Hz. This value determines the "height" of each frequency bin and, therefore, the window length to be used (sampling rate / frequency resolution).
#' @param cutoff Numeric. The amplitude threshold (in dBFS) for removing low-amplitude values in the spectrogram. Default is `-50`.
#' @param click.height Numeric. The minimum height (in frequency bins) for a detected click to be kept. Default is `10`.
#' @param difference Numeric. The maximum difference in amplitude between adjacent frequency bins to be considered part of a single 'click'. Default is `20`.
#' @param gap.allowance Numeric. The size of gaps (in frequency bins) allowed between contiguous parts of a click. Default is `2`. Gaps larger than this value will split clicks.
#' @param plot Logical. Should a spectrogram with highlighted clicks be plotted? Default is `TRUE`.
#' @param dark.plot Logical. Should the plot use a dark theme (black background)? Default is `FALSE`.
#' @param plot.title Character. The title for the plot, if `plot` is `TRUE`. Default is `NULL`.
#' @param n.cores The number of cores to use for parallel processing. Use `n.cores = -1` to use all but one core. Default is NULL (single-core processing).
#' @param verbose Logical. If TRUE, details of dynamic range will be printed on the console.

#' @return A tibble containing the following columns:
#'   - `index`: The name of the index. Useful later when merging data with other indices.
#'   - `value`: The number of clicks detected in the recording.
#'   - `mean`: The mean click height (in frequency bins).
#'   - `variance`: The variance of the click height.
#'   - `sd`: The standard deviation of the click height.

#' @export
#'
#' @import tuneR
#' @import seewave
#' @import patchwork
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom tibble tibble
#'
#'
#' @examples
#' files.list <- list_waves(path/to/folder)
#' bbai_results <- bbai_list(files.list)
bbai_list <- function(audio.list,
                        channel = 'each',
                        hpf = 0,
                        rm.offset = TRUE,
                        freq.res = 50,
                        cutoff = -60,
                        click.length = 10,
                        difference = 10,
                        gap.allowance = 2,
                        verbose = FALSE,
                        output.csv = "bbai_results.csv",
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

  # setwd(folder)
  # files <- list.files(path=folder, pattern = ".wav|.WAV")

  filename <- tibble(filename = audio.list)
  nFiles <- length(audio.list)


  # Evaluate the duration of the analysis
  # Measure processing time for a single file
  startTime <- Sys.time()

  sound1 <- readWave(audio.list[1])
  type <- ifelse(sound1@stereo, "stereo", "mono")

  bbai1 <- quiet(bbai(sound1, channel = channel))
  tibble::tibble(file_name = "filename") %>% bind_cols(bbai1)

  # Assess how long it takes to parse 1 file
  timePerFile <-  Sys.time() - startTime
  # Add overhead per file
  timePerFile <- timePerFile + as.numeric(seconds(2.2))

  rm(sound1)
  rm(bbai1)


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
  results <- foreach(file = audio.list, .packages = c("tuneR", "seewave", "tibble")) %dopar% {

    filename <- basename(file)  # Get file name without path

    # Read audio file
    audio <- readWave(file)

    # Initialize an empty tibble for the results
    result_list <- list()

    # # Handle different channel selections
    # if (channel == "each") {
    #   # Process each channel separately
    #   bbai_left <- bbai(audio, plot.title = filename, channel = "left")$summary
    #   bbai_right <- bbai(audio, plot.title = filename, channel = "right")$summary
    #
    #   bbai <- bbai_left %>%
    #     dplyr::select(-c(value))
    #
    #   # Combine the results into one row with 'value_l' and 'value_r' columns
    #   result_list <- list(tibble(
    #     file_name = filename,
    #     value_l = bbai_left$value,
    #     value_r = bbai_right$value,
    #     value_avg = round((bbai_left$value + bbai_right$value) / 2, 3),
    #     nbai
    #   ))


    # } else {
    # Process selected channel ("left", "right", or "mix")
    bbai <- bbai(audio,
                 channel = channel,
                 hpf = hpf,
                 rm.offset = rm.offset,
                 freq.res = freq.res,
                 cutoff = cutoff,
                 click.length = click.length,
                 difference = difference,
                 gap.allowance = gap.allowance,
                 verbose = FALSE
    )

    result_list <- list(
      tibble(file_name = filename, channel = channel, bbai)
    )
    # }

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



