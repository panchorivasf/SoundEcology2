#' Check the Integrity of a Wave File Dataset
#'
#' @param folder Character. The path to the folder containing the wave files.
#' @param type Character. Type of plot to use in the geom call when ggplot is TRUE. Options are "lines", "bars", and "points".
#' @param ggplot Logical. If TRUE, the plot is produced with ggplot2. If FALSE, base plot is used.
#' @param plot.file Character. A name for the output plot. Include the extension in the string.
#' @param log.file Character. A name for the log file (text file) output.
#' @param n.cores Numeric. Number of cores to use in parallel processing. Defults to NULL (no parallel processing).
#' @param dump.folder Logical. If TRUE, a 'dump' folder is created in the working directory and all the DUMP files are moved to that folder.
#' @param tail.folder Logical. If TRUE, a 'tail' folder is created in the working directory and all files belonging to dates after the last date with complete recordings are moved to that folder.
#'
#' @return A tibble with summary data on file size, including the last day with complete recordings and days with corrupted files (i.e., smaller than most).
#' @description This function parses all the WAVE files in a folder and checks their size, reporting corrupted files.
#' @noRd
#'
#' @importFrom dplyr mutate group_by summarize filter pull
#' @importFrom ggplot2 ggplot aes geom_line geom_bar geom_point scale_x_date labs theme element_text ggsave
#' @importFrom writexl write_xlsx
#'
#' @examples integrity <- wave_integrity("pathToFolder", ggplot = FALSE)

wave_integrity <- function(folder,
                           type = "lines",
                           ggplot = TRUE,
                           plot.file = "wave_integrity_plot.png",
                           log.file = "wave_integrity_log.txt",
                           n.cores = -1,
                           dump.folder = TRUE,
                           tail.folder = TRUE) {

  cat("Analyzing the folder... please wait...\n")

  # List .wav files and extract size and modification date
  wav_files <- list.files(folder, pattern = "\\.wav$", full.names = TRUE)

  if (length(wav_files) == 0) {
    stop("No .wav files found in the specified folder.")
  }

  # Create "dump" folder only if there are "dump" files to move
  if (dump.folder) {
    # Detect "dump" files
    dump_files <- list.files(folder, pattern = "dump$", full.names = TRUE)

    if (length(dump_files) > 0) {
      # Create the "dump" folder in the working directory if it doesn't exist
      dump_folder <- file.path(getwd(), "dump")
      if (!dir.exists(dump_folder)) {
        dir.create(dump_folder)
        cat("Created 'dump' folder in the working directory.\n")
      }

      # Move "dump" files to the "dump" folder
      file.rename(dump_files, file.path(dump_folder, basename(dump_files)))
      cat(length(dump_files), "Dump files have been moved to the 'dump' folder.\n")
    } else {
      cat("No dump files were found.\n")
    }
  }

  # Set up parallel processing
  if (is.null(n.cores)) {
    num_cores <- 1
  } else if (n.cores == -1) {
    num_cores <- parallel::detectCores() - 1
  } else {
    num_cores <- n.cores
  }

  if (num_cores > 1) {
    cat("Using parallel processing with", num_cores, "cores\n")
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  } else {
    cat("Using single core (no parallel processing)\n")
  }

  # Parallel processing for gathering file info
  if (num_cores > 1) {
    wav_info <- foreach::foreach(file = wav_files, .packages = c("tuneR", "dplyr")) %dopar% {
      size <- file.info(file)$size / (1024 * 1024)  # size in MB
      date <- as.POSIXct(file.info(file)$mtime)
      list(file = basename(file), size = size, date = date)
    }
    wav_info <- do.call(rbind, lapply(wav_info, function(x) as.data.frame(x)))
  } else {
    wav_info <- data.frame(
      file = basename(wav_files),
      size = file.info(wav_files)$size / (1024*1024),  # size in MB
      date = as.POSIXct(file.info(wav_files)$mtime)
    )
  }

  # Summarize by day
  summary_data <- wav_info %>%
    mutate(date = as.Date(date)) %>%
    group_by(date) %>%
    summarize(mean_size = round(mean(size)), .groups = 'drop')

  # Calculate the rounded median size of all wav files
  median_size <- round(median(wav_info$size, na.rm = TRUE))

  # Identify corrupted dates where mean size is less than the median size
  corrupted_dates <- summary_data %>%
    filter(mean_size < median_size) %>%
    pull(date)

  # Write log file with summary information
  writeLines(c(
    paste("Median file size: ", round(median_size, 2), " MB"),
    "Dates with corrupted files:",
    paste(as.character(corrupted_dates), collapse = "\n")
  ), con = log.file)

  # Print diagnostics to the console
  cat("Median file size: ", round(median_size, 2), " MB\n")
  cat("Dates with corrupted files: \n")
  print(as.character(corrupted_dates))

  # Create "tail" folder and move small files if tail.folder is TRUE
  if (tail.folder) {
    tail_folder <- file.path(getwd(), "tail")
    if (!dir.exists(tail_folder)) {
      dir.create(tail_folder)
      cat("Created 'tail' folder in the working directory.\n")
    }

    # Move files smaller than the rounded median size and within corrupted dates
    small_files <- wav_files[wav_info$size < median_size & as.Date(wav_info$date) %in% corrupted_dates]

    if (length(small_files) > 0) {
      if (num_cores > 1) {
        foreach::foreach(file = small_files) %dopar% {
          file.rename(file, file.path(tail_folder, basename(file)))
        }
      } else {
        lapply(small_files, function(file) {
          file.rename(file, file.path(tail_folder, basename(file)))
        })
      }
      cat(length(small_files), "files from corrupted dates have been moved to the 'tail' folder.\n")
    } else {
      cat("No corrupted fileswere found.\n")
    }
  }

  # Plot the results
  if (ggplot) {
    p <- ggplot(summary_data, aes(x = date, y = mean_size)) +
      scale_x_date(date_breaks = "1 day", date_labels = "%Y-%m-%d") +
      labs(x = "Date", y = "Mean Size (MB)", title = "WAV files integrity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (type == "lines") {
      p <- p + geom_line()
    } else if (type == "bars") {
      p <- p + geom_bar(stat = "identity", width = 0.5)
    } else if (type == "points") {
      p <- p + geom_point()
    }

    # Highlight corrupted files in red
    p <- p + geom_point(data = summary_data[summary_data$mean_size < median_size, ],
                        aes(x = date, y = mean_size), color = 'red', size = 3)

    print(p)

    # Save the plot
    ggsave(plot.file, plot = p, width = 8, height = 6)
    cat("Plot saved to:", plot.file, "\n")
  } else {
    # Base R plot
    png(filename = plot.file, width = 800, height = 600)
    plot(summary_data$date, summary_data$mean_size, type = "n", xlab = "Date", ylab = "Mean Size (MB)",
         main = "WAV file size over time", xaxt = "n")

    if (type == "lines") {
      lines(summary_data$date, summary_data$mean_size)
    } else if (type == "bars") {
      barplot(summary_data$mean_size, names.arg = summary_data$date, xlab = "Date", ylab = "Mean Size (MB)")
    } else if (type == "points") {
      points(summary_data$date, summary_data$mean_size, pch = 19)
    }

    axis(1, at = summary_data$date, labels = format(summary_data$date, "%Y-%m-%d"), las = 2)

    # Highlight corrupted files in red
    points(summary_data$date[summary_data$mean_size < median_size],
           summary_data$mean_size[summary_data$mean_size < median_size], col = "red", pch = 19)

    dev.off()
    cat("Plot saved to:", plot.file, "\n")
  }

  # Stop the cluster if parallel processing was used
  if (num_cores > 1) {
    parallel::stopCluster(cl)
  }

  invisible(summary_data)
}






