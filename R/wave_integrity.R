#' Check the Integrity of a Wave File Dataset
#'
#' @param folder Character. The path to the folder containing the wave files.
#' @param type Character. Type of plot to use in the geom call when ggplot is TRUE. Options are "lines", "bars", and "points".
#' @param ggplot Logical. If TRUE, the plot is produced with ggplot2. If FALSE, base plot is used.
#' @param output.file Character. A name for the output .xlsx file. Include the extension in the string.
#' @param plot.file Character. A name for the output plot. Include the extension in the string.
#'
#' @return A tibble with summary data on file size, including the last day with complete recordings and days with corrupted files (i.e., smaller than most).
#' @description This function parses all the WAVE files in a folder and checks their size, reporting corrupted files.
#' @export
#'
#' @importFrom dplyr mutate group_by summarize filter pull
#' @importFrom ggplot2 ggplot aes geom_line geom_bar geom_point scale_x_date labs theme element_text ggsave
#' @importFrom writexl write_xlsx
#'
#' @examples integrity <- wave_integrity("pathToFolder", ggplot = FALSE)

wave_integrity <- function(folder, type = "lines", ggplot = TRUE, output.file = "wave_integrity_summary.xlsx", plot.file = "wave_integrity_plot.png") {
  # List .wav files and extract size and modification date
  wav_files <- list.files(folder, pattern = "\\.wav$", full.names = TRUE)

  if (length(wav_files) == 0) {
    stop("No .wav files found in the specified folder.")
  }

  wav_info <- data.frame(
    file = basename(wav_files),
    size = file.info(wav_files)$size / 1024, # size in KB
    date = as.POSIXct(file.info(wav_files)$mtime)
  )

  # Summarize by day
  summary_data <- wav_info %>%
    mutate(date = as.Date(date)) %>%
    group_by(date) %>%
    summarize(mean_size = mean(size), .groups = 'drop')

  # Calculate the first day's mean size
  first_day_mean_size <- summary_data$mean_size[1]

  # Detect the last complete date where the mean size is equal to the first day's mean size
  last_complete_date <- summary_data %>%
    filter(mean_size == first_day_mean_size) %>%
    summarize(last_date = max(date)) %>%
    pull(last_date)

  # Detect dates with corrupted files (outliers)
  corrupted_dates <- summary_data %>%
    filter(mean_size < first_day_mean_size) %>%
    pull(date)

  # Print diagnostics
  cat("Last complete date: ", format(last_complete_date, "%Y-%m-%d"), "\n")

  cat("Dates with corrupted files: \n")
  print(as.character(corrupted_dates))

  # Save the summary_data to an Excel file
  write_xlsx(summary_data, output.file)
  cat("Summary data saved to:", output.file, "\n")

  # Plot based on ggplot or base R
  if (ggplot) {
    # GGplot2
    p <- ggplot(summary_data, aes(x = date, y = mean_size)) +
      scale_x_date(date_breaks = "1 day", date_labels = "%Y-%m-%d") +  # Ensure daily ticks and labels
      labs(x = "Date", y = "Mean Size (KB)", title = "WAV files integrity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

    if (type == "lines") {
      p <- p + geom_line()
    } else if (type == "bars") {
      p <- p + geom_bar(stat = "identity", width = 0.5)
    } else if (type == "points") {
      p <- p + geom_point()
    }

    # Add outliers as red points
    p <- p + geom_point(data = summary_data[summary_data$mean_size < first_day_mean_size, ],
                        aes(x = date, y = mean_size), color = 'red', size = 3)

    print(p)

    # Save the plot as a PNG
    ggsave(plot.file, plot = p, width = 8, height = 6)
    cat("Plot saved to:", plot.file, "\n")

  } else {
    # Base R plot
    png(filename = plot.file, width = 800, height = 600)  # Open PNG device

    plot(summary_data$date, summary_data$mean_size, type = "n", xlab = "Date", ylab = "Mean Size (KB)",
         main = "Size of WAV files over time", xaxt = "n")

    if (type == "lines") {
      lines(summary_data$date, summary_data$mean_size)
    } else if (type == "bars") {
      barplot(summary_data$mean_size, names.arg = summary_data$date, xlab = "Date", ylab = "Mean Size (KB)")
    } else if (type == "points") {
      points(summary_data$date, summary_data$mean_size, pch = 19)
    }

    # Customize x-axis to show daily ticks and labels
    axis(1, at = summary_data$date, labels = format(summary_data$date, "%Y-%m-%d"), las = 2)  # Rotate labels

    # Add outliers in red
    points(summary_data$date[summary_data$mean_size < first_day_mean_size],
           summary_data$mean_size[summary_data$mean_size < first_day_mean_size], col = "red", pch = 19)

    dev.off()  # Close PNG device
    cat("Plot saved to:", plot.file, "\n")
  }

  invisible(summary_data)
}

