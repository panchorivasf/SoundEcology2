#' Generate and Save Summary Spectrograms
#'
#' This function generates and saves spectrogram images for various summary categories (mean, max, min, median, mode) for each recording in the input data frame. Spectrograms are created for specific wave files as indicated in the data, with options to normalize and remove DC offset.
#'
#' @param df A data frame containing summary statistics of recordings. It should include columns `sensor_id`, `channel`, `index`, and columns named `closest_to_<stat>` for categories `mean`, `max`, `min`, `median`, and `mode`, each referring to the closest wave file for each statistic.
#' @param parent_dir A character string specifying the parent directory containing the `.wav` files. Defaults to `NULL`, which sets `parent_dir` to the current working directory.
#' @param output_dir A character string specifying the directory to save generated spectrogram images. Defaults to `"summary_spectrograms_se2"`.
#' @param rmdcoff A logical value indicating whether to remove the DC offset from wave files before generating spectrograms. Defaults to `TRUE`.
#' @param norm A logical value indicating whether to normalize the spectrogram. Defaults to `FALSE`.
#'
#' @return None. Spectrogram images are saved in the specified output directory.
#' @export
#'
#' @examples
#' # Example usage with a data frame `df` that includes necessary summary columns
#' summary_spectrograms(df, parent_dir = "path/to/wav/files", output_dir = "spectrograms", rmdcoff = TRUE, norm = FALSE)
summary_spectrograms <- function(df,
                                 parent_dir = NULL,
                                 output_dir = "summary_spectrograms_se2",
                                 rmdcoff = TRUE,
                                 norm = FALSE) {

  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (is.null(parent_dir)) {
    par_dir <- getwd()
  } else {
    par_dir <- parent_dir
  }

  # Get all wave files in subdirectories
  wave_files <- list.files(par_dir, pattern = "\\.wav$", recursive = TRUE, full.names = TRUE)

  # Helper function to process each row of the data
  process_row <- function(row) {
    sensor_id <- row$sensor_id
    channel <- row$channel
    index <- row$index

    # Create separate output folder for each sensor_id and channel combination
    unit_channel_dir <- file.path(output_dir, paste0(sensor_id, "_", channel))
    if (!dir.exists(unit_channel_dir)) {
      dir.create(unit_channel_dir, recursive = TRUE)
    }

    # Define the list of file categories to process (mean, max, min, median, mode)
    categories <- c("mean", "max", "min", "median", "mode")

    # Loop through each category and create a spectrogram
    for (cat in categories) {
      # Construct the file name from the dataframe
      file_name <- row[[paste0("closest_to_", cat)]]
      if (!is.na(file_name)) {
        # Find the full path to the wave file (search in subdirectories)
        file_path <- wave_files[grep(file_name, wave_files)]

        if (length(file_path) == 1 && file.exists(file_path)) {
          # Generate the output file name and path
          output_file <- file.path(unit_channel_dir, paste0(sensor_id, "_", channel, "_", index, "_", cat, ".png"))

          # Read the wave file
          wave <- tuneR::readWave(file_path)

          if (rmdcoff) {
            wave <- seewave::rmoffset(wave, output = "Wave")
          }

          # Generate the spectrogram
          png(output_file, width = 2400, height = 1200, res = 300) # High resolution, landscape format

          seewave::spectro(wave,
                           norm = norm,
                           scale = FALSE)
          # collevels = seq(-60, 0, by = 1)) # Customize frequency limits as needed

          grDevices::dev.off()
        } else {
          message(paste("File not found or multiple matches:", file_name))
        }
      }
    }
  }

  # Use future.apply for parallel processing to process each row of the dataframe
  future::plan(future::multisession) # Set up parallel backend (works on Windows)
  future.apply::future_lapply(seq_len(nrow(df)), function(i) process_row(df[i, ]))
}
