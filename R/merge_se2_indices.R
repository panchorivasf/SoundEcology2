#' Merge CSV outputs of several indices calculated with SoundEcology2
#'
#' @param folder_path character. The path to the folder containing the csv files with the "results.csv" suffix
#' @param recursive logical. If TRUE, it looks for csv files in the subfolders.
#'
#' @return A tibble with the merged tables for different indices.
#' @export
#'
#' @examples
#' \dontrun{
#' indices_all <- merge_se2_indices("path_to_folder", recursive = TRUE)
#' }
merge_se2_indices <- function(folder_path = NULL, recursive = TRUE) {

  # Set the folder path to current working directory if none is provided
  folder <- ifelse(is.null(folder_path), getwd(), folder_path)

  # Get list of all files with pattern "*results.csv" (recursively if specified)
  files <- list.files(path = folder,
                      pattern = ".*results.csv$",
                      full.names = TRUE,
                      recursive = recursive)
  cat("Detected files: \n")
  print(files)

  # Read, filter, and bind rows
  merged_data <- files |>
    purrr::map_dfr(~ {
      # Read the CSV file
      data <- utils::read.csv(.x)

      # Check if required columns ("file_name", "value_l", "value_r") are present
      if (all(c("file_name", "value_l", "value_r", "value_avg") %in% names(data))) {

        # Select columns from "file_name" to "value_r" and filter out any "channel" column
        selected_data <- data |>
          dplyr::select(file_name:value_avg) |>
          dplyr::select(-matches("^channel$"))  # Remove original "channel" column if present

        return(selected_data)

      } else {
        # If required columns are missing, skip this file and return NULL
        return(NULL)
      }
    }, .id = NULL)  # This helps to safely skip files

  # Pivot data longer for "value_l" and "value_r"
  # pivoted_data <- merged_data |>
  #   tidyr::pivot_longer(cols = c("value_l", "value_r", "value_avg"),
  #                names_to = "channel",     # Create a "channel" column 
  #                values_to = "value") #|>
    # rename channels
    # dplyr::mutate(channel = case_when(
    #   channel == "value_l" ~ "left",
    #   channel == "value_r" ~ "right",
    #   channel == "value_avg" ~ "mean"
    # ))

  return(merged_data)
}

