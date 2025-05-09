#' Merge Index Results from CSV Files
#'
#' This function takes a folder path, reads all CSV files in the folder that match the pattern "*results.csv", 
#' and merges them into a single data frame. It keeps only the relevant columns: `file_name`, `sensor_id`, `datetime`, 
#' `date`, `hour`, `index`, and `value_avg`.
#'
#' @param folder_path A character string specifying the folder containing the CSV files to be merged.
#'
#' @return A data frame with the merged contents of the CSV files, containing only the specified columns.
#' 
#' @export
#' 
#' @importFrom dplyr select bind_rows
#'
#' @examples
#' \dontrun{
#'   # Merge CSV files in the specified folder
#'   merged_data <- merge_index_results("/path/to/csv/folder")
#' }
merge_index_results <- function(folder_path) {
  # List all CSV files in the folder that match the pattern "*results.csv"
  file_list <- list.files(path = folder_path, pattern = "*results.csv", full.names = TRUE)
  
  # Import and merge the CSV files, keeping only the desired columns
  merged_data <- file_list  |> 
    lapply(read.csv)  |> 
    lapply(function(df) select(df, file_name, sensor_id, datetime, date, hour, index, value_avg)) %>%
    bind_rows()
  
  return(merged_data)
}

