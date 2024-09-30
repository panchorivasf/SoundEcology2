#' Import and Merge CSV Files Matching a Pattern
#'
#' This function reads all CSV files in the specified folder (or the current working directory, if no folder is provided) that match a given pattern. It selects a predefined set of columns from each CSV and merges the files into a single data frame.
#'
#' @param folder A string specifying the folder path where the CSV files are located. If `NULL`, the function uses the current working directory.
#' @param pattern A string representing the file name pattern to search for (default is "_results.csv"). This should be a regular expression that matches the files to be imported.
#'
#' @return A data frame containing the merged data from all CSV files, with the following columns: `file_name`, `sensor_id`, `datetime`, `date`, `hour`, `index`, `value_l`, `value_r`, and `value_avg`.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Import all CSV files matching the pattern "_results.csv" from a folder
#' folder_path <- "path/to/folder"
#' combined_data <- import_indices(folder = folder_path, pattern = "_results.csv")
#' 
#' # Example: Import from the current working directory with default pattern
#' combined_data <- import_indices()
#' }
import_indices <- function(folder = NULL, pattern = "_results.csv") {
  
  if(is.null(folder)){
    folder <- getwd()
  }
  # List all files in the folder matching the pattern
  files <- list.files(path = folder, pattern = pattern, full.names = TRUE)
  
  # Check if any files are found
  if (length(files) == 0) {
    stop("No files matching the pattern were found.")
  }
  
  # Define the columns to be selected
  selected_columns <- c("sensor_id", "datetime", "date", "hour",
                        "index", "value_l", "value_r", "value_avg")
  
  # Read each CSV, select specific columns, and store them in a list
  df_list <- lapply(files, function(file) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Ensure only the selected columns are kept and return the trimmed data frame
    df <- df[, selected_columns, drop = FALSE]
    
    return(df)
  })
  
  # Combine all data frames into one using rbind
  merged_df <- do.call(rbind, df_list)
  
  return(merged_df)
}
