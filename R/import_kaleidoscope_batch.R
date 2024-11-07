#' Batch Import and Process Kaleidoscope Files in Parallel
#'
#' This function imports and processes multiple Excel files (.xlsx or .xls) from the specified directory or file list using the `import_kaleidoscope` function. The function operates in parallel to speed up the process by utilizing multiple cores. The processed data from each file is row-bound into a single data frame, and a new column (specified by the user) is added with a constant value.
#'
#' @param file_list A list of Excel files to be imported. If NULL, the function will search for Excel files in the specified `path` (default is the working directory).
#' @param new.col.name The name of the new column to be added to the output data frame.
#' @param new.col.value The value that will be assigned to the new column for all rows.
#' @param path The path to the directory containing the Excel files. Default is the working directory.
#' @param cores The number of CPU cores to use for parallel processing. If NULL, the function will use all but one available core.
#'
#' @return A data frame combining the imported and processed data from all Excel files, with the new column added.
#' @export
#'
#' @examples
#' # Assuming there are Excel files in the working directory:
#' result <- import_kaleidoscope_batch(new.col.name = "treatment", new.col.value = "test")
#' 
#' # Specifying a list of files and a custom path:
#' file_list <- c("file1.xlsx", "file2.xlsx")
#' result <- import_kaleidoscope_batch(file_list = file_list, new.col.name = "batch", new.col.value = "batch_1", path = "path/to/files")

import_kaleidoscope_batch <- function(file_list = NULL, new.col.name, new.col.value, path = getwd(), cores = NULL) {
  
  # Load required packages
  library(parallel)
  library(doParallel)
  library(foreach)
  
  # If no file list is provided, find all Excel files in the working directory or specified path
  if (is.null(file_list)) {
    file_list <- list.files(path = path, pattern = "\\.(xlsx|xls)$", full.names = TRUE)
  }
  
  # Check if there are any Excel files
  if (length(file_list) == 0) {
    stop("No Excel files found in the specified directory.")
  }
  
  # Determine the number of cores to use
  if (is.null(cores)) {
    cores <- detectCores() - 1  # Use all but one core
  }
  
  # Set up parallel backend
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Parallel processing using foreach
  new.df <- foreach(file = file_list, .combine = dplyr::bind_rows, .packages = c("dplyr", "tidyr", "readxl", "lubridate", "tools", "SoundEcology2")) %dopar% {
    
    # Import the data using the import_kaleidoscope function
    df <- import_kaleidoscope(file)
    
    # Add the new column with the specified name and value
    df[[new.col.name]] <- new.col.value
    
    return(df)
  }
  
  # Stop the cluster after processing
  stopCluster(cl)
  
  return(new.df)
}
