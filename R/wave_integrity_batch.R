#' Compute file integrity diagnostics for several folders containing WAV files
#'
#' @param parent_folder Character. The path to the parent folder containing folders with WAV files.
#' @param n.cores Numeric. Number of cores to be used in the parallel processing. Defaults to -1, which means all but one core.
#'
#' @return A tibble with summary data on file size, including the last day with complete recordings and days with corrupted files (i.e., smaller than most).
#' @export
#'
#' @examples wave_integrity_batch(parentFolder)
wave_integrity_batch <- function(parent_folder, n.cores = -1) {
  # List all subfolders in the parent folder
  subfolders <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)
  
  # Check if there are subfolders
  if (length(subfolders) == 0) {
    stop("No subfolders found in the parent folder.")
  }
  
  # Loop through each subfolder and process with wave_integrity
  for (subfolder in subfolders) {
    cat("\nProcessing subfolder:", subfolder, "\n")
    
    tryCatch({
      # Call the wave_integrity function for each subfolder
      wave_integrity(
        folder = subfolder, 
        type = "lines", 
        ggplot = TRUE, 
        plot.file = paste0("wave_integrity_plot_", basename(subfolder), ".png"), 
        log.file = paste0("wave_integrity_log_", basename(subfolder), ".txt"), 
        n.cores = n.cores
      )
    }, error = function(e) {
      cat("Error processing subfolder:", subfolder, "\n", conditionMessage(e), "\n")
    })
  }
  
  cat("All subfolders processed successfully.\n")
}
