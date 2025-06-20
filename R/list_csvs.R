#' List the CSV files in a directory
#'
#' @param folder Directory to search. If NULL (default), the current working directory will be used.
#' @param recursive Logical. Whether to search in subfolders. Default is TRUE.
#' @return A list with the csv files in the chosen directory. 
#' @export
#'
#' @examples 
#' \dontrun{
#' list_csvs()
#' }
list_csvs <- function(folder = NULL, recursive = TRUE){
  
  print_list <- function(list){
    for (i in seq_along(list)) {
      cat(i, list[[i]], "\n")
    }
  }
  
  if (is.null(folder)) {
    setwd(getwd())
  } else {
    setwd(folder)
  }
  
  list <- list.files(pattern = "*.csv$", recursive = recursive)
  
  print_list(list)
  
  invisible(list)
  
}
