#' List the CSV files in a directory
#'
#' @param folder Directory to search. If NULL (default), the current working directory will be used.
#'
#' @return A list with the csv files in the chosen directory. 
#' @export
#'
#' @examples 
#' \dontrun{
#' list_csvs()
#' }
list_csvs <- function(folder = NULL){
  
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
  
  print_list(list)
  
  invisible(list)
  
}