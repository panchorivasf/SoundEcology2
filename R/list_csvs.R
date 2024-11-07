#' List the CSV files in a directory
#'
#' @param path Directory to search. If NULL (default), the current working directory will be used.
#'
#' @return A list with the csv files in the chosen directory. 
#' @export
#'
#' @examples 
#' \dontrun{
#' list_csvs()
#' }
list_csvs <- function(path = NULL){
  
  if (is.null(path)) {
    setwd(getwd())
  } else {
    setwd(path)
  }
  
  list <- list.files(pattern = "*.csv$")
  
  return(list)
  
}