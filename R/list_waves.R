#' List wave files in a directory
#'
#' @param folder Character. The path to the folder where the WAV files are stored.
#'
#' @return A list of wave files
#' @export
#'
#' @examples 
#' \dontrun{
#' list_waves(folder)
#' }
list_waves <- function(folder = NULL, recursive = FALSE){
  
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
  

  list <- list.files(pattern = "*.wav$", recursive = recursive)
  
  print_list(list)

  invisible(list)

}
