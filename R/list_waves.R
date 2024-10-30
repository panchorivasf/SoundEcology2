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
list_waves <- function(folder = NULL){
  
  if(is.null(folder)){
    folder <- getwd()
  }

  setwd(folder)

  list <- list.files(pattern = "*.wav$")

  return(list)

}
