#' List wave files in a directory
#'
#' @param folder Character. The path to the folder where the WAV files are stored.
#'
#' @return A list of wave files
#' @export
#'
#' @examples list_waves(folder)
list_waves <- function(folder){

  setwd(folder)

  list <- list.files(pattern = "*.wav$")

  return(list)

}
