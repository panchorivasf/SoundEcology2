#' Harmonize "old" soundecology2 output with new SoundEcology2 table format
#'
#' @param data.frame A data frame containing the results from soundecology2
#' @param original.index Character. A name for the 'old' index, e.g., "adi_original"
#'
#' @return A new data frame with the harmonized columns.
#' @export
#' @importFrom dplyr mutate
#' 
#' @examples harmonize_index(ndsi_original, "ndsi_original")
harmonize_index <- function(data.frame, original.index = NULL){
  data.frame <- data.frame |>
    rename(file_name = FILENAME,
           index = INDEX,
           value_l = LEFT_CHANNEL,
           value_r = RIGHT_CHANNEL,
    ) |>
    mutate(value_avg = round(((data.frame$LEFT_CHANNEL + data.frame$RIGHT_CHANNEL)/2),3))
  
  data.frame <- data.frame |>
    select(file_name, index, value_l, value_r, value_avg)
  
  data.frame <- addMetadata(data.frame)
  
  data.frame$index <- original.index
  
  write.csv(data.frame, file = paste0(original.index, "_results.csv"))
  
  return(data.frame)
}
