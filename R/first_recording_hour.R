#' Subsample a dataset with only the first recording in each hour
#'
#' @param data a data frame or tibble
#'
#' @return a data frame with the subsetted observations (1 per hour).
#' @export
#' @importFrom dplyr group_by slice_min ungroup
#'
#' @examples first_recording_hour(dataframe)
first_recording_hour <- function(data) {
  data  |> 
    dplyr::group_by(sensor_id, date, hour)  |> 
    dplyr::slice_min(order_by = datetime, with_ties = FALSE)  |> 
    dplyr::ungroup()
}