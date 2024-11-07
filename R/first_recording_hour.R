#' Subsample a dataset with only the first recording in each hour
#'
#' @param data a data frame or tibble
#'
#' @return a data frame with the subsetted observations (1 per hour).
#' @export
#' @import tidyverse
#'
#' @examples first_recording_hour(dataframe)
first_recording_hour <- function(data) {
  data %>%
    group_by(sensor_id, date, hour) %>%
    slice_min(order_by = datetime, with_ties = FALSE) %>%
    ungroup()
}