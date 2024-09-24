#' Add Metadata to a data frame or tibble containing a file_name column for Wildlife Acoustics' SM recorders.
#'
#' @param df a data frame obtained with one of the acoustic index functions for lists of folders.
#'
#' @return a new data frame with added metadata including date-time and its division and the unit id
#'
#' @details
#' The input data frame must contain a column called "file_name"
#'
#' @import lubridate
#' @import tidyverse
#'
#' @noRd
#' @examples addMetadata(df)
addMetadata <- function(df){

  df$datetime <- sapply(strsplit(df$file_name, "[_]"), function(x) paste(x[2], gsub(".wav", "", x[3]), sep = " "))

  # # delete the '.wav' extension
  # df$datetime<-gsub(".wav","",as.character(df$datetime))
  df$datetime <- gsub("\\.wav$", "", df$datetime)

  # transform to date time format
  df$datetime <- lubridate::as_datetime(df$datetime, , format = "%Y%m%d %H%M%S")
  # add a date column
  df$date <- lubridate::date(df$datetime)
  # add an hour column
  df$hour <- lubridate::hour(df$datetime)

  # extract the sensor id
  df$sensor_id <- sapply(strsplit(df$file_name,"[_]"), function(x){x[1]})

  df <- df %>%
    select(c(file_name, sensor_id, datetime, date, hour, everything()))

  return(df)

}
