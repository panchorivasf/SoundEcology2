
addMetadata <- function(df){

  require(lubridate)
  require(tidyverse)

  # Extract date and time and create new columns for them
  df$datetime <- sapply(strsplit(df$file_name,"[_]"),function(x){x[2]}) %>%
    paste(sapply(strsplit(df$file_name,"[_]"),function(x){x[3]}))
  # delete the '.wav' extension
  df$datetime<-gsub(".wav","",as.character(df$datetime))
  # transform to date time format
  df$datetime <- as_datetime(df$datetime)
  # add a date column
  df$date <- date(df$datetime)
  # add an hour column
  df$hour <- hour(df$datetime)

  # extract the sensor id
  df$unit_id <- sapply(strsplit(df$file_name,"[_]"), function(x){x[1]})

  df <- df %>%
    select(c(file_name, unit_id, datetime, date, hour, everything()))

  return(df)

}
