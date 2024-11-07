#' Import tidy dataset from Kaleidoscope
#'
#' @param data An excel
#'
#' @return a data frame with the following columns 'file_name', 'sensor_id', 'datetime', 'date', 'hour', 'value_l', 'value_r', 'value_avg'
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples import_kaleidoscope(data)
import_kaleidoscope <- function(file) {
  
  file_extension <- tools::file_ext(file)
  
  if (file_extension %in% "csv"){
    
    df <- read.csv(file)
  } else if (file_extension %in% c("xlsx", "xls")){
    
    df <- read_excel(file)
  } else {
    stop("The file extension must be .csv, .xls, or .xlsx.\n")
  }
  
  # Pivot longer on NDSI:CENT columns
  df <- df %>%
    pivot_longer(cols = NDSI:CENT, names_to = "INDEX", values_to = "VALUE") %>%
    select(`IN FILE`, CHANNEL, INDEX, VALUE) %>%
    group_by(`IN FILE`, CHANNEL, INDEX) %>%
    # Ensure unique values for each combination of file_name, CHANNEL, and INDEX
    summarise(VALUE = mean(VALUE, na.rm = TRUE), .groups = "drop") %>%
    # Spread the values into VALUE_L and VALUE_R based on the CHANNEL column
    pivot_wider(names_from = CHANNEL, values_from = VALUE, names_prefix = "VALUE_") %>%
    # Rename columns for consistency
    rename(FILENAME = `IN FILE`, LEFT_CHANNEL = VALUE_0, RIGHT_CHANNEL = VALUE_1)
  
  df <- df %>%
    rename(file_name = FILENAME,
           index = INDEX,
           value_l = LEFT_CHANNEL,
           value_r = RIGHT_CHANNEL,
    )
  
  df <- df %>%
    mutate(value_avg = round(((df$value_l + df$value_r)/2),3))
  
  df <- df %>%
    select(file_name, index, value_l, value_r, value_avg)
  
  
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
    select(sensor_id, datetime, date, hour, index, value_l, value_r, value_avg)
  
  return(df)
}
