#' Summary of Acoustic Indices
#'
#' This function processes and summarizes acoustic indices from a data frame. It formats data based on its origin, 
#' either `"kaleidoscope"` or `"soundecology2"`, and computes summary statistics such as mean, median, mode, 
#' minimum, and maximum values for each acoustic index. Additionally, it identifies the recording files that 
#' are closest to each of these summary statistics.
#'
  #' @param df A data frame containing acoustic indices. Expected column names depend on the `origin` parameter. 
  #' For `"kaleidoscope"` data, expected columns include `IN.FILE`, `CHANNEL`, and several acoustic indices 
  #' (`ndsi`, `cent`, etc.). For `"soundecology2"` data, columns should include `file_name`, `sensor_id`, 
  #' and `datetime`.
#' @param origin A character string specifying the source of the data. Options are `"kaleidoscope"` or 
#' `"soundecology2"`. Default is `"kaleidoscope"`.
#'
#' @return A data frame containing summary statistics (`mean`, `median`, `mode`, `max`, `min`) for each 
#' acoustic index, with the associated recording file closest to each statistic.
#' @export
#'
#' @examples
#' # Example with kaleidoscope data
#' df <- data.frame(in.file = c("unit1_20220101_123000.wav", "unit2_20220101_124000.wav"),
#'                  channel = c(0, 1), ndsi = c(0.3, 0.7), cent = c(200, 300))
#' indices_summary(df, origin = "kaleidoscope")
#'
#' # Example with soundecology2 data
#' df2 <- data.frame(file_name = c("unit1_2022-01-01 12:30:00", "unit2_2022-01-01 12:40:00"),
#'                   sensor_id = c("unit1", "unit2"), datetime = c("2022-01-01 12:30:00", "2022-01-01 12:40:00"))
#' indices_summary(df2, origin = "soundecology2")
indices_summary <- function(df, origin = "soundecology2") {
  
  if (origin == "kaleidoscope") {
    # Reformat data frame
    colnames(df) <- tolower(colnames(df))

    df <- dplyr::select(df, -offset)

    df <- df |>
      dplyr::mutate(channel = dplyr::case_when(
        channel == 0 ~ 'left',
        channel == 1 ~ 'right'
      )) |>
      tidyr::pivot_longer(cols = c(ndsi:cent), names_to = "index", values_to = "value")

    df$datetime <- sapply(strsplit(df$in.file, "[_]"), function(x) paste(x[2], gsub(".wav", "", x[3]), sep = " "))
    df$datetime <- lubridate::as_datetime(df$datetime, format = "%Y%m%d %H%M%S")
    df$date <- lubridate::date(df$datetime)
    df$hour <- lubridate::hour(df$datetime)
    df$sensor_id <- sapply(strsplit(df$in.file, "[_]"), function(x) { x[1] })

    df <- dplyr::select(df, c(in.file, sensor_id, datetime, date, hour, dplyr::everything()))
    
    df <- df |>
      dplyr::rename(file_name = in.file)

  } else if (origin == "soundecology2") {
    # Adjust to the structure of `indices_soundec2`
    df <- df |>
      dplyr::mutate(datetime = lubridate::as_datetime(datetime, format = "%Y-%m-%d %H:%M:%S"),
                    date = lubridate::date(datetime),
                    hour = lubridate::hour(datetime))
  } else {
    stop("Origin should be either 'kaleidoscope' or 'soundecology2'")
  }

  # Get summary statistics and find recordings closest to each statistic
  summary <- df |>
    dplyr::group_by(sensor_id, channel, index) |>
    dplyr::reframe(
      mean_value = mean(value, na.rm = TRUE),
      closest_to_mean = file_name[which.min(abs(value - mean(value, na.rm = TRUE)))],

      max_value = max(value, na.rm = TRUE),
      closest_to_max = file_name[which.min(abs(value - max(value, na.rm = TRUE)))],

      min_value = min(value, na.rm = TRUE),
      closest_to_min = file_name[which.min(abs(value - min(value, na.rm = TRUE)))],

      median_value = median(value, na.rm = TRUE),
      closest_to_median = file_name[which.min(abs(value - median(value, na.rm = TRUE)))],

      mode_value = pracma::Mode(round(value)),
      closest_to_mode = file_name[which.min(abs(value - pracma::Mode(round(value))))]
    )
  
  return(summary)
}
