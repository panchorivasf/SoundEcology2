#' Create time series analysis and visualization for acoustic index data
#'
#' This function processes acoustic index data, summarizes it by day/week/month,
#' and creates interactive visualizations with optional smoothing. It handles datetime
#' parsing from filenames, data summarization at different temporal resolutions,
#' and generates plots with range ribbons and average values.
#'
#' @param data A data frame containing acoustic index data. Must include columns:
#'   'index', 'file_name', and 'value_avg'.
#' @param index_name The name of the index to analyze (default: "adi"). Used for plot labels.
#' @param extra_title Additional text to include in the plot title (default: "").
#' @param start Optional start date for filtering data (format: "YYYY-MM-DD").
#' @param end Optional end date for filtering data (format: "YYYY-MM-DD").
#' @param span Span parameter for LOESS smoothing (default: 0.3).
#' @param summarize_by Temporal resolution for the main plot: "day", "week", or "month" (default: "day").
#'
#' @return A list containing:
#' \itemize{
#'   \item daily_tibble: Data summarized by day
#'   \item weekly_tibble: Data summarized by week
#'   \item monthly_tibble: Data summarized by month
#'   \item plot: Interactive plotly visualization of the time series
#' }
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Extracts datetime information from filenames
#'   \item Parses datetime using multiple possible formats
#'   \item Summarizes data at daily, weekly, and monthly levels
#'   \item Filters data based on start/end dates if provided
#'   \item Creates an interactive plot with range ribbons and average line
#'   \item Optionally adds LOESS smoothing when enough data points exist
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- index_ts(acoustic_data)
#' 
#' # With custom parameters
#' result <- index_ts(acoustic_data, 
#'                   index_name = "ndsi",
#'                   extra_title = "Forest Site A",
#'                   start = "2023-01-01",
#'                   end = "2023-06-30",
#'                   summarize_by = "week")
#' }
#'
#' @importFrom dplyr mutate group_by summarise ungroup filter rename case_when
#' @importFrom lubridate parse_date_time ymd floor_date
#' @importFrom plotly plot_ly add_ribbons add_trace layout
#' @importFrom stats predict loess median na.exclude
#' @export

index_ts <- function(data, 
                     index_name = "adi",
                     extra_title = "",
                     start = NULL, 
                     end = NULL, 
                     span = 0.3, 
                     summarize_by = "day") {
  
  # Validate summarize_by argument
  if (!summarize_by %in% c("day", "week", "month")) {
    stop("summarize_by must be one of: 'day', 'week', or 'month'")
  }
  
  # Check if index exists in data
  if (!"index" %in% colnames(data)) {
    stop("Data must contain an 'index' column")
  }
  
  data$datetime <- sapply(strsplit(data$file_name, "[_]"), function(x) paste(x[2], gsub(".wav", "", x[3]), sep = " "))
  
  
  # Filter data for the specified index
  data_filtered <- data
  
  if (nrow(data_filtered) == 0) {
    stop(paste("No data found for index:", index_name))
  }
  
  # Improved datetime conversion with multiple format attempts
  if (is.character(data_filtered$datetime)) {
    data_filtered <- data_filtered |> 
      mutate(
        datetime = case_when(
          # Try lubridate's automatic parsing first
          !is.na(lubridate::parse_date_time(datetime, orders = c("ymd HMS", "mdy HMS", "dmy HMS"))) ~ 
            lubridate::parse_date_time(datetime, orders = c("ymd HMS", "mdy HMS", "dmy HMS")),
          # Fallback to specific common formats
          !is.na(as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) ~ 
            as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
          !is.na(as.POSIXct(datetime, format = "%m/%d/%Y %H:%M")) ~ 
            as.POSIXct(datetime, format = "%m/%d/%Y %H:%M"),
          !is.na(as.POSIXct(datetime, format = "%d-%m-%Y %H:%M")) ~ 
            as.POSIXct(datetime, format = "%d-%m-%Y %H:%M"),
          TRUE ~ as.POSIXct(NA)  # If all parsing fails
        )
      )
    
    # Check if any dates failed to parse
    if (any(is.na(data_filtered$datetime))) {
      warning("Some datetime values could not be parsed. Check your datetime format.")
    }
  }
  
  # Rest of your function remains the same...
  # Create daily tibble (base for all other summaries)
  daily_tibble <- data_filtered |>
    mutate(date = as.Date(datetime)) |>
    group_by(date) |>
    summarise(
      index_min = min(value_avg, na.rm = TRUE),
      index_max = max(value_avg, na.rm = TRUE),
      index_avg = mean(value_avg, na.rm = TRUE),
      index_median = median(value_avg, na.rm = TRUE),
      n_obs = n()
    ) |>
    ungroup()
  
  # Create weekly tibble
  weekly_tibble <- daily_tibble |>
    mutate(week_start = floor_date(date, "week")) |>
    group_by(week_start) |>
    summarise(
      index_min = min(index_avg, na.rm = TRUE),
      index_max = max(index_avg, na.rm = TRUE),
      index_avg = mean(index_avg, na.rm = TRUE),
      index_median = median(index_avg, na.rm = TRUE),
      n_obs = sum(n_obs, na.rm = TRUE)
    ) |>
    rename(date = week_start) |>
    ungroup()
  
  # Create monthly tibble
  monthly_tibble <- daily_tibble |>
    mutate(month_start = floor_date(date, "month")) |>
    group_by(month_start) |>
    summarise(
      index_min = min(index_avg, na.rm = TRUE),
      index_max = max(index_avg, na.rm = TRUE),
      index_avg = mean(index_avg, na.rm = TRUE),
      index_median = median(index_avg, na.rm = TRUE),
      n_obs = sum(n_obs, na.rm = TRUE)
    ) |>
    rename(date = month_start) |>
    ungroup()
  
  # Subset data based on start and end dates if provided
  if (!is.null(start)) {
    start_date <- ymd(start)
    daily_tibble <- daily_tibble |> filter(date >= start_date)
    weekly_tibble <- weekly_tibble |> filter(date >= start_date)
    monthly_tibble <- monthly_tibble |> filter(date >= start_date)
  }
  
  if (!is.null(end)) {
    end_date <- ymd(end)
    daily_tibble <- daily_tibble |> filter(date <= end_date)
    weekly_tibble <- weekly_tibble |> filter(date <= end_date)
    monthly_tibble <- monthly_tibble |> filter(date <= end_date)
  }
  
  # Select which tibble to use for plotting based on summarize_by
  plot_tibble <- switch(summarize_by,
                        "day" = daily_tibble,
                        "week" = weekly_tibble,
                        "month" = monthly_tibble)
  
  # Create the plot
  index_plot <- plot_ly(plot_tibble) |>
    add_ribbons(x = ~date, 
                ymin = ~index_min, 
                ymax = ~index_max, 
                name = paste(toupper(index_name), "Range"),
                fillcolor = 'rgba(100, 149, 237, 0.2)',
                line = list(color = 'rgba(100, 149, 237, 0.1)'),
                hoverinfo = 'text',
                text = ~paste("Date:", date, 
                              "<br>Min:", round(index_min, 3),
                              "<br>Max:", round(index_max, 3))) |>
    add_trace(x = ~date, y = ~index_avg, name = paste("Avg", toupper(index_name)), 
              type = 'scatter', mode = 'lines',
              line = list(color = 'rgba(70, 130, 180, 0.7)', width = 2),
              hoverinfo = 'text',
              text = ~paste("<br>Mean:", round(index_avg, 3))) |>
    layout(yaxis = list(title = paste(toupper(index_name), "Index Value")),
           showlegend = FALSE,
           title = paste(extra_title, unique(data$unit_id), toupper(index_name), "by", summarize_by),
           margin = list(t = 60),
           xaxis = list(title = "Date"),
           hovermode = "x unified")
  
  # Add LOESS smoothing only if there's enough data
  if (nrow(plot_tibble) > 3) {
    index_plot <- index_plot |>
      add_trace(x = ~date, 
                y = ~predict(loess(index_avg ~ as.numeric(date), 
                                   span = span, na.action = na.exclude)),
                name = "Smoothed", 
                type = 'scatter', mode = 'lines',
                line = list(color = 'darkblue', dash = 'dash', width = 3),
                hoverinfo = 'none')
  }
  # Return list with all tibbles and the plot
  list(daily_tibble = daily_tibble,
       weekly_tibble = weekly_tibble,
       monthly_tibble = monthly_tibble,
       plot = index_plot)
}
