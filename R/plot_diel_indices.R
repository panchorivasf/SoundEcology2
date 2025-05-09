#' Plot Diel Indices with Optional LOESS Smoothing and Scaling Transformations
#'
#' This function takes a dataset containing datetime and index values and creates 
#' an interactive plot of the scaled indices over time. Optionally, LOESS smoothing
#' can be applied to the line plots, replacing the original lines. Users can apply
#' various transformations to the data for better comparison.
#'
#' @param data A data frame containing `datetime`, `index`, and `value_avg` columns.
#' @param sensor.id Character. Sensor identifier as expressed in the 'sensor_id' column. If not NULL, only the data for the requester sensor will be plotted.
#' @param plot.title A character string specifying the title of the plot. Defaults to "Diel Index Patterns".
#' @param loess A logical value indicating whether to apply LOESS smoothing to the lines. Defaults to `TRUE`.
#' @param span A numeric value indicating the smoothing parameter for LOESS. Defaults to 0.3.
#' @param scaling A character string specifying the type of scaling/transformation to apply. 
#'   Options are "none" (default), "minmax", "maxabs", "zscore", "log", and "robust".
#'
#' @return A `plotly` object showing a time series of scaled indices for each index group.
#' @export
#' 
#' @importFrom dplyr mutate group_by filter ungroup
#' @importFrom plotly plot_ly layout add_trace
#' @importFrom viridis viridis
#' @importFrom stats loess predict
#'
#' @examples
#' \dontrun{
#'   # Assuming `data` is a data frame with the necessary columns
#'   plot_diel_indices(data, loess = TRUE, scaling = "log")
#' }

plot_diel_indices <- function(data, 
                              sensor_id = NULL, 
                              plot.title = "Circadian Indices", 
                              loess = TRUE, 
                              span = 0.3, 
                              scaling = "none") {
  
  if (!scaling %in% c("none", "minmax", "maxabs", "zscore", "log", "robust")) {
    stop("'scaling' can only be 'none', 'minmax', 'maxabs', 'zscore', 'log', or 'robust'")
  }
  
  # Ensure the datetime column is properly parsed as POSIXct
  data <- data |>
    mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
           hour = format(datetime, "%H")) # Extract hour for grouping
  
  # Check for missing or infinite values in the `value_avg` column
  if (any(is.na(data$value_avg))) {
    cat("Warning: NA values found in value_avg\n")
  }
  
  if (any(is.infinite(data$value_avg))) {
    cat("Warning: Infinite values found in value_avg\n")
  }
  
  # Convert index to uppercase and ensure no indices are missing
  data <- data |>
    mutate(index = toupper(index))  # Convert to uppercase
  
  # Apply the selected scaling/transformation to the data
  if (scaling == "minmax") {
    data <- data |>
      group_by(index) |>
      mutate(value_avg = (value_avg - min(value_avg)) / (max(value_avg) - min(value_avg))) |>
      ungroup()
  } else if (scaling == "maxabs") {
    data <- data |>
      group_by(index) |>
      mutate(value_avg = value_avg / max(abs(value_avg))) |>
      ungroup()
  } else if (scaling == "zscore") {
    data <- data |>
      group_by(index) |>
      mutate(value_avg = (value_avg - mean(value_avg)) / sd(value_avg)) |>
      ungroup()
  } else if (scaling == "log") {
    data <- data |>
      mutate(value_avg = log(value_avg + 1))  # Apply log transformation (with +1 to avoid log(0))
  } else if (scaling == "robust") {
    data <- data |>
      group_by(index) |>
      mutate(value_avg = (value_avg - median(value_avg)) / IQR(value_avg, na.rm = TRUE)) |>
      ungroup()
  }
  
  # If a specific sensor_id is provided, filter for that sensor
  if (!is.null(sensor_id)) {
    data <- data |>
      filter(sensor_id == sensor_id)
    
    if (nrow(data) == 0) {
      stop("No data found for the specified sensor_id.")
    }
    
    # Group by hour and index only for the selected sensor
    data_grouped <- data |>
      group_by(hour, index) |>
      summarise(value_avg = mean(value_avg, na.rm = TRUE)) |>
      ungroup()
  } else {
    # Group by hour and index across all sensors
    data_grouped <- data |>
      group_by(hour, index) |>
      summarise(value_avg = mean(value_avg, na.rm = TRUE)) |>
      ungroup()
  }
  
  # Create color palette that can handle more than 8 colors
  unique_indices <- length(unique(data_grouped$index))
  colors <- viridis(unique_indices, option = "D")
  
  # Initialize plotly object
  plot <- plot_ly()
  
  # Loop through each unique index and add either raw or LOESS-smoothed lines
  for (idx in unique(data_grouped$index)) {
    data_idx <- data_grouped |> filter(index == idx)
    
    if (loess) {
      # Apply LOESS smoothing with interpolation to match the original 'hour' values
      loess_model <- loess(value_avg ~ as.numeric(hour), data = data_idx, span = span)
      
      smoothed_values <- predict(loess_model, newdata = as.numeric(data_idx$hour))
      
      plot <- plot |>
        add_trace(data = data_idx, x = ~hour, y = smoothed_values, 
                  type = 'scatter', mode = 'lines', 
                  name = idx, line = list(width = 3, shape = 'spline'),
                  showlegend = TRUE,
                  hovertemplate = paste('Time: %{x}<br>',
                                        'Value: %{y:.2f}<extra></extra>'))
    } else {
      # Plot raw lines
      plot <- plot |>
        add_trace(data = data_idx, x = ~hour, y = ~value_avg, 
                  type = 'scatter', mode = 'lines', 
                  name = idx, line = list(width = 2,shape = 'spline'),
                  showlegend = TRUE,
                  hovertemplate = paste('Time: %{x}<br>',
                                        idx,': %{y:.2f}<extra></extra>'))
    }
  }
  
  
  # Customize layout
  plot <- plot |>
    layout(title = plot.title,
           xaxis = list(title = "Hour of the Day", tickangle = 90),
           yaxis = list(title = "Scaled Indices (0-1)"))
  
  return(plot)
}
