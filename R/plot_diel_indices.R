#' Plot Diel Indices with Optional LOESS Smoothing and Scaling Transformations
#'
#' This function takes a dataset containing datetime and index values and creates 
#' an interactive plot of the scaled indices over time. Optionally, LOESS smoothing
#' can be applied to the line plots, replacing the original lines. Users can apply
#' various transformations to the data for better comparison.
#'
#' @param data A data frame containing `datetime`, `index`, and `value_avg` columns.
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

plot_diel_indices <- function(data, plot.title = "Index Diel Plot", 
                              loess = TRUE, span = 0.3, 
                              scaling = "none") {
  
  # Ensure the datetime column is properly parsed as POSIXct
  data <- data %>%
    mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"))
  
  # Check for missing or infinite values in the `value_avg` column
  if (any(is.na(data$value_avg))) {
    cat("Warning: NA values found in value_avg\n")
  }
  
  if (any(is.infinite(data$value_avg))) {
    cat("Warning: Infinite values found in value_avg\n")
  }
  
  # Convert index to uppercase and ensure no indices are missing
  data <- data %>%
    mutate(index = toupper(index))  # Convert to uppercase
  
  # Apply the selected scaling/transformation to the data
  if (scaling == "minmax") {
    data <- data %>%
      group_by(index) %>%
      mutate(value_avg = (value_avg - min(value_avg)) / (max(value_avg) - min(value_avg))) %>%
      ungroup()
  } else if (scaling == "maxabs") {
    data <- data %>%
      group_by(index) %>%
      mutate(value_avg = value_avg / max(abs(value_avg))) %>%
      ungroup()
  } else if (scaling == "zscore") {
    data <- data %>%
      group_by(index) %>%
      mutate(value_avg = (value_avg - mean(value_avg)) / sd(value_avg)) %>%
      ungroup()
  } else if (scaling == "log") {
    data <- data %>%
      mutate(value_avg = log(value_avg + 1))  # Apply log transformation (with +1 to avoid log(0))
  } else if (scaling == "robust") {
    data <- data %>%
      group_by(index) %>%
      mutate(value_avg = (value_avg - median(value_avg)) / IQR(value_avg, na.rm = TRUE)) %>%
      ungroup()
  }
  
  # Group by index and check for any groups with no valid data
  data_grouped <- data %>%
    group_by(index) %>%
    filter(!is.na(value_avg) & !is.infinite(value_avg))  # Remove invalid entries within each group
  
  # Create color palette that can handle more than 8 colors
  unique_indices <- length(unique(data$index))
  colors <- viridis(unique_indices, option = "D")
  
  # Initialize plotly object
  plot <- plot_ly()
  
  # Loop through each unique index and add either raw or LOESS-smoothed lines
  for (idx in unique(data_grouped$index)) {
    data_idx <- data_grouped %>% filter(index == idx)
    
    if (loess) {
      # Apply LOESS smoothing
      loess_model <- loess(value_avg ~ as.numeric(datetime), data = data_idx, span = span)
      smoothed_values <- predict(loess_model)
      plot <- plot %>%
        add_trace(data = data_idx, x = ~datetime, y = smoothed_values, 
                  type = 'scatter', mode = 'lines', 
                  name = idx, line = list(width = 3),
                  showlegend = TRUE)
    } else {
      # Plot raw lines
      plot <- plot %>%
        add_trace(data = data_idx, x = ~datetime, y = ~value_avg, 
                  type = 'scatter', mode = 'lines', 
                  name = idx, line = list(width = 2),
                  showlegend = TRUE)
    }
  }
  
  # Customize layout
  plot <- plot %>%
    layout(title = plot.title,
           xaxis = list(title = "Time of the Day", tickangle = 45, tickformat = "%H:%M"),
           yaxis = list(title = "Scaled Indices (0-1)"))
  
  return(plot)
}
