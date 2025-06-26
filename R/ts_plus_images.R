#' Create Interactive Time Series Plot with Image Popups
#'
#' Creates an interactive time series plot using Chart.js where users can click
#' on data points to view corresponding images in a modal popup. Images are
#' matched to data points based on dates extracted from filenames. Includes
#' zoom functionality on the horizontal axis.
#'
#' @param data A data.frame containing the time series data. Must include a 'date' column.
#' @param image_folder Character. Path to folder containing images. Images should have
#'   dates in YYYYMMDD format in their filenames (e.g., "sensor_data_20231215.png").
#'   If NULL, no images will be displayed. Default is NULL.
#' @param output_file Character. Path for the output HTML file. If provided, creates
#'   a standalone HTML file with accompanying image folder. If NULL, creates a
#'   temporary file and opens in browser. Default is NULL.
#' @param y_column Character. Name of the column in data to use for the y-axis values.
#'   Default is "mean_centroid".
#' @param variable_name Character. Descriptive name for the y-variable to use in
#'   plot labels and title (e.g., "Temperature", "BBAI Centroid"). If NULL,
#'   uses the y_column name. Default is NULL.
#' @param identifier Character. An identifier to include in the plot title
#'   (e.g., sensor ID, location name). If provided, will be shown as
#'   "Interactive Time Series - [variable_name] ([identifier])". Default is NULL.
#'
#' @return Invisibly returns the HTML content as a character string.
#'
#' @details
#' The function creates an interactive time series plot where:
#' \itemize{
#'   \item Data points with corresponding images appear as red circles
#'   \item Data points without images appear as blue circles
#'   \item Clicking on data points shows detailed information
#'   \item Clicking on points with images opens a modal with the full-size image
#'   \item The plot includes tooltips and keyboard support (Esc to close modals)
#'   \item Zoom controls for the horizontal axis (drag to zoom, double-click to reset)
#' }
#'
#' Image files should contain dates in YYYYMMDD format in their filenames.
#' Supported image formats: PNG, JPG, JPEG, GIF, BMP.
#'  
#' When output_file is specified, the function creates:
#' \itemize{
#'   \item A standalone HTML file at the specified location
#'   \item A "ts_images" subfolder containing copies of all relevant images
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' ts_plus_images(my_data)
#'
#' # With images and custom labels
#' ts_plus_images(my_data, 
#'                image_folder = "path/to/images",
#'                y_column = "temperature",
#'                variable_name = "Temperature (°C)",
#'                identifier = "Sensor-001")
#'
#' # Save to HTML file
#' ts_plus_images(my_data,
#'                image_folder = "path/to/images",
#'                output_file = "temperature_analysis.html",
#'                y_column = "temp",
#'                variable_name = "Daily Temperature",
#'                identifier = "Weather Station A")
#' }
#'
#' @importFrom lubridate ymd
#' @importFrom jsonlite toJSON
#'
#' @export
ts_plus_images <- function(data, 
                           image_folder = NULL, 
                           output_file = NULL, 
                           var_column = "mean_centroid",
                           variable_name = NULL,
                           identifier = NULL) {
  
  # Load required libraries
  if (!require(lubridate)) stop("lubridate package is required")
  if (!require(jsonlite)) stop("jsonlite package is required")
  
  # Validate inputs
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (!"date" %in% names(data)) stop("data must contain a 'date' column")
  if (!var_column %in% names(data)) stop(paste("Column", var_column, "not found in data"))
  
  # Set variable name if not provided
  if (is.null(variable_name)) {
    variable_name <- var_column
  }
  
  # Create plot title
  plot_title <- if (!is.null(identifier)) {
    paste0("Interactive Time Series - ", variable_name, " (", identifier, ")")
  } else {
    paste0("Interactive Time Series - ", variable_name)
  }
  
  # Convert date to proper format
  data$date <- lubridate::ymd(data$date)
  
  # Remove rows with NA dates or y values
  data <- data[!is.na(data$date) & !is.na(data[[var_column]]), ]
  
  # Sort by date
  data <- data[order(data$date), ]
  
  # Function to extract date from filename
  extract_date_from_filename <- function(filename) {
    date_pattern <- "\\d{8}"
    date_match <- regmatches(filename, regexpr(date_pattern, filename))
    if (length(date_match) > 0) {
      parsed_date <- tryCatch({
        lubridate::ymd(date_match)
      }, error = function(e) {
        return(NA)
      })
      if (!is.na(parsed_date) && inherits(parsed_date, "Date")) {
        return(parsed_date)
      }
    }
    return(as.Date(NA))
  }
  
  # Prepare image data if image_folder is provided
  data$image_path <- NA
  if (!is.null(image_folder) && dir.exists(image_folder)) {
    # Get all image files
    image_files <- list.files(image_folder, 
                              pattern = "\\.(png|jpg|jpeg|gif|bmp)$", 
                              ignore.case = TRUE, 
                              full.names = TRUE)
    
    if (length(image_files) > 0) {
      # Create a data frame of images with their extracted dates
      image_data <- data.frame(
        filepath = image_files,
        filename = basename(image_files),
        stringsAsFactors = FALSE
      )
      
      # Extract dates from filenames
      image_data$image_date <- sapply(image_data$filename, extract_date_from_filename)
      if (!inherits(image_data$image_date, "Date")) {
        image_data$image_date <- as.Date(image_data$image_date, origin = "1970-01-01")
      }
      
      # Remove images where date extraction failed
      image_data <- image_data[!is.na(image_data$image_date), ]
      
      if (nrow(image_data) > 0) {
        # Match images to data points
        for (i in 1:nrow(data)) {
          matching_image <- image_data[image_data$image_date == data$date[i], ]
          if (nrow(matching_image) > 0) {
            data$image_path[i] <- matching_image$filepath[1]
          }
        }
        # Convert file paths to relative paths for web display
        # Copy images to a temporary directory for web access
        temp_dir <- file.path(tempdir(), "ts_images")
        if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
        
        data$web_image_path <- NA
        for (i in 1:nrow(data)) {
          if (!is.na(data$image_path[i])) {
            # Copy image to temp directory
            img_name <- basename(data$image_path[i])
            new_path <- file.path(temp_dir, img_name)
            file.copy(data$image_path[i], new_path, overwrite = TRUE)
            data$web_image_path[i] <- file.path("ts_images", img_name)
          }
        }
      }
    }
  }
  
  # Prepare data for JSON
  chart_data <- data.frame(
    date = as.character(data$date),
    value = data[[var_column]],
    image = ifelse(is.na(data$web_image_path), "", data$web_image_path),
    stringsAsFactors = FALSE
  )
  
  # Convert to JSON
  json_data <- jsonlite::toJSON(chart_data, auto_unbox = FALSE)
  
  # Create HTML content
  html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>', plot_title, '</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/3.9.1/chart.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-zoom@2.0.1"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 2500px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .chart-container {
            position: relative;
            height: 500px;
            margin-bottom: 20px;
        }
        .info-panel {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.8);
        }
        .modal-content {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            max-width: 90%;
            max-height: 90%;
            overflow: auto;
        }
        .modal-image {
            max-width: 100%;
            max-height: 70vh;
            object-fit: contain;
        }
        .close {
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            cursor: pointer;
            line-height: 1;
        }
        .close:hover {
            color: black;
        }
        .point-info {
            margin-top: 15px;
            padding: 10px;
            background-color: #e9ecef;
            border-radius: 4px;
        }
        .instructions {
            color: #666;
            font-style: italic;
            margin-bottom: 15px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>', plot_title, '</h1>
        <div class="instructions">
            Click on any data point to view the corresponding image (if available).<br>
            Use the mouse wheel to zoom or drag to select zoom area.
        </div>
        <div class="info-panel">
            <div id="point-info">Click on a data point to see details</div>
        </div>
        <div class="chart-container">
            <canvas id="timeSeriesChart"></canvas>
        </div>
    </div>

    <!-- Modal for image display -->
    <div id="imageModal" class="modal">
        <div class="modal-content">
            <span class="close">&times;</span>
            <div id="modal-body">
                <img id="modal-image" class="modal-image" src="" alt="">
                <div id="modal-info" class="point-info"></div>
            </div>
        </div>
    </div>

    <script>
        // Data from R
        const chartData = ', json_data, ';
        
        // Prepare data for Chart.js
        const labels = chartData.map(d => d.date);
        const values = chartData.map(d => d.value);
        
        // Chart configuration
        const config = {
            type: "line",
            data: {
                labels: labels,
                datasets: [{
                    label: "', variable_name, '",
                    data: values,
                    borderColor: "rgb(30, 144, 255)",
                    backgroundColor: "rgba(30, 144, 255, 0.1)",
                    borderWidth: 2,
                    pointRadius: 6,
                    pointHoverRadius: 8,
                    pointBackgroundColor: function(context) {
                        const index = context.dataIndex;
                        return chartData[index].image ? "rgb(255, 99, 132)" : "rgb(30, 144, 255)";
                    },
                    pointBorderColor: function(context) {
                        const index = context.dataIndex;
                        return chartData[index].image ? "rgb(255, 99, 132)" : "rgb(30, 144, 255)";
                    }
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,

                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        callbacks: {
                            afterLabel: function(context) {
                                const index = context.dataIndex;
                                return chartData[index].image ? "📷 Click to view image" : "No image available";
                            }
                        }
                    },
                    zoom: {
                        limits: {
                            x: {min: "original", max: "original"}
                        },
                        pan: {
                            enabled: true,
                            mode: "x"
                        },
                        zoom: {
                            wheel: {
                                enabled: true
                            },
                            drag: {
                                enabled: true,
                                mode: "x",
                                backgroundColor: "rgba(30, 144, 255, 0.3)",
                                borderColor: "rgb(30, 144, 255)",
                                borderWidth: 1
                            },
                            pinch: {
                                enabled: true,
                                mode: "x"
                            },
                            mode: "x"
                        }
                    }
                },
                scales: {
                    x: {
                        display: true,
                        title: {
                            display: true,
                            text: "Date"
                        }
                    },
                    y: {
                        display: true,
                        title: {
                            display: true,
                            text: "', variable_name, '"
                        }
                    }
                },
                onClick: function(event, elements) {
                    if (elements.length > 0) {
                        const index = elements[0].index;
                        const dataPoint = chartData[index];
                        
                        // Update info panel
                        document.getElementById("point-info").innerHTML = 
                            `<strong>Date:</strong> ${dataPoint.date}<br>
                             <strong>', variable_name, ':</strong> ${dataPoint.value.toFixed(4)}<br>
                             <strong>Image:</strong> ${dataPoint.image ? "Available (click point to view)" : "Not available"}`;
                        
                        // Show image if available
                        if (dataPoint.image) {
                            showImageModal(dataPoint);
                        }
                    }
                },
                onDoubleClick: function(event, elements) {
                    // Reset zoom on double click
                    chart.resetZoom();
                }
            }
        };
        
        // Register the zoom plugin
        Chart.register(ChartZoom);
        
        // Create chart
        const ctx = document.getElementById("timeSeriesChart").getContext("2d");
        const chart = new Chart(ctx, config);
        
        // Modal functionality
        const modal = document.getElementById("imageModal");
        const modalImage = document.getElementById("modal-image");
        const modalInfo = document.getElementById("modal-info");
        const closeBtn = document.getElementsByClassName("close")[0];

        function showImageModal(dataPoint) {
            modalImage.src = dataPoint.image;
            modalInfo.innerHTML = 
                `<strong>Date:</strong> ${dataPoint.date}<br>
                 <strong>', variable_name, ':</strong> ${dataPoint.value.toFixed(4)}`;
            modal.style.display = "block";
        }

        closeBtn.onclick = function() {
            modal.style.display = "none";
        }

        window.onclick = function(event) {
            if (event.target == modal) {
                modal.style.display = "none";
            }
        }

        // Keyboard support
        document.addEventListener("keydown", function(event) {
            if (event.key === "Escape") {
                modal.style.display = "none";
            }
        });
    </script>
</body>
</html>')
  
  # Save to file if specified
  if (!is.null(output_file)) {
    # Ensure the file has .html extension
    if (!grepl("\\.html$", output_file, ignore.case = TRUE)) {
      output_file <- paste0(output_file, ".html")
    }
    
    # Get the directory of the output file
    output_dir <- dirname(normalizePath(output_file, mustWork = FALSE))
    
    # Check if images are already accessible
    if (dirname(output_file) == dirname(image_folder)) {
      # Use relative paths directly, no copying needed
      data$web_image_path <- file.path(basename(image_folder), basename(data$image_path))
    } else {
      # Copy images to the same directory as the HTML file
      if (!is.null(image_folder) && any(!is.na(data$image_path))) {
        img_output_dir <- file.path(output_dir, "ts_images")
        if (!dir.exists(img_output_dir)) dir.create(img_output_dir, recursive = TRUE)
        
        for (i in 1:nrow(data)) {
          if (!is.na(data$image_path[i])) {
            img_name <- basename(data$image_path[i])
            file.copy(data$image_path[i], file.path(img_output_dir, img_name), overwrite = TRUE)
          }
        }
      }
    }
    
    # Write HTML file
    writeLines(html_content, output_file)
    message(paste("Interactive plot saved to:", output_file))
    
    if (!is.null(image_folder)) {
      message(paste("Images copied to:", file.path(output_dir, "ts_images")))
    }
  } else {
    # Create temporary file for viewing
    temp_file <- file.path(tempdir(), "ts_plus_images.html")
    writeLines(html_content, temp_file)
    
    # Try to open in browser
    if (interactive()) {
      browseURL(paste0("file://", temp_file))
    }
    
    message(paste("Temporary plot created at:", temp_file))
  }
  
  return(invisible(html_content))
}