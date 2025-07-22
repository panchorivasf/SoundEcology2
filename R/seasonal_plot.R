#' Stack multiple diel ribbon plots vertically in chronological order
#'
#' @description 
#' This function takes individual diel ribbon PNG files and combines them into 
#' a single vertical stack with properly positioned date and time axis labels.
#'
#' @param ribbon_folder Path to folder containing ribbon PNG files (default: "./diel_plots")
#' @param output_file Output filename (default: "stacked_ribbons.png")
#' @param output_width Width of output image in inches (default: 10)
#' @param ribbon_height Height of each ribbon in inches (default: 0.3)
#' @param time_axis_breaks Positions of time axis ticks (default: seq(0, 23, by = 2))
#' @param bg_color Background color (default: "white")
#' @param label_size Font size for labels (default: 12)
#' 
#' @return Invisibly returns the ggplot object
#' @export
#' 
#' @importFrom magick image_read image_info image_ggplot
#' @importFrom lubridate ymd
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot ggsave aes annotation_raster scale_x_continuous 
#' @importFrom ggplot2 scale_y_continuous labs theme_minimal theme element_text element_blank element_rect
#' @importFrom grid unit
seasonal_plot <- function(ribbon_folder = "./diel_plots",
                          output_file = "./diel_plots/seasonal.png",
                          output_width = 10,
                          ribbon_height = 0.3,
                          time_axis_breaks = seq(0, 23, by = 2),
                          bg_color = "white",
                          label_size = 12) {
  
  # List and sort ribbon files by date
  ribbon_files <- list.files(ribbon_folder, 
                             pattern = "^ribbon_.*\\.png$", 
                             full.names = TRUE)
  
  if (length(ribbon_files) == 0) stop("No ribbon files found")
  
  # Extract and sort dates
  dates <- lubridate::ymd(str_extract(basename(ribbon_files), "([0-9]{8})"))
  ribbon_files <- ribbon_files[order(dates)]
  dates <- format(dates[order(dates)], "%Y-%m-%d")
  
  # Read and process all ribbon images
  ribbons <- lapply(ribbon_files, function(f) {
    img <- magick::image_read(f)
    # Remove alpha channel if exists
    # if (magick::image_info(img)$colorspace == "sRGBAlpha") {
    #   img <- magick::image_flatten(img, "white")
    # }
    # # Convert to raster array
    # as.raster(img)
  })
  
  # Calculate plot dimensions
  n_ribbons <- length(ribbon_files)
  total_height <- n_ribbons * ribbon_height
  
  # Create plot
  p <- ggplot() +
    # Add ribbons with exact pixel alignment
    lapply(1:n_ribbons, function(i) {
      annotation_raster(
        ribbons[[i]],
        xmin = -0.5, xmax = 23.5,  # Slightly extended to prevent gaps
        ymin = (i-1)*ribbon_height, 
        ymax = i*ribbon_height,
        interpolate = FALSE  # Prevents blurring that can cause artifacts
      )
    }) +
    # Axes and labels
    scale_x_continuous(
      "Time (hours)",
      breaks = time_axis_breaks,
      expand = c(0, 0),
      limits = c(0, 23)
    ) +
    scale_y_continuous(
      "Date",
      breaks = (1:n_ribbons)*ribbon_height - ribbon_height/2,
      labels = dates,
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = bg_color, colour = NA),
      # plot.background = element_rect(fill = bg_color, colour = NA),
      axis.text = element_text(size = label_size),
      axis.title = element_text(size = label_size + 2),
      # Explicitly remove all grid lines
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # Save with exact dimensions
  ggsave(
    output_file,
    plot = p,
    width = output_width,
    height = total_height + 1,  # Extra space for labels
    dpi = 300,
    bg = bg_color
  )
  
  invisible(p)
}