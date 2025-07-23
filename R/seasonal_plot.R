#' Stack multiple diel ribbon plots vertically in chronological order
#'
#' @description 
#' This function takes individual diel ribbon PNG files and combines them into 
#' a single vertical stack with properly positioned date and time axis labels,
#' including small ticks on both axes.
#'
#' @param ribbon_folder Path to folder containing ribbon PNG files (default: "./ribbons")
#' @param output_file Output filename (default: "./seasonal.png")
#' @param output_width Width of output image in inches (default: 10)
#' @param ribbon_height Height of each ribbon in inches (default: 0.3)
#' @param time_axis_breaks Positions of time axis ticks (default: seq(0, 24, by = 1))
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
#' @importFrom ggplot2 scale_y_continuous labs theme_minimal theme element_text 
#' @importFrom ggplot2 element_blank element_rect element_line
#' @importFrom grid unit
seasonal_plot <- function(ribbon_folder = "./ribbons",
                          output_file = "./seasonal.png",
                          output_width = 10,
                          ribbon_height = 0.05,
                          time_axis_breaks = seq(0, 24, by = 1),
                          bg_color = "white",
                          label_size = 10) {
  
  # List and sort ribbon files by date
  ribbon_files <- list.files(ribbon_folder, 
                             pattern = "^ribbon_.*\\.png$", 
                             full.names = TRUE)
  
  if (length(ribbon_files) == 0) stop("No ribbon files found")
  
  # Extract and sort dates
  dates <- lubridate::ymd(stringr::str_extract(basename(ribbon_files), "([0-9]{8})"))
  ribbon_files <- ribbon_files[order(dates)]
  dates <- format(dates[order(dates, decreasing = TRUE)], "%Y-%m-%d")
  
  # Read all ribbon images
  ribbons <- lapply(ribbon_files, magick::image_read)
  
  # Calculate plot dimensions
  n_ribbons <- length(ribbon_files)
  total_height <- n_ribbons * ribbon_height
  
  # Create plot with proper axis labels and ticks
  p <- ggplot2::ggplot() +
    # Add ribbons
    lapply(1:n_ribbons, function(i) {
      ggplot2::annotation_raster(
        ribbons[[i]],
        xmin = 0, xmax = 24,
        ymin = total_height - i*ribbon_height,  # Changed this line
        ymax = total_height - (i-1)*ribbon_height  # And this line
        # ymin = (i-1)*ribbon_height, 
        # ymax = i*ribbon_height
      )
    }) +
    # X-axis (time) with minor ticks
    ggplot2::scale_x_continuous(
      "Time (hours)",
      breaks = time_axis_breaks,
      minor_breaks = 0:23,  
      expand = c(0, 0),
      limits = c(0, 24)
    ) +
    ggplot2::scale_y_continuous(
      breaks = (1:n_ribbons)*ribbon_height - ribbon_height/2,
      minor_breaks = (0:n_ribbons)*ribbon_height,  
      labels = dates,
      expand = c(0, 0),
      limits = c(0, total_height)
    ) +
    ggplot2::labs(x = "Time (hours)", y = NULL) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = bg_color, colour = NA),
      plot.background = ggplot2::element_rect(fill = bg_color, colour = NA),
      axis.text = ggplot2::element_text(size = label_size),
      axis.title = ggplot2::element_text(size = label_size + 2),
      axis.text.x = ggplot2::element_text(vjust = 1),
      axis.ticks.y = element_blank()
    )
  
  # Save with exact dimensions
  ggplot2::ggsave(
    output_file,
    plot = p,
    width = output_width,
    height = total_height + 1,  
    dpi = 300,
    bg = bg_color
  )
  
  invisible(p)
}