#' Generate a Diel Spectrogram from multiple Audio Recordings
#'
#' Creates a spectrogram visualization of daily acoustic patterns using 
#' multiple recordings, with options for noise floor clipping and parallel 
#' processing. Specifically, it takes the variance of acoustic energy in each
#' frequency bin for each recording and merges the resulting vectors into a 
#' summary spectrogram. 
#'
#' @param folder Path to the folder containing WAV files (default: current 
#' working directory).
#' @param recursive Logical. Whether to search files in subfolders. Default: TRUE.
#' @param list A list with WAV files to parse. The waves must be in the working 
#' directory (default: NULL).
#' @param title Character string. Plot title (default: empty). 
#' @param freq_res Frequency resolution in Hz for FFT (default: 50)
#' @param n.cores Number of cores for parallel processing (-1 for all available 
#' cores)
#' @param cutoff Minimum dBFS value (noise floor). Values below this will be 
#' clipped (default: -100)
#' @param max_amp Maximum dBFS value to be represented by the color scale. 
#' @param plot Whether to generate the plot (default: TRUE)
#' @param save_plot Logical. Whether to save the plot
#' @param save_to Character. If save_plot = TRUE, provide a path to the folder 
#' where the file should be stored. If no path is provided, a new folder
#' "diel_plots" is created in the current working directory (default) to save 
#' the exported files.
#' @param dc_on Number of minutes where the recorded was "on". Used to identify
#' expected number of files per day.
#' @param dc_off Number of minutes where the recorded was "off". Used to identify
#' expected number of files per day.
#' @param target_dates Either a single date or multiple dates in the format 
#' "YYYY-MM-DD" to be analyzed. If NULL (default), all the dates in the folder
#' will be parsed. To provide a range, see Details. 
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - The ggplot2 spectrogram object
#'   \item spectral_data - Wide-format spectral data (freq bins × recordings)
#'   \item plot_data - Long-format data for plotting
#'   \item file_metadata - File information with extracted metadata
#'   \item overall_range - The dBFS range of all data
#' }
#' @export
#'
#' @importFrom tuneR readWave mono
#' @importFrom seewave meanspec
#' @importFrom stringr str_extract
#' @importFrom lubridate as_datetime hour date
#' @importFrom dplyr mutate select everything case_when
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn scale_x_continuous scale_y_continuous labs theme_minimal
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' 
#' @details
#' This function parses files in a folder or list and groups them by date, 
#' creating a single diel spectrogram for each. To parse a subset range of dates 
#' contained in a folder, use seq().
#' Example: target_dates = seq(as.Date("2022-10-20"), as.Date("2022-10-25), "days")
diel_spectrogram <- function(folder = NULL,
                             recursive = TRUE,
                             list = NULL,
                             title = "",
                             freq_res = 50,
                             n.cores = -1,
                             cutoff = -100,
                             max_amp = -10,
                             plot = TRUE,
                             save_plot = TRUE,
                             save_to = "./diel_plots",
                             dc_on = NULL,
                             dc_off = NULL,
                             target_dates = NULL) {
  
  # Validate input parameters
  if (!is.null(folder) && !is.null(list)) {
    stop("Please specify only 'folder' OR 'list', not both.")
  }
  
  if (is.null(folder) && is.null(list)) {
    folder <- getwd()
    message("No input specified, using current working directory")
  }
  
  # Validate duty cycle parameters
  if (!is.null(dc_on) && !is.null(dc_off)) {
    if (!is.numeric(dc_on) || !is.numeric(dc_off) || dc_on <= 0 || dc_off <= 0) {
      stop("'dc_on' and 'dc_off' must be positive numbers")
    }
  }
  
  # Validate target_dates parameter
  if (!is.null(target_dates)) {
    tryCatch({
      target_dates <- as.Date(target_dates)
    }, error = function(e) {
      stop("'target_dates' must be in format 'YYYY-MM-DD' or a vector of such dates")
    })
  }
  
  # Create file dataframe based on input type
  if (!is.null(list)) {
    # Validate list input
    if (!is.character(list) || !all(file.exists(list))) {
      stop("'list' must be a character vector of valid file paths")
    }
    
    file_df <- tibble::tibble(
      file_name = basename(list),
      file_path = list
    )
  } else {
    file_paths <- list.files(folder, 
                             pattern = "\\.wav$", 
                             full.names = TRUE, 
                             ignore.case = TRUE,
                             recursive = recursive)
    
    file_df <- tibble::tibble(
      file_name = basename(file_paths),  # This extracts just the filename without path
      file_path = file_paths
    )
    
    if (nrow(file_df) == 0) stop("No WAV files found in folder")
  }
  
  # Add metadata using internal function
  file_df <- addMetadata(file_df)
  
  # Filter by target_dates if specified
  if (!is.null(target_dates)) {
    file_df <- file_df |>
      dplyr::mutate(date = as.Date(datetime)) |>
      dplyr::filter(date %in% target_dates)
    
    if (nrow(file_df) == 0) {
      stop("No files found for the specified target_dates")
    }
    
    message("Filtered to ", length(target_dates), " target date(s): ", 
            paste(target_dates, collapse = ", "))
  }
  
  # If duty cycle parameters provided, check for expected files
  if (!is.null(dc_on) && !is.null(dc_off)) {
    dc_interval <- dc_on + dc_off  # Total duty cycle interval in minutes
    
    # Split files by date and check each day's sequence
    file_df <- file_df |>
      dplyr::mutate(date = as.Date(datetime)) |>
      dplyr::group_by(date) |>
      dplyr::group_modify(~ {
        day_files <- .x
        expected_times <- seq(
          from = min(day_files$datetime),
          to = max(day_files$datetime),
          by = dc_interval * 60  # Convert minutes to seconds
        )
        
        # Create complete sequence of expected files
        all_times <- data.frame(datetime = expected_times)
        day_files <- dplyr::full_join(day_files, all_times, by = "datetime") |>
          dplyr::arrange(datetime) |>
          dplyr::mutate(
            file_exists = !is.na(file_path),
            time_diff = as.numeric(difftime(datetime, dplyr::lag(datetime), units = "mins")),
            gap_detected = ifelse(is.na(time_diff), FALSE, time_diff > dc_interval * 1.1)  # 10% tolerance
          )
        
        return(day_files)
      }) |>
      dplyr::ungroup()
  } else {
    file_df <- file_df |>
      dplyr::mutate(
        date = as.Date(datetime),
        file_exists = TRUE,
        gap_detected = FALSE
      )
  }
  
  # Split files by date
  file_dfs <- file_df |>
    dplyr::group_split(date)
  
  # Process each date separately
  results <- lapply(file_dfs, function(date_files) {
    current_date <- unique(as.Date(date_files$datetime))
    message("\nProcessing date ", current_date, "...")
    
    # Set up parallel processing
    if (n.cores == -1) n.cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    
    # Export necessary variables to parallel workers
    parallel::clusterExport(cl, varlist = c("freq_res", "cutoff"), envir = environment())
    
    # Get frequency bins from first valid file
    first_valid <- which(date_files$file_exists)[1]
    if (is.na(first_valid)) {
      warning("No valid files found for date ", current_date)
      return(NULL)
    }
    
    wave <- tryCatch(
      tuneR::readWave(date_files$file_path[first_valid]),
      error = function(e) NULL
    )
    
    if (is.null(wave)) {
      warning("Cannot determine frequency bins - first file is corrupted for date ", current_date)
      return(NULL)
    }
    
    wl <- round(wave@samp.rate / freq_res)
    if (wl %% 2 == 1) wl <- wl + 1
    freq_bins <- seewave::meanspec(wave, wl = wl, plot = FALSE)[, 1]
    na_vector <- rep(NA, length(freq_bins))
    
    message("Found ", sum(date_files$file_exists), " valid files out of ", nrow(date_files), " expected")
    if (any(date_files$gap_detected, na.rm = TRUE)) {
      message("Gaps detected in the recording sequence")
    }
    
    # Process files in parallel (including NA placeholders for missing files)
    spec_results <- foreach::foreach(i = 1:nrow(date_files), .packages = c("tuneR", "seewave")) %dopar% {
      if (!date_files$file_exists[i]) {
        return(list(
          freq = freq_bins,
          amp = na_vector,
          datetime = date_files$datetime[i],
          is_missing = TRUE
        ))
      }
      
      tryCatch({
        wave <- tuneR::readWave(date_files$file_path[i])
        
        if (wave@stereo) wave <- tuneR::mono(wave, "both")
        
        amp_max <- switch(as.character(wave@bit),
                          "16" = 32768,
                          "24" = 8388607,
                          "32" = 2147483647,
                          stop("Unsupported bit depth"))
        
        spec <- seewave::meanspec(wave, 
                                  FUN = var, 
                                  norm = FALSE, 
                                  correction = "amplitude",
                                  wl = wl, 
                                  plot = FALSE)
        
        # Convert to dBFS with clipping at cutoff and 0 dBFS
        spec_dBFS <- 10 * log10(spec[, 2] / (amp_max^2))
        spec_dBFS <- pmin(pmax(spec_dBFS, cutoff), 0)
        
        list(
          freq = freq_bins,
          amp = spec_dBFS,
          datetime = date_files$datetime[i],
          is_missing = FALSE
        )
      }, error = function(e) {
        warning("Error processing file ", date_files$file_path[i], ": ", e$message)
        list(
          freq = freq_bins,
          amp = na_vector,
          datetime = date_files$datetime[i],
          is_missing = TRUE
        )
      })
    }
    
    parallel::stopCluster(cl)
    
    # Create frequency bin dataframe
    spec_df <- tibble::tibble(Frequency = freq_bins)
    
    # Add each recording as a column
    for (i in seq_along(spec_results)) {
      col_name <- format(spec_results[[i]]$datetime, "%Y%m%d_%H%M%S")
      spec_df[[col_name]] <- spec_results[[i]]$amp
    }
    
    # Create plotting dataframe
    plot_df <- spec_df |>
      tidyr::pivot_longer(
        cols = -Frequency,
        names_to = "datetime",
        values_to = "Amplitude"
      ) |>
      dplyr::mutate(
        datetime = as.POSIXct(datetime, format = "%Y%m%d_%H%M%S"),
        Time = as.numeric(difftime(datetime, trunc(min(datetime), "days"), 
                                   units = "hours")),
        is_missing = ifelse(all(is.na(Amplitude)), TRUE, FALSE)
      )
    
    # Calculate and report range (excluding missing data)
    valid_amps <- plot_df$Amplitude[!plot_df$is_missing]
    overall_range <- if (length(valid_amps) > 0) range(valid_amps, 
                                                       na.rm = TRUE) else c(NA, NA)
    message("Overall dBFS range for date ", current_date, ": [", 
            round(overall_range[1], 2), ", ", 
            round(overall_range[2], 2), "]")
    
    # Generate plot if requested
    p <- NULL
    if (plot) {
      
      # Create base subtitle without missing data note
      plot_subtitle <- paste("Sensor:", unique(date_files$sensor_id),
                             "\nDate:", current_date)
      
      # Add missing data note if needed
      if (any(plot_df$is_missing)) {
        plot_subtitle <- paste0(plot_subtitle, "  ! Missing data (gray)")
      }
      
      
      # Create base plot
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Time, y = Frequency, 
                                                 fill = Amplitude)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(
          colors = c("black", "blue", "yellow", "red"),
          name = "Variance\n(dBFS)",
          limits = c(cutoff, max_amp),
          na.value = "gray70"  
        ) +
        ggplot2::scale_x_continuous(
          expand = c(0, 0),
          breaks = seq(0, 24, 2),
          labels = seq(0, 24, 2),
          limits = c(0, 24)
        ) +
        ggplot2::scale_y_continuous(
          labels = function(x) x,
          expand = c(0, 0)
        ) +
        ggplot2::labs(
          x = "Time (hours)",
          y = "Frequency (kHz)",
          title = title,
          subtitle = plot_subtitle
        ) +
        ggplot2::theme_classic()
      
      
      
      # Add rectangles for missing data periods
      if (any(plot_df$is_missing)) {
        missing_ranges <- plot_df |>
          dplyr::filter(is_missing) |>
          dplyr::group_by(grp = cumsum(c(1, diff(Time) > 1))) |>
          dplyr::summarize(
            xmin = min(Time),
            xmax = max(Time),
            .groups = "drop"
          )
        
        for (i in 1:nrow(missing_ranges)) {
          p <- p + ggplot2::annotate(
            "rect",
            xmin = missing_ranges$xmin[i],
            xmax = missing_ranges$xmax[i],
            ymin = min(freq_bins),
            ymax = max(freq_bins),
            fill = "gray70",
            alpha = 0.7
          )
        }
      }
      
      print(p)
      
      if (save_plot) {
        if (is.null(save_to)) {
          destination <- getwd()
        } else {
          destination <- save_to
          
          # Create directory if it doesn't exist
          if (!dir.exists(destination)) {
            dir.create(destination, recursive = TRUE)
            message("Created directory: ", destination)
          }
          
        }
        
        sensor <- unique(date_files$sensor_id)
        date <- current_date
        
        ggplot2::ggsave(
          filename = paste0("diel_",sensor, "_", gsub("-", "", date), ".png"),
          path = destination,
          plot = p,
          width = 10,
          height = 5.6,
          dpi = 300,
          units = "in"
        )
      }
    }
    
    list(
      plot = p,
      spectral_data = spec_df,
      plot_data = plot_df,
      file_metadata = date_files,
      overall_range = overall_range,
      missing_data = any(plot_df$is_missing)
    )
  })
  
  # Remove NULL results (dates with no valid files)
  results <- purrr::compact(results)
  
  # If only one date, return that directly, otherwise return list of results
  if (length(results) == 1) {
    invisible(results[[1]])
  } else {
    invisible(results)
  }
}
