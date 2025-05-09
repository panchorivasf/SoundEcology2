#' ADI Heatmap
#' @description
#' Produces a heatmap with the proportions used to calculate ADI for a day of data. The current version supports only data frames describing stereo files with a duty cycle of 1 on; 9 off.
#'
#' @param df A data frame with ADI data for 1 day
#' @param plot.title Character; a title for the plot
#'
#' @return A heatmap plot created with ggplot2.
#' @export
#'
#' @import tidyverse
#' @import ggplot2
#' @import lubridate
#' @import patchwork
#'
#' @examples
#' adi_heatap(adi20230903, "ADI proportions")
adi_heatmap <- function(df, plot.title) {

  if (!any(names(df) == "noise_red")) {
    df$noise_red <- "none"
  }
  df <- df |> relocate(noise_red, .after = "norm")

  # Extract settings
  settings <- df |> select(cutoff, norm, noise_red, prop_den, w_len, w_fun) |> distinct()
  settings$cutoff <- paste0("cutoff.", settings$cutoff)
  settings$norm <- paste0("norm.", settings$norm)
  settings$noise_red <- paste0("noise_red.", settings$noise_red)
  settings$prop_den <- paste0("prop.", settings$prop_den)
  settings$w_len <- paste0("wl.", settings$w_len)
  settings$w_fun <- paste0("wf.", settings$w_fun)

  # Select relevant columns dynamically
  relevant_cols <- c("file_name", "datetime", "value_l", "value_r", grep("left_|right_", names(df), value = TRUE))
  df <- df[, relevant_cols]

  # mean ADI
  df <- df |>
    rowwise() |>
    mutate(adi = mean(c(value_l, value_r))) |>
    select(-value_l, -value_r)


  # Convert wide format to long format, dynamically identifying frequency bands
  df_long <- df |>
    pivot_longer(cols = starts_with("left_") | starts_with("right_"),
                 names_to = c(".value", "frequency_band"),
                 names_pattern = "(left|right)_(.*)") |>
    rowwise() |>
    mutate(proportion = mean(c(left, right))) |>
    select(-left, -right)

  # Ordering frequency bands based on their appearance in the data frame
  frequency_order <- unique(df_long$frequency_band)
  df_long$frequency_band <- factor(df_long$frequency_band, levels = frequency_order)

  colors <- c("white", "yellow", "orange", "red")

  # Heatmap plot
  heatmap <- ggplot(df_long, aes(x = datetime, y = frequency_band, fill = proportion)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors) +
    labs(title = plot.title, x = "", y = "Frequency Band (kHz)",
         subtitle = paste("Settings:", paste(settings, collapse = " / "))) +
    theme_classic() +
    scale_x_datetime(date_breaks = "1 hours", date_minor_breaks = "30 min", date_labels = "%H",
                     expand = c(0.01,0.01)) +
    theme(plot.margin = unit(c(1,1,0,1), "lines"), axis.ticks.x = element_line(size = 1))


  # ADI plot
  adi_plot <- ggplot(df_long, aes(x = datetime, y = adi)) +
    geom_line(color = "black", lwd = 0.5) +
    geom_smooth(color = "darkgreen", lwd = 1, alpha = .5, se = TRUE, method = 'loess', span = 0.1) +
    labs(x = "Time", y = "ADI") +
    theme_classic() +
    scale_x_datetime(date_breaks = "1 hours", date_minor_breaks = "30 min", date_labels = "%H",
                     expand = c(0.01,0.01)) +
    theme(plot.margin = unit(c(1,1,0,1), "lines"), axis.ticks.x = element_line(size = 1))

  # Combine the two plots
  combined_plot <- heatmap / adi_plot + plot_layout(heights = c(3, 1))


  return(combined_plot)
}


