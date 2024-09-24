#' Plot a Mean Spectrum with Flipped Axes
#'
#' This function is a wrapper around the `meanspec` function from the `seewave` package. It calculates 
#' the mean spectrum of a wave object, flips the axes so that amplitude is on the x-axis and frequency 
#' is on the y-axis, and fills the area under the curve with a black color.
#'
#' @param wave A wave object from the `tuneR` or `seewave` package, representing the audio signal to be analyzed.
#' @param freq.res A numeric value representing the desired frequency resolution. This controls the window length (`wl`) used in the `meanspec` function. Default is 50.
#' @param normalize Logical. Should the spectrum be normalized (i.e., expressed in decibels)? Passed to the `norm` argument of `meanspec`. Default is `FALSE`.
#' @param PSD Logical. If `TRUE`, compute the power spectral density (PSD). Default is `FALSE`.
#' @param PMF Logical. If `TRUE`, compute the probability mass function (PMF) of the spectrum. Default is `FALSE`.
#' @param plot.title A string for the plot's title. Default is an empty string (`""`).
#' @param verbose Logical. If `TRUE`, print information about the window length and the number of frequency bins. Default is `TRUE`.
#' 
#' @return The function returns a plot of the mean spectrum with flipped axes (amplitude on the x-axis and frequency on the y-axis) and the area under the curve filled with black.
#' 
#' @details The function computes the mean spectrum using the `meanspec` function from the `seewave` package. It calculates the window length (`wl`) based on the provided frequency resolution, ensuring it is an even number. The function then plots the amplitude on the x-axis and the frequency on the y-axis, with a filled black area under the curve.
#' 
#' @examples
#' \dontrun{
#'   # Load example sound data
#'   data(tico)
#'   # Plot the mean spectrum with flipped axes
#'   mean_spectrum_flip(tico, freq.res = 100, normalize = TRUE, plot.title = "Flipped Mean Spectrum")
#' }
#' 
#' @importFrom seewave meanspec
mean_spectrum_flip <- function(wave, freq.res = 50, normalize = FALSE, PSD=FALSE, PMF=FALSE, plot.title = "", verbose = TRUE){
  
  
  wl = wave@samp.rate / freq.res
  
  # Force it to be even
  if (wl %% 2 == 1) { wl <- wl + 1 }
  


  spec_data <- meanspec(wave, main = "Mean Spectrum", norm = normalize, PSD=PSD, PMF=PMF)
  
  if(verbose){
    cat("Window length:", wl, "samples \n",
        "Number of frequency bins:", nrow(spec_data), "\n")
  }
  
  # Extract frequency (x-axis) and amplitude (y-axis)
  frequency <- spec_data[, 1]
  amplitude <- spec_data[, 2]
  
  # Create the base plot (no line yet, just the axes)
  plot(amplitude, frequency, type = "n", xlab = "Amplitude", ylab = "Frequency (kHz)", 
       xlim = c(0, max(amplitude)), ylim = c(0, max(frequency)),
       xaxs = "i", yaxs = "i",
       main = plot.title)
  
  # Use polygon to fill the area under the curve
  polygon(c(amplitude, rep(0, length(amplitude))), c(frequency, rev(frequency)), col = "black", border = NA)
  

}
