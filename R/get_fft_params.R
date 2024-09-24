#' Calculate FFT Parameters
#'
#' This function calculates the time and frequency resolution based on the provided sampling rate and either window length, time resolution, or frequency resolution. Depending on which parameters are given, it calculates and returns the appropriate values for the others.
#'
#' @param samp.rate Numeric. The sampling rate of the audio signal in Hz.
#' @param win.length Numeric. The window length in samples. Can be `NULL` if time or frequency resolution is provided.
#' @param time.res Numeric. The time resolution in seconds. Can be `NULL` if window length or frequency resolution is provided.
#' @param freq.res Numeric. The frequency resolution in Hz. Can be `NULL` if window length or time resolution is provided.
#'
#' @return Prints the calculated time and frequency resolution, or window length and either time or frequency resolution depending on the provided inputs.
#' @export
#'
#' @examples
#' get_fft_params(samp.rate = 44100, win.length = 1024)
#' get_fft_params(samp.rate = 44100, time.res = 0.023)
#' get_fft_params(samp.rate = 44100, freq.res = 43)

get_fft_params <- function(samp.rate=NULL, win.length=NULL, time.res=NULL, freq.res=NULL){

  if(is.null(time.res) && is.null(freq.res)){
    time.res = win.length/samp.rate
    freq.res = samp.rate/win.length

    cat(" Time resolution:", round(time.res, 3), "s\n",
        "Frequency resolution:", round(freq.res,3), "Hz\n")

  }else if(is.null(win.length) && is.null(freq.res)){
    win.length = samp.rate * time.res
    freq.res = samp.rate/win.length

    cat(" Window length:", win.length, "samples\n",
        "Frequency resolution:", round(freq.res,3), "Hz\n")

  }else if(is.null(win.length) && is.null(time.res)){
    win.length = samp.rate/freq.res
    time.res = win.length/samp.rate

    cat(" Window length:", win.length, "samples\n",
        "Time resolution:", round(time.res, 3), "s\n")

  }

}


