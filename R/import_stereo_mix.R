#' Import a Stereo Wave File as Mono and Optionally Remove DC Offset
#'
#' This function reads a stereo wave file, converts it to mono by averaging both channels, 
#' and optionally removes the DC offset from the audio signal. The removal of the DC offset 
#' is controlled by the \code{remove.dc.off} parameter. This function uses the \pkg{tuneR} 
#' package to read and convert the wave file and the \pkg{seewave} package to remove the DC offset.
#'
#' @param audio.file A string specifying the path to the wave file.
#' @param remove.dc.off A logical parameter; if TRUE, the DC offset is removed from the wave object.
#'
#' @return Returns a wave object processed based on the input parameters.
#' @export
#' 
#' @import tuneR
#' @import seewave
#' 
#' @examples
#' # Import a wave file as mono with DC offset removal
#' wave <- import.as.mono("path/to/your/audiofile.wav")
#'
import_stereo_mix <- function(audio.file, remove.dc.off = TRUE){
  # audio.file should be in the working directory
  wave.object <- tuneR::readWave(audio.file)
  wave.object <- tuneR::mono(wave.object, "both")
  if(remove.dc.off){
    wave.object <- seewave::rmoffset(wave.object, output = "Wave")
  }
  return(wave.object)
}
